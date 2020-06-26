#!/usr/bin/env python
import argparse
import os
import sys
from glob import glob
from os.path import basename

import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as n2r
from rpy2.rinterface import RRuntimeWarning
from scipy.sparse import csr_matrix, load_npz, save_npz, vstack
from scipy.stats import mode
from sklearn.linear_model import LogisticRegression
from tqdm import tqdm
n2r.activate()

r = ro.r
r.library('glmpath')

def get_pe(pe_pair_file):
    ep_dis = dict()
    ep_pair = dict()
    with open(pe_pair_file, 'r') as pairs:
        for pair in pairs:
            p, e, d = pair.strip("\n").split("\t")
            ep_dis[p + e] = int(d)
            if p in ep_pair:
                ep_pair[p].append(e)
            else:
                ep_pair[p] = []
                ep_pair[p].append(e)
    return ep_dis, ep_pair


def get_pg(pg_pair_file):
    pg_pair = dict()
    glist = []
    with open(pg_pair_file, 'r') as f:
        for line in f:
            loc, gid = line.strip('\n').split('\t')
            glist.append(gid)
            if gid in pg_pair:
                pg_pair[gid].append(loc)
            else:
                pg_pair[gid] = []
                pg_pair[gid].append(loc)
    return pg_pair, np.unique(np.array(glist))


def get_weight(distance, breaks, weight):
    diff = breaks.values - abs(distance)
    print(diff)
    print(np.where(diff > 0)[0])
    print(diff[diff > 0].argmin())
    if sum(diff > 0) > 0:
        return weight[np.where(diff > 0)[0][diff[diff > 0].argmin()]]
    else:
        return 0


def get_distance(p, e, ep_dis):
    if (p + e) in ep_dis:
        return ep_dis[p + e]
    else:
        return 1000000 - 1


def get_gene(pg_pair, glist, ep_pair, sample_count, motif_count):
    xs = []
    for i, g in enumerate(glist):
        x = np.zeros((sample_count, motif_count))
        promoters = np.array(pg_pair[g])
        enhancers = [
            np.array(ep_pair[p]) for p in promoters if (p in ep_pair)
        ]
        if len(enhancers) > 0:
            enhancers = np.unique(np.concatenate(enhancers))
        else:
            continue
        pids = np.array([
            np.where(promoter_results['regions'] == promoter)[0][0]
            for promoter in promoters
        ])
        for enhancer in enhancers:
            eid = np.where(enhancer_results['regions'] == enhancer)[0][0]
            weight = get_weight(
                get_distance(promoters[0], enhancer, ep_dis), distance_breaks,
                distance_weight)
            x += weight * enhancer_profiles[eid * sample_count:
                                            (eid + 1) * sample_count, :].todense()
        x += distance_weight[0] * sum([
            promoter_profiles[pid * sample_count:
                              (pid + 1) * sample_count, :].todense() for pid in pids
        ])
        x = x / float(len(promoters) + len(enhancers))
        xs.append(x)

    return xs


def normalize(x):
    x = np.array(x - x.mean(0))
    x_std = np.array(x.std(0))
    x[:, x_std != 0] = x[:, x_std != 0] / x_std[x_std != 0]
    return x


def perm_test(x, y, covariates):
    x = normalize(x)
    if x.std() == 0:
        return np.repeat(-1, args.size_multiplier), np.repeat(0, args.size_multiplier)
    stats_j = np.zeros(args.size_multiplier)
    sels_j = np.zeros(args.size_multiplier)
    logistic_alt = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
    logistic_null = LogisticRegression(
        C=1e42, penalty='l2', solver='liblinear')

    covariates = normalize(covariates)
    for i in range(args.size_multiplier):
        y_perm = np.random.permutation(y)
        nullmodel = logistic_null.fit(covariates, y_perm)
        Lnullmodel = np.log(nullmodel.predict_proba(covariates))

        try:
            reg = r['glmpath'](
                x, y_perm, **{'min.lambda': 1e-2, 'max.steps': 10, 'max.vars': 10})
            coefs_ = np.asanyarray(reg.rx2('b.predictor'))
            coefs = coefs_[-1][1:]
            sel = np.where(coefs != 0)[0]
            if sel.shape[0] == 0:
                stat = -1
                nsel = 0
            else:
                xalt = np.concatenate((x[:, sel], covariates),
                                      axis=1).astype(float)
                altmodel = logistic_alt.fit(xalt, y_perm)
                _Laltmodel = altmodel.predict_proba(xalt)
                if (_Laltmodel <= 0).any():
                    stat = -1
                    nsel = len(sel)
                else:
                    Laltmodel = np.log(_Laltmodel)
                    _y = np.vstack([1 - y_perm, y_perm])
                    stat = 2 * (np.matmul(_y, Laltmodel - Lnullmodel).trace())
                    nsel = len(sel)
                    if np.isnan(stat):
                        stat = -1
                        nsel = len(sel)

        except rpy2.rinterface_lib.embedded.RRuntimeError:
            stat = -1
            nsel = 0

        stats_j[i] = stat
        sels_j[i] = nsel

    return stats_j, sels_j


def test(x, y, covariates):
    x = normalize(x)
    if x.std() == 0:
        return -1, None, None
    covariates = normalize(covariates)
    logistic_alt = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
    logistic_null = LogisticRegression(
        C=1e42, penalty='l2', solver='liblinear')
    nullmodel = logistic_null.fit(covariates, y)
    Lnullmodel = nullmodel.predict_log_proba(covariates)
    try:
        reg = r['glmpath'](
            x, y, **{'min.lambda': 1e-2, 'max.steps': 10, 'max.vars': 10})
        coefs_ = np.asanyarray(reg.rx2('b.predictor'))
        coefs = coefs_[-1][1:]
        sel = np.where(coefs != 0)[0]
        if sel.shape[0] == 0:
            # not selected
            sel = None
            coefs = None
            stat = -1
        else:
            coefs = coefs[sel]
            xalt = np.concatenate(
                (x[:, sel], covariates), axis=1).astype(float)
            altmodel = logistic_alt.fit(xalt, y)
            _Laltmodel = altmodel.predict_proba(xalt)
            if (_Laltmodel <= 0).any():
                stat = -1
            else:
                Laltmodel = np.log(_Laltmodel)
                _y = np.vstack([1 - y, y])
                stat = 2 * (np.matmul(_y, Laltmodel - Lnullmodel).trace())
                if np.isnan(stat):
                    stat = -1

    except rpy2.rinterface_lib.embedded.RRuntimeError:

        sel = None
        coefs = None
        stat = -1

    return stat, sel, coefs


# argument
parser = argparse.ArgumentParser()
parser.add_argument(
    '-w', '--weights', type=str, default='enhancer_promoter_distance_weight.csv',
    help='CSV file that specifying enhancer-promoter distance weights. (default: enhancer_promoter_distance_weight.csv)')

parser.add_argument(
    '-r', '--results', type=str, default=None, required=True,
    help='The results folder of region-based test. (default: None)')

parser.add_argument(
    '-e', '--promoter-enhancer-pair', type=str, default='promoter_enhancer_pair.txt',
    help='File specifying promoter enhancer relationship (default: Promoter_enhancer_pair.txt)')

parser.add_argument(
    '-g', '--promoter-gene-pair', type=str, default='promoter_gene_pair.txt',
    help='File specifying promoter gene relationship (default: promoter_gene_pair.txt)')

parser.add_argument(
    '-m', '--motifs', type=str, default=None, required=True,
    help='Path to motif PWMs (default: None)')

parser.add_argument(
    '-p', '--phenocov', type=str, default='pheno_covar.txt',
    help='Phenotype and covariate file (default: pheno_covar.txt)')

parser.add_argument('-n', '--size-multiplier', action='store',  type=int, help='Specify the number of permutation per region (default: 10)', default=10)

parser.add_argument('-o', "--output-file", default="gene",
                    type=str, help="output directory")
args = parser.parse_args()

enhancer_results = np.load(args.results + '/enhancer_results.npz', allow_pickle=True)
promoter_results = np.load(args.results + '/promoter_results.npz', allow_pickle=True)
enhancer_profiles = load_npz(args.results + '/enhancer_profiles.npz')
promoter_profiles = load_npz(args.results + '/promoter_profiles.npz')
ep_dis, ep_pair = get_pe(args.promoter_enhancer_pair)
pg_pair, glist = get_pg(args.promoter_gene_pair)
distance_breaks = pd.read_csv(args.weights)['breaks']
distance_weight = pd.read_csv(args.weights)['weights']

phenocov = pd.read_csv(args.phenocov, sep="\t")
phenocov = phenocov.sort_values("samples")
sample_count = phenocov['samples'].shape[0]
motif_count = enhancer_profiles.shape[1]
matrix_names = np.array([basename(x)[0:-20]
                         for x in glob(args.motifs+'/*.pwm')])
y = phenocov.pheno.values
covariates = phenocov.iloc[:, 3:].values
profiles = []
regions = []
stats = []
sels = []
coefs = []
perm_regions = []
perm_stats = []
perm_sels = []

for g in tqdm(glist):
    regions.append(g)
    perm_regions.append(g)
    x = np.zeros((sample_count, motif_count))
    promoters = np.array(pg_pair[g])
    enhancers = [
        np.array(ep_pair[p]) for p in promoters if p in ep_pair
    ]
    if len(enhancers) > 0:
        enhancers = np.unique(np.concatenate(enhancers))
    else:
        continue

    print(promoter_results['regions'], promoters)
    pids = np.array([
        np.where(promoter_results['regions'] == promoter)[0][0]
        for promoter in promoters
    ])
    for enhancer in enhancers:
        eid = np.where(enhancer_results['regions'] == enhancer)[0][0]
        weight = get_weight(
            get_distance(promoters[0], enhancer, ep_dis), distance_breaks,
            distance_weight)
        x += weight * enhancer_profiles[eid * sample_count:
                                        (eid + 1) * sample_count, :].todense()
    x += distance_weight[0] * sum([
        promoter_profiles[pid * sample_count:
                          (pid + 1) * sample_count, :].todense() for pid in pids
    ])
    x = x / float(len(promoters) + len(enhancers))
    x = csr_matrix(x)
    profiles.append(x)
    stat, sel, coef = test(x, y, covariates)
    stats.append(stat)
    sels.append(sel)
    coefs.append(coef)
    perm_stat, perm_sel = perm_test(x, y, covariates)
    perm_stats.append(perm_stat)
    perm_sels.append(perm_sel)

stats = np.array(stats)
sels = np.array(sels)
coefs = np.array(coefs)
perm_stats = np.array(perm_stats)
perm_sels = np.array(perm_sels)
np.savez(file=args.output_file+"_real_results.npz",
         stats=stats, sels=sels, coefs=coefs, regions=regions)
np.savez(file=args.output_file+"_perm_results.npz",
         stats=perm_stats, sels=perm_sels, regions=perm_regions)
save_npz(args.output_file+"_profiles.npz", vstack(profiles))
