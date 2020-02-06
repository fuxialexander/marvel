#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
from glob import glob

import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as n2r
from rpy2.rinterface import RRuntimeWarning
from scipy.sparse import csr_matrix, save_npz, vstack
from scipy.stats import mode
from sklearn.linear_model import LogisticRegression
from tqdm import tqdm

n2r.activate()

r = ro.r
r.library('glmpath')


parser = argparse.ArgumentParser(description="Scan motif profiles")
parser.add_argument("-v", "--verbosity", action="count", default=1, help='verbose (-vv, -vvv for more)')

# input file
input_group = parser.add_argument_group("input files (at least one matrix and sequence file required)")
input_group.add_argument('-P','--phenocov', metavar='phenocov', type=str, action='store', dest='phenocov_file', help='phenotype and covariate file')
# option
option_group = parser.add_argument_group("Test and permutation behaviour (optional)")
option_group.add_argument('--permutation-size-multiplier', metavar='x', action='store', dest='size_multiplier', type=int, help='specify the number of permutation per region', default=10)
# output file
out_group = parser.add_argument_group("output")
out_group.add_argument('-o','--outfile', metavar='outfile', action='store', dest='output_file', help='output filename key')

args = parser.parse_args()

def normalize(x):
    x = np.array(x - x.mean(0))
    x_std = np.array(x.std(0))
    x[:, x_std != 0 ] = x[:, x_std != 0] / x_std[x_std != 0]
    return x

def perm_test(x, y, covariates):
    x = normalize(x)
    if x.std()==0:
        return np.repeat(-1, args.size_multiplier), np.repeat(0, args.size_multiplier)
    stats_j = np.zeros(args.size_multiplier)
    sels_j = np.zeros(args.size_multiplier)
    logistic_alt = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
    logistic_null = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')

    covariates = normalize(covariates)
    for i in range(args.size_multiplier):
        y_perm = np.random.permutation(y)
        nullmodel = logistic_null.fit(covariates, y_perm)
        Lnullmodel = np.log(nullmodel.predict_proba(covariates))

        try:
            reg = r['glmpath'](x, y_perm, **{'min.lambda':1e-2, 'max.steps':10, 'max.vars':10})
            coefs_ = np.asanyarray(reg.rx2('b.predictor'))
            coefs = coefs_[-1][1:]
            sel = np.where(coefs != 0)[0]
            if sel.shape[0] == 0:
                stat = -1
                nsel=0
            else:
                xalt = np.concatenate((x[:, sel], covariates),
                                        axis=1).astype(float)
                altmodel = logistic_alt.fit(xalt, y_perm)
                _Laltmodel = altmodel.predict_proba(xalt)
                if (_Laltmodel<=0).any():
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
            stat=-1
            nsel = 0

        stats_j[i] = stat
        sels_j[i] = nsel

    return stats_j, sels_j

def test(x, y, covariates):
    x = normalize(x)
    if x.std()==0:
        return -1, None, None
    covariates = normalize(covariates)
    logistic_alt = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
    logistic_null = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
    nullmodel = logistic_null.fit(covariates, y)
    Lnullmodel = nullmodel.predict_log_proba(covariates)
    try:
        reg = r['glmpath'](x, y, **{'min.lambda':1e-2, 'max.steps':10, 'max.vars':10})
        coefs_ = np.asanyarray(reg.rx2('b.predictor'))
        coefs = coefs_[-1][1:]
        sel = np.where(coefs != 0)[0]
        if sel.shape[0] == 0 :
            # not selected
            sel = None
            coefs = None
            stat = -1
        else:
            coefs = coefs[sel]
            xalt = np.concatenate((x[:, sel], covariates), axis=1).astype(float)
            altmodel = logistic_alt.fit(xalt, y)
            _Laltmodel = altmodel.predict_proba(xalt)
            if (_Laltmodel<=0).any():
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

    phenocov = pd.read_csv(args.phenocov_file, sep="\t")
    phenocov = phenocov.sort_values("samples")

    print("{}: scanning each regions".format(os.path.basename(__file__)), file=sys.stderr)


    y = phenocov.pheno.values
    covariates = phenocov.iloc[:,3:].values
    profiles = []
    stats = []
    sels = []
    coefs = []
    perm_stats = []
    perm_sels = []

    for profile in tqdm(profiles):
        stat, sel, coef = test(profile, y, covariates)
        stats.append(stat)
        sels.append(sel)
        coefs.append(coef)

        perm_stat, perm_sel = perm_test(profile, y, covariates)
        perm_stats.append(perm_stat)
        perm_sels.append(perm_sel)

    stats = np.array(stats)
    sels = np.array(sels)
    coefs = np.array(coefs)
    perm_stats = np.array(perm_stats)
    perm_sels = np.array(perm_sels)
    np.savez(file=args.output_file+"_real_results.npz", stats=stats, sels=sels, coefs=coefs)
    np.savez(file=args.output_file+"_perm_results.npz", stats=perm_stats, sels=perm_sels)
