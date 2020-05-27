#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
from glob import glob

import MOODS.parsers
import MOODS.scan
import MOODS.tools
import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as n2r
from Bio.SeqIO import parse as seqparse
from rpy2.rinterface import RRuntimeWarning
from scipy.sparse import csr_matrix, save_npz, vstack
from scipy.stats import mode
from sklearn.linear_model import LogisticRegression
from tqdm import tqdm

n2r.activate()

r = ro.r
r.library('glmpath')

### Scan and test and save motif profiles and stats for both permutated and real enhancers

# --- Argument parsing ---

parser = argparse.ArgumentParser(description="matches position weight matrices against DNA sequences")
parser.add_argument("-v", "--verbosity", action="count", default=1, help='verbose (-vv, -vvv for more)')

# input files
input_group = parser.add_argument_group("input files (at least one matrix and sequence file required)")
input_group.add_argument('-m','--matrices', metavar='M', nargs='+', action='store', dest='matrix_files', help='matrix files (count/frequency, will be converted to PWM before matching)', default = [])
input_group.add_argument('-S','--score-matrices', metavar='M', nargs='+', action='store', dest='lo_matrix_files', help='matrix files (PWM/other score matrix, will be matched directly)', default = [])
input_group.add_argument('-s','--sequences', metavar='S', type=str, action='store', dest='sequence_file', help='sequence file')
input_group.add_argument('-P','--phenocov', metavar='phenocov', type=str, action='store', dest='phenocov_file', help='phenotype and covariate file')

# thresholds
th_group = parser.add_argument_group("threshold selection (exactly one required)")
th_group.add_argument('-p','--p-value', metavar='p', action='store', dest='p_val', type=float, help='compute threshold from p-value')
th_group.add_argument('-t','--threshold', metavar='T', action='store', dest='t', type=float, help='use specified absolute threshold')
th_group.add_argument('-B','--best-hits', metavar='n', action='store', dest='max_hits', type=int, help='return at least the specified amount of best matches')

# output
out_group = parser.add_argument_group("output")
out_group.add_argument('-o','--outfile', metavar='outfile', action='store', dest='output_file', help='output filename key')


# behaviour
option_group = parser.add_argument_group("search and model behaviour (optional)")
option_group.add_argument('-R','--no-rc', dest='no_rc', action="store_true", help='disable matching versus reverse complement strand')
option_group.add_argument('--batch', dest='batch', action="store_true", help='do not recompute thresholds from p-values for each sequence separately (recommended when dealing with lots of short sequences)')
option_group.add_argument('--bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='bg', default = [0.25,0.25,0.25,0.25], help='background distribution for computing thresholds from p-value with --batch (default is 0.25 for all alleles)')
option_group.add_argument('--ps', metavar='p', action='store', dest='ps', type=float, help='total pseudocount added to each matrix column in log-odds conversion (default = 0.01)', default = 0.01)# bg
option_group.add_argument('--log-base', metavar='x', action='store', dest='log_base', type=float, help='logarithm base for log-odds conversion (default natural logarithm)')
option_group.add_argument('--lo-bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='lo_bg', default = [0.25,0.25,0.25,0.25], help='background distribution for log-odds conversion (default is 0.25 for all alleles)')
option_group.add_argument('--threshold-precision', metavar='x', action='store', dest='threshold_precision', type=float, help='specify the precision used for computing the thresholds from p-values (default = {})'.format(MOODS.tools.DEFAULT_DP_PRECISION), default = MOODS.tools.DEFAULT_DP_PRECISION)
option_group.add_argument('--permutation-size-multiplier', metavar='x', action='store', dest='size_multiplier', type=int, help='specify the number of permutation per region', default=10)


args = parser.parse_args()

# sanity check for parameters, apply effects

if len(args.matrix_files) == 0 and len(args.lo_matrix_files) == 0:
    print("{}: error: no matrix files given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
    sys.exit(1)
if args.sequence_file is None:
    print("{}: error: no sequence files given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
    sys.exit(1)
if (args.p_val is None) and (args.t is None) and (args.max_hits is None):
    print("{}: error: no threshold given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
    sys.exit(1)
if (args.p_val is not None and args.t is not None) or (args.t is not None and args.max_hits is not None) or (args.p_val is not None and args.max_hits is not None):
    print("{}: error: only one threshold specification allowed (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
    sys.exit(1)
if args.max_hits is not None and args.batch:
    print("{}: warning: ignoring --batch when used with -B".format(os.path.basename(__file__)), file=sys.stderr)
if args.t is not None and args.batch:
    print("{}: warning: --batch is redundant when used with -t".format(os.path.basename(__file__)), file=sys.stderr)


if args.log_base is not None and (args.log_base <= 1):
        print("{}: error: --log-base has to be > 1".format(os.path.basename(__file__)), file=sys.stderr)
        sys.exit(1)

def pfm_to_log_odds(filename):
    if args.log_base is not None:
        mat = MOODS.parsers.pfm_to_log_odds(filename, args.lo_bg, args.ps,
                                            args.log_base)
    else:
        mat = MOODS.parsers.pfm_to_log_odds(filename, args.lo_bg, args.ps)
    if len(mat) != 4:
        return False, mat
    else:
        return True, mat

def adm_to_log_odds(filename):
    if args.log_base is not None:
        mat = MOODS.parsers.adm_to_log_odds(filename, args.lo_bg, args.ps, 4, args.log_base)
    else:
        mat = MOODS.parsers.adm_to_log_odds(filename, args.lo_bg, args.ps, 4)
    if len(mat) != 16:
        return False, mat
    else:
        return True, mat


def pfm(filename):
    mat = MOODS.parsers.pfm(filename)
    if len(mat) != 4:
        return False, mat
    else:
        return True, mat

def adm(filename):
    mat = MOODS.parsers.adm_1o_terms(filename, 4)
    if len(mat) != 16:
        return False, mat
    else:
        return True, mat



def aggregate_results(matrix_names, results):
    # Get forward strand matching according to number of motifs
    # fr = results
    fr = results[:len(matrix_names)]
    # Get reverse strand matching according to number of moti
    rr = results[len(matrix_names):]

    # Calculate sum of motif matching score
    # results is a list of sums
    # print([r.score for r in fr[1]])
    results = [
        sum([0] + [r.score for r in fr[i]] + [r.score for r in rr[i]])
        for i in range(len(matrix_names))
    ]
    return np.array(results)

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

# --- Load matrices ---
if args.verbosity >= 1:
    print("{}: reading matrices ({} matrix files)".format(os.path.basename(__file__), len(args.matrix_files)+len(args.lo_matrix_files)), file=sys.stderr)
matrix_names = [os.path.basename(n) for n in args.matrix_files + args.lo_matrix_files]
matrices = []
matrices_rc = []
for filename in args.matrix_files:
    valid = False
    if filename[-4:] != '.adm': # let's see if it's pfm
        valid, matrix = pfm_to_log_odds(filename)
        if valid and args.verbosity >= 3:
            print("{}: read matrix file {} as pfm (converted to PWM)".format(os.path.basename(__file__), filename), file=sys.stderr)
    if not valid: # well maybe it is adm
        valid, matrix = adm_to_log_odds(filename)
        if valid and args.verbosity >= 3:
            print("{}: read matrix file {} as adm (converted to high-order PWM)".format(os.path.basename(__file__), filename), file=sys.stderr)
    if not valid:
        print("{}: error: could not parse matrix file {}".format(os.path.basename(__file__), filename), file=sys.stderr)
        sys.exit(1)
    else:
        matrices.append(matrix)
        matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))
for filename in args.lo_matrix_files:
    valid = False
    if filename[-4:] != '.adm': # let's see if it's pfm
        valid, matrix = pfm(filename)
        if valid and args.verbosity >= 3:
            print("{}: read matrix file {} as pfm".format(os.path.basename(__file__), filename), file=sys.stderr)
    if not valid: # well maybe it is adm
        valid, matrix = adm(filename)
        if valid and args.verbosity >= 3:
            print("{}: read matrix file {} as adm".format(os.path.basename(__file__), filename), file=sys.stderr)
    if not valid:
        print("{}: error: could not parse matrix file {}".format(os.path.basename(__file__), filename), file=sys.stderr)
        sys.exit(1)
    else:
        matrices.append(matrix)
        matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))

if args.no_rc:
    matrices_all = matrices
else:
    matrices_all = matrices + matrices_rc

# --- Scanning ---

if args.p_val is not None or args.t is not None:
    scanner = None
    # build scanner if using fixed threshold
    if args.t is not None:
        thresholds = [args.t for m in matrices_all]
        scanner = MOODS.scan.Scanner(7)
        bg = args.bg # used for search optimisation only
        scanner.set_motifs(matrices_all, bg, thresholds)
    # build scanner if using p-value and --batch
    if args.p_val is not None and args.batch:
        bg = args.bg
        if args.verbosity >= 1:
            print("{}: computing thresholds from p-value".format(os.path.basename(__file__)), file=sys.stderr)
        thresholds = [MOODS.tools.threshold_from_p_with_precision(m,bg,args.p_val,args.threshold_precision,4) for m in tqdm(matrices_all)]
        if args.verbosity >= 3:
            for (m, t) in zip(matrix_names, thresholds):
                print("{}: threshold for {} is {}".format(os.path.basename(__file__), m, t), file=sys.stderr)
        if args.verbosity >= 1:
            print("{}: preprocessing matrices for scanning".format(os.path.basename(__file__)), file=sys.stderr)
        scanner = MOODS.scan.Scanner(7)
        bg = args.bg
        scanner.set_motifs(matrices_all, bg, thresholds)


    phenocov = pd.read_csv(args.phenocov_file, sep="\t")
    phenocov = phenocov.sort_values("samples")
    seq = seqparse(args.sequence_file, "fasta")
    seq_table = pd.DataFrame([(i.id, i.description.split(' ')[1], str(i.seq)) for i in seq], columns=['samples', 'region', 'seq'])
    seq_table = seq_table[seq_table.samples.isin(phenocov.samples)].drop_duplicates()
    region_list = pd.unique(seq_table.region.sort_values())
    print("{}: scanning each regions".format(os.path.basename(__file__)), file=sys.stderr)


    y = phenocov.pheno.values
    covariates = phenocov.iloc[:,3:].values
    profiles = []
    regions = []
    stats = []
    sels = []
    coefs = []
    perm_regions = []
    perm_stats = []
    perm_sels = []

    for region in tqdm(region_list):
        profile = [aggregate_results(matrix_names, scanner.scan(s)) for s in seq_table[seq_table['region']==region].sort_values("samples").seq]
        profile = np.vstack(profile)
        profile = csr_matrix(profile - mode(profile, 0)[0])
        profiles.append(profile)

        stat, sel, coef = test(profile, y, covariates)
        stats.append(stat)
        sels.append(sel)
        coefs.append(coef)
        regions.append(region)

        perm_stat, perm_sel = perm_test(profile, y, covariates)
        perm_stats.append(perm_stat)
        perm_sels.append(perm_sel)
        perm_regions.append(region)

    stats = np.array(stats)
    sels = np.array(sels)
    coefs = np.array(coefs)
    regions = np.array(regions)
    perm_stats = np.array(perm_stats)
    perm_sels = np.array(perm_sels)
    perm_regions = np.array(perm_regions)
    np.savez(file=args.output_file+"_real_results.npz", stats=stats, sels=sels, coefs=coefs, regions=regions)
    np.savez(file=args.output_file+"_perm_results.npz", stats=perm_stats, sels=perm_sels, regions=perm_regions)
    save_npz(args.output_file+"_profiles.npz", vstack(profiles))