#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from glob import glob

import numpy as np
from scipy.sparse import load_npz, save_npz, vstack
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

# argument
parser = argparse.ArgumentParser()
parser.add_argument("region_type", default="",type=str,
                    help="Type of region, e.g. enhancer or promoter")
parser.add_argument("result_path", default="",type=str,
                    help="path of chunk testing results")
parser.add_argument("output_path", default="./", type=str, help="output directory")
args = parser.parse_args()

# load profiles
profiles = []

for i in sorted(glob(args.result_path + "*_profiles.npz")):
    profiles.append(load_npz(i))

profiles = vstack(profiles)

# load results
regions = []
stats = []
sels = []
coefs = []

for i in sorted(glob(args.result_path + "*real_results.npz")):
    regions.append(np.load(i, 'r', True)['regions'])
    stats.append(np.load(i, 'r', True)['stats'])
    # print([np.array(i) for i in np.load(i, 'r', True,)['sels']])
    sels.append([np.array(i) for i in np.load(i, 'r', True,)['sels']])
    coefs.append([np.array(i) for i in np.load(i, 'r', True,)['coefs']])

regions = np.concatenate(regions)
stats = np.concatenate(stats)
sels = np.array(sum(sels, []))
coefs = np.array(sum(coefs, []))

# load permutation stats
perm_stats = []
for i in sorted(glob(args.result_path + "*perm_results.npz")):
    perm_stats.append(np.load(i, 'r', True)['stats'])

perm_stats = np.concatenate(perm_stats)

# compute P-values
all_stats = np.concatenate((stats, perm_stats.reshape(-1)))
all_stats[all_stats == -1] = 0
all_stats = np.sort(all_stats)
pvals = []
for i in stats:
    pvals.append(1 - np.searchsorted(all_stats, i,
                                     'right') / float(len(all_stats)))
pvals = np.array(pvals)
pvals[pvals == 0] = 1/float(len(all_stats))

_, fdrs = fdrcorrection0(pvals)

# save results
np.savez(args.output_path + args.region_type + '_results.npz',
         regions=regions,
         pvals=pvals,
         fdrs=fdrs,
         stats=stats,
         perm_stats=perm_stats,
         sels=sels,
         coefs=coefs
         )

# save profiles
save_npz(args.output_path + args.region_type + '_profiles.npz', profiles)
