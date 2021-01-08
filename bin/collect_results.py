#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from glob import glob
import pandas as pd
import numpy as np
from scipy.sparse import load_npz, save_npz, vstack
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

# argument
parser = argparse.ArgumentParser()
parser.add_argument('-n', "--region-name", default="", type=str, required=True,
                    help="Type of region, e.g. enhancer or promoter")
parser.add_argument('-r', "--result-path", default="", type=str, required=True,
                    help="path of chunk testing results")
parser.add_argument('-a', "--region_annotation", default="", type=str, required=False, help="path of region annotation file")
parser.add_argument('-o', "--output_path", default="./", required=True,type=str, help="output directory")
args = parser.parse_args()

if args.region_annotation:
    regions_annotations = pd.read_csv(args.region_annotation, header=None)
# load profiles
profiles = []
if args.region_name=="gene":
    name_prefix = args.result_path + "promoter_gene_pair"
else:
    name_prefix = args.result_path + args.region_name

for i in sorted(glob(name_prefix + "*_profiles.npz")):
    profiles.append(load_npz(i))

profiles = vstack(profiles)

# load results
regions = []
stats = []
sels = []
coefs = []

for i in sorted(glob(name_prefix + "*real_results.npz")):
    regions.append(np.load(i, 'r', True)['regions'])
    stats.append(np.load(i, 'r', True)['stats'])
    sels.append([np.array(i) for i in np.load(i, 'r', True,)['sels']])
    coefs.append([np.array(i) for i in np.load(i, 'r', True,)['coefs']])

regions = np.concatenate(regions)
stats = np.concatenate(stats)
sels = np.array(sum(sels, []))
coefs = np.array(sum(coefs, []))

# load permutation stats
perm_stats = []
for i in sorted(glob(name_prefix + "*perm_results.npz")):
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
np.savez(args.output_path + args.region_name + '_results.collected.npz',
         regions=regions,
         pvals=pvals,
         fdrs=fdrs,
         stats=stats,
         perm_stats=perm_stats,
         sels=sels,
         coefs=coefs
         )

# save profiles
save_npz(args.output_path + args.region_name + '_profiles.collected.npz', profiles)
