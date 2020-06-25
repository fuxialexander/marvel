#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse
from glob import glob
from os.path import basename
from sys import argv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from assocplots import qqplot
from scipy.sparse import load_npz
from scipy.stats import beta, binom, chi2, linregress, norm
from scipy.stats.mstats import mquantiles
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from tqdm import tqdm

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.figsize'] = 3, 3
# %%
# argument
# parser = argparse.ArgumentParser()
# parser.add_argument("results", type=str, help="results npz file")
# parser.add_argument("profiles", default="", type=str,
#                     help="profiles npz file")
# parser.add_argument("phenocov", default="", type=str,
#                     help="Phenotype and covariates file")
# parser.add_argument("motif_dir", default="", type=str,
#                     help="Path to motif PWMs")
# args = parser.parse_args()

# Load results
profiles = load_npz("../results/test/enhancer_profiles.npz")
results = np.load("../results/test/enhancer_results.npz", mmap_mode='r', allow_pickle=True)

pvals = results['pvals']
fdrs = results['fdrs']
sels = results['sels']
coefs = results['coefs']
regions = results['regions']
# %%
# Load phenotype and covariates
pheno = pd.read_table("../test_input/samples/pheno_covar.txt", header=None)
y = pheno.values[:, 1]
covariates = pheno.values[:, 2:]
N_SAMPLE = pheno.shape[0]
# %%
# Load motif and annotation
matrix_names = [basename(x)[0:-20] for x in glob("../results/motifs")]
matrix_names_arr = np.array(matrix_names)
tfs = pd.read_csv("motif_annotation.txt", sep="\t")
# tfs_id = np.unique(np.concatenate(
    # [np.where(matrix_names_arr == m)[0] for m in tfs['Motif'].values]))
# tfs_conversion = dict(zip(tfs['Motif'].values, tfs['Symbol'].values))
# %%
# Preparation
loga = LogisticRegression(C=1e42, penalty='l2', solver='liblinear')
# %%
def get_tfs(ids):
    if ids is not None:
        return ','.join([tfs_conversion[matrix_names_arr[i]] if matrix_names_arr[i] in tfs_conversion else matrix_names_arr[i] for i in ids])

def get_tfs_coefs(i, coefs):
    if coefs[i] is not None:
        return ','.join(map(str, coefs[i].round(3)))
    else:
        return ''

def get_predict_from_id(id, sels, profiles, y):
    if sels[id] is not None:
        x_sel = profiles[id*N_SAMPLE: (id+1)*N_SAMPLE][:, sels[id]].todense()
        xalt = np.concatenate((x_sel, covariates), axis=1).astype(float)
    else:
        xalt = covariates.astype(float)

    return loga.fit(xalt, y).predict_proba(xalt)

def get_region_annots(regions, annots_file, sep='\t'):
    # translate regions location to annotations (gene names, hg38 coordinates, etc)
    annots = {}
    if annots_file:
        with open(annots_file, 'r') as f:
            for line in f:
                key, val = line.strip().split(sep)
                annots[key] = val
        annots = {}
    else:
        for r in regions:
            annots[r] = r

    return annots

def print_regions_table(fname, ids, regions, regions_annot, pvals, fdrs, sels, coefs, profiles):
    locs = [regions[i].decode('UTF-8') for i in ids]  # Region coordinates
    # Annotation for each region, e.g. gene name for promoters
    annots = [regions_annot[n] for n in locs]
    motifs = [get_tfs(sels[i]) for i in ids]  # Altered TF motifs
    coeffs = [get_tfs_coefs(i, coefs)
              for i in ids]  # Coefficient of corresponding TFs
    prediction = [get_predict_from_id(i, sels, profiles, y) for i in tqdm(ids)]
    aurocs = [roc_auc_score(y, i[:, 1]) for i in prediction]
    table = {
        'locations': locs,
        'annotations': annots,
        'p-vals': [pvals[i] for i in ids],
        'fdr': [fdrs[i] for i in ids],
        'motifs': motifs,
        'coefs': coeffs,
        'auroc': aurocs
    }
    pdtable = pd.DataFrame.from_dict(table)
    pdtable.to_csv(fname+'.table.tsv', sep='\t')
    pdtable.to_excel(fname+'.table.xlsx')
    return pdtable

def myqqplot(data, labels, n_quantiles=100, alpha=0.95, error_type='theoretical', distribution='binomial', log10conv=True, color=['k', 'r', 'b'], fill_dens=[0.1, 0.1, 0.1], type='uniform', title='title'):
    '''
    Function for plotting Quantile Quantile (QQ) plots with confidence interval (CI)
    :param data: NumPy 1D array with data
    :param labels:
    :param type: type of the plot
    :param n_quantiles: number of quntiles to plot
    :param alpha: confidence interval
    :param distribution: beta/normal/binomial -- type of the error estimation. Most common in the literature is 'beta'.
    :param log10conv: conversion to -log10(p) for the figure
    :return: nothing
    '''
    xmax = 0
    ymax = 0
    if type == 'uniform':
        # we expect distribution from 0 to 1
        for j in range(len(data)):
            # define quantiles positions:
            q_pos = np.concatenate([np.arange(
                99.)/len(data[j]), np.logspace(-np.log10(len(data[j]))+2, 0, n_quantiles)])
            print(q_pos.shape)
            print(q_pos)
            # define quantiles in data
            # linear interpolation
            q_data = mquantiles(data[j], prob=q_pos,
                                alphap=0, betap=1, limit=(0, 1))
            # define theoretical predictions
            q_th = q_pos.copy()
            # evaluate errors
            q_err = np.zeros([len(q_pos), 2])
            if np.sum(alpha) > 0:
                for i in range(0, len(q_pos)):
                    if distribution == 'beta':
                        q_err[i, :] = beta.interval(alpha, len(
                            data[j])*q_pos[i], len(data[j]) - len(data[j])*q_pos[i])

                        print(q_err.shape)
                    elif distribution == 'binomial':
                        q_err[i, :] = binom.interval(
                            alpha=alpha, n=len(data[j]), p=q_pos[i])
                    elif distribution == 'normal':
                        q_err[i, :] = norm.interval(alpha, len(
                            data[j])*q_pos[i], np.sqrt(len(data[j])*q_pos[i]*(1.-q_pos[i])))
                    else:
                        print('Distribution is not defined!')
                    q_err[i, q_err[i, :] < 0] = 1e-15
                if (distribution == 'binomial') | (distribution == 'normal'):
                    q_err /= 1.0*len(data[j])
                    for i in range(0, 100):
                        q_err[i, :] += 1e-15
            # slope, intercept, r_value, p_value, std_err = linregress(
            #     q_th, q_data)
            plt.plot(-np.log10(q_th[n_quantiles-1:]), -
                     np.log10(q_data[n_quantiles-1:]), '.', color=color[j])
            plt.plot(-np.log10(q_th[:n_quantiles]), -np.log10(
                q_data[:n_quantiles]), '.', color=color[j], label=labels[j])
            xmax = np.max([xmax, - np.log10(q_th[1])])
            ymax = np.max([ymax, - np.log10(q_data[0])])
            # print(- np.log10(q_th[:]))
            if np.sum(alpha) > 0:
                if error_type == 'experimental':
                    plt.fill_between(-np.log10(q_th), -np.log10(q_data/q_th*q_err[:, 0]), -np.log10(
                        q_data/q_th*q_err[:, 1]), color=color[j], alpha=fill_dens[j], label='%1.3f CI' % alpha)
        if np.sum(alpha) > 0:
            if error_type == 'theoretical':
                plt.fill_between(-np.log10(q_th), -np.log10(q_err[:, 0]), -np.log10(
                    q_err[:, 1]), color='#fee8c8', label='%1.3f CI' % alpha)
    plt.legend(loc=4)
    plt.xlabel('Expected -log10 P')
    plt.ylabel('Observed -log10 P')
    plt.plot([0, 100], [0, 100], '--k')
    plt.xlim([0, np.ceil(xmax)])
    plt.ylim([0, np.ceil(ymax*1.05)])
    plt.title(title)
    plt.tight_layout()
    return plt
    # return q_data, q_th, q_err

def p_plot_qqplot(pvals, filename, distribution='beta', error_type='theoretical'):
    mp = myqqplot([pvals], [''],
                  color=['#e34a33'],
                  fill_dens=[0.1],
                  error_type=error_type,
                  distribution=distribution,
                  title=None)
    if filename:
        mp.savefig(filename, dpi=600)

# p_plot_qqplot(pvals, argv[1]+"_qqplot.pdf")  # "fantom_qqplot.pdf"


# %%
myqqplot([pvals], ['test'], error_type='theoretical', distribution='beta')

# %%
