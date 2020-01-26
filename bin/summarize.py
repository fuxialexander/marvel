from glob import glob
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

stats = []
coefs = []
sels = []

perm_stats = []
perm_sels = []

for i in tqdm(glob(argv[1])):  # 'fantom_real/*.npz'
    with np.load(i, allow_pickle=True) as data:
        stats.append(data['stats'])
        coefs.append(data['coefs'])
        sels.append(data['sels'])

stats = np.concatenate(stats)
sels = np.concatenate(sels)
coefs = np.concatenate(coefs)

for i in tqdm(glob(argv[2])):  # 'fantom_perm/*.npz'
    with np.load(i) as data:
        perm_stats.append(data['stats'])
        perm_sels.append(data['sels'])

perm_stats = np.concatenate(perm_stats)
perm_stats = perm_stats.reshape(-1)

perm_stats_non_zero = np.concatenate((stats, perm_stats))
perm_stats_non_zero[perm_stats_non_zero == -1] = 0
perm_stats_non_zero = np.sort(perm_stats_non_zero)

p_stats = []
for i in stats:
    p_stats.append(1 - np.searchsorted(perm_stats_non_zero, i,
                                       'right') / float(len(perm_stats_non_zero)))
p_stats = np.array(p_stats)
p_stats[p_stats == 0] = 1/float(len(perm_stats_non_zero))


def get_lambda(p, definition='median'):
    '''
    Evaluates Lambda value
    :param p: distribution of p-values
    :param definition: definition of lambda
    :return:
    '''
    if definition == 'median':
        pm = np.median(p)
        Chi = chi2.ppf(1. - pm, 1)
        return Chi / chi2.ppf(0.5, 1)
    else:
        raise Exception(
            "Only 'median' definition of lambda is implemented at this moment.")


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
            slope, intercept, r_value, p_value, std_err = linregress(
                q_th, q_data)
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
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['savefig.dpi'] = 600
    plt.rcParams['figure.figsize'] = 3, 3
    mp = myqqplot([pvals], [''],
                  color=['#e34a33'],
                  fill_dens=[0.1],
                  error_type=error_type,
                  distribution=distribution,
                  title=None)
    if filename:
        mp.savefig(filename, dpi=600)


p_plot_qqplot(p_stats, argv[3])  # "fantom_qqplot_ex1.pdf"
