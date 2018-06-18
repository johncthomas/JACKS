import os
#import subprocess
#from subprocess import run
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


"""
for(i in 1:dim(metadata(expt)$jacks_w)[2]){ # for each experiment
        # calculate nominal p-values from z-scores
        m$jacks_neg_pval[,i] = stats::pnorm(m$jacks_w[,i] / m$jacks_sdw[,i])
        m$jacks_pos_pval[,i] = 1. - m$jacks_neg_pval[,i]
        # calculate fdr corrected p-values
        m$jacks_neg_fdr[,i] = stats::p.adjust(m$jacks_neg_pval[,i], 'fdr')
        m$jacks_pos_fdr[,i] = stats::p.adjust(m$jacks_pos_pval[,i], 'fdr')"""


def tabulate_jacksres(jacks_results, groups, grnas):
    """
    All data by group and by guide.
    Args:
        tabulate_jacksres(*pickle.load(<results_file>))
        groups as defined in the repmap file.

    Returns:
        dict of DF, keyed by groups.

    This is super fragile, any change to jacks.py formatting will mess it up
    need to review all changes when updating jacks
    """
    cols = ['guide', 'gene', 'condition', 'essentiality', 'ess_sd', 'fold_change',
            'tau', 'efficacy', 'eff_sd']
    ess_cols = ['gene', 'essentiality', 'ess_sd', 'p_pos', 'p_neg', 'fdr_pos', 'fdr_neg']
    jktabs = {}
    ess_tables = {}
    nguides = sum([len(gds) for gds in grnas.values()])
    ngenes = len(grnas)

    for group_i, group in enumerate(groups):
        ps = np.empty(ngenes, dtype=np.float64)
        rowi = 0  #
        guidetab = {k: np.empty(nguides, dtype=object, ) for k in cols}
        esstab = {k: np.empty(len(jacks_results.keys()), dtype=object) for k in ess_cols}
        for gene_i, (gene, res) in enumerate(jacks_results.items()):
            # non-ess values listed first by guide, then by condition
            fold_change, tau, eff, x2, ess, w2 = res
            # convert to standard deviation
            ess_sd = (w2 - ess ** 2) ** 0.5
            eff_sd = (x2 - eff ** 2) ** 0.5

            for guide_i, guide in enumerate(grnas[gene]):
                guidetab['guide'][rowi] = guide
                # dat['seq'][rowi]   = seq
                guidetab['fold_change'][rowi] = fold_change[guide_i][group_i]
                guidetab['tau'][rowi] = tau[guide_i][group_i]
                guidetab['efficacy'][rowi] = eff[guide_i]
                guidetab['eff_sd'][rowi] = eff_sd[guide_i]
                guidetab['essentiality'][rowi] = ess[group_i]
                guidetab['ess_sd'][rowi] = ess_sd[group_i]
                guidetab['condition'][rowi] = group
                guidetab['gene'][rowi] = gene
                # dat['p'][rowi] = p = norm.cdf(ess[group_i]/ess_sd[group_i])
                # ps[guide_i] = p
                rowi += 1

            # p of gene essentiality, for FDR later
            ps[gene_i] = p = norm.cdf(ess[group_i] / ess_sd[group_i])

            esstab['p_neg'][gene_i] = p
            esstab['p_pos'][gene_i] = 1 - p
            esstab['gene'][gene_i] = gene
            esstab['essentiality'][gene_i] = ess[group_i]
            esstab['ess_sd'][gene_i] = ess_sd[group_i]

        # get false discovery rate
        passed, neg_adj_p, _, _ = multipletests(esstab['p_neg'], method='fdr_bh')
        esstab['fdr_neg'] = neg_adj_p
        esstab['fdr_pos'] = multipletests(esstab['p_pos'], method='fdr_bh')[1]

        # convert dicts to pd.DF
        datdf = pd.DataFrame(guidetab, columns=cols)
        jktabs[group] = datdf.set_index('guide')
        essdf = pd.DataFrame(esstab, columns=ess_cols)
        ess_tables[group] = essdf.set_index('gene')

    return ess_tables, jktabs


def plot_jkstat(table, sampx, sampy, statx, staty=None):
    if not staty:
        staty = statx
    datx, daty = table[sampx][statx], table[sampy][staty]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.scatter(datx, daty, facecolors='none', edgecolors='b')
    axismax = int(max(max(datx), max(daty)))
    axismin = int(min(min(datx), min(daty)))

    # diag
    plt.plot(
        range(axismin, axismax),
        range(axismin, axismax),
        'k--'
    )

    plt.plot(
        np.zeros(-axismin + axismax),
        range(axismin, axismax),
        'g--'
    )

    plt.plot(
        range(axismin, axismax),
        np.zeros(-axismin + axismax),
        'g--'
    )

    plt.xlabel(sampx + ' ' + statx)
    plt.ylabel(sampy + ' ' + staty)
    plt.show()


def plot_volcano(esstab):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    pos = esstab['essentiality'] > 0
    neg = ~pos

    faces = dict(facecolors='none', edgecolors='b')
    plt.scatter(esstab.loc[pos, 'essentiality'], esstab.loc[pos, 'fdr_pos'], **faces)
    plt.scatter(esstab.loc[neg, 'essentiality'], esstab.loc[neg, 'fdr_neg'], **faces)
    plt.yscale('log')

    # get lowest not zero and set zeroes to 1/10 that value
    min_pos = min(esstab.loc[esstab['fdr_pos'] > 0, 'fdr_pos'])
    min_neg = min(esstab.loc[esstab['fdr_neg'] > 0, 'fdr_neg'])

    for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
        esstab.loc[esstab[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

    plt.ylim(min(min_pos, min_neg) / 10)
    plt.gca().invert_yaxis()
    plt.show()

#if __name__ == '__main__':
if True:
    os.chdir('/Users/johnc.thomas/becca1')
    scripts_dir = '/Users/johnc.thomas/Applications/guide_trimming_scripts'

    jkfn = './jacks/48hrWRef2_JACKS_results_full.pickle'

    # check efficacy file guides vs becs data
    guide_eff = pd.read_table('/Users/johnc.thomas/Dropbox/crispr/gRNA_efficacy_w_X2.txt', sep='\t')

    with open(jkfn, 'rb') as f:
        jkres = jacks_results, cell_lines, grnas = pickle.load(f)

    with open('./jacks/all_counts.tsv') as f:
        counts = pd.read_table(f, index_col=0)
    counts.head()

    ess_tables, guide_tables = tabulate_jacksres(*jkres)
    # for group, tab in ess_tables.items():
    #     print(tab.isna().sum())