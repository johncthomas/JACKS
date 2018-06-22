import os
#import subprocess
#from subprocess import run
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def tabulate_jacksres(prefix):
    """Return a dict of tables giving the essentiality, fdr_pos and fdr_neg
    for each gene. Pulls the data from jacks .txt results files, so prefix
    should contain the path to the files if they aren't in os.getcwd()"""

    kwtab = dict(sep='\t', index_col=0)

    # othertab = pd.DataFrame(columns=("IC10","IC90","D14"), index=essen['D14'].index)
    # Tables produced by jacks have columns that are the groups
    genes = pd.read_table(prefix + '_gene_JACKS_results.txt', sep='\t', index_col=0)
    genes_index = sorted(genes.index)
    genes.reindex(genes_index)
    genesstd = pd.read_table(prefix + '_genestd_JACKS_results.txt', sep='\t', index_col=0)
    ps = genes / genesstd
    ps = ps.apply(norm.cdf)

    # multiindex DF for each experiment giving results
    sig_cols = pd.MultiIndex.from_product((ps.columns, ['essentiality', 'fdr_pos', 'fdr_neg']),
                                      names=('exp', 'stat'))
    sig_df = pd.DataFrame(index=genes_index, columns=sig_cols)

    for exp in ps.columns:
        sig_df.loc[:, (exp, 'fdr_neg')] = multipletests(ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'fdr_pos')] = multipletests(1 - ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'essentiality')] = genes[exp]
        # tab.loc[:, 'essentiality'] = genes[exp]
        # tab.loc[:, 'fdr_pos'] = fdr_pos
        # tab.loc[:, 'fdr_neg'] = fdr_neg
        # sig_tables[exp] = tab

    # get guide data, foldchange and efficacies
    guide_cols = pd.MultiIndex.from_product((ps.columns, ['foldchange','fold_std', 'eff', 'eff_std']),
                                           names=['exp','stat'])
    fchange_df = pd.DataFrame(columns=guide_cols)
    foldchange = pd.read_table(prefix + '_logfoldchange_means.txt', **kwtab)
    foldstd = pd.read_table(prefix + '_logfoldchange_std.txt', **kwtab)
    eff_tab = pd.read_table(prefix + '_grna_JACKS_results.txt', **kwtab)
    for exp in ps.columns:

        fchange_df.loc[:, (exp, 'foldchange')] = foldchange[exp]
        fchange_df.loc[:, (exp, 'fold_std')] = foldstd[exp]

    efficacies = pd.DataFrame(columns=('eff', 'eff_std'))
    efficacies.loc[:, 'eff'] = eff_tab['X1']
    efficacies.loc[:, 'eff_std'] = (eff_tab['X2'] - eff_tab['X1']**2)**0.5

    return sig_df, efficacies, fchange_df


def plot_volcano(esstab, savefn=None):
    """scatter plot of essentiality vs"""
    faces = dict(facecolors='none', edgecolors='b')
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    pos = esstab['essentiality'] > 0
    neg = ~pos

    faces = dict(facecolors='none', edgecolors='b')
    plt.scatter(esstab.loc[pos, 'essentiality'], esstab.loc[pos, 'fdr_pos'], **faces)
    plt.scatter(esstab.loc[neg, 'essentiality'], esstab.loc[neg, 'fdr_neg'], **faces)
    plt.yscale('log')

    plt.plot([min(esstab['essentiality']), max(esstab['essentiality'])],
             [0.05, 0.05], 'k--')

    # get lowest not zero and set zeroes to 1/10 that value
    min_pos = min(esstab.loc[esstab['fdr_pos'] > 0, 'fdr_pos'])
    min_neg = min(esstab.loc[esstab['fdr_neg'] > 0, 'fdr_neg'])
    print(min_neg, min_pos)
    for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
        esstab.loc[esstab[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

    plt.xlabel('Essentiality')
    plt.ylabel('FDR')

    plt.ylim(min(min_pos, min_neg) / 10)
    plt.gca().invert_yaxis()
    if not savefn:
        plt.show()
    else:
        plt.savefig(savefn)
    return fig, ax

if __name__ == '__main__':
    pass
    #
    # os.chdir('/Users/johnc.thomas/becca1')
    # scripts_dir = '/Users/johnc.thomas/Applications/guide_trimming_scripts'
    #
    # jkfn = './jacks/48hrWRef2_JACKS_results_full.pickle'
    #
    # # check efficacy file guides vs becs data
    # guide_eff = pd.read_table('/Users/johnc.thomas/Dropbox/crispr/gRNA_efficacy_w_X2.txt', sep='\t')
    #
    # with open(jkfn, 'rb') as f:
    #     jkres = jacks_results, cell_lines, grnas = pickle.load(f)
    #
    # with open('./jacks/all_counts.tsv') as f:
    #     counts = pd.read_table(f, index_col=0)
    # counts.head()
    #
    # ess_tables, guide_tables = tabulate_jacksres(*jkres)
    # # for group, tab in ess_tables.items():
    # #     print(tab.isna().sum())