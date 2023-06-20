import pandas as pd
import scanpy as sc
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


def process_gene_ranks(gene_ranks, key, shift_by = 1e-200):
    """ Process output from scanpy.tl.rank_genes_groups for volcano plotting 
        
    Keyword arguments:
    gene_ranks: dictionary of np.ndarrays stored in ad.uns['rank_genes_groups']
    key: string of the name of the comparison group or groups
    shift_by: float representing a very small floating point number to prevent log10p values smaller than minimum float representation. (default: 1e-200)
    
    Returns:
    DataFrame containing the p values, adjusted p-values, log2 fold change and -log10(p) values for each gene in the indicated comparison.
    """
    df = pd.DataFrame(index = gene_ranks['names'][key]) 
    df['pvals'] = pd.Series(data = np.asarray(gene_ranks['pvals'][key], dtype = float), index = df.index)
    df['pvals_adj'] =  pd.Series(data = np.asarray(gene_ranks['pvals_adj'][key], dtype = float), index = df.index)
    df['logfoldchanges'] = pd.Series(data = np.asarray(gene_ranks['logfoldchanges'][key], dtype = float), index = df.index)
    shift = pd.Series(shift_by, index = df.index)
    df['log10p'] = np.negative(np.log10(np.asarray(df['pvals_adj'].add(shift))))
    return df
def ad_binner(ad, limits, key = "pseudotime", axis = 0):
    """ Assigns pseudotime bin to each cell in an anndata based on the provided pseudotime.
    
    Keyword arguments:
    ad: an AnnData object containing calculated pseudotime in .obs
    limits: a list of floats indicating the upper limits for each bin, not including the maximum pseudotime along the axis
    key: a string indicating the key in .obs to use as pseudotime
    axis: an integer indicating which axis to bin along (0: obs, 1: var).
    """
    if axis == 0: 
        bins = pd.Series(index = ad.obs_names, dtype = 'float64')
        pseudo_axis = ad.obs[key]
        max_pt = pseudo_axis.max()
        lims = limits.copy()
       	for lim in list(lims.keys()):	
            if lim >= max_pt:	
                lims.pop(lim)	
        lim_list = [0] + list(lims.keys())+[max_pt]	
        lims[max_pt] = lims[max(lims.keys())] + 1
        bin_cut = pd.cut(ad.obs[key], lim_list, duplicates = 'drop', include_lowest = True)
        bin_dict = dict(zip(bin_cut.cat.categories, lims.values()))
        bins = pd.Series([bin_dict[x] for x in bin_cut], index = ad.obs_names)
        cat_list = list(bins.unique())
        cat_list.sort()
        cat_type = pd.CategoricalDtype(categories=cat_list, ordered=True)
        bins = bins.astype(cat_type)
        ad.obs[f'{key}_bin'] = bins.astype(cat_type)
        return bins
    if axis == 1:
        pseudo_axis = ad.var[key]
        max_pt = pseudo_axis.max()
        lims = limits.copy()
        for lim in list(lims.keys()):
            if lim >= max_pt:
                lims.pop(lim)
        lim_list = [0] + list(lims.keys())+[max_pt]
        lims[max_pt] = limits[max(lims.keys())] + 1
        bin_cut = pd.cut(ad.var[key], lim_list, duplicates = 'drop', include_lowest = True)
        bin_dict = dict(zip(bin_cut.cat.categories, lims.values()))
        bins = pd.Series([bin_dict[x] for x in bin_cut], index = ad.var_names)
        cat_list = list(bins.unique())
        cat_list.sort()
        cat_type = pd.CategoricalDtype(categories=cat_list, ordered=True)
        bins = bins.astype(cat_type)
        ad.var[f'{key}_bin'] = bins.astype(cat_type)
        return bins
    
def binner(df, limits, axis = 0, key = 'rna_pseudotime'):
    """ Assigns pseudotime bin to each cell in an dataframe where the column axis labels are pseudotime values.
    
    Keyword arguments:
    df: a DataFrame with the index representing genes and the columns representing pseudotime values
    limits: a list of floats indicating the upper limits for each bin, not including the maximum pseudotime along the axis
    axis: axis along which to bin. 'O' bins along df.index, '1' along df.columns
    Returns:
    bins: a Series with index = df.index, containing the bin where each gene is maximally expressed
    """
def binner(df, limits, axis = 0, key = 'rna_pseudotime'):
    """ Assigns pseudotime bin to each cell in an dataframe where the column axis labels are pseudotime values.
    
    Keyword arguments:
    df: a DataFrame with the index representing genes and the columns representing pseudotime values
    limits: a list of floats indicating the upper limits for each bin, not including the maximum pseudotime along the axis
    axis: axis along which to bin. 'O' bins along df.index, '1' along df.columns
    Returns:
    bins: a Series with index = df.index, containing the bin where each gene is maximally expressed
    """
    if axis == 0:
        maxes = df.idxmax(axis = 1)
        max_pt = max(maxes)
        lims = limits.copy()
        for lim in list(lims.keys()):
            if lim >= max_pt:
                lims.pop(lim)
        lim_list = [0] + list(lims.keys())+[max_pt]
        lims[max_pt] = lims[max(lims.keys())] + 1
        bin_cut = pd.cut(pd.Series(maxes, index = df.index), lim_list, duplicates = 'drop', include_lowest = True)
        bin_dict = dict(zip(bin_cut.cat.categories, lims.values()))
        bins = pd.Series([bin_dict[x] for x in bin_cut], index = df.index)
        cat_list = list(bins.unique())
        cat_list.sort()
        cat_type = pd.CategoricalDtype(categories=cat_list, ordered=True)
        bins = bins.astype(cat_type)
    elif axis ==1:
        max_pt = max(df.columns)
        lims = limits.copy()
        for lim in list(lims.keys()):
            if lim >= max_pt:
                lims.pop(lim)
        lim_list = [0] + list(lims.keys())+[max_pt]
        lims[max_pt] = lims[max(lims.keys())] + 1
        bin_cut = pd.cut(pd.Series(df.columns, index = df.columns), lim_list, duplicates = 'drop', include_lowest = True)
        bin_dict = dict(zip(bin_cut.cat.categories, lims.values()))
        bins = pd.Series([bin_dict[x] for x in bin_cut], index = df.columns)
        cat_list = list(bins.unique())
        cat_list.sort()
        cat_type = pd.CategoricalDtype(categories=cat_list, ordered=True)
        bins = bins.astype(cat_type)
    return bins

def hv_gene_list(ad):
    """ Returns a list of genes marked highly variable in the given annData.
    
    Keyword arguments:
    ad: an AnnData object with ad.var['highly_variable'] assigned with scanpy.pp.highly_variable_genes
    
    Returns:
    a list of strings containing gene names
    """
    hv_names = ad.var['highly_variable']
    hv_loc = pd.Series(data = range(0, len(hv_names)), index = hv_names.index)
    cols = {'highly_variable': hv_names, 'loc' : hv_loc}
    hv = pd.DataFrame(cols)
    hv= hv[hv['highly_variable'] == True]
    return hv.index.tolist()

def plot_volcano(ax, df, title, top = 100, cutoff = 2, minchange = 2, labels = True):
    """ Plots a volcano plot of differential expression analysis results using the log2fold change and -log10(adj_p) values. Returns a list of a user-specified number of top genes sorted by adjusted-pvalue and abs(log fold change).
    
    Keyword arguments:
    ax: a Matplotlib axes object for plotting
    df: a dataframe with an index of gene names and columns 'logfoldchanges' and 'pvals_adj'
    title: a string containing the title of the plot
    top: int or None, specifies the number of top genes to be returned as a list. None returns all genes above the p-value and fold change cutoffs.
    cutoff: int, the -log10 adjusted p-value cutoff for significance
    minchange: int, the minimum fold change to be reported as a significant hit
    labels: boolean, plots x and y axis labels if True and no axis labels if Talse
    """
    colors = list()
    for gene in df.index:
        log = df.loc[gene, 'logfoldchanges']
        p = df.loc[gene,'pvals_adj']
        if p > 0.05:
            colors.append('grey')
        elif log >= 1: 
            colors.append('green')
        elif log <= -1:
            colors.append('red')
        else:
            colors.append('grey')
    ax.scatter(y = df['log10p'], x = df['logfoldchanges'], s = 0.5, c = colors  )
    if labels == True:
        ax.set_xlabel('log2(fold change)')
        ax.set_ylabel('-log10(adjusted p-value)')
    ax.hlines(y = cutoff, xmin =-30, xmax = 30, linestyle = '--', color = 'lightgray')
    ax.vlines(x = minchange, ymin = -10, ymax = 300, linestyle = '--', color = 'lightgray')
    ax.vlines(x = -minchange, ymin = -10, ymax = 300, linestyle = '--', color = 'lightgray')
    sort = df
    sort['abs'] = abs(sort['logfoldchanges'])
    sort = sort[sort['abs']>=minchange].sort_values(['log10p', 'abs'], ascending = [False, False])
    sort = sort[sort['log10p'] > cutoff].sort_values(['log10p', 'abs'], ascending = [False, False])
    if top == None:
        top = len(sort)
    for gene in sort.iloc[0:top].index:
        ax.annotate(gene, (df.loc[gene, 'logfoldchanges']+0.1,df.loc[gene, 'log10p']+0.1), fontsize = 3.0 )
    ax.set_title(title)
    return sort.iloc[0:top]

