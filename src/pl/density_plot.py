import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import density_analysis as da
from scipy import signal


def sort_gene_trends(label,lin, trend_set, gene_sets, limits, layers = False,bin_label = "gex_bins"):
    """
    A function to sort gene expression trends by the pseudotime bin where each gene is most highly expressed.
    
    Keyword arguments:
    trend_dict: a dictionary with genes or loci as keys and dataframes of predicted expression and associated statistics as values (output of getspan.span.calc_reg)
    gene_set: a list of genes/keys which is a subset of trend_dict.keys(). if None, use all genes
    limits: a list of floats giving the upper bounds of each pseudotime bin, not inclusive of the maximum value along the axis.
    """
    trend_types = list(trend_set[label][lin].keys())
    keys = list(trend_set[label][lin][trend_types[0]].keys())
    bin_index = trend_set[label][lin][trend_types[0]][keys[0]].index
    gene_index = pd.Index(gene_sets[label][lin])
    df = pd.DataFrame(data = np.zeros((len(bin_index), len(gene_index))),columns = gene_index, index = bin_index)
    for gene in gene_index:
        for i in bin_index:
            df[gene].loc[i] = trend_set[label][lin][trend_types[0]][gene]['expression'][i]
    sorted_df = df.sort_values(axis = 1, by = df.first_valid_index(), ascending = False)
    sorted_df.index = trend_set[label][lin][trend_types[0]][keys[0]]['pseudo_axis']
    T_df = sorted_df.transpose()
    T_df[45.0] = da.binner(T_df, limits, axis = 0) #not exactly good code, but this is a temporary column with an index that will not convert to string
    T_df.sort_values(axis = 0, by = 45.0, ascending = True, inplace = True)
    bins = T_df[45.0]
    T_df = T_df.drop(45.0, axis = 1)
    ad = anndata.AnnData(T_df)
    ad.obs[bin_label] =  pd.Series(bins, index = ad.obs_names) 
    #print(trend_types[0])
    if layers == True:
        for i in range(len(trend_types)):
            _add_trend_layer(ad, bin_index, gene_index, trend_types[i], trend_set[label][lin], gene_sets, limits)
    return ad


def sort_trends(trends, gene_set, limits, layers = True, bin_label = 'gex_bins'):
    """
    A function to sort gene expression trends by the pseudotime bin where each gene is most highly expressed.
    
    Keyword arguments:
    trend_dict: a dictionary with genes or loci as keys and dataframes of predicted expression and associated statistics as values (output of getspan.span.calc_reg)
    gene_set: a list of genes/keys which is a subset of trend_dict.keys(). if None, use all genes
    limits: a list of floats giving the upper bounds of each pseudotime bin, not inclusive of the maximum value along the axis.
    """
    trends = trends
    trend_types = list(trends.keys())#list(trends.layers.keys())
    keys = list(trends[trend_types[0]].keys())
    bin_index = trends[trend_types[0]][keys[0]].index
    if gene_set == None:
        gene_set = list(trends.keys())
    gene_index = pd.Index(gene_set)
    df = pd.DataFrame(data = np.zeros((len(bin_index), len(gene_index))),columns = gene_index, index = bin_index)
    for gene in gene_index:
        df[gene].loc[i] = trends[trend_types[0]][gene]['expression'][bin_index]
    sorted_df = df.sort_values(axis = 1, by = df.first_valid_index(), ascending = False)
    sorted_df.index = trends[trend_types[0]][keys[0]]['pseudo_axis']
    T_df = sorted_df.transpose()
    bins = da.binner(T_df, limits, axis = 0)
    T_df[3.0] = bins #temporary column name which will not cause conversion to string
    T_df.sort_values(axis = 0, by = 3.0, ascending = True, inplace = True)
    bins = T_df[3.0]
    T_df = T_df.drop(3.0, axis = 1)
    T_df.columns = T_df.columns.astype(str)
    #T_df.columns = pd.NumericIndex(T_df.index()) ## will be needed for compatibility with pandas 2.x
    ad = anndata.AnnData(T_df, dtype = np.float64)
    ad.obs[bin_label] = bins
    if layers == True:
        for i in range(len(trend_types)):
            _add_trend_layer(ad , bin_index, gene_index, trend_types[i], trends, gene_set, limits)
    return ad

def _add_trend_layer(ad, bin_index, gene_index, trend_type, trend_set, gene_set, limits):
    """
    helper function for sort_gene_trends. Adds a data layer to the annData object created by the parent function for each trend type in the given input.
    
    keyword arguments:
    ad: the annData object passed from sort_gene_trends
    bin_index: the index for the columns of the dataframe, passed in from sort_gene_trends
    gene_index: the index for the rows of the dataframe, passed in from sort_gene_trends
    trend_type: a key in trends.keys(), from parent function
    trend_set: the dictionary of predicted trends generated by getspan.span.calc_reg
    gene_set: the gene set selected in the parent function (user-specified, or all genes if None)
    limits: a list of floats giving the upper bounds of each pseudotime bin, not inclusive of the maximum value along the axis.
    """
    #print(trend_type)
    #keys = list(trend_set[label][lin][trend_type].keys())
    bin_index = bin_index#trend_set[label][lin][trend_type][keys[0]].index
    gene_index = gene_index#pd.Index(gene_set[label][lin])
    df = pd.DataFrame(data = np.zeros( (len(gene_index), len(bin_index))), index = ad.obs_names, columns = bin_index)
    for gene in ad.obs_names:
        for i in bin_index:
            df.loc[gene, i] = trend_set[trend_type][gene]['expression'][i]
    ad.layers[trend_type] = df
    return ad


def plot_gene_trends(ad, lin, modality, col_map= None, row_map = None, row_labels = False, col_labels = False, bin_max = None, bin_min = None):
    """
    A function which plots the Z-scored predicted expression for each gene in the given trend matrix.
    
    Keyword arguments:
    ad: an AnnData object containing predicted expression trends and the associated metadata.
    lin: the celltype lineage corresponding to the given expression trends
    modality: default "rna", the layer of the trend matrix which should be accessed for plotting
    col_map: default None, a dictionary mapping category labels to colors - should correspond to colors in ad.obs['density_color']
    row_map: default None, a dictionary mapping category labels to colors - should correspond to colors in ad.var["bin_color"]
    row_labels: default False, a boolean value which labels rows with gene names if True
    col_labels: default False, a boolean value which labels columns with column values if True
    bin_max: int or float (default None), indicates the maximum bin level to be plotted
    bin_max: int or float (default None), indicates the minimum bin level to be plotted. if bin_min == bin_max only that bin will be plotted.
    """
    if bin_max == None:
        bin_max = max(ad.var['density_bins'].values.unique())
    series = pd.Series(ad.var['density_bins']<=bin_max)
    if bin_min != None and bin_min <= bin_max:
        series = pd.Series(ad.var['density_bins']>=bin_min)
    df =  pd.DataFrame(ad[:, series].layers[modality], index = ad[:, series].obs_names, columns = ad[:, series].var_names)
    if row_map == None:
        tab20 = dict()
        for num in range(0, 20):
            tab20[num] = matplotlib.colors.to_hex(matplotlib.cm.tab20(num))
        row_map = {int(x):tab20[x%20] for x in ad[ :, series].obs['gex_bins'].unique()}
    row_colors = ad[:, series].obs['bin_color']
    col_colors = ad[ :, series].var['density_color']
    sort = ad[ :, series].obs['gex_bins'].sort_values(ascending = True)
    df = df.loc[sort.index, :]
    g = sns.clustermap(df ,cmap = 'RdBu_r', 
                   center = 0, col_cluster=False, 
                   row_cluster= False, yticklabels = row_labels, xticklabels = col_labels,
                   z_score = 0, vmin = -2, vmax = 2, col_colors = col_colors, row_colors = row_colors[df.index])
    if col_map != None:
        for label in col_map.keys():
            g.ax_col_dendrogram.bar(0, 0, color=col_map[label],
                                    label=label, linewidth=0)
    if row_map != None:
        for label in row_map.keys():
            g.ax_row_dendrogram.bar(0, 0, color=row_map[label],
                                    label=label, linewidth=0)
    #g.ax_row_dendrogram.set_visible(False)
    g.ax_row_dendrogram.legend(loc = "center", ncol = 2, title = 'gene expression bins' )
    g.ax_col_dendrogram.legend(loc="center", ncol=6, title = f"{lin} density level")
    plt.ylabel('Z-score')
    return g
    #plt.savefig('figures/{lin}_trends.png')
    
    
def side_by_side_plots(path1, path2, out_path, figsize = (20,10), titles = ['RNA', 'ATAC'] ):
    """
    plot any two images side by side, with the titles given as a list.
    
    Keyword arguments: 
    path1: a string indicating the path to the image to be plotted on the left.
    path2: a string indicating the path to the image to be plotted on the right.
    out_path: a string indicating the path to save the side-by-side plot to.
    figsize: default (20, 10), a tuple of integers setting the width and height of the plot
    titles: a list of length 2 containing strings for the titles of the left and right plots.
    """ 
    fig, axs = plt.subplots(1,2, figsize =figsize)
    img = matplotlib.image.imread(path1)
    axs[0].imshow(img)
    axs[0].set_title(titles[0])
    img = matplotlib.image.imread(path2)
    axs[1].imshow(img)
    axs[1].set_title(titles[1])
    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[0].set_yticks([])
    axs[1].set_yticks([])
    sns.despine(bottom = True, left = True)
    plt.savefig(out_path)
    
def cross_correlate_derivatives(gene, return_lags = False):
    lineage = 'primed'
    lag_ticks = list(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'].iloc[::-1].drop(0).multiply(-1)) + list(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'])
    primed_trend = np.gradient(np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten())
    n_primed_trend = (primed_trend - np.mean(primed_trend)) /(np.std(primed_trend) * len(primed_trend))
    rna_trend = np.gradient(np.array(atac_trends['rna']['Bcells'][gene]['expression']).flatten())
    n_rna_trend = (rna_trend - np.mean(rna_trend)) / np.std(primed_trend)
    primed_corr = np.correlate(n_primed_trend, n_rna_trend, mode = 'full')
    primed_lags = signal.correlation_lags(len(np.array(atac_trend_ads[lineage]['Bcells'][gene].X).flatten()), 
                                   len(np.array(atac_trend_ads['rna']['Bcells'][gene].X).flatten()), mode = 'full')
    print(lag_ticks[np.argmax(primed_corr)])
    lineage = 'lineage_specific'
    lineage_trend = np.gradient(np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten())
    n_lineage_trend = (lineage_trend - np.mean(lineage_trend)) /(np.std(lineage_trend) * len(lineage_trend))
    lin_corr = np.correlate(n_lineage_trend, n_rna_trend, mode = 'full')
    lin_lags = signal.correlation_lags(len(np.array(atac_trend_ads[lineage]['Bcells'][gene].X).flatten()), 
                                   len(np.array(atac_trend_ads['rna']['Bcells'][gene].X).flatten()), mode = 'full' )
    print(lag_ticks[np.argmax(lin_corr)])
    fig, (ax_primed, ax_corr) = plt.subplots(2, 1, figsize = (10, 8))
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten()), color = 'cornflowerblue')
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends['primed']['Bcells'][gene]['expression']).flatten()), color =  'indianred')
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends['rna']['Bcells'][gene]['expression']).flatten()), color = 'green')
    max_y = np.max([np.max(lin_corr), np.max(primed_corr)])
    min_y = np.min([np.min(lin_corr), np.min(primed_corr)])
    ax_corr.vlines(x = lag_ticks[np.argmax(primed_corr)], ymin = min_y, ymax = max_y, linestyle = '--', color = 'gray')
    ax_corr.vlines(x = lag_ticks[np.argmax(lin_corr)], ymin = min_y, ymax = max_y, linestyle = '--', color = 'gray')
    ax_corr.plot(lag_ticks, lin_corr,  color = 'cornflowerblue')
    ax_corr.plot(lag_ticks, primed_corr,  color = 'indianred')
    ax_corr.plot(lag_ticks[np.argmax(lin_corr)], np.max(lin_corr), marker = 'o', color =  'cornflowerblue')
    ax_corr.plot(lag_ticks[np.argmax(primed_corr)], np.max(primed_corr), marker = 'o', color = 'indianred')
    ax_primed.set_title(f'{gene} trends')
    ax_primed.set_ylabel('Z-score')
    ax_primed.set_xlabel('pseudotime')
    ax_corr.set_ylabel('zero normalized\n cross-correlation')
    ax_corr.set_xlabel('lag')
    patches = list()
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='indianred', label=f'Primed'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='cornflowerblue', label=f'Lineage specific'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='green', label='RNA'))
    ax_primed.legend(handles=patches, fontsize = 'xx-small', frameon = False )
    patches = list()
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='indianred', label=f'Primed vs. RNA'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='cornflowerblue', label=f'Lineage specific vs. RNA'))
    ax_corr.legend(handles=patches, fontsize = 'xx-small', frameon = False )
    #ax_corr.set_xticks(lag_ticks)
    #ax_corr.set_title(f'{gene} trend cross-correlation lags')
    sns.despine()
    lags = [lag_ticks[np.argmax(primed_corr)], lag_ticks[np.argmax(lin_corr)]]
    if return_lags == True:
        return fig, lags
    else:
        return fig
    
def cross_correlate(gene, return_lags = False):
    lineage = 'primed'
    lag_ticks = list(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'].iloc[::-1].drop(0).multiply(-1)) + list(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'])
    primed_trend = np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten()
    n_primed_trend = (primed_trend - np.mean(primed_trend)) /(np.std(primed_trend) * len(primed_trend))
    rna_trend = np.array(atac_trends['rna']['Bcells'][gene]['expression']).flatten()
    n_rna_trend = (rna_trend - np.mean(rna_trend)) / np.std(primed_trend)
    primed_corr = np.correlate(n_primed_trend, n_rna_trend, mode = 'full')
    primed_lags = signal.correlation_lags(len(np.array(atac_trend_ads[lineage]['Bcells'][gene].X).flatten()), 
                                   len(np.array(atac_trend_ads['rna']['Bcells'][gene].X).flatten()), mode = 'full')
    print(lag_ticks[np.argmax(primed_corr)])
    lineage = 'lineage_specific'
    lineage_trend = np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten()
    n_lineage_trend = (lineage_trend - np.mean(lineage_trend)) /(np.std(lineage_trend) * len(lineage_trend))
    lin_corr = np.correlate(n_lineage_trend, n_rna_trend, mode = 'full')
    lin_lags = signal.correlation_lags(len(np.array(atac_trend_ads[lineage]['Bcells'][gene].X).flatten()), 
                                   len(np.array(atac_trend_ads['rna']['Bcells'][gene].X).flatten()), mode = 'full' )
    print(lag_ticks[np.argmax(lin_corr)])
    fig, (ax_primed, ax_corr) = plt.subplots(2, 1, figsize = (10, 8))
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends[lineage]['Bcells'][gene]['expression']).flatten()), color = 'cornflowerblue')
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends['primed']['Bcells'][gene]['expression']).flatten()), color =  'indianred')
    ax_primed.plot(atac_trends[lineage]['Bcells'][gene]['pseudo_axis'], stats.zscore(np.array(atac_trends['rna']['Bcells'][gene]['expression']).flatten()), color = 'green')
    max_y = np.max([np.max(lin_corr), np.max(primed_corr)])
    min_y = np.min([np.min(lin_corr), np.min(primed_corr)])
    ax_corr.vlines(x = lag_ticks[np.argmax(primed_corr)], ymin = min_y, ymax = max_y, linestyle = '--', color = 'gray')
    ax_corr.vlines(x = lag_ticks[np.argmax(lin_corr)], ymin = min_y, ymax = max_y, linestyle = '--', color = 'gray')
    ax_corr.plot(lag_ticks, lin_corr,  color = 'cornflowerblue')
    ax_corr.plot(lag_ticks, primed_corr,  color = 'indianred')
    ax_corr.plot(lag_ticks[np.argmax(lin_corr)], np.max(lin_corr), marker = 'o', color =  'cornflowerblue')
    ax_corr.plot(lag_ticks[np.argmax(primed_corr)], np.max(primed_corr), marker = 'o', color = 'indianred')
    ax_primed.set_title(f'{gene} trends')
    ax_primed.set_ylabel('Z-score')
    ax_primed.set_xlabel('pseudotime')
    ax_corr.set_ylabel('zero normalized\n cross-correlation')
    ax_corr.set_xlabel('lag')
    patches = list()
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='indianred', label=f'Primed'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='cornflowerblue', label=f'Lineage specific'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='green', label='RNA'))
    ax_primed.legend(handles=patches, fontsize = 'xx-small', frameon = False )
    patches = list()
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='indianred', label=f'Primed vs. RNA'))
    patches.append(plt.Line2D([],[], color = 'white', marker = 'o', markerfacecolor='cornflowerblue', label=f'Lineage specific vs. RNA'))
    ax_corr.legend(handles=patches, fontsize = 'xx-small', frameon = False )
    #ax_corr.set_xticks(lag_ticks)
    #ax_corr.set_title(f'{gene} trend cross-correlation lags')
    sns.despine()
    lags = [lag_ticks[np.argmax(primed_corr)], lag_ticks[np.argmax(lin_corr)]]
    if return_lags == True:
        return fig, lags
    else:
        return fig