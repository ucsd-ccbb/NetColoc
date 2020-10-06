"""
-------------------------------------------
Author: Brin Rosenthal
Date: 7/3/19
-------------------------------------------
"""

# common packages, most likely already installed
import scipy
import math
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys
from scipy.special import ndtr

# uncommon packages required for this analysis
import seaborn as sns # pip install seaborn

# -------------------- LOCALIZATION ---------------------------------#

def localization(Gint, focal_genes, num_reps = 10, sample_frac = 0.8, method = 'numedges', plot = True, print_counter = False,
                background_list=None):
    
    """
        Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
        Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
        localization analysis. Calculates by sampling sub-sections of the focal genes/random set. Percentage to sample
        is set by sample_frac. Option to plot the distributions of random and focal gene localizaiton.
        
        Args:
            Gint: Networkx Graph, background network to randomly sample from
            focal_genes: List, set of genes to calculate localization of
            num_reps: Int, number of times to randomly sample
            sample_frac: Float, percent of sampled genes
            method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
            plot: Bool, whether to plot the distributions in the output jupyter notebook cell
            print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                           Useful when the num_reps is very high.
            background_list: list of background genes to sample from. If none, jsut use all interactome genes
            
        Returns: 
            numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
                Empty if method = 'LLC'. 
            numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
                similar degree in the background network. Empty if method = 'LLC'.
            LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
                Empty if method = 'numedges'. 
            LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
                similar degree in the background network. Empty if method = 'numedges'. 
    """

    # Create degree bins to sample from
    bins = get_degree_binning(Gint, 10)
    min_degree, max_degree, genes_binned = zip(*bins)
    bin_df = pd.DataFrame({'min_degree':min_degree, 'max_degree':max_degree, 'genes_binned':genes_binned})
    
    # create a lookup table for degree and index
    actual_degree_to_bin_df_idx = {}
    for i in range(0, bin_df['max_degree'].max() + 1):
        idx_temp = bin_df[ (bin_df['min_degree'].lt(i + 1)) & (bin_df['max_degree'].gt(i - 1)) ].index.tolist()

        if len(idx_temp) > 0: # there are some degrees which aren't represented in the graph
            actual_degree_to_bin_df_idx[i] = idx_temp[0]
    
    focal_genes = list(np.intersect1d(focal_genes, Gint.nodes())) # only use focal_genes which are in Gint

    numedges_list = []
    numedges_rand = []
    LCC_list = []
    LCC_rand = []
    
    if background_list==None:
        background_list=Gint.nodes()
    
    for r in range(num_reps):

        if print_counter == True:
            # so user knows how far along the process is
            if (r % 25) == 0:
                print(r)

        focal_80 = focal_genes
        np.random.shuffle(focal_80)
        focal_80 = focal_80[:int(len(focal_80)*sample_frac)]

        # find genes with similar degrees to focal gene degree
        seed_random = []
        for g in focal_80:
            degree_temp = nx.degree(Gint,g)
            genes_temp = bin_df.loc[actual_degree_to_bin_df_idx[degree_temp]]['genes_binned'] # use the lookup table for speed
            np.random.shuffle(genes_temp) # shuffle them
            while (genes_temp[0] in seed_random) or (genes_temp[0] not in background_list): # make sure the gene isn't already in the list, but is in the background_list
                np.random.shuffle(genes_temp) # shuffle them
            
            seed_random.append(genes_temp[0]) # build the seed_D1_random list
            
        #print(len(focal_80))
        #print(len(seed_random))
        #print(len(np.unique(seed_random)))

        if (method == 'numedges') or (method == 'both'):
            
            # number edges calc on focal set
            numedges_temp = len(nx.subgraph(Gint,focal_80).edges())
            numedges_list.append(numedges_temp)
            
            # number edges calc on random sample
            numedges_temp_rand = len(nx.subgraph(Gint,seed_random).edges())
            numedges_rand.append(numedges_temp_rand)
            
        if (method == 'LCC') or (method == 'both'):
            
            # LLC calc on focal set
            G_sub_temp = nx.Graph(nx.subgraph(Gint, focal_80))
            G_sub_temp = max(nx.connected_component_subgraphs(G_sub_temp), key = len)
            LCC_list.append(len(G_sub_temp.nodes()))
            
            # LLC calc on random sample
            G_sub_temp = nx.Graph(nx.subgraph(Gint, seed_random))
            G_sub_temp = max(nx.connected_component_subgraphs(G_sub_temp), key=len)
            LCC_rand.append(len(G_sub_temp.nodes()))
    
    if plot == True:
        if (method == 'numedges') or (method == 'both'):
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(numedges_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(numedges_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('number of edges', fontsize = 16)
            plt.title('Number of Edges Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
                
        if (method == 'LCC') or (method == 'both'):
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(LCC_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(LCC_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('largest connected component size', fontsize = 16)
            plt.title('Largest Connected Component Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
            
    return numedges_list, numedges_rand, LCC_list, LCC_rand
    

def localization_full(Gint, focal_genes, 
                       num_reps = 200,
                       method = 'LCC', 
                       print_counter = False,
                       label = 'focal genes',
                       line_height = 0.1,
                       legend_loc = 'upper left'):
                       
    """
        Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
        Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
        localization analysis. DOes no sub-sampling. Plots the distribution of random gene localizaiton, and 
        marks the focal set localization on distribution. Includes p-value of focal set localization.
        
        Args:
            Gint: Networkx Graph, background network to randomly sample from
            focal_genes: List, set of genes to calculate localization of
            num_reps: Int, number of times to randomly sample
            method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
            print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                           Useful when the num_reps is very high.
            label: String, label for focal genes in graph legend
            line_height: Float, the height of the red line that marks the focal gene localization
            legend_loc: String, relative position of legend in graph. Something similar to 'upper left'.
            
        Returns: 
            numedges_list: List, the number of edges calculated for each rep, over focal genes. 
                Empty if method = 'LLC'. 
            numedges_rand: List, the number of edges calculated for each rep, over random genes of 
                similar degree in the background network. Empty if method = 'LLC'.
            LCC_list: List, the size of the largest connected component, calculated for each rep, over focal genes. 
                Empty if method = 'numedges'. 
            LCC_rand: List, the size of the largest connected component, calculated for each rep, over random genes of 
                similar degree in the background network. Empty if method = 'numedges'. 
    """
    
    numedges_list, numedges_rand, LCC_list, LCC_rand = localization(Gint, focal_genes, num_reps, 
                                                                    sample_frac = 1, 
                                                                    method = method, 
                                                                    plot = False,
                                                                    print_counter = print_counter)
    if method == 'numedges':
        analysis_list = numedges_list
        analysis_rand = numedges_rand
        title = 'number of edges'
    else:
        analysis_list = LCC_list
        analysis_rand = LCC_rand
        title = 'largest connected component'

    # plot distributions for non-sampled case
    fig, ax = plt.subplots(figsize = (12, 7))
    sns.set_style('white') 
    plt.vlines(np.mean(analysis_list), ymin = 0, ymax = line_height, color = 'r', lw = 2, label = label)
    sns.kdeplot(analysis_rand, ax = ax, color = 'k', lw = 2, alpha = 0.5, shade = True, label = 'random')
    plt.legend(loc = legend_loc, fontsize = 12)
    plt.ylabel('frequency', fontsize = 16)
    plt.xlabel(title, fontsize = 16)

    # print the z-score and fdr
    analysis_z = (np.mean(analysis_list) - np.mean(analysis_rand))/float(np.std(analysis_rand))

    print(1 - ndtr(analysis_z))

    plt.title('permutation p = ' + str(1 - ndtr(analysis_z)))
    
    return numedges_list, numedges_rand, LCC_list, LCC_rand
    
    
def get_degree_binning(g, bin_size, lengths = None):
    
    """
        Helper function for localization(). This function comes from network_utilities.py of emregtoobox. https://github.com/emreg00/toolbox 
    """
    
    degree_to_nodes = {}
    
    if sys.version_info >= (3, 0):
        for node, degree in dict(g.degree()).items():
            if lengths is not None and node not in lengths:
                continue
            degree_to_nodes.setdefault(degree, []).append(node)
    
    else:
        for node, degree in dict(g.degree()).iteritems():
            if lengths is not None and node not in lengths:
                continue
            degree_to_nodes.setdefault(degree, []).append(node)
        
    values = list(degree_to_nodes.keys())
    values.sort()
    bins = []
    i = 0
    while i < len(values):
        
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
            
        if i == len(values):
            i -= 1
            
        high = values[i]
        i += 1 
        
        #print low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
            
    return bins
            