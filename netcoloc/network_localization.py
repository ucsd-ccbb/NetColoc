# -*- coding: utf-8 -*-

'''Functions for performing network localization
'''

# External library imports
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import ndtr

# Internal module convenience imports
from .netcoloc_utils import *

def __init__(self):
    pass

# -------------------- LOCALIZATION ---------------------------------#

# TODO: add netprop localization functionality (outline here- needs to be tested)
def netprop_localization(z_scores, random_final_heats, seed_genes, z_score_threshold=3, plot=True):
    '''
    Function to calculate size of netowrk proximity, and evaluate if larger than expected by chance
    z_scores: netprop zscores output from netprop_zscore.py
    random_final_heats: random network propagation scores seeded from degree-matched input sets (output from netprop_zscore.py)
    seed_genes: list of seed genes for z_scores, random_final_heats
    z_threshold: threshold to call significant subgraph
    plot: if true plot the distribution
    '''
    
    #Precalculate mean and standard deviation (should really remove focal column, but takes a long time to do this)
    random_final_heats_mean = np.nanmean(np.log(random_final_heats), axis=0)
    random_final_heats_std = np.nanstd(np.log(random_final_heats), axis=0)

    random_proximal_network_sizes = []
    for r in range(len(random_final_heats)):
        if r % 100 == 0:
            print(r)
        #Average over all colums except focal r 
        random_z_score = (np.log(random_final_heats[r]) - random_final_heats_mean) / random_final_heats_std
        random_proximal_network_sizes.append(sum(random_z_score > z_score_threshold) - len(seed_genes)) #Don't count seeds in proximal network
      
    #Calculate the size of the true subgraph
    proximal_network_size = sum(z_scores > z_score_threshold) - len(seed_genes) #Don't count seeds in proximal network

    if plot == True:
        plt.figure(figsize=(5, 4))
        dfig=sns.distplot(random_proximal_network_sizes, label='random', kde=False)
        plt.vlines(proximal_network_size, ymin=0, ymax=dfig.dataLim.bounds[3], color='r', label='true gene set')
        plt.xlabel('size of proximal network, z > ' + str(z_score_threshold), fontsize=16)
        plt.ylabel('count', fontsize=16) 
        plt.legend(loc='upper left')

    return random_proximal_network_sizes, proximal_network_size




# TODO: add function to call significantly proximal subgraph (outline here)
# def proximal_subgraph(Gint,netprop_z,zthresh=3):
#    proximal_genes = netprop_z[netprop_z['z']>zthresh].index.tolist()
#    G_sub = nx.subgraph(Gint, proximal_genes)    
#    return G_sub



# TODO: need to rename the old localization functions, may need to simplify them to avoid confusion with netprop_localization
def localization(interactome, focal_genes, background_list=None, num_reps=10, method='both', plot=True, 
                 print_counter=False, sample_fraction=1):
    '''Calculates localization of an input set of focal genes on a background 
    interactome network.
        
    There is an option to compute localization using the number of edges between
    focal genes or using the largest connected component created by focal genes.
    Localization is compared to the localization of random sets of genes of 
    similar degree to the focal genes.
        
    Args:
        interactome (NetworkX graph): Background network to calculate 
            localization in.
        focal_genes (list): Set of genes to use to calculate localization.
        background_list (list): List of interactome genes to sample from when
            forming random gene sets. If this is set to None, all interactome
            genes will be available for sampling. (Default: None)
        num_reps (int): Number of times to calculate localization using random
            focal genes. (Default: 10)
        method (string): Either 'num_edges', 'LCC', or 'both'. The method that
            should be used to calculate localization. If 'num_edges' is chosen,
            localization will be calculated using the number of edges between
            members of the focal gene set. If 'LCC' is chosen, localization will
            be calculated using the size of the largest connected component
            formed by genes from the focal gene set. If 'both' is chosen, both
            methods will be used.
        plot (bool): If this is True, then the localization distributions will
            be outputted as a plot. If this is False, no plot will be drawn.
        print_counter (bool): If this is True, a counter will print the current
            repetition number every 25 repetitions. Useful when num_reps is a
            large number.
            
    Returns: 
        
          numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
              Empty if method = 'LLC'. 
          numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
              similar degree in the background network. Empty if method = 'LLC'.
          LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
              Empty if method = 'numedges'. 
          LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
              similar degree in the background network. Empty if method = 'numedges'. 
    '''
    degrees_dict = dict(interactome.degree)
    # Create degree bins to sample from
    bins, actual_degree_to_bin_index = get_degree_binning(degrees_dict, 10)
    
    focal_genes = list(np.intersect1d(focal_genes, interactome.nodes())) # only use focal_genes which are in Gint
    numedges_list = []
    numedges_rand = []
    LCC_list = []
    LCC_rand = []
    
    if background_list==None:
        background_list=interactome.nodes()
    
    for r in range(num_reps):

        if print_counter == True:
            # so user knows how far along the process is
            if (r % 25) == 0:
                print(r)

        focal_80 = focal_genes
        np.random.shuffle(focal_80)
        focal_80 = focal_80[:int(len(focal_80)*sample_fraction)]

        # find genes with similar degrees to focal gene degree
        seed_random = []
        for g in focal_80:
            degree_temp = degrees_dict[g]
            #genes_temp = bin_df.loc[actual_degree_to_bin_df_idx[degree_temp]]['genes_binned'] # use the lookup table for speed
            genes_temp = bins[actual_degree_to_bin_index[degree_temp]]
            np.random.shuffle(genes_temp) # shuffle them
            while (genes_temp[0] in seed_random) or (genes_temp[0] not in background_list): # make sure the gene isn't already in the list, but is in the background_list
                np.random.shuffle(genes_temp) # shuffle them
            
            seed_random.append(genes_temp[0]) # build the seed_D1_random list

        if (method == 'numedges') or (method == 'both'):
            
            # number edges calc on focal set
            numedges_temp = len(nx.subgraph(interactome,focal_80).edges())
            numedges_list.append(numedges_temp)
            
            # number edges calc on random sample
            numedges_temp_rand = len(nx.subgraph(interactome,seed_random).edges())
            numedges_rand.append(numedges_temp_rand)
            
        if (method == 'LCC') or (method == 'both'):
            
            # LLC calc on focal set
            G_sub_temp = nx.Graph(nx.subgraph(interactome, focal_80))
            G_sub_temp = max([G_sub_temp.subgraph(c) for c in nx.connected_components(G_sub_temp)], key = len)
            LCC_list.append(len(G_sub_temp.nodes()))
            
            # LLC calc on random sample
            G_sub_temp = nx.Graph(nx.subgraph(interactome, seed_random))
            G_sub_temp = max([G_sub_temp.subgraph(c) for c in nx.connected_components(G_sub_temp)], key = len)
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
            plt.ylim([0, 1])
                
        if (method == 'LCC') or (method == 'both'):
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(LCC_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(LCC_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('largest connected component size', fontsize = 16)
            plt.title('Largest Connected Component Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
            plt.ylim([0, 1])
            
    return numedges_list, numedges_rand, LCC_list, LCC_rand
    
def localization_fractional(interactome, focal_genes, background_list=None, num_reps=10, sample_fraction=1, method='both',
                 plot=True, print_counter=False):
    '''Calculates localization of an input set of focal genes on a background 
    interactome network.
        
    There is an option to compute localization using the number of edges between
    focal genes or using the largest connected component created by focal genes.
    Localization is compared to the localization of random sets of genes of 
    similar degree to the focal genes. Localization is calculated by sampling
    sub-sections of the focal genes or random gene set.
        
    Args:
        interactome (NetworkX graph): Background network to calculate 
            localization in.
        focal_genes (list): Set of genes to use to calculate localization.
        background_list (list): List of interactome genes to sample from when
            forming random gene sets. If this is set to None, all interactome
            genes will be available for sampling. (Default: None)
        num_reps (int): Number of times to calculate localization using random
            focal genes. (Default: 10)
        sample_fraction (float): A number between 0 and 1. Percent of gene set
            to sample. (Default: 1)
        method (string): Either 'num_edges', 'LCC', or 'both'. The method that
            should be used to calculate localization. If 'num_edges' is chosen,
            localization will be calculated using the number of edges between
            members of the focal gene set. If 'LCC' is chosen, localization will
            be calculated using the size of the largest connected component
            formed by genes from the focal gene set. If 'both' is chosen, both
            methods will be used.
        plot (bool): If this is True, then the localization distributions will
            be outputted as a plot. If this is False, no plot will be drawn.
        print_counter (bool): If this is True, a counter will print the current
            repetition number every 25 repetitions. Useful when num_reps is a
            large number.
            
    Returns: 
          numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
              Empty if method = 'LLC'. 
          numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
              similar degree in the background network. Empty if method = 'LLC'.
          LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
              Empty if method = 'numedges'. 
          LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
              similar degree in the background network. Empty if method = 'numedges'. 
    '''
    degrees_dict = dict(interactome.degree)
    # Create degree bins to sample from
    bins, actual_degree_to_bin_index = get_degree_binning(degrees_dict, 10)
    
    focal_genes = list(np.intersect1d(focal_genes, interactome.nodes())) # only use focal_genes which are in Gint
    numedges_list = []
    numedges_rand = []
    LCC_list = []
    LCC_rand = []
    
    if background_list==None:
        background_list=interactome.nodes()
    
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
            degree_temp = degrees_dict[g]
            #genes_temp = bin_df.loc[actual_degree_to_bin_df_idx[degree_temp]]['genes_binned'] # use the lookup table for speed
            genes_temp = bins[actual_degree_to_bin_index[degree_temp]]
            np.random.shuffle(genes_temp) # shuffle them
            while (genes_temp[0] in seed_random) or (genes_temp[0] not in background_list): # make sure the gene isn't already in the list, but is in the background_list
                np.random.shuffle(genes_temp) # shuffle them
            
            seed_random.append(genes_temp[0]) # build the seed_D1_random list

        if (method == 'numedges') or (method == 'both'):
            
            # number edges calc on focal set
            numedges_temp = len(nx.subgraph(interactome,focal_80).edges())
            numedges_list.append(numedges_temp)
            
            # number edges calc on random sample
            numedges_temp_rand = len(nx.subgraph(interactome,seed_random).edges())
            numedges_rand.append(numedges_temp_rand)
            
        if (method == 'LCC') or (method == 'both'):
            
            # LLC calc on focal set
            G_sub_temp = nx.Graph(nx.subgraph(interactome, focal_80))
            G_sub_temp = max([G_sub_temp.subgraph(c) for c in nx.connected_components(G_sub_temp)], key = len)
            LCC_list.append(len(G_sub_temp.nodes()))
            
            # LLC calc on random sample
            G_sub_temp = nx.Graph(nx.subgraph(interactome, seed_random))
            G_sub_temp = max([G_sub_temp.subgraph(c) for c in nx.connected_components(G_sub_temp)], key = len)
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
            plt.ylim([0, 1])
                
        if (method == 'LCC') or (method == 'both'):
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(LCC_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(LCC_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('largest connected component size', fontsize = 16)
            plt.title('Largest Connected Component Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
            plt.ylim([0, 1])
            
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
