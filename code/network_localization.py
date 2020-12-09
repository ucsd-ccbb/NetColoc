"""
-------------------------------------------
Author: Brin Rosenthal, Sophie Liu
Date: 10/7/20
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
import netprop_zscore

# uncommon packages required for this analysis
import seaborn as sns # pip install seaborn

# -------------------- LOCALIZATION ---------------------------------#

# TODO: add netprop localization functionality (outline here- needs to be tested)
def netprop_localization(netprop_z,Fnew_rand,seed_genes,zthresh=3,plot=True):
    '''
    
    Function to calculate size of netowrk proximity, and evaluate if larger than expected by chance
    netprop_z: netprop zscores output from netprop_zscore.py
    Fnew_rand: random network propagation scores seeded from degree-matched input sets (output from netprop_zscore.py)
    seed_genes: list of seed genes for z, Fnew_rand
    zthresh: threshold to call significant subgraph
    plot: if true plot the distribution
    
    '''
    
    # create randomized z-scores
    z_rand = (np.log(Fnew_rand[0])-np.nanmean(np.log(Fnew_rand.T),axis=0))/np.nanstd(np.log(Fnew_rand.T),axis=0)
    z_rand = pd.DataFrame(z_rand)
    z_rand.columns=['z_rand']
    
    size_rand_pnet=[]
    # precalculate mean and sd (should really remove focal column, but takes a long time to do this)
    Fnew_rand_mean = np.nanmean(np.log(Fnew_rand),axis=1)
    Fnew_rand_std = np.nanstd(np.log(Fnew_rand),axis=1)
    for r in np.arange(len(Fnew_rand.columns)):
        if (r%100)==0:
            print(r)
        zrand_temp = (np.log(Fnew_rand[r])-Fnew_rand_mean)/Fnew_rand_std
        size_rand_pnet.append(sum(zrand_temp>zthresh)-len(seed_genes)) # don't count seeds in proximal network

    # calculate the size of the true subgraph
    focal_size_pnet = sum(netprop_z['z']>zthresh)-len(seed_genes) # don't count seeds in proximal network
    if plot==True:
        plt.figure(figsize=(5,4))
        dfig=sns.distplot(size_rand_pnet,label='random',kde=False)
        plt.vlines(focal_size_pnet,ymin=0,ymax=dfig.dataLim.bounds[3],color='r',label='true gene set')
        plt.xlabel('size of proximal network, z>'+str(zthresh),fontsize=16)
        plt.ylabel('count',fontsize=16)
        plt.legend(loc='upper left')

    
    return size_rand_pnet, focal_size_pnet




# TODO: add function to call significantly proximal subgraph (outline here)
# def proximal_subgraph(Gint,netprop_z,zthresh=3):
#    proximal_genes = netprop_z[netprop_z['z']>zthresh].index.tolist()
#    G_sub = nx.subgraph(Gint, proximal_genes)    
#    return G_sub



# TODO: need to rename the old localization functions, may need to simplify them to avoid confusion with netprop_localization
def localization(interactome, focal_genes, num_reps=10, sample_fraction=0.8, method = 'numedges', plot=True, print_counter=False, background_list=None):
    
    """
        Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
        Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LCC') 
        localization analysis. Calculates by sampling sub-sections of the focal genes/random set. Percentage to sample
        is set by sample_frac. Option to plot the distributions of random and focal gene localizaiton.
        
        Args:
            Gint: Networkx Graph, background network to randomly sample from
            focal_genes: List, set of genes to calculate localization of
            num_reps: Int, number of times to randomly sample
            sample_frac: Float, percent of sampled genes
            method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LCC', or 'both'.
            plot: Bool, whether to plot the distributions in the output jupyter notebook cell
            print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                           Useful when the num_reps is very high.
            background_list: list of background genes to sample from. If none, jsut use all interactome genes
            
        Returns: 
            numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
                Empty if method = 'LCC'. 
            numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
                similar degree in the background network. Empty if method = 'LCC'.
            LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
                Empty if method = 'numedges'. 
            LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
                similar degree in the background network. Empty if method = 'numedges'. 
    """

    nodes = list(interactome.nodes)
    degrees = dict(interactome.degree)

    # Create degree bins to sample from
    bins, actual_degree_to_bin_index = netprop_zscore.get_degree_binning(degrees, 10)
    
    focal_genes = list(np.intersect1d(focal_genes, nodes)) # only use focal_genes which are in Gint

    numedges_list = []
    numedges_rand = []
    LCC_list = []
    LCC_rand = []
    
    if background_list is None:
        background_list = nodes

    sample_number = int(len(focal_genes) * sample_fraction)
    if sample_number == len(focal_genes):
        #Only perform determinate procedures one time
        subgraph = nx.subgraph(interactome, focal_genes)
        if method == 'numedges' or method == 'both':
            numedges_list = [len(subgraph.edges())] * num_reps
        if method == 'LCC' or method == 'both':
            largest_connected_component = max(nx.connected_component_subgraphs(subgraph), key=len)
            LCC_list = [max(nx.connected_component_subgraphs(subgraph), key=len)] * num_reps
        #Random genes
        for r in range(num_reps):
            random_genes = get_random_genes_of_similar_degree(focal_genes, bins, actual_degree_to_bin_index, degrees, background_list)
            #Find subgraph of random genes
            subgraph = nx.subgraph(interactome, random_genes)
            if method == 'numedges' or method == 'both':
                numedges_rand.append(len(subgraph.edges()))
            if method == 'LCC' or method == 'both':
                largest_connected_component = max(nx.connected_component_subgraphs(subgraph), key=len)
                LCC_rand.append(len(largest_connected_component.nodes))
    else:
        for r in range(num_reps):
            #Shuffle focal genes
            np.random.shuffle(focal_genes)
            focal_genes_subset = focal_genes[:sample_number]
            #Get subgraph
            subgraph = nx.subgraph(interactome, focal_genes_subset)
            if method == 'numedges' or method == 'both':
                numedges_list.append(len(subgraph.edges()))
            if method == 'LCC' or method == 'both':
                largest_connected_component = max(nx.connected_component_subgraphs(subgraph), key=len)
                LCC_list.append(len(largest_connected_component.nodes))
            #Get random genes
            random_genes = get_random_genes_of_similar_degree(focal_genes_subset, bins, actual_degree_to_bin_index, degrees, background_list)
            #Find subgraph of random genes
            subgraph = nx.subgraph(interactome, random_genes)
            if method == 'numedges' or method == 'both':
                numedges_rand.append(len(subgraph.edges()))
            if method == 'LCC' or method == 'both':
                largest_connected_component = max(nx.connected_component_subgraphs(subgraph), key=len)
                LCC_rand.append(len(largest_connected_component.nodes))
    
    if plot == True:
        if method == 'numedges' or method == 'both':
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(numedges_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(numedges_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('number of edges', fontsize = 16)
            plt.title('Number of Edges Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
                
        if method == 'LCC' or method == 'both':
            fig, ax = plt.subplots(figsize = (12, 7))
            sns.distplot(LCC_list, ax = ax, hist = True, label = 'focal genes')
            sns.distplot(LCC_rand, ax = ax, hist = True, label = 'random set')
            plt.ylabel('frequency', fontsize = 16)
            plt.xlabel('largest connected component size', fontsize = 16)
            plt.title('Largest Connected Component Localization', fontsize = 18)
            plt.legend(loc = 'upper right', fontsize = 14)
            
    return numedges_list, numedges_rand, LCC_list, LCC_rand

def get_random_genes_of_similar_degree(focal_genes, bins, actual_degree_to_bin_index, degrees, background_list):
    random_genes = []
    for gene in focal_genes:
        genes_of_similar_degree = bins[actual_degree_to_bin_index[degrees[gene]]]
        np.random.shuffle(genes_of_similar_degree)
        index = 0
        while genes_of_similar_degree[index] in random_genes or genes_of_similar_degree[index] not in background_list:
            index += 1
        random_genes.append(genes_of_similar_degree[index])
    return random_genes
    

def localization_full(Gint, focal_genes, 
                       num_reps = 200,
                       method = 'LCC', 
                       print_counter = False,
                       label = 'focal genes',
                       line_height = 0.1,
                       legend_loc = 'upper left'):
                       
    """
        Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
        Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LCC') 
        localization analysis. DOes no sub-sampling. Plots the distribution of random gene localizaiton, and 
        marks the focal set localization on distribution. Includes p-value of focal set localization.
        
        Args:
            Gint: Networkx Graph, background network to randomly sample from
            focal_genes: List, set of genes to calculate localization of
            num_reps: Int, number of times to randomly sample
            method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LCC', or 'both'.
            print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                           Useful when the num_reps is very high.
            label: String, label for focal genes in graph legend
            line_height: Float, the height of the red line that marks the focal gene localization
            legend_loc: String, relative position of legend in graph. Something similar to 'upper left'.
            
        Returns: 
            numedges_list: List, the number of edges calculated for each rep, over focal genes. 
                Empty if method = 'LCC'. 
            numedges_rand: List, the number of edges calculated for each rep, over random genes of 
                similar degree in the background network. Empty if method = 'LCC'.
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
            