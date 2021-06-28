# -*- coding: utf-8 -*-

'''Functions for performing network colocalization
'''

# External library imports
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import random
import seaborn as sns


def __init__(self):
    pass

def calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=3,
                              z1_threshold=1.5,z2_threshold=1.5):
    '''Function to determine which genes overlap. Returns a list of the 
    overlapping genes.
    
    Args:
    z_scores_1 (pandas.Series): Pandas Series resulting from the 
        netprop_zscore.netprop_zscore or netprop_zscore.calc_zscore_heat
        methods, containing the z-scores of each gene following network
        propagation. The index consists of gene names.
    z_scores_2 (pandas.Series): Similar to z_scores_1. The two pandas Series
        must contain the same genes (ie. come from the same interactome
        network).
    z_score_threshold (float): The threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded. (Default: 3)
    z1_threshold (float): The individual z1-score threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with z1-scores
        below this threshold will be discarded. (Default: 1.5)
    z2_threshold (float): The individual z2-score threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with z2-scores
        below this threshold will be discarded. (Default: 1.5)
    

    Returns:
        list: List of genes in the network overlap (genes with high combined
            z-scores).
    '''
    z_scores_1 = z_scores_1.to_frame(name='z_scores_1')
    z_scores_2 = z_scores_2.to_frame(name='z_scores_2')
    z_scores_joined = z_scores_1.join(z_scores_2)
    z_scores_combined = (z_scores_joined['z_scores_1'] 
                        * z_scores_joined['z_scores_2'] 
                        * (z_scores_joined['z_scores_1'] > 0) 
                        * (z_scores_joined['z_scores_2'] > 0))
    # get rid of unlikely genes which have low scores in either z1 or z2
    high_z_score_genes = z_scores_combined[
        (z_scores_combined >= z_score_threshold) 
         & (z_scores_joined['z_scores_1'] > z1_threshold) 
         & (z_scores_joined['z_scores_2'] > z2_threshold)
    ].index.tolist()
    
    return high_z_score_genes

def calculate_network_overlap_subgraph(interactome, z_scores_1, z_scores_2, z_score_threshold=3,
                                      z1_threshold=1.5,z2_threshold=1.5):
    '''Function to return subgraph of network intersection. 
    
    Code to create subgraph is from NetworkX documentation:
    https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.subgraph.html

    Args:
    interactome (NetworkX graph): The network whose subgraph will be returned.
    z_scores_1 (pandas.Series): Pandas Series resulting from the 
        netprop_zscore.netprop_zscore or netprop_zscore.calc_zscore_heat
        methods, containing the z-scores of each gene following network
        propagation. The index consists of gene names.
    z_scores_2 (pandas.Series): Similar to z_scores_1. The two pandas Series
        must contain the same genes (ie. come from the same interactome
        network).
    z_score_threshold (float): The threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded. (Default: 3)
    z1_threshold (float): The individual z1-score threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with z1-scores
        below this threshold will be discarded. (Default: 1.5)
    z2_threshold (float): The individual z2-score threshold to determine whether a gene is 
        a part of the network overlap or not. Genes with z2-scores
        below this threshold will be discarded. (Default: 1.5)

    Returns:
        NetworkX graph: Subgraph of the interactome containing only genes that
            are in the network intersection (genes with high combined z-scores).
    '''
    network_overlap = calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=z_score_threshold,
                                               z1_threshold=z1_threshold,z2_threshold=z1_threshold)
    
    # Create subgraph that has the same type as original graph
    network_overlap_subgraph = interactome.__class__()
    network_overlap_subgraph.add_nodes_from((node, interactome.nodes[node]) for node in network_overlap)
    if network_overlap_subgraph.is_multigraph():
        network_overlap_subgraph.add_edges_from((node, neighbor, key, dictionary)
            for node, neighbors in interactome.adj.items() if node in network_overlap
            for neighbor, item in neighbors.items() if neighbor in network_overlap
            for key, dictionary in item.items())
    else:
        network_overlap_subgraph.add_edges_from((node, neighbor, dictionary)
            for node, neighbors in interactome.adj.items() if node in network_overlap
            for neighbor, dictionary in neighbors.items() if neighbor in network_overlap)
    network_overlap_subgraph.graph.update(interactome.graph)
    
    return network_overlap_subgraph

def calculate_expected_overlap(z_scores_1, z_scores_2, gene_set_name_1='Gene Set 1', gene_set_name_2='Gene Set 2', 
                               z_score_threshold=3, z1_threshold=1.5,z2_threshold=1.5,
                               num_reps=1000, save_random_network_overlap=False, plot=False):    
    '''Function to determine size of expected network overlap by randomly
    shuffling gene names.

    Args:
        z_scores_1 (pandas.Series): Pandas Series resulting from the 
            netprop_zscore.netprop_zscore or netprop_zscore.calc_zscore_heat
            methods, containing the z-scores of each gene following network
            propagation. The index consists of gene names.
        z_scores_2 (pandas.Series): Similar to z_scores_1. The two pandas Series
            must contain the same genes (ie. come from the same interactome
            network).
        z_score_threshold (float): The threshold to determine whether a gene is 
            a part of the network overlap or not. Genes with combined z-scores
            below this threshold will be discarded. (Default: 3)
        z1_threshold (float): The individual z1-score threshold to determine whether a gene is 
            a part of the network overlap or not. Genes with z1-scores
            below this threshold will be discarded. (Default: 1.5)
        z2_threshold (float): The individual z2-score threshold to determine whether a gene is 
            a part of the network overlap or not. Genes with z2-scores
            below this threshold will be discarded. (Default: 1.5)
            num_reps (int): The number of times that gene names will be shuffled.
        plot (bool): If True, the distribution will be plotted. If False, it
            will not be plotted. (Default: False)

    Returns:
        float: 

    '''
    # Build a distribution of expected network overlap sizes by shuffling node names
    random_network_overlap_sizes = []
    z_scores_1_copy = z_scores_1.copy()
    z_scores_2_copy = z_scores_2.copy()
    gene_set_1 = z_scores_1.index.tolist()
    gene_set_2 = z_scores_2.index.tolist()
    for _ in range(num_reps):
        # Shuffle gene name labels
        np.random.shuffle(gene_set_1)
        z_scores_1_copy.index = gene_set_1

        np.random.shuffle(gene_set_2)
        z_scores_2_copy.index = gene_set_2

        random_size = len(calculate_network_overlap(z_scores_1_copy, z_scores_2_copy, z_score_threshold=z_score_threshold,
                                                   z1_threshold=z1_threshold,z2_threshold=z2_threshold))
        random_network_overlap_sizes.append(random_size)
    
    network_overlap_size = len(calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=z_score_threshold,
                                                        z1_threshold=z1_threshold,z2_threshold=z2_threshold))

    if plot:
        plt.figure(figsize=(5, 4))
        dfig = sns.histplot(random_network_overlap_sizes, label='Expected network intersection size')
        plt.vlines(network_overlap_size, ymin=0, ymax=dfig.dataLim.bounds[3], color='r', label='Observed network intersection size')
        plt.xlabel('Size of proximal subgraph, z > ' + str(z_score_threshold), fontsize=16)
        plt.legend(fontsize=12)
    
    return network_overlap_size, random_network_overlap_sizes

