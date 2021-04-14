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

def calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=3):
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
         & (z_scores_joined['z_scores_1'] > 1.5) 
         & (z_scores_joined['z_scores_2'] > 1.5)
    ].index.tolist()
    
    return high_z_score_genes

def calculate_network_overlap_subgraph(interactome, z_scores_1, z_scores_2, z_score_threshold=3):
    '''Function to return subgraph of network intersection.

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

    Returns:
        NetworkX graph: Subgraph of the interactome containing only genes that
            are in the network intersection (genes with high combined z-scores).
    '''
    network_overlap = calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=z_score_threshold)
    network_overlap_subgraph = nx.subgraph(interactome, network_overlap)
    
    return(network_overlap_subgraph)

def calculate_expected_overlap(z_scores_1, z_scores_2, gene_set_name_1='Gene Set 1', gene_set_name_2='Gene Set 2', 
                               z_score_threshold=3, num_reps=1000, plot=False):    
    '''Function to determine size of expected network overlap.

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
    num_reps (int):
    plot (bool): 
    '''

    high_z_rand = []
    z_scores_1_copy = z_scores_1.copy()
    z_scores_2_copy = z_scores_2.copy()
    gene_set_1 = z_scores_1.index.tolist()
    gene_set_2 = z_scores_2.index.tolist()
    for _ in range(num_reps):
        # use permutation shuffling method instead of Fnew comparison
        np.random.shuffle(gene_set_1)
        z_scores_1_copy.index = gene_set_1

        np.random.shuffle(gene_set_2)
        z_scores_2_copy.index = gene_set_2

        high_z_temp = len(calculate_network_overlap(z_scores_1_copy, z_scores_2_copy, z_score_threshold=z_score_threshold))
        high_z_rand.append(high_z_temp)
    
    network_overlap_size =len(calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=z_score_threshold))

    if plot==True:
        sns.distplot(high_z_rand, label='expected network intersection size')
        plt.plot([network_overlap_size, network_overlap_size], [0, 0.015], label='observed ' + gene_set_name_1 + '-' + gene_set_name_2 + ' network intersection size')
        plt.xlabel('size of proximal subgraph, z >= ' + str(z_score_threshold), fontsize=16)
        plt.legend(fontsize=12)
    
    return network_overlap_size, high_z_rand

