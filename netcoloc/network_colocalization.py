# -*- coding: utf-8 -*-

'''Functions for performing network colocalization
'''


import warnings

# External library imports
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance


# gprofiler prelim annotation
from gprofiler import GProfiler


gp = GProfiler("MyToolName/0.1")


def __init__(self):
    pass


def calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=3,
                              z1_threshold=1.5, z2_threshold=1.5):
    """
    Function to determine which genes overlap. Returns a list of the
    overlapping genes

    :param z_scores_1: Result from :py:func:`~netcoloc.netprop_zscore.netprop_zscore`
                       or :py:func:`~netcoloc.netprop_zscore.calculate_heat_zscores`
                       containing the z-scores of each gene following network
                       propagation. The index consists of gene names
    :type z_scores_1: :py:class:`pandas.Series`
    :param z_scores_2: Similar to **z_scores_1**. This and **z_scores_1**
                       must contain the same genes (ie. come from the same
                       interactome network)
    :type z_scores_2: :py:class:`pandas.Series`
    :param z_score_threshold: threshold to determine whether a gene is
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded
    :type z_score_threshold: float
    :param z1_threshold: individual z1-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z1-scores
        below this threshold will be discarded
    :type z1_threshold: float
    :param z2_threshold: individual z2-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z2-scores
        below this threshold will be discarded
    :type z2_threshold: float
    :return: genes in the network overlap (genes with high combined
            z-scores)
    :rtype: list
    """
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


def calculate_network_overlap_subgraph(interactome, z_scores_1,
                                       z_scores_2, z_score_threshold=3,
                                       z1_threshold=1.5,z2_threshold=1.5):
    """
    Function to return subgraph of network intersection.

    Code to create subgraph is from
    `NetworkX subgraph documentation <https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.subgraph.html>`__

    :param interactome: network whose subgraph will be returned
    :type interactome: :py:class:`networkx.Graph`
    :param z_scores_1: Result from :py:func:`~netcoloc.netprop_zscore.netprop_zscore`
                       or :py:func:`~netcoloc.netprop_zscore.calculate_heat_zscores`
                       containing the z-scores of each gene following network
                       propagation. The index consists of gene names
    :type z_scores_1: :py:class:`pandas.Series`
    :param z_scores_2: Similar to **z_scores_1**. This and **z_scores_1**
                       must contain the same genes (ie. come from the same
                       interactome network)
    :type z_scores_2: :py:class:`pandas.Series`
    :param z_score_threshold: threshold to determine whether a gene is
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded
    :type z_score_threshold: float
    :param z1_threshold: individual z1-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z1-scores
        below this threshold will be discarded
    :type z1_threshold: float
    :param z2_threshold: individual z2-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z2-scores
        below this threshold will be discarded
    :type z2_threshold: float
    :return: Subgraph of the interactome containing only genes that
            are in the network intersection (genes with high combined z-scores)
    :rtype: :py:class:`networkx.Graph`
    """
    network_overlap = calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold=z_score_threshold,
                                               z1_threshold=z1_threshold,z2_threshold=z2_threshold)

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


def calculate_expected_overlap(z_scores_1, z_scores_2,
                               z_score_threshold=3, z1_threshold=1.5,
                               z2_threshold=1.5,
                               num_reps=1000, plot=False):
    """
    Determines size of expected network overlap by randomly
    shuffling gene names

    :param z_scores_1: Result from :py:func:`~netcoloc.netprop_zscore.netprop_zscore`
                       or :py:func:`~netcoloc.netprop_zscore.calculate_heat_zscores`
                       containing the z-scores of each gene following network
                       propagation. The index consists of gene names
    :type z_scores_1: :py:class:`pandas.Series`
    :param z_scores_2: Similar to **z_scores_1**. This and **z_scores_1**
                       must contain the same genes (ie. come from the same
                       interactome network)
    :type z_scores_2: :py:class:`pandas.Series`
    :param z_score_threshold: threshold to determine whether a gene is
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded
    :type z_score_threshold: float
    :param z1_threshold: individual z1-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z1-scores
        below this threshold will be discarded
    :type z1_threshold: float
    :param z2_threshold: individual z2-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z2-scores
        below this threshold will be discarded
    :type z2_threshold: float
    :param num_reps:
    :param plot: If ``True``, distribution will be plotted
    :type plot: bool
    :return:
    :rtype: float
    """
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

        random_size = len(calculate_network_overlap(z_scores_1_copy, z_scores_2_copy,
                                                    z_score_threshold=z_score_threshold,
                                                    z1_threshold=z1_threshold,
                                                    z2_threshold=z2_threshold))
        random_network_overlap_sizes.append(random_size)

    network_overlap_size = len(calculate_network_overlap(z_scores_1, z_scores_2,
                                                         z_score_threshold=z_score_threshold,
                                                         z1_threshold=z1_threshold,
                                                         z2_threshold=z2_threshold))

    if plot:
        plt.figure(figsize=(5, 4))
        dfig = sns.histplot(random_network_overlap_sizes,
                            label='Expected network intersection size')
        plt.vlines(network_overlap_size, ymin=0, ymax=dfig.dataLim.bounds[3], color='r',
                   label='Observed network intersection size')
        plt.xlabel('Size of proximal subgraph, z > ' + str(z_score_threshold),
                   fontsize=16)
        plt.legend(fontsize=12)

    return network_overlap_size, random_network_overlap_sizes


def transform_edges(G,method='cosine_sim', edge_weight_threshold=0.95):
    """
    Transforms binary edges using selected method (currently only cosine similarity is implemented).
    Cosine similarity measures the similarity between neighbors of node pairs in the input network

    :param G: network whose edges will be transformed
    :type G: :py:class:`networkx.Graph`
    :param method: Method to use, only ``cosine_sim`` supported. Any other value will
                   cause this method to output a warning and immediately return
    :type method: str
    :param edge_weight_threshold: Transformed edges will be returned which have values greater than this
    :type edge_weight_threshold: float
    :return: Graph with nodes identical to input G, but with transformed edges (values > edge_weight_threshold)
    :rtype: :py:class:`networkx.Graph`
    """

    if not method in ['cosine_sim']: # update this if we add more methods
        warnings.warn('Error: ' + method + ' method not yet implemented')
        return

    # compute the adjacency matrix
    print('computing the adjacency matrix...')
    adj_temp = pd.DataFrame(nx.to_numpy_matrix(G))
    adj_temp.index = G.nodes()
    adj_temp.columns = G.nodes()

    nodelist = list(G.nodes())

    # compute the cosine similarity
    print('computing the cosine similarity...')
    cos_pc = pd.DataFrame(np.zeros((len(nodelist),len(nodelist))),index=nodelist)
    cos_pc.columns = nodelist

    counter = -1
    for i in np.arange(len(nodelist)-1):
        n1 = nodelist[i]
        neigh1 = list(nx.neighbors(G,n1))
        counter += 1
        for j in np.arange(i+1,len(nodelist)):
            n2 = nodelist[j]
            neigh2 = list(nx.neighbors(G,n2))
            # make sure they have some neighbors
            if len(np.union1d(neigh1,neigh2))>0:
                cos_sim_temp = distance.cosine(adj_temp[n1],adj_temp[n2])
            else:
                cos_sim_temp = 1

            cos_pc.loc[n1][n2] = cos_sim_temp
            cos_pc.loc[n2][n1] = cos_sim_temp

    # Rank transform 1-cos distance
    print('rank transforming...')
    m1cos = 1-cos_pc
    m1cos = m1cos.replace(np.nan,0)
    sim_names = m1cos.index.tolist()
    sim_rank = m1cos.rank(0) / (m1cos.shape[0] - 1)

    sim_rank = pd.DataFrame((sim_rank.values + sim_rank.values.T) / 2.0, columns=sim_names, index=sim_names)

    # remove self edges
    sim_rank.values[[np.arange(sim_rank.shape[0])]*2] = 0

    sim_rank['gene_temp'] = sim_rank.index.tolist()
    sim_rank_EL = sim_rank.melt(id_vars=['gene_temp'])
    sim_rank_EL.columns = ['node1','node2','sim']
    sim_rank_EL = sim_rank_EL[sim_rank_EL['sim'] > edge_weight_threshold]

    G_transf = nx.Graph()
    G_transf.add_nodes_from(G)
    G_transf.add_weighted_edges_from(zip(sim_rank_EL['node1'],
                                         sim_rank_EL['node2'],
                                         sim_rank_EL['sim']))

    print('number of transformed edges returned = ')
    print(len(G_transf.edges()))

    return G_transf

