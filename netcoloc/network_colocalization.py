# -*- coding: utf-8 -*-

'''Functions for performing network colocalization
'''


import warnings
import logging

# External library imports
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance

# NEW:
import ipycytoscape
import ipywidgets as widgets
from scipy.stats import hypergeom
from scipy.stats import norm


from netcoloc.netprop_zscore import *


def __init__(self):
    pass


logger = logging.getLogger(__name__)


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


def transform_edges(G, method='cosine_sim', edge_weight_threshold=0.95):
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

    if method not in ['cosine_sim']:  # update this if we add more methods
        warnings.warn('Error: ' + method + ' method not yet implemented')
        return

    # compute the adjacency matrix
    logging.info('computing the adjacency matrix...')

    nodelist = list(G.nodes())
    # get graph as adjacency matrix
    graph_as_adj = nx.to_numpy_array(G)

    # add transpose to matrix to remove edge direction
    adj_temp = graph_as_adj + np.transpose(graph_as_adj)

    # compute the cosine similarity
    logging.info('computing the cosine similarity...')
    cos_pc = pd.DataFrame(np.zeros((adj_temp.shape[0],
                                   adj_temp.shape[1])),
                          index=nodelist)
    cos_pc.columns = nodelist

    for i in np.arange(0, len(nodelist)-1):
        n1 = nodelist[i]
        # this node has no neighbors so set
        # cosine distance to maximum distance aka 1.0
        # in cos_pc and continue
        if max(adj_temp[i]) < 1.0:
            for j in np.arange(i + 1, len(nodelist)):
                n2 = nodelist[j]
                cos_pc.loc[n1][n2] = 1.0
                cos_pc.loc[n2][n1] = 1.0
            continue
        for j in np.arange(i+1, len(nodelist)):
            n2 = nodelist[j]

            # make sure they have some neighbors
            if max(adj_temp[j]) > 0.0:
                cosine_distance = distance.cosine(adj_temp[i],
                                                  adj_temp[j])
            else:
                # no neighbors in neigh2 so set
                # cosine distance to maximum distance aka 1.0
                cosine_distance = 1.0

            cos_pc.loc[n1][n2] = cosine_distance
            cos_pc.loc[n2][n1] = cosine_distance

    # Rank transform 1-cos distance
    logger.info('rank transforming...')
    m1cos = 1-cos_pc
    m1cos = m1cos.replace(np.nan, 0)
    sim_names = m1cos.index.tolist()
    sim_rank = m1cos.rank(0) / (m1cos.shape[0] - 1)
    sim_rank = np.matrix(sim_rank)
    np.fill_diagonal(sim_rank,0) # remove self edges
    sim_rank = pd.DataFrame(sim_rank)
    sim_rank = pd.DataFrame((sim_rank.values + sim_rank.values.T) / 2.0,
                            columns=sim_names, index=sim_names)

    # remove self edges
    #sim_rank.values[[np.arange(sim_rank.shape[0])]*2] = 0


    sim_rank['gene_temp'] = sim_rank.index.tolist()
    sim_rank_EL = sim_rank.melt(id_vars=['gene_temp'])
    sim_rank_EL.columns = ['node1', 'node2', 'sim']
    sim_rank_EL = sim_rank_EL[sim_rank_EL['sim'] > edge_weight_threshold]
    G_transf = nx.Graph()
    G_transf.add_nodes_from(G)
    G_transf.add_weighted_edges_from(zip(sim_rank_EL['node1'],
                                         sim_rank_EL['node2'],
                                         sim_rank_EL['sim']))

    logger.info('number of transformed edges returned = ' +
                str(len(G_transf.edges())))

    return G_transf

def view_G_hier(G_hier,layout='cose'):
    
    """
    In-notebook visualization of NetColoc hierarchy, using ipycytoscape.

    :param G_hier: network to visualize. Expects output of  cdapsutil.CommunityDetection(), transformed to networkx format. 'CD_MemberList_LogSize' is a required field of the network to map to the node size.
    :type G: :py:class:`networkx.Graph`
    :param layout: Layout method to use, any layout supported by cytoscape.js is supported. Suggest 'cose' or 'breadthfirst'.
    :type layout: str
    :param edge_weight_threshold: Transformed edges will be returned which have values greater than this
    :type edge_weight_threshold: float
    :return: Nothing
    """
    
    # add a new field to map to node size
    new_vals = nx.get_node_attributes(G_hier,'CD_MemberList_LogSize')
    new_vals_keys = list(new_vals.keys())
    new_vals_vals = [np.float64(v)*10 for v in list(new_vals.values())]
    new_vals = dict(zip(new_vals_keys,[str(v) for v in new_vals_vals]))
    new_vals
    nx.set_node_attributes(G_hier,new_vals,name='CD_MemberList_LogSize_viz')
    
    ipyviewer = ipycytoscape.CytoscapeWidget()
    ipyviewer.graph.add_graph_from_networkx(G_hier)
    ipyviewer.set_style([{
                            'selector': 'node',
                            'css': {
                                'content': 'data(name)',
                                'text-valign': 'center',
                                'color': 'white',
                                'text-outline-width': 2,
                                'text-outline-color': 'green',
                                'background-color': 'green',
                                'width': 'data(CD_MemberList_LogSize_viz)',
                                'height': 'data(CD_MemberList_LogSize_viz)'
                            }
                            }
                            ])

    #ipyviewer.set_layout(name='breadthfirst',directed='true')
    ipyviewer.set_layout(name='cose')
    display(ipyviewer)
    
def sweep_input_pvals(D1_df,D2_df,
                      individual_heats_matrix,nodes,degrees,
                      cutoff_list = [.01,.02,.03,.04,.05,.1],
                      gene_column='gene',
                      score_column='pval',
                      cutoff_max=True,
                     num_reps=100,
                     verbose=True,
                     z_score_threshold=3,
                     z12_threshold=1.5):
    
    """
    Evaluate NetColoc enrichment for a range of thresholds on input gene lists.
    
    :param D1_df: DataFrame containing gene names and scores for the first gene set
    :type D1_df: :py:class:`pandas.DataFrame`
    :param D2_df: DataFrame containing gene names and scores for the second gene set
    :type D2_df: :py:class:`pandas.DataFrame`
    :param individual_heats_matrix: output of the
            netprop.get_individual_heats_matrix. A square matrix containing the
            final heat contributions of each gene
    :type individual_heats_matrix: :py:class:`numpy.ndarray`
    :param nodes: nodes, in the order in which they were supplied to
            the :py:func:`~netcoloc.netprop.get_normalized_adjacency_matrix` method
            which returns the precursor to the individual_heats_matrix
    :type nodes: list
    :param degrees: Mapping of node names to node degrees
    :type degrees: dict
    :param cutoff_list: list of values to threshold the input gene sets by
    :type cutoff_list: list
    :param gene_column: name of column containing genes in D1_df and D2_df
    :type gene_column: string
    :param score_column: name of column containing scores (usually p-value or log fold change) in D1_df and D2_df
    :type score_column: string
    :param cutoff_max: if True, genes will be selected which have scores less than the cutoff value. If false, genes will be selected which have scores greater than the cutoff value.
    :type cutoff_max: Boolean
    :param num_reps: Number of times the network propagation algorithm should
        be run using random seed genes in order to build the null model
    :type num_reps: int
    :param verbose: if True, print out some diagnostics
    :type verbose: Boolean
    :param z_score_threshold: threshold to determine whether a gene is
        a part of the network overlap or not. Genes with combined z-scores
        below this threshold will be discarded
    :type z_score_threshold: float
    :param z12_threshold: individual z1/z2-score threshold to determine whether a gene is
        a part of the network overlap or not. Genes with z1/z2-scores
        below this threshold will be discarded
    :type z12_threshold: float

    :return: netcoloc_pval_df: DataFrame containing NetColoc enrichment results
    :rtype: :py:class:`pandas.DataFrame`
    """
    
    D1_num_genes, D2_num_genes = [],[]
    num_shared_genes = []
    plist = []
    obs_overlap_list=[]
    network_exp_mean_overlap_list=[]
    network_exp_std_overlap_list=[]
    
    for c in cutoff_list:

        if cutoff_max==True:
            D1_genes = D1_df[D1_df[score_column]<c][gene_column].tolist()
            D2_genes = D2_df[D2_df[score_column]<c][gene_column].tolist()
        else:
            D1_genes = D1_df[D1_df[score_column]>c][gene_column].tolist()
            D2_genes = D2_df[D2_df[score_column]>c][gene_column].tolist()
            
        D1_num_genes.append(len(D1_genes))
        D2_num_genes.append(len(D2_genes))
        num_shared_genes.append(len(np.intersect1d(D1_genes,D2_genes)))
            
        # check that the gene sizes meet the input criteria of >5 and <500
        if (5<=len(D1_genes)<=500) & (5<=len(D2_genes)<=500):
            
            # D1 variant network propagation
            print('\nCalculating D1 variant z-scores: ')
            z_D1, Fnew_D1, Fnew_rand_D1 = calculate_heat_zscores(individual_heats_matrix, nodes, 
                                                                                degrees, 
                                                                                D1_genes, num_reps=num_reps)

            z_D1 = pd.DataFrame({'z':z_D1})
            z_D1.sort_values('z',ascending=False).head()

            # D2 variant network propagation
            print('\nCalculating D2 variant z-scores: ')
            z_D2, Fnew_D2, Fnew_rand_D2 = calculate_heat_zscores(individual_heats_matrix, nodes, 
                                                                                degrees, 
                                                                                D2_genes, num_reps=num_reps)

            z_D2 = pd.DataFrame({'z':z_D2})

            z_d1d2_size, high_z_rand = calculate_expected_overlap(
                z_D1['z'],
                z_D2['z'],
                plot=False,
                num_reps=num_reps,
                z_score_threshold=z_score_threshold,
                z1_threshold=z12_threshold,
                z2_threshold=z12_threshold
            )

            ztemp = (z_d1d2_size - np.mean(high_z_rand)) / np.std(high_z_rand)
            ptemp = norm.sf(ztemp)

            # save the number of overlapping genes and overlap p-value
            network_num_overlap = z_d1d2_size
            network_pval_overlap = ptemp
            
            network_exp_mean_overlap = np.mean(high_z_rand)
            network_exp_std_overlap = np.std(high_z_rand)

            
            obs_exp_temp = float(z_d1d2_size) / np.mean(high_z_rand)
            
        else:
            ptemp = np.nan
            network_num_overlap=np.nan
            network_exp_mean_overlap=np.nan
            network_exp_std_overlap=np.nan
            z_d1d2_size=np.nan
            obs_exp_temp=np.nan
            
        plist.append(ptemp)
        obs_overlap_list.append(network_num_overlap)
        network_exp_mean_overlap_list.append(network_exp_mean_overlap)
        network_exp_std_overlap_list.append(network_exp_std_overlap)
            
        if verbose:
            print('\ncutoff = '+str(c))
            print('number of D1 genes = '+str(len(D1_genes)))
            print('number of D2 genes = '+str(len(D2_genes)))
            print('number of shared genes = '+str(len(np.intersect1d(D1_genes,D2_genes))))
            
            
            print('size of network intersection = ' + str(z_d1d2_size))
            print('observed size/ expected size = ' + str(obs_exp_temp))
            print('p = ' + str(ptemp))
            



    netcoloc_pval_df = pd.DataFrame({'score_cutoff':cutoff_list,
                                           'D1_num_genes':D1_num_genes,
                                           'D2_num_genes':D2_num_genes,
                                            'num_shared_genes':num_shared_genes,
                                          'observed_overlap':obs_overlap_list,
                                          'expected_overlap_mean':network_exp_mean_overlap_list,
                                          'expected_overlap_std':network_exp_std_overlap_list,
                                           'empirical_p':plist
                                          })
    netcoloc_pval_df['obs_exp']=netcoloc_pval_df['observed_overlap']/netcoloc_pval_df['expected_overlap_mean']

            
    
    return netcoloc_pval_df



def calculate_network_enrichment(z_D1,z_D2,zthresh_list = list(np.arange(1,15)),z12thresh_list=[1,1.5,2],
                                 verbose=True):

    """
    Evaluate NetColoc enrichment for a range of thresholds on network proximity z-scores.
    
    :param z_D1: DataFrame containing gene names and network proximity z-scores for the first gene set
    :type z_D1: :py:class:`pandas.DataFrame`
    :param z_D2: DataFrame containing gene names and network proximity z-scores for the second gene set
    :type z_D2: :py:class:`pandas.DataFrame`
    :param zthresh_list: list of combined z-score thresholds to iterate over
    :type zthresh_list: :list
    :param z12thresh_list: list of individual z-score thresholds to iterate over
    :type z12thresh_list: :list
    :param verbose: if True, print out some diagnostics
    :type verbose: Boolean
    
    :return: netcoloc_enrichment_df: DataFrame containing NetColoc enrichment results
    :rtype: :py:class:`pandas.DataFrame`
    """
    
    plist = []
    obs_overlap_list=[]
    network_exp_mean_overlap_list=[]
    network_exp_std_overlap_list=[]
    zlist=[]
    z12list=[]

    for zthresh in zthresh_list:
        for z12 in z12thresh_list:

            z_d1d2_size, high_z_rand = calculate_expected_overlap(
                z_D1['z'],
                z_D2['z'],
                plot=False,
                num_reps=100,
                z_score_threshold=zthresh,
                z1_threshold=z12,
                z2_threshold=z12
            )

            ztemp = (z_d1d2_size - np.mean(high_z_rand)) / np.std(high_z_rand)
            ptemp = norm.sf(ztemp)
            plist.append(ptemp)

            # save the number of overlapping genes and overlap p-value
            network_num_overlap = z_d1d2_size
            obs_overlap_list.append(network_num_overlap)
            network_pval_overlap = ptemp
            obs_exp_temp = float(z_d1d2_size) / np.mean(high_z_rand)
            network_obs_exp = obs_exp_temp
            network_exp_mean_overlap = np.mean(high_z_rand)
            network_exp_mean_overlap_list.append(network_exp_mean_overlap)
            network_exp_std_overlap = np.std(high_z_rand)
            network_exp_std_overlap_list.append(network_exp_std_overlap)

            zlist.append(zthresh)
            z12list.append(z12)

            if verbose:
                print('z = '+str(zthresh))
                print('z12 = '+str(z12))
                print('size of network intersection = ' + str(z_d1d2_size))
                print('observed size/ expected size = ' + str(obs_exp_temp))
                print('p = ' + str(ptemp))



    netcoloc_enrichment_df = pd.DataFrame({'z_comb':zlist,
                                           'z_12':z12list,
                                          'observed_overlap':obs_overlap_list,
                                          'expected_overlap_mean':network_exp_mean_overlap_list,
                                          'expected_overlap_std':network_exp_std_overlap_list,
                                           'empirical_p':plist
                                          })
    netcoloc_enrichment_df['obs_exp']=netcoloc_enrichment_df['observed_overlap']/netcoloc_enrichment_df['expected_overlap_mean']
    
    return netcoloc_enrichment_df

    

