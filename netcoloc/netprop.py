# -*- coding: utf-8 -*-

'''Functions for performing network propagation
'''

# External library imports
import networkx as nx
import numpy as np
import pandas as pd


def __init__(self):
    pass


def get_normalized_adjacency_matrix(graph, conserve_heat=True, weighted=False):
    """
    Returns normalized adjacency matrix (W'), as detailed in:

    Vanunu, Oron, et al. 'Associating genes and protein complexes with disease
    via network propagation.'

    :param graph: Interactome from which to calculate normalized
            adjacency matrix.
    :type graph: :py:class:`networkx.Graph`
    :param conserve_heat: If ``True``, heat will be conserved
            (ie. the sum of the heat vector will be equal to 1),
            and the graph will be asymmetric. Otherwise, heat will
            not be conserved, and the graph will be symmetric.
    :type conserve_heat: bool
    :param weighted: If ``True``, then the graph's edge weights
            will be taken into account. Otherwise, all edge weights
            will be set to 1.
    :type weighted: bool
    :return: Square normalized adjacency matrix
    :rtype: :py:class:`numpy.ndarray`
    """
    # Create graph
    if conserve_heat:
        # If conserving heat, make G_weighted a di-graph (not symmetric)
        graph_weighted= nx.DiGraph()
    else:
        # If not conserving heat, make G_weighted a simple graph (symmetric)
        graph_weighted = nx.Graph()

    #Create edge weights
    edge_weights = []
    node_to_degree_dict = dict(graph.degree)
    for e in graph.edges(data=True):
        v1 = e[0]
        v2 = e[1]
        deg1 = node_to_degree_dict[v1]
        deg2 = node_to_degree_dict[v2]
        
        if weighted:
            weight = e[2]['weight']
        else:
            weight = 1
        
        if conserve_heat:
            edge_weights.append((v1, v2, weight / float(deg2)))
            edge_weights.append((v2, v1, weight / float(deg1)))
        else:
            edge_weights.append((v1, v2, weight / np.sqrt(deg1 * deg2)))
    
    # Apply edge weights to graph
    graph_weighted.add_weighted_edges_from(edge_weights)
    
    # Transform graph to adjacency matrix
    w_prime = nx.to_numpy_matrix(graph_weighted, nodelist=graph.nodes())
    w_prime = np.array(w_prime)
    
    return w_prime


def get_individual_heats_matrix(normalized_adjacency_matrix, alpha=0.5):
    """
    Returns the pre-calculated contributions of each individual gene in the
    interactome to the final heat of each other gene in the interactome after
    propagation.

    :param normalized_adjacency_matrix: square normalized
            adjacency matrix from :py:func:`~netcoloc.netprop.get_normalized_adjacency_matrix`
    :type normalized_adjacency_matrix: :py:class:`numpy.ndarray`
    :param alpha: heat dissapation coefficient between 1 and 0. The
           contribution of the heat propagated from adjacent nodes in
           determining the final heat of a node, as opposed to the contribution
           from being a part of the gene set initially
    :type alpha: float
    :return: square individual heats matrix
    :rtype: :py:class:`numpy.ndarray`
    """
    return np.linalg.inv(
        np.identity(normalized_adjacency_matrix.shape[0]) 
        - alpha * normalized_adjacency_matrix
    ) * (1 - alpha)


def network_propagation(individual_heats_matrix, nodes, seed_genes):
    """
    Implements network propagation, as detailed in:

    Vanunu, Oron, et al. 'Associating genes and protein complexes with
    disease via network propagation.'

    Using this function, the final heat of the network is calculated directly,
    instead of iteratively. This method is faster when many different
    propagations need to be performed on the same network (with different seed
    gene sets). It is slower than
    :py:func:`~netcoloc.netprop.iterative_network_propagation` for a
    single propagation.

    :param individual_heats_matrix: Square matrix that is the
                                    output of :py:func:`~netcoloc.netprop.get_individual_heats_matrix`
    :type individual_heats_matrix: :py:class:`numpy.ndarray`
    :param nodes: List of nodes in the network represented by the
            individual_heats_matrix, in the same order in which they were
            supplied to :py:func:`~netcoloc.netprop.get_individual_heats_matrix`
    :type nodes: list
    :param seed_genes:
    :return: Final heat of each node after propagation, with the name
             of the nodes as the index
    :rtype: :py:class:`pandas.Series`
    """
    # Remove genes that are not in network
    seed_genes = list(np.intersect1d(nodes, seed_genes))

    # Initialize results vector
    F = np.zeros(len(nodes))

    # Add up resulting heats from each gene in seed genes set
    for gene in seed_genes:
        F += individual_heats_matrix[:,nodes.index(gene)]

    # Normalize results by number of seed genes
    F /= len(seed_genes)

    #Return as pandas series
    return pd.Series(F, index=nodes)
