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
    '''Returns normalized adjacency matrix (W'), as detailed in:
    Vanunu, Oron, et al. 'Associating genes and protein complexes with disease
    via network propagation.'

    Args:
        graph (NetworkX graph): Interactome from which to calculate normalized
            adjacency matrix.
        conserve_heat (bool): If this is set to True, heat will be conserved
            (ie. the sum of the heat vector will be equal to 1), and the graph
            will be asymmetric. Otherwise, heat will not be conserved, and the
            graph will be symmetric. (Default: True)
        weighted (bool): If this is set to true, then the graph's edge weights
            will be taken into account. Otherwise, all edge weights will be set 
            to 1. (Default: False)

    Returns:
        numpy.ndarray: A square normalized adjacency matrix
    '''
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
    '''Returns the pre-calculated contributions of each individual gene in the
    interactome to the final heat of each other gene in the interactome after
    propagation.

    Args:
        normalized_adjacency_matrix (numpy.ndarray): A square normalized 
            adjacency matrix from netprop.get_normalized_adjacency_matrix()
        alpha (float): A heat dissapation coefficient between 1 and 0. The 
           contribution of the heat propagated from adjacent nodes in 
           determining the final heat of a node, as opposed to the contribution 
           from being a part of the gene set initially. (Default: 0.5)

    Returns:
        numpy.ndarray: A square individual heats matrix.
    '''

    return np.linalg.inv(
        np.identity(normalized_adjacency_matrix.shape[0]) 
        - alpha * normalized_adjacency_matrix
    ) * (1 - alpha)

def network_propagation(individual_heats_matrix, nodes, seed_genes):
    '''Implements network propagation, as detailed in: Vanunu, Oron, et al. 
    'Associating genes and protein complexes with disease via network 
    propagation.'

    Using this function, the final heat of the network is calculated directly, 
    instead of iteratively. This method is faster when many different 
    propagations need to be performed on the same network (with different seed 
    gene sets). It is slower than netprop.iterative_network_propagation() for a 
    single propagation.

    Args:
        individual_heats_matrix (numpy.ndarray): Square matrix that is the 
            output of netprop.get_individual_heats_matrix().
        nodes (list): List of nodes in the network represented by the
            individual_heats_matrix, in the same order in which they were 
            supplied to netprop.get_individual_heats_matrix().

    Returns: 
        pandas.Series: Final heat of each node after propagation, with the name
            of the nodes as the index.
    '''

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