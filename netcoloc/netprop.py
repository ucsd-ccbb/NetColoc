# -*- coding: utf-8 -*-

'''Functions for performing network propagation
'''

import networkx as nx
import numpy as np
import pandas as pd
import warnings


def __init__(self):
    pass


def get_normalized_adjacency_matrix(graph, conserve_heat=True, weighted=False):
    """
    Returns normalized adjacency matrix (W'), as detailed in:

    Vanunu, Oron, et al. 'Associating genes and protein complexes with disease
    via network propagation.'


    With version `0.1.6` and newer, the :py:class:`networkx.Graph`
    can be directly passed into
    :py:func:`~netcoloc.netprop.get_individual_heats_matrix` and
    this method will be invoked to create the normalized adjacency matrix

    .. note::
        Resulting matrix from this function can be saved to a file with :py:func:`numpy.save`
        and loaded later with :py:func:`numpy.load`, but resulting file can be several gigabytes
        and take a minute or more to save/load.

        .. code-block:: python

            numpy.save('nam.npy', adjacency_matrix)
            adjacency_matrix = numpy.load('nam.npy')


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
    assert 0 not in dict(graph.degree).values(), "Graph cannot have nodes with degree=zero"

    # Create graph
    if conserve_heat:
        # If conserving heat, make G_weighted a di-graph (not symmetric)
        graph_weighted = nx.DiGraph()
    else:
        # If not conserving heat, make G_weighted a simple graph (symmetric)
        graph_weighted = nx.Graph()

    # Create edge weights
    edge_weights = []
    node_to_degree_dict = dict(graph.degree)
    if weighted and not nx.is_weighted(G=graph):
        warnings.warn("Input graph is not weighted. All edge weights will be set to 1.")

    for e in graph.edges(data=True):
        v1 = e[0]
        v2 = e[1]
        deg1 = node_to_degree_dict[v1]
        deg2 = node_to_degree_dict[v2]

        if weighted and nx.is_weighted(G=graph):
            weight = e[2]['weight']
        else:
            weight = 1

        if conserve_heat:
            # created asymmetrically weighted edges - each directed edge u->v normalized by the degree of v
            edge_weights.append((v1, v2, weight / float(deg1)))
            edge_weights.append((v2, v1, weight / float(deg2)))
        else:
            # normalize single undirected edge by the degree of both endpoints as per Vanunu, Oron, et al. 2010
            edge_weights.append((v1, v2, weight / np.sqrt(deg1 * deg2)))

    # Apply edge weights to graph
    graph_weighted.add_weighted_edges_from(edge_weights)

    # Transform graph to adjacency array
    if len(graph.nodes) != len(graph_weighted):
        raise ValueError("Input graph has nodes with zero degrees. Please remove these nodes.")

    w_prime = nx.to_numpy_array(graph_weighted, nodelist=graph.nodes())

    return w_prime


def get_individual_heats_matrix(nam_or_graph, alpha=0.5,
                                conserve_heat=True, weighted=False):
    """
    Returns the pre-calculated contributions of each individual gene in the
    interactome to the final heat of each other gene in the interactome after
    propagation.

    .. versionchanged:: 0.1.6
        In addition, to a normalized adjacency matrix, this function
        now also supports :py:class:`networkx.Graph` network as input


    If a :py:class:`networkx.Graph` network is passed in as the **nam_or_graph**
    parameter, the function :py:func:`~netcoloc.netprop.get_normalized_adjacency_matrix`
    is called to generate the normalized adjacency matrix using **conserve_heat** and
    **weighted** parameters

    .. note::
        Resulting matrix from this function can be saved to a file with :py:func:`numpy.save`
        and loaded later with :py:func:`numpy.load`, but resulting file can be several gigabytes
        and take a minute or more to save/load.

        .. code-block:: python

            numpy.save('heats_matrix.npy', w_double_prime)
            w_double_prime = numpy.load('heats_matrix.npy')


    :param nam_or_graph: square normalized
            adjacency matrix or network
    :type nam_or_graph: :py:class:`numpy.ndarray` or :py:class:`networkx.Graph`
    :param alpha: heat dissipation coefficient between 1 and 0. The
           contribution of the heat propagated from adjacent nodes in
           determining the final heat of a node, as opposed to the contribution
           from being a part of the gene set initially
    :type alpha: float
    :param conserve_heat: If ``True``, heat will be conserved
            (ie. the sum of the heat vector will be equal to 1),
            and the graph will be asymmetric. Otherwise, heat will
            not be conserved, and the graph will be symmetric.
            **NOTE:** Only applies if **nam_or_graph** is :py:class:`networkx.Graph`
    :type conserve_heat: bool
    :param weighted: If ``True``, then the graph's edge weights
            will be taken into account. Otherwise, all edge weights
            will be set to 1.
            **NOTE:** Only applies if **nam_or_graph** is :py:class:`networkx.Graph`
    :type weighted: bool
    :return: square individual heats matrix
    :rtype: :py:class:`numpy.ndarray`
    """
    assert 1 >= alpha >= 0, "Alpha must be between 0 and 1"

    nam = nam_or_graph
    if isinstance(nam_or_graph, nx.Graph):
        nam = get_normalized_adjacency_matrix(nam_or_graph,
                                              conserve_heat=conserve_heat,
                                              weighted=weighted)
    nam = np.transpose(nam)

    d_name = np.linalg.inv(np.identity(nam.shape[0]) - alpha * nam) * (1 - alpha)

    return d_name


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
    :param seed_genes: # TODO
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
        # TODO check that this is the correct orientation
        F += individual_heats_matrix[:,nodes.index(gene)]

    # Normalize results by number of seed genes
    F /= len(seed_genes)

    #Return as pandas series
    # TODO does this need to be a pandas series?
    return pd.Series(F, index=nodes)
