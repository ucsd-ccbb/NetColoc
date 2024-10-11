#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_netcoloc
----------------------------------

Tests for `netcoloc` module.
"""


import sys
import unittest
import os
import ndex2
import networkx as nx
import numpy as np

from netcoloc import netprop


class TestNetcolocNetProp(unittest.TestCase):

    def setUp(self):
        """
        Create directed and undirected toy example networks.
        :return:
        """
        net = ndex2.create_nice_cx_from_file(self._get_5node_network())
        self.g = net.to_networkx(mode='default')
        graph_undirected = nx.Graph()
        graph_undirected.add_nodes_from([1, 2, 3])
        graph_undirected.add_edges_from([(1, 2), (2, 3)])
        self.undirected = graph_undirected
        graph_directed = nx.DiGraph()
        graph_directed.add_nodes_from([1, 2, 3])
        graph_directed.add_weighted_edges_from([(1, 2, 1), (2, 1, 0.5), (2, 3, 0.5), (3, 2, 1)])
        self.directed = graph_directed

    def tearDown(self):
        pass

    def _get_5node_network(self):
        """
        Gets testppi.cx from data directory
        under tests
        :return:
        """
        testdir = os.path.dirname(__file__)
        return os.path.join(testdir, 'data',
                            '5node.cx')

    def _get_testppi_network(self):
        """
        Gets testppi.cx from data directory
        under tests
        :return:
        """
        testdir = os.path.dirname(__file__)
        return os.path.join(testdir, 'data',
                            'testppi.cx')
        
    ### Test get_normalized_adjacency_matrix ###

    # Handle graphs with no zero-degree nodes correctly
    def test_handle_graphs_with_no_zero_degree_nodes(self):
        # Create a graph with no zero-degree nodes
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])

        # Call the function under test
        result = netprop.get_normalized_adjacency_matrix(graph)

        # Check if the result is not empty
        self.assertTrue(result.any())
    
    def test_pass_adjacency_matrix(self):
        # Create a test graph
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3])
        G.add_edges_from([(1, 2), (2, 3)])
        
        # Calculate the normalized adjacency matrix
        result = netprop.get_normalized_adjacency_matrix(nx.to_numpy_array(G), conserve_heat=True)
        
        # Define the expected result
        expected_result = np.array([[0.0, 1, 0.0],
                                    [0.5, 0.0, 0.5],
                                    [0.0, 1, 0.0]])
        np.testing.assert_array_almost_equal(result, expected_result)
        
        
        # Calculate normalized adjacency matrix for a weighted graph with conserved heat
    def test_normalized_adjacency_matrix_weighted_with_conserve_heat(self):
        # Create a test graph
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3])
        G.add_weighted_edges_from([(1, 2, 0.5), (2, 3, 0.8)])

        # Calculate the normalized adjacency matrix
        result = netprop.get_normalized_adjacency_matrix(G, conserve_heat=True, weighted=True)

        # Define the expected result
        expected_result = np.array([[0.0, 0.5, 0.0],
                                    [0.25, 0.0, 0.4],
                                    [0.0, 0.8, 0.0]])

        # Check if the result matches the expected result
        np.testing.assert_array_almost_equal(result, expected_result)
        
    def test_normalized_adjacency_matrix_weighted_without_conserve_heat(self):
        # Create a test graph
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3])
        G.add_weighted_edges_from([(1, 2, 0.5), (2, 3, 0.8)])

        # Calculate the normalized adjacency matrix
        result = netprop.get_normalized_adjacency_matrix(G, conserve_heat=False, weighted=True)

        # Define the expected result
        expected_result = np.array([[0.0, 0.5/np.sqrt(2), 0.0],
                                    [0.5/np.sqrt(2), 0.0, 0.8/np.sqrt(2)],
                                    [0.0, 0.8/np.sqrt(2), 0.0]])

        # Check if the result matches the expected result
        np.testing.assert_array_almost_equal(result, expected_result)
        
    def test_normalization_conservation(self):
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])

        # Call the function under test
        result = netprop.get_normalized_adjacency_matrix(graph, conserve_heat=True)
        # assert row sums are 1
        self.assertTrue(np.allclose(result.sum(axis=1), np.ones(3)))
        
        result2 = netprop.get_normalized_adjacency_matrix(graph, conserve_heat=False)
        # assert row sums are equal to degrees
        self.assertTrue(np.allclose(result2.sum(axis=1), np.array([1, 2, 1])/np.sqrt(2)))

    def test_normalization_warns_weighted_unweighted_input(self):
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])

        # Call the function under test
        with self.assertWarns(UserWarning):
            netprop.get_normalized_adjacency_matrix(graph, weighted=True)

    def test_normalization_warns_degree_zero(self):
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])
        graph.add_node(4)
        self.assertRaises(AssertionError, netprop.get_normalized_adjacency_matrix, graph)
        
    def test_normalization_multigraph(self):
        G = nx.MultiGraph()
        self.assertRaises(ValueError, netprop.get_normalized_adjacency_matrix, G)

    def test_normalization_directed(self):
        G = nx.DiGraph()
        self.assertRaises(ValueError, netprop.get_normalized_adjacency_matrix, G)
        
    ### Test network propagation ###
    
    def test_heats_catch_bad_alpha(self):
        self.assertRaises(AssertionError, netprop.get_individual_heats_matrix, np.array([0]), alpha=10)
        self.assertRaises(AssertionError, netprop.get_individual_heats_matrix, np.array([0]), alpha=-1)


    def test_netprop_undirected(self):
        """
        Test that method returns correct result for undirected graph i.e. symmetrical adjacency matrix
        :return:
        """
        adj = nx.to_numpy_array(self.undirected)
        heats = netprop.get_individual_heats_matrix(adj, alpha=0.5)
        F = netprop.network_propagation(heats, list(self.undirected.nodes), [1])
        expected = [0.75, 0.5, 0.25]
        for i, node in enumerate(F.keys()):
            self.assertAlmostEqual(expected[i], F[node])

    def test_netprop_directed_unidirectional(self):
        """
        Test that model can handle directed networks having nodes with uneven in/out degree
        :return:
        """
        adj = nx.to_numpy_array(self.g)
        heats = netprop.get_individual_heats_matrix(adj, alpha=0.5)
        F1 = netprop.network_propagation(heats, list(self.g.nodes), [185])
        # 1 - test no propagation out of seed with out-degree=0
        expected1 = [0, 0, 0.5, 0, 0]
        for i, node in enumerate(F1.keys()):
            self.assertAlmostEqual(expected1[i], F1[node])
        # 2 - test propagation out of seeds with out-degree > 0
        F2 = netprop.network_propagation(heats, list(self.g.nodes), [180, 190])
        expected2 = [0.25, 0.25, 0.125, 0.25, 0.25]
        for i, node in enumerate(F2.keys()):
            self.assertAlmostEqual(expected2[i], F2[node])

    def test_netprop_directed_bidirectional(self):
        """
        Test that method returns correct result for directed network with all nodes have equal in/out degree,
        but asymmetric adjacency matrix.
        :return:
        """
        adj = nx.to_numpy_array(self.directed)
        heats = netprop.get_individual_heats_matrix(adj, alpha=0.5)
        F = netprop.network_propagation(heats, list(self.directed.nodes), [1])
        expected = [7/12, 2/6, 1/12]
        for i, node in enumerate(F.keys()):
            self.assertAlmostEqual(expected[i], F[node], places=5)

    def test_netprop_unconnected(self):
        """
        Test that method appropriately deals with networks that contain unconnected components
        :return:
        """
        g32 = self.g.copy()
        g32.remove_edge(180, 195, 0)
        g32.remove_edge(190, 175, 0)
        # network now has two components (195, 190) and (185, 180, 175)
        adj = nx.to_numpy_array(g32)
        heats = netprop.get_individual_heats_matrix(adj, alpha=0.5)
        F = netprop.network_propagation(heats, list(g32.nodes), [180])
        # test that nodes with no path to seed get no heat
        self.assertEqual(0, F[195])
        self.assertEqual(0, F[190])
        # test that nodes with path to seed do get heat
        for node in [185, 180, 175]:
            self.assertGreater(F[node], 0)

    def test_netprop_alpha(self):
        """
        Test that the resulting heats change appropriately when alpha is modified from default.
        :return:
        """
        adj = nx.to_numpy_array(self.undirected)
        heats = netprop.get_individual_heats_matrix(adj, alpha=0.1)
        F = netprop.network_propagation(heats, list(self.undirected.nodes), [1])
        expected = [891/980, 9/98, 9/980]
        for i, node in enumerate(F.keys()):
            self.assertAlmostEqual(expected[i], F[node])


if __name__ == '__main__':
    sys.exit(unittest.main())
