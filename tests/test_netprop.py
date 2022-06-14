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

    def test_normalization_multigraph(self):
        #TODO
        pass

    def test_normalization_undirected(self):
        #TODO
        pass

    def test_normalization_directed(self):
        #TODO
        pass

    def test_normalization_conservation(self):
        #TODO
        pass

    def test_normalization_weighting(self):
        #TODO
        pass

    def test_normalization_warns_weighted_unweighted_input(self):
        #TODO
        pass

    def test_normalization_warns_degree_zero(self):
        g41 = self.g.copy()
        g41.remove_edge(180, 195, 0)
        g41.remove_edge(190, 195, 0)
        self.assertRaises(AssertionError, netprop.get_normalized_adjacency_matrix, g41)

    def test_heats_symmetric(self):
        #TODO
        pass

    def test_heats_asymmetric(self):
        #TODO
        pass

    def test_heats_catch_bad_alpha(self):
        self.assertRaises(AssertionError, netprop.get_individual_heats_matrix, np.array([0]), alpha=10)
        self.assertRaises(AssertionError, netprop.get_individual_heats_matrix, np.array([0]), alpha=-1)


if __name__ == '__main__':
    sys.exit(unittest.main())
