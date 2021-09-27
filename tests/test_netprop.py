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
from contextlib import contextmanager
from click.testing import CliRunner

#from netcoloc import netcoloc
from netcoloc import cli
from netcoloc import netprop
from netcoloc import netprop_zscore
from netcoloc import network_localization
from netcoloc import network_colocalization
from netcoloc import netcoloc_utils

class TestNetcoloc(unittest.TestCase):

    def setUp(self):
        pass

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

    def test_000_something(self):
        assert 'one' == 'one'

    def test_normalize_adj_unweighted(self):
        """
        Tests for correct normalization of an unweighted input graph under all parameter options
        :return:
        """
        net = ndex2.create_nice_cx_from_file(self._get_5node_network())
        g = net.to_networkx(mode='default')
        # 1 - non-weighted graph, weighted=True, conserve=False
        adj = netprop.get_normalized_adjacency_matrix(g, conserve_heat=False, weighted=True)
        self.assertEqual([adj[0, 3], adj[3, 0], adj[4, 3], adj[3, 4]], [0.5, 0.5, 1/np.sqrt(6), 1/np.sqrt(6)])
        # 2 - non-weighted graph, weighted=False, conserve = False
        adj = netprop.get_normalized_adjacency_matrix(g, conserve_heat=False, weighted=False)
        self.assertEqual([adj[0, 3], adj[3, 0], adj[4, 3], adj[3, 4]], [0.5, 0.5, 1/np.sqrt(6), 1/np.sqrt(6)])
        # 3 - non-weighted graph, weighted=False, conserve = True
        adj = netprop.get_normalized_adjacency_matrix(g, conserve_heat=True, weighted=False)
        self.assertEqual([adj[0, 3], adj[3, 0], adj[4, 3], adj[3, 4]], [0.5, 0.5, 0.5, 1/3])

    def test_normalize_adj_weighted(self):
        """
        Weights the test graph and tests for correct normalization of an weighted input graph
        under all parameter options
        :return:
        """
        net = ndex2.create_nice_cx_from_file(self._get_5node_network())
        g = net.to_networkx(mode='default')
        weights = {e: [1, 5, 0.1, 0.9, 2][i] for i, e in enumerate(g.edges)}
        g_weighted = g.copy()
        nx.set_edge_attributes(g_weighted, values=weights, name="weight")
        self.assertTrue(nx.is_weighted(g_weighted))
        # 1 - weighted input graph, weighted=False, conserve=True
        adj = netprop.get_normalized_adjacency_matrix(g_weighted, conserve_heat=True, weighted=False)
        self.assertEqual([adj[0, 3], adj[3, 0], adj[1, 4], adj[4, 1]], [0.5, 0.5, 1/3, 1/2])
        # 2 - weighted input graph, weighted=True, conserve= True
        adj = netprop.get_normalized_adjacency_matrix(g_weighted, conserve_heat=True, weighted=True)
        self.assertEqual([adj[0, 1], adj[1, 0], adj[3, 4], adj[4, 3]], [5/2, 5/2, 0.9/3, 0.9/2])
        # 3 - weighted input graph, weighted=True, conserve=False
        adj = netprop.get_normalized_adjacency_matrix(g_weighted, conserve_heat=False, weighted=True)
        self.assertEqual([adj[0, 3], adj[3, 0], adj[4, 2], adj[2, 4]], [0.1/2, 0.1/2, 2/np.sqrt(3), 2/np.sqrt(3)])


    def test_individual_heats_directed(self):
        #TODO
        pass

    def test_individual_heats_undirected(self):
        #TODO
        pass

    def test_network_propagation_directed(self):
        #TODO
        pass

    def test_network_propagation_undirected(self):
        #TODO
        pass
    #def test_command_line_interface(self):
     #   runner = CliRunner()
      #  result = runner.invoke(cli.main)
       # assert result.exit_code == 0
        #assert 'netcoloc.cli.main' in result.output
        #help_result = runner.invoke(cli.main, ['--help'])
        #assert help_result.exit_code == 0
        #assert '--help  Show this message and exit.' in help_result.output


if __name__ == '__main__':
    sys.exit(unittest.main())
