#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_netcoloc
----------------------------------

Tests for `netcoloc` module.
"""


import os
import sys
import unittest
import json
import ndex2
import networkx as nx
import pandas as pd

from netcoloc import network_colocalization


class TestNetworkColocalization(unittest.TestCase):

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

    def test_transform_edges_5node_network(self):
        net = ndex2.create_nice_cx_from_file(self._get_5node_network())

        g = net.to_networkx(mode='default')

        name_dict = {}
        for entry in g.nodes.data():
            name_dict[entry[1]['name']] = entry[0]

        res = network_colocalization.transform_edges(g)

        self.assertEqual(2, len(list(res.edges())))
        self.assertEqual(5, len(list(res.nodes())))
        self.assertEqual(1.125, res.get_edge_data(name_dict['B'],
                                                  name_dict['D'])['weight'])
        self.assertEqual(1.0, res.get_edge_data(name_dict['A'],
                                                name_dict['E'])['weight'])

    @unittest.skip('Tests overlap network from notebook. Too slow...')
    def test_transform_edges_testppi_network(self):
        net = ndex2.create_nice_cx_from_file(self._get_testppi_network())

        g = net.to_networkx(mode='default')

        name_dict = {}
        for entry in g.nodes.data():
            name_dict[entry[1]['name']] = entry[0]

        res = network_colocalization.transform_edges(g)
        self.assertEqual(7862, len(list(res.edges())))
        self.assertEqual(773, len(list(res.nodes())))


if __name__ == '__main__':
    sys.exit(unittest.main())
