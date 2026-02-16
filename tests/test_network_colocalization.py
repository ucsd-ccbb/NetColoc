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
import unittest.mock
import ndex2
import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt

from netcoloc.network_colocalization import *


class TestNetworkColocalization(unittest.TestCase):

    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.data_path = os.path.join(self.dir_path, 'data')
        self.mock_interactome_file = os.path.join(self.data_path, 'mock_interactome.gpickle')
        with open(self.mock_interactome_file, 'rb') as f:
            self.mock_interactome = pickle.load(f)
        self.zscores1 = pd.read_csv(os.path.join(self.data_path, 'zscores1.txt'), sep='\t', index_col=0)
        self.zscores2 = pd.read_csv(os.path.join(self.data_path, 'zscores2.txt'), sep='\t', index_col=0)

    def tearDown(self):
        plt.close()

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

        res = transform_edges(g)

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

        res = transform_edges(g)
        self.assertEqual(7862, len(list(res.edges())))
        self.assertEqual(773, len(list(res.nodes())))
        
    def test_calculate_network_overlap_series(self):
        z1 = self.zscores1.z
        z2 = self.zscores2.z
        zz = calculate_network_overlap(z1, z2)
        self.assertListEqual(zz, [22, 32, 84, 41, 6, 30, 75, 47, 76, 27, 28])
    
    def test_calculate_network_overlap_df(self):
        zz = calculate_network_overlap(self.zscores1, self.zscores2)
        self.assertListEqual(zz, [22, 32, 84, 41, 6, 30, 75, 47, 76, 27, 28])
    
    def test_calculate_network_overlap_array(self):
        z1 = self.zscores1.z.values
        z2 = self.zscores2.z.values
        zz = calculate_network_overlap(z1, z2)
        self.assertEqual(len(zz), 11)
        
    def test_calculate_network_overlap_threshold(self):
        z1 = self.zscores1.z
        z2 = self.zscores2.z
        zz = calculate_network_overlap(z1, z2, z_score_threshold=10)
        self.assertListEqual(zz, [])
        
    def test_calculate_network_overlap_unequal_length(self):
        z1 = self.zscores1.z
        z2 = self.zscores2.z[:10]
        with self.assertRaises(AssertionError):
            calculate_network_overlap(z1, z2)
            
    def test_calculate_network_overlap_raises_warnings(self):
        z1 = self.zscores1.z
        z2 = self.zscores2.z
        # create mock warning function
        with unittest.mock.patch('warnings.warn') as mock_warn:
            calculate_network_overlap(z1, z2, z_score_threshold=0, z1_threshold=2, z2_threshold=2)
            mock_warn.assert_called()

        with unittest.mock.patch('warnings.warn') as mock_warn:
            calculate_network_overlap(z1, z2, z_score_threshold=3, z1_threshold=1, z2_threshold=2)
            mock_warn.assert_called()


    def test_calculate_network_overlap_subgraph(self):
        z1 = self.zscores1.z
        z2 = self.zscores2.z
        Gsub = calculate_network_overlap_subgraph(self.mock_interactome, z1, z2)
        self.assertIsInstance(Gsub, nx.Graph)
        self.assertEqual(len(Gsub.nodes()), 11)
        self.assertEqual(len(Gsub.edges()), 20)
        Gsub2 = calculate_network_overlap_subgraph(self.mock_interactome, z1, z2, z_score_threshold=10)
        self.assertIsInstance(Gsub2, nx.Graph)
        self.assertEqual(len(Gsub2.nodes()), 0)
        
    def test_calculate_network_overlap_subgraph_raises_node_error(self):
        z1 = self.zscores1.z
        z1.index = [x+100 for x in range(100)]
        z2 = self.zscores2.z
        z2.index = [x+100 for x in range(100)]
        with self.assertRaises(AssertionError):
            calculate_network_overlap_subgraph(self.mock_interactome, z1, z2)
        
    def test_calculate_expected_overlap_no_overlap_control(self):
        overlap_size, null_sizes = calculate_expected_overlap(self.zscores1, self.zscores2, num_reps=100, overlap_control=None, 
                                                              seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertEqual(overlap_size, 11)
        self.assertEqual(len(null_sizes), 100)
    
    def test_calculate_expected_overlap_bin_overlapping_genes(self):
        overlap_size, null_sizes = calculate_expected_overlap(self.zscores1, self.zscores2, num_reps=100, overlap_control='bin', 
                                                              seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertEqual(overlap_size, 11)
        self.assertEqual(len(null_sizes), 100)
    
    def test_calculate_expected_overlap_remove_overlapping_genes(self):
        overlap_size, null_sizes = calculate_expected_overlap(self.zscores1, self.zscores2, num_reps=100, overlap_control='remove', 
                                                              seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertEqual(overlap_size, 9) # two of the three overlapping genes meet the thresholds
        self.assertEqual(len(null_sizes), 100)
        
    def test_calculate_expected_overlap_plot(self):
        with unittest.mock.patch('matplotlib.pyplot.legend') as mock_legend:
            _, _ = calculate_expected_overlap(self.zscores1, self.zscores2, num_reps=100, plot=True)
            mock_legend.assert_called()
    
    def test_calculate_mean_z_score_distribution_no_control(self):
        mu, random_mu = calculate_mean_z_score_distribution(self.zscores1, self.zscores2, num_reps=100, overlap_control=None,
                                                             zero_double_negatives=False,
                                                            seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertAlmostEqual(mu, 0.85409, places=3)
        self.assertEqual(len(random_mu), 100)
        
    def test_calculate_mean_z_score_distribution_bin_overlapping(self):
        mu, random_mu = calculate_mean_z_score_distribution(self.zscores1, self.zscores2, num_reps=100, overlap_control='bin',
                                                            zero_double_negatives=False,
                                                            seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertAlmostEqual(mu, 0.85409, places=3)
        self.assertEqual(len(random_mu), 100)
        
    def test_calculate_mean_z_score_distribution_remove_overlapping(self):
        mu, random_mu = calculate_mean_z_score_distribution(self.zscores1, self.zscores2, num_reps=100, overlap_control='remove',
                                                            zero_double_negatives=False,
                                                            seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertAlmostEqual(mu, 0.78398, places=3)
        self.assertEqual(len(random_mu), 100)
        
    def test_calculate_mean_z_score_distribution_zero_double_negatives(self):
        mu, random_mu = calculate_mean_z_score_distribution(self.zscores1, self.zscores2, num_reps=100, zero_double_negatives=True, overlap_control=None,
                                                            seed1=[22, 7, 84, 6, 35, 90], seed2=[22, 7, 84, 71, 3, 89])
        self.assertAlmostEqual(mu, 0.82479, places=3)
        self.assertEqual(len(random_mu), 100)

    def test_get_p_from_permutation_results(self):
        obs = 0.97
        permuted = [np.random.random() for _ in range(100)]
        p_value = get_p_from_permutation_results(obs, permuted, alternative='greater')
        self.assertIsInstance(p_value, float)
        p_value2 = get_p_from_permutation_results(obs, permuted, alternative='two-sided')
        self.assertAlmostEqual(p_value2, 2*p_value, places=3)
        

    def test_get_p_from_permutation_results_prints_warning(self):
        obs = 100
        permuted = [np.nan for _ in range(100)]
        with unittest.mock.patch('warnings.warn') as mock_warn:
            _ = get_p_from_permutation_results(obs, permuted)
            mock_warn.assert_called()
        
    def test_get_p_from_permutation_results_invalid_alternative(self):
        obs = 0.97
        permuted = [np.random.random() for _ in range(100)]
        with self.assertRaises(AssertionError):
            _ = get_p_from_permutation_results(obs, permuted, alternative='invalid')



if __name__ == '__main__':
    sys.exit(unittest.main())
