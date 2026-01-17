#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_netcoloc
----------------------------------

Tests for `netcoloc` module.
"""


import sys
import unittest
from contextlib import contextmanager
from click.testing import CliRunner
import pickle
import os
from unittest import mock
import pickle
import os
from unittest import mock

#from netcoloc import netcoloc
from netcoloc.netprop_zscore import *

class TestNetcolocZscore(unittest.TestCase):

    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.data_path = os.path.join(self.dir_path, 'data')
        # load mock_interactome from pickle
        self.mock_interactome_file = os.path.join(self.data_path, 'mock_interactome.gpickle')
        with open(self.mock_interactome_file, 'rb') as f:
            self.mock_interactome = pickle.load(f)
            
        self.mock_seed_file = os.path.join(self.data_path, 'mock_seed_file.txt')
        with open(self.mock_seed_file, 'r') as f:
            self.mock_seeds = [int(x) for x in f.read().splitlines()]
            
        self.mock_degree = dict(self.mock_interactome.degree())
        self.mock_nodes = list(self.mock_interactome.nodes())
        self.mock_heat_file = os.path.join(self.data_path, 'mock_heat_matrix.pkl')
        with open(self.mock_heat_file, 'rb') as f:
            self.heat_matrix = pickle.load(f)
            

    def tearDown(self):
        file_list = []
        for f in file_list:
            if os.path.exists(os.path.join(self.data_path, f)):
                os.remove(os.path.join(self.data_path, f))

    # Successfully calculates z-scores when given valid seed genes and interactome
    def test_load_interactome_calls(self):
        # Mocking the networkx and ndex2 functions
        with mock.patch('pickle.load', return_value=nx.Graph()) as mock_read_gpickle:
            with mock.patch('ndex2.create_nice_cx_from_server', return_value=mock.Mock(to_networkx=lambda: nx.Graph())) as mock_create_nice_cx:
                # Test with interactome_file
                z_scores, random_final_heats = netprop_zscore(
                    seed_gene_file=self.mock_seed_file,
                    interactome_file=self.mock_interactome_file
                )
                self.assertIsInstance(z_scores, pd.Series)
                self.assertIsInstance(random_final_heats, np.ndarray)

                # Test with interactome_uuid
                z_scores, random_final_heats = netprop_zscore(
                    seed_gene_file=self.mock_seed_file,
                    interactome_uuid='f93f402c-86d4-11e7-a10d-0ac135e8bacf'
                )
                self.assertIsInstance(z_scores, pd.Series)
                self.assertIsInstance(random_final_heats, np.ndarray)

        # Ensure the mocks were called
        mock_read_gpickle.assert_called()
        mock_create_nice_cx.assert_called()
        
    def test_netprop_zscores_save_z_scores(self):
        # mock to_csv
        with mock.patch('pandas.Series.to_csv') as mock_to_csv:
            z_scores, random_final_heats = netprop_zscore(
                seed_gene_file=self.mock_seed_file,
                interactome_file=self.mock_interactome_file,
                save_z_scores=True, verbose=False
            )
            # Check if to_csv was called
            mock_to_csv.assert_called_once()
            # Check if the output is a pandas Series
            self.assertIsInstance(z_scores, pd.Series)
    
    def test_netprop_zscores_save_final_heat(self):
        # mock to_csv
        with mock.patch('pandas.DataFrame.to_csv') as mock_to_csv:
            z_scores, random_final_heats = netprop_zscore(
                seed_gene_file=self.mock_seed_file,
                interactome_file=self.mock_interactome_file,
                save_final_heat=True,verbose=False)
            # Check if to_csv was called
            mock_to_csv.assert_called_once()
    
    def test_netprop_zscores_save_random_final_heats(self):
        # mock to_csv
        with mock.patch('pandas.DataFrame.to_csv') as mock_to_csv:
            z_scores, random_final_heats = netprop_zscore(
                seed_gene_file=self.mock_seed_file,
                interactome_file=self.mock_interactome_file,
                save_random_final_heats=True, verbose=False
            )
            # Check if to_csv was called
            mock_to_csv.assert_called_once()
            # Check if the output is a numpy array
            self.assertIsInstance(random_final_heats, np.ndarray)
    
    def test_netprop_zscore_success(self):
        num_reps=10
        z_scores, random_final_heats = netprop_zscore(
            seed_gene_file=self.mock_seed_file,
            interactome_file=self.mock_interactome_file,
            verbose=False, num_reps=num_reps
        )
        self.assertIsInstance(z_scores, pd.Series)
        self.assertIsInstance(random_final_heats, np.ndarray)
        self.assertEqual(len(z_scores), len(self.mock_nodes))
        self.assertEqual(random_final_heats.shape[0], num_reps)
        self.assertEqual(random_final_heats.shape[1], len(self.mock_nodes))
        
    # Raises TypeError if neither interactome_file nor interactome_uuid is provided
    def test_netprop_zscore_missing_interactome_source_raises_typeerror(self):
        with self.assertRaises(TypeError):
            netprop_zscore(seed_gene_file='dummy_seed_file.txt', interactome_uuid=None)
            
    def test_calculate_heat_zscores_raises_no_seed_gene_assertion_error(self):
        with self.assertRaises(AssertionError):
            calculate_heat_zscores(
                individual_heats_matrix=self.heat_matrix,
                nodes=self.mock_nodes,
                degrees=self.mock_degree,
                seed_genes=[523]
            )
            
    def test_calculate_heat_zscores_raises_minimum_bin_size_assertion_error(self):
        with self.assertRaises(AssertionError):
            calculate_heat_zscores(
                individual_heats_matrix=self.heat_matrix,
                nodes=self.mock_nodes,
                degrees=self.mock_degree,
                seed_genes=self.mock_seeds,
                minimum_bin_size=1000  # Larger than the number of nodes
            )

    def test_calculate_heat_zscores_success(self):
        num_reps=10
        z_scores, final_heat, random_final_heats = calculate_heat_zscores(
            individual_heats_matrix=self.heat_matrix,
            nodes=self.mock_nodes,
            degrees=self.mock_degree,
            seed_genes=self.mock_seeds,
            num_reps=num_reps
        )
        self.assertIsInstance(z_scores, pd.Series)
        self.assertIsInstance(random_final_heats, np.ndarray)
        self.assertEqual(len(z_scores), len(self.mock_nodes))
        self.assertEqual(len(final_heat), len(self.mock_nodes))
        self.assertEqual(random_final_heats.shape[0], num_reps)
        self.assertEqual(pd.DataFrame(random_final_heats).isna().sum(axis=1).max(), len(self.mock_seeds))

    def test_netprop_zscore(self):
        pass
    
    



if __name__ == '__main__':
    sys.exit(unittest.main())