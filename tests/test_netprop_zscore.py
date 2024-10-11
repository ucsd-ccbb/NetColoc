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

    def tearDown(self):
        file_list = []
        for f in file_list:
            if os.path.exists(os.path.join(self.data_path, f)):
                os.remove(os.path.join(self.data_path, f))

    def test_000_something(self):
        pass

    # Successfully calculates z-scores when given valid seed genes and interactome
    def test_load_interactome_success(self):
        # Mocking the networkx and ndex2 functions
        with mock.patch('networkx.read_gpickle', return_value=nx.Graph()) as mock_read_gpickle:
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
        
        # Raises TypeError if neither interactome_file nor interactome_uuid is provided
    def test_missing_interactome_source_raises_typeerror(self):
        with self.assertRaises(TypeError):
            netprop_zscore(seed_gene_file='dummy_seed_file.txt')
            


if __name__ == '__main__':
    sys.exit(unittest.main())
