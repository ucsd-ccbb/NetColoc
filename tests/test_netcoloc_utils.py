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
from netcoloc import netcoloc_utils


class TestNetcolocUtil(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_degree_binning(self):
        pass

    # Bins nodes correctly based on their degree when all nodes are included
    def test_bins_nodes_correctly_with_all_nodes(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 2
        expected_bins = [['A', 'B', 'C'], ['D', 'E', 'F']]
        expected_degree_to_bin_index = {1: 0, 2: 0, 3: 1}
    
        bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)
    
        self.assertEqual(bins, expected_bins)
        self.assertEqual(degree_to_bin_index, expected_degree_to_bin_index)
        
    # Handles empty `node_to_degree_dict` gracefully
    def test_handles_empty_node_to_degree_dict(self):
        node_to_degree_dict = {}
        min_bin_size = 2
        
        with self.assertRaises(AssertionError):
            netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)
    
        min_bin_size = 0
        bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)
        self.assertEqual(bins, [])
        self.assertEqual(degree_to_bin_index, {})
    
    # Handles cases where `lengths` parameter is provided, binning only specified nodes
    def test_bins_nodes_correctly_with_specified_nodes(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 2
        lengths = ['A', 'B', 'D', 'E']
        expected_bins = [['A', 'B'], ['D', 'E']]
        expected_degree_to_bin_index = {1: 0, 2: 0, 3: 1}

        bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size, lengths)

        self.assertEqual(bins, expected_bins)
        self.assertEqual(degree_to_bin_index, expected_degree_to_bin_index)
        
    # Combines the last bin with the second last bin when the last bin is too small
    def test_combine_bins_when_last_bin_small(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 4
        expected_bins = [['A', 'B', 'C', 'D', 'E', 'F']]
        expected_degree_to_bin_index = {1: 0, 2: 0, 3: 0}

        bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)

        self.assertEqual(bins, expected_bins)
        self.assertEqual(degree_to_bin_index, expected_degree_to_bin_index)
    
        # Asserts an error is raised when min_bin_size is larger than the total number of nodes
    def test_assert_error_min_bin_size_too_high(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 10

        with self.assertRaises(AssertionError):
            netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)
    
        # Test that the function raises a warning when there are missing values in lengths
    def test_missing_values_in_lengths(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 2
        lengths = ['G', 'H']  # Values not present in node_to_degree_dict

        with self.assertWarns(UserWarning):
            bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size, lengths)

        self.assertEqual(bins, [])
        self.assertEqual(degree_to_bin_index, {})
        
    # Test that a warning is raised when lengths parameter is empty
    def test_warning_raised_empty_lengths_parameter(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 2
        lengths = []

        with self.assertWarns(UserWarning):
            bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size, lengths)

        self.assertEqual(bins, [])
        self.assertEqual(degree_to_bin_index, {})
    
    # test that min_bin_size equal to the number of nodes is allowed
    def test_min_bin_equal_same_as_number_of_nodes(self):
        node_to_degree_dict = {'A': 1, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3}
        min_bin_size = 6
    
        expected_bins = [['A', 'B', 'C', 'D', 'E', 'F']]
        expected_degree_to_bin_index = {1: 0, 2: 0, 3: 0}

        bins, degree_to_bin_index = netcoloc_utils.get_degree_binning(node_to_degree_dict, min_bin_size)

        self.assertEqual(bins, expected_bins)
        self.assertEqual(degree_to_bin_index, expected_degree_to_bin_index)



if __name__ == '__main__':
    sys.exit(unittest.main())
