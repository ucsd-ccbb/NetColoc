
"""
test_netcoloc.validation
----------------------------------

Tests for `validation` module.
"""


import os
import sys
import unittest
import json
import unittest.mock
import ndex2
import networkx as nx
import pandas as pd
import requests
from unittest import mock
from netcoloc.validation import *


class TestNetworkColocalization(unittest.TestCase):

    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.data_path = os.path.join(self.dir_path, 'data')
        # make a file called 'MPheno_OBO.ontology' in the current directory
        mpo_file = os.path.join(self.data_path, 'MPheno_OBO.ontology.exists')
        with open(mpo_file, 'w') as f:
            f.write('')
        mgi_file = os.path.join(self.data_path, 'MGI_PhenoGenoMP.rpt.exists')
        with open(mgi_file, 'w') as f:
            f.write('MGI:123456\tTP53\tinvolves\tMP:000001\tPMID:123456\tMGI:123456\n'
                    'MGI:789012\tSCN2A\tinvolves\tMP:000002\tPMID:789012\tMGI:789012\n'
                   )
        self.mock_mgi = os.path.join(self.data_path , 'mock_mgi.txt')

    def tearDown(self):
        # remove the file created in setUp
        mpo_file = os.path.join(self.data_path, 'MPheno_OBO.ontology.exists')
        if os.path.exists(mpo_file):
            os.remove(mpo_file)
        mgi_file = os.path.join(self.data_path, 'MGI_PhenoGenoMP.rpt.exists')
        if os.path.exists(mgi_file):
            os.remove(mgi_file)

    def test_check_MGI_download_url(self):
        url = 'http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt'
        response = requests.head(url, allow_redirects=True)
        self.assertEqual(response.status_code, 200)
        
    def test_MGI_already_exists(self):
        with unittest.mock.patch('requests.get') as mock_get:
            _ = load_MGI_mouseKO_data(url = 'http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt.exists',
                                      update=False, data_loc=self.data_path)
            mock_get.assert_not_called()
    
    def test_update_MGI(self):
        mock_response = unittest.mock.Mock()
        mock_response.content = b'MGI:123456\tTP53\tinvolves\t'
        with unittest.mock.patch('requests.get', return_value=mock_response) as mock_get:
            _ = load_MGI_mouseKO_data(url = 'http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt.exists',
                                      update=True, data_loc=self.data_path)
            mock_get.assert_called_once_with('http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt.exists', 
                                             allow_redirects=True)
    
    def test_load_MGI_mouseKO_data(self):
        mgi = load_MGI_mouseKO_data(url = 'mock_mgi.txt',
                                     update=False, data_loc=self.data_path)
        self.assertIn('human_ortholog', mgi.columns)
        
    def test_check_MPO_download_url(self):
        url = 'http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology'
        response = requests.head(url, allow_redirects=True)
        self.assertEqual(response.status_code, 200)
        
    def test_MPO_already_exists(self):
        mpo_file = 'MPheno_OBO.ontology.exists'
        # mock the obo.read_obo function
        with unittest.mock.patch('obonet.read_obo') as mock_read_obo:
            _ = load_MPO(url = 'http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology.exists',
                           update=False, data_loc=self.data_path)
            mock_read_obo.assert_called_once_with(os.path.join(self.data_path, mpo_file))

    def test_update_MPO(self):
        mpo_file = 'MPheno_OBO.ontology.exists'
        # mock the obo.read_obo function
        with unittest.mock.patch('obonet.read_obo') as mock_read_obo:
            _ = load_MPO(url = 'http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology.exists',
                           update=True, data_loc=self.data_path)
            mock_read_obo.assert_called_once_with('http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology.exists')


if __name__ == '__main__':
    sys.exit(unittest.main())
