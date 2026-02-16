
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
        self.hier_df = pd.DataFrame({'CD_MemberList':['1 2 3', '2 3', '4 5 6 7 1 2 3', '7 3 4 3 0 3 2 1', '6 9 7'], 'name':['C1', 'C2', 'C3', 'C4', 'C5']}, index=['C1', 'C2', 'C3', 'C4', 'C5'])
        self.hier_df['CD_MemberList_Size'] = self.hier_df['CD_MemberList'].apply(lambda x: len(x.split(' ')))
        self.mock_mpo_file = os.path.join(self.data_path, 'mock_mpo.obo')
        self.mock_assoc = pd.read_csv(os.path.join(self.data_path, 'propagated_associations.tsv'), sep='\t')
        self.mock_assoc.human_ortholog = self.mock_assoc.human_ortholog.astype(str)
    
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
        # mock the obo.read_obo function
        with unittest.mock.patch('obonet.read_obo') as mock_read_obo:
            _ = load_MPO(url = 'http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology.exists',
                           update=True, data_loc=self.data_path)
            mock_read_obo.assert_called_once_with('http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology.exists')
            
    def test_check_MRK_url(self):
        url = 'http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt'
        response = requests.head(url, allow_redirects=True)
        self.assertEqual(response.status_code, 200)
        
    def test_check_HMD_url(self):
        url = 'http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt'
        response = requests.head(url, allow_redirects=True)
        self.assertEqual(response.status_code, 200)
        
    def test_load_MPO(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        self.assertIsInstance(mpo, nx.MultiDiGraph)
        self.assertEqual(len(mpo.nodes), 16)
        
    def test_hypergeometric_test(self):
        mgi_genes = ['1', '2', '3', '4']
        G_int = nx.Graph()
        G_int.add_nodes_from([str(x) for x in range(0, 10)])
        p1, x1, l1 = perform_hypergeometric_test(self.hier_df, 'C1', G_int, mgi_genes)
        self.assertAlmostEqual(p1, 0.033, places=2)
        self.assertIsInstance(x1, int)
        self.assertIsInstance(l1, str)
        p2, x2, l2 = perform_hypergeometric_test(self.hier_df, 'C4', G_int, mgi_genes)
        self.assertAlmostEqual(p2, 0.33, places=2)
        self.assertEqual(x2, 4)
        self.assertEqual(l2, '1 2 3 4')
        #test with no overlap
        p3, x3, l3 = perform_hypergeometric_test(self.hier_df, 'C1', G_int, ['4', '5', '6'])
        self.assertTrue(np.isnan(p3))
        self.assertIsInstance(x3, int)
        self.assertIsInstance(l3, str)
        
    def test_log_odds_ratio(self):
        mgi_genes = ['1', '2', '3', '4']
        G_int = nx.Graph()
        G_int.add_nodes_from([str(x) for x in range(0, 10)])
        p1, x1, l1 = perform_log_odds_z_test(mgi_genes, self.hier_df.at['C1', 'CD_MemberList'].split(' '), list(G_int.nodes()), verbose=False, name='')
        self.assertAlmostEqual(p1, 0.05543, places=2)
        self.assertEqual(len(x1), 2)
        self.assertIsInstance(l1, float)
        p2, _, _ = perform_log_odds_z_test([], self.hier_df.at['C4', 'CD_MemberList'].split(' '), list(G_int.nodes()), verbose=False, name='')
        self.assertTrue(np.isnan(p2))
        p3, _, _ = perform_log_odds_z_test(mgi_genes, [], list(G_int.nodes()), verbose=False, name='')
        self.assertTrue(np.isnan(p3))
        
    def test_focus_ontology(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        test_term = 'MP:0000007'
        focal_terms = focus_ontology(mpo, test_term, include_children=True, include_parents=False)
        self.assertIn('MP:0000007', focal_terms)
        self.assertNotIn('MP:0000001', focal_terms)  # MP:0000001 is a parent of MP:0000007
        focal_terms = focus_ontology(mpo, test_term, include_children=False, include_parents=True)
        self.assertIn('MP:0000007', focal_terms)
        self.assertIn('MP:0000001', focal_terms)

    def test_get_focal_terms(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        test_term = 'MP:0000007'
        terms, term_name = get_focal_terms(mpo, test_term, use_ddot=False, include_children=True, include_parents=False)
        self.assertIn('MP:0000007', terms)
        self.assertNotIn('MP:0000001', terms)
        self.assertEqual(term_name, 'abnormal learning/memory')
        terms, term_name = get_focal_terms(mpo, 'MP:0000016',include_children=True, include_parents=False )
        self.assertIn('MP:0000016', terms)
        self.assertEqual(len(terms), 1)
                
    def test_test_single_MPO_term(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        focal_terms, name = get_focal_terms(mpo, 'MP:0000007', use_ddot=False, include_children=True, include_parents=False)
        G_int = nx.Graph()
        G_int.add_nodes_from([str(x) for x in range(0, 10)])
        p, ci, odds, term_genes = test_single_MPO_term(focal_terms, self.hier_df, 'C1', self.mock_assoc, list(G_int.nodes()), verbose=False,
                         min_genes=1, max_genes=2000, name=name)
        self.assertIsInstance(p, float)
        self.assertEqual(len(ci), 2)
        self.assertIsInstance(odds, float)
        self.assertEqual(len(term_genes), 5)

    def test_MPO_enrichment_root(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        focal_terms, name = get_focal_terms(mpo, 'MP:0000007', use_ddot=False, include_children=True, include_parents=False)
        G_int = nx.Graph()
        G_int.add_nodes_from([str(x) for x in range(0, 10)])
        
        res_df = MPO_enrichment_root(self.hier_df,mpo,self.mock_assoc,focal_terms,G_int,verbose=False, use_ddot=False,
                        min_genes=1, max_genes=2000)
        self.assertIsInstance(res_df, pd.DataFrame)
        self.assertListEqual(res_df.columns.tolist(), ['OR_p', 'log_OR', 'log_OR_CI_lower', 'log_OR_CI_upper','num_genes_in_term', 'MP_description'])
        self.assertEqual(len(res_df), 4)
        res_df2 = MPO_enrichment_root(self.hier_df,mpo,self.mock_assoc,focal_terms,G_int,verbose=False, use_ddot=False,
                        min_genes=2, max_genes=2000)
        self.assertEqual(len(res_df2), 3)
        
    def test_MPO_enrichment_full(self):
        mpo = load_MPO(url=self.mock_mpo_file, use_ddot=False, update=False, data_loc=self.data_path)
        focal_terms, name = get_focal_terms(mpo, 'MP:0000007', use_ddot=False, include_children=True, include_parents=False)
        G_int = nx.Graph()
        G_int.add_nodes_from([str(x) for x in range(0, 10)])
        
        res_df = MPO_enrichment_full(self.hier_df,mpo,self.mock_assoc,focal_terms,G_int,verbose=False, use_ddot=False,
                        min_genes=1, max_genes=2000)
        self.assertIsInstance(res_df, pd.DataFrame)
        self.assertEqual(res_df.shape, (5, 20))
        res_df2 = MPO_enrichment_full(self.hier_df,mpo,self.mock_assoc,focal_terms,G_int,verbose=False, use_ddot=False,
                        min_genes=2, max_genes=2000)
        self.assertEqual(res_df2.shape, (5, 20))

if __name__ == '__main__':
    sys.exit(unittest.main())
