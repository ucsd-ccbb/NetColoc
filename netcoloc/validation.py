# -*- coding: utf-8 -*-

'''Functions for performing validation of NetColoc subgraph
'''

import warnings
import numpy as np
import pandas as pd
import pickle
import obonet as obo
import os

import requests
from scipy.stats import hypergeom
from statsmodels.stats import contingency_tables
import networkx as nx
from collections import defaultdict

try:
    # if forcing use of ddot
    import ddot
    from ddot import Ontology
    DDOT_LOADED = True
    warnings.warn('Use of ddot will be deprecated in netcoloc>=1.0.1. By default, ddot will no longer be used, to force use of ddot set use_ddot=True for all validation functions.', DeprecationWarning)
except ImportError as ie:
    DDOT_LOADED = False

# find human orthologs of mouse genes
import mygene
mg = mygene.MyGeneInfo()


def focus_ontology(ont, start_node, include_children=True, include_parents=False):
    """
    Focus the ontology on a specific node and its neighbors.
    :param ont: The ontology to focus on
    :type ont: :py:class:`networkx.MultiDiGraph`
    :param start_node: The node to focus on
    :type start_node: str
    :param include_children: Whether to include children of the start node
    :type include_children: bool
    :param include_parents: Whether to include parents of the start node. Note that this will include the root node.
    :type include_parents: bool
    :return: List of nodes in the focused ontology
    :rtype: list
    """
    successors = nx.dfs_successors(ont, start_node)
    predecessors = nx.dfs_successors(nx.reverse(ont, copy=True), start_node)
    focus_list = set([start_node])
    if include_parents:
        for k,v in successors.items():
            focus_list = focus_list.union(set(v))
    if include_children:
        for k,v in predecessors.items():
            focus_list = focus_list.union(set(v))
    return list(focus_list)

def load_MGI_mouseKO_data(url='http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt',
                          update=False, data_loc=None, map_using='mygeneinfo', verbose=False):
    """
    Function to parse and load mouse knockout data from MGI.

    :param url: location of MGI knockout data
    :type url: str
    :param data_loc: location to save the downloaded file, if None, saves in current directory
    :type data_loc: str
    :param update: whether to update the data if it already exists
    :type update: bool
    :return: parsed MGI knockout dataframe, including column for human orthologs
    :rtype: :py:class:`pandas.DataFrame`
    
    """
    if data_loc is not None:
        rpt_file_target = os.path.join(data_loc, url.split('/')[-1])
    else:
        rpt_file_target = url.split('/')[-1]
    
    if update or (not os.path.exists(rpt_file_target)):
        r = requests.get(url,allow_redirects=True)
        with open(rpt_file_target, 'wb') as f:
            f.write(r.content)
    
    mgi_df = pd.read_csv(rpt_file_target, sep='\t',
                        names=['MGI_Allele_Accession_ID',
                               'Allele symbol', 'involves',
                               'MP', 'PMID', 'MGI_marker_accession_ID'])

    mapping = map_mgi_to_human_orthologs(mgi_df, map_using=map_using, verbose=verbose, data_loc=data_loc, update=update)

    mgi_df['human_ortholog']=mgi_df['gene_name'].map(mapping)
    return mgi_df

def map_mgi_to_human_orthologs(mgi_df, map_using='mygeneinfo', verbose=False, data_loc=None, update=False):
    # TODO needs testing for formatting.
    assert map_using in ['mygeneinfo', 'mgi'], 'map_using must be one of mygeneinfo, mgi'
    if map_using == 'mygeneinfo':
        # extract gene names
        gene_name = [a.split('<')[0] for a in mgi_df['Allele symbol'].tolist()]
        mgi_df['gene_name']=gene_name
        mgi_df.index=mgi_df['gene_name']

        # map mouse genes to human orthologs
        mouse_genes = list(np.unique(mgi_df['gene_name']))
        mg_mapped = mg.querymany(mouse_genes,
                                as_dataframe=True, species=['mouse','human'],
                                scopes='symbol', fields='symbol')
        # map mouse genes to human orthologs
        mouse_genes = list(np.unique(mgi_df['gene_name']))
        mg_mapped = mg.querymany(mouse_genes,
                                as_dataframe=True, species=['mouse','human'],
                                scopes='symbol', fields='symbol')

        # drop genes with no human ortholog
        if verbose:
            print('Raw', len(mg_mapped))
        mg_mapped = mg_mapped.dropna(subset=['symbol'])
        if verbose:
            print('With symbol', len(mg_mapped))
        # drop duplicates
        mg_mapped = mg_mapped[~mg_mapped.index.duplicated(keep='first')]
        if verbose:
            print('Deduplicated', len(mg_mapped))
            mg_mapped.head()
        return dict(mg_mapped['symbol'])
    elif map_using == 'mgi':
        if data_loc is not None:
            mrk_target = os.path.join(data_loc, 'MRK_List2.rpt')
            hmd_target = os.path.join(data_loc, 'HMD_HumanPhenotype.rpt')
        else:
            mrk_target = 'MRK_List2.rpt'
            hmd_target = 'HMD_HumanPhenotype.rpt'
        
        if not os.path.exists(mrk_target) or update:
            keep_url = "http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
            r_map = requests.get(keep_url, allow_redirects=True)
            open(mrk_target, 'wb').write(r_map.content)
        keep = pd.read_csv(mrk_target, sep="\t", usecols=["MGI Accession ID", "Marker Symbol",
                            "Feature Type", "Marker Name"])
        keep = keep.loc[keep["Feature Type"].isin(["protein coding gene"])].reset_index(drop=True)
        mgi_df["MGI"] = mgi_df.MGI_marker_accession_ID.apply(lambda x: x.split("|"))
        mgi_df = mgi_df.explode("MGI", ignore_index=True)
        mgi_df["MGI"] = [mg if type(mg) is str else mg[0] for mg in mgi_df.MGI]
        mgi_df = mgi_df.loc[mgi_df["MGI"].isin(keep["MGI Accession ID"])]
        mgi_df = mgi_df.merge(keep.loc[:, ("MGI Accession ID", "Marker Symbol")], left_on="MGI",
                            right_on="MGI Accession ID", how="left")

        if not os.path.exists(hmd_target) or update:
            map_url = "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
            r_map = requests.get(map_url, allow_redirects=True)
            open(hmd_target, 'wb').write(r_map.content)
        mapping = pd.read_csv(hmd_target, sep="\t", header=None, usecols=[0, 2, 3],
                            index_col=False, names=["symbol", "gene_name", "MGI"])
        mapping = mapping.loc[mapping["MGI"].isin(keep["MGI Accession ID"])]

        mg_mapped = mgi_df.merge(mapping, on="MGI", how="left")
        mg_mapped.loc[mg_mapped.symbol.isna(), "gene_name"] = mg_mapped.loc[mg_mapped.symbol.isna(), "Marker Symbol"]
        mg_mapped = mg_mapped.drop_duplicates()
        mg_mapped.rename(columns={"symbol": 'human_ortholog'}, inplace=True)
        
        return mg_mapped.set_index('gene_name')['human_ortholog'].to_dict()


def load_MPO(url='http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology', 
            use_ddot=False, update=False, data_loc=None):
    """
    Function to parse and load mouse phenotype ontology, using DDOT's ontology module

    :param url: URL containing MPO ontology file
    :type url: str
    :param use_ddot: Use the deprecated DDOT package to load the MGI ontology. Default False to use obonet
    :type use_ddot: boolean
    :return: MPO parsed using DDOT
    :rtype: :py:class:`ddot.Ontology`
    :raises ImportError: If DDOT package is not found
    """
    if data_loc is not None:
        obo_file_target = os.path.join(data_loc, url.split('/')[-1])
    else:
        obo_file_target = url.split('/')[-1]
    
    if use_ddot:
        # add deprecation warning
        warnings.warn('Use of ddot will be deprecated from version 0.1.9', DeprecationWarning)
        # download the mammalian phenotype ontology, parse with ddot
        if update or (not os.path.exists(obo_file_target)):
            r = requests.get(url,allow_redirects=True)
            open('MPheno_OBO.ontology','wb').write(r.content)
            if DDOT_LOADED is False:
                raise ImportError('ddot package is required to use this method')
        
        ddot.parse_obo('MPheno_OBO.ontology',
                       'parsed_mp.txt',
                      'id2name_mp.txt',
                      'id2namespace_mp.txt',
                      'altID_mp.txt')


        MP2desc = pd.read_csv('id2name_mp.txt',sep='\t',
                              names=['MP','description'],index_col='MP')

        MP2desc=MP2desc.loc[MP2desc.index.dropna()] # drop NAN from index
        print(len(MP2desc))

        hierarchy = pd.read_table('parsed_mp.txt',
                                  sep='\t',
                                  header=None,
                                  names=['Parent', 'Child', 'Relation', 'Namespace'])
        
        MPO = Ontology.from_table(
            table=hierarchy,
            parent='Parent',
            child='Child',
            add_root_name='MP:00SUPER',
            ignore_orphan_terms=True)

        # add description to node attribute
        terms_keep = list(np.unique(hierarchy['Parent'].tolist()+hierarchy['Child'].tolist()))
        MPO.node_attr=MP2desc.loc[terms_keep]
        
    else:
        if update or (not os.path.exists(obo_file_target)):
            MPO = obo.read_obo(url)
        else:
            MPO = obo.read_obo(obo_file_target)
    return MPO


def format_mapping(mapping, gene_col='human_ortholog', term_col='MP'):
        formatted_mapping = mapping.loc[:, (gene_col, term_col)].dropna().reset_index()
        formatted_mapping = formatted_mapping.loc[:, (gene_col, term_col)]
        formatted_mapping.columns = ["Gene", "Term"]
        return formatted_mapping


def map_genes_to_MPO(MPO, mapping, restrict_to=None, map_col='human_ortholog', MP_col='MP'):
    if ('Gene' not in mapping.columns) or ('Term' not in mapping.columns):
        mapping = format_mapping(mapping, gene_col=map_col, term_col=MP_col)
    MPO_mapped = MPO.copy()
    if restrict_to is not None:
        mapping = mapping.loc[mapping.Gene.isin(restrict_to)]
    mapping_dict = defaultdict(set)
    for row in mapping.iterrows():
        mapping_dict[row[1]['Term']].add(row[1]['Gene'])
    # add any terms not present in mapping with empty sets to ensure attribute exists for all terms in the ontology
    missing_terms = [t for t in MPO.nodes() if t not in mapping_dict]
    mapping_dict = {**mapping_dict, **{t:set() for t in missing_terms}}
    
    # set node attributes will ignore additional terms, attribute will not exist for term not in mapping.    
    nx.set_node_attributes(MPO_mapped, mapping_dict, 'term2genes')
    return MPO_mapped

def MPO_enrichment_root(hier_df,MPO,mgi_df,MP_focal_list,G_int,verbose=False, use_ddot=False,
                        min_genes=10, max_genes=2000):
    """
    Function to test for enrichment of genes resulting in selected phenotypes
    when knocked out in root node of NetColoc hierarchy.

    The returned :py:class:`pandas.DataFrame` will have the following columns:

    * **OR_p** - Odds ratio p-vlaue
    * **log_OR** - natural log odds ratio
    * **log_OR_CI_lower** - lower 95% confidence interval on log_OR
    * **log_OR_CI_upper** - upper 95% confidence interval on log_OR
    * **num_genes_in_term** - number of genes in MPO term
    * **MP_description** - description of MPO phenotype

    :param hier_df: NetColoc systems map (processed output from cdaps_util)
    :type hier_df: :py:class:`pandas.DataFrame`
    :param MPO: DDOT ontology containing the parsed mammalian phenotype ontology
    :type MPO: :py:class:`ddot.Ontology`
    :param mgi_df: parsed MGI knockout dataframe
    :type mgi_df: :py:class:`pandas.DataFrame`
    :param MP_focal_list: List of MPO phenotypes to check for enrichment against
    :type MP_focal_list: list
    :param G_int: Background interactome
    :type G_int: :py:class:`networkx.Graph`
    :param verbose: If true, print out some progress
    :type verbose: bool
    :param use_ddot: Use the deprecated DDOT package to load the MGI ontology. Default False to use obonet
    :type use_ddot: boolean
    :return: Dataframe containing enrichment results
    :rtype: :py:class:`pandas.DataFrame`
    """
    # test for enrichment in root node
    OR_p_list, OR_CI_list,log_OR_list = [], [], []
    num_genes_in_term_list = []

    MP_keep_list = []

    # root node is the largest node
    hier_df.index=hier_df['name']
    root_node = hier_df['CD_MemberList_Size'].sort_values(ascending=False).head(1).index.tolist()[0]

    # get list of node names from G_int
    G_int_nodes = list(G_int.nodes())

    for MP_focal in MP_focal_list:

        focal_terms, MP_desc_focal = get_focal_terms(MPO, MP_focal, use_ddot=use_ddot, include_parents=False, include_children=True)
        OR_p_temp, OR_CI_temp, log_OR_temp, mgi_genes = test_single_MPO_term(focal_terms, hier_df, root_node, mgi_df, G_int_nodes, verbose=verbose,
                                                                             min_genes=min_genes, max_genes=max_genes, name=MP_desc_focal)
        # test uf OR_p_temp is np.nan:
        if OR_p_temp is not None and not np.isnan(OR_p_temp):
            OR_p_list.append(OR_p_temp)
            OR_CI_list.append(OR_CI_temp)
            log_OR_list.append(log_OR_temp)
            num_genes_in_term_list.append(len(mgi_genes))

            MP_keep_list.append(MP_focal)
        else:
            if verbose:
                print(f'Number of intersecting genes for {MP_focal} in root node {root_node}: {len(mgi_genes)}')
                print(f'No genes found for {MP_focal} in root node {root_node}. Skipping enrichment test.')

    OR_CI_lower, OR_CI_upper = zip(*OR_CI_list)

    root_KO_df = pd.DataFrame({'OR_p':OR_p_list,'log_OR':log_OR_list,
                               'log_OR_CI_lower':OR_CI_lower,'log_OR_CI_upper':OR_CI_upper,
                              'num_genes_in_term':num_genes_in_term_list},
                              index=MP_keep_list)
    
    if use_ddot:
        warnings.warn('Use of ddot will be deprecated from version 0.1.9', DeprecationWarning)
        root_KO_df['MP_description'] = root_KO_df.index.map(dict(MPO.node_attr['description']))
    else:
        term_description_dict = {t: MPO.nodes[t]['name'] for t in MPO.nodes}
        root_KO_df['MP_description'] = root_KO_df.index.map(term_description_dict)

    return root_KO_df


def get_focal_terms(MPO, MP_focal, use_ddot=False, include_children=True, include_parents=False):
    """
    Function to get the focal terms for a given MPO term.
    :param MPO: DDOT ontology containing the parsed mammalian phenotype ontology
    :param MP_focal: Focal MPO term to get the description for
    :param use_ddot: Whether to use the DDOT package
    :param include_children: Whether to include child terms
    :param include_parents: Whether to include parent terms
    :return: Tuple containing the focal terms and their description
    """
    if use_ddot:
        warnings.warn('Use of ddot will be deprecated from version 0.1.9', DeprecationWarning)
        assert DDOT_LOADED, 'DDOT not successfully loaded. Please check installation or set use_ddot=False'

        MP_desc_focal = dict(MPO.node_attr['description'])[MP_focal]
        # focus the hierarchy on one branch, and look up all terms within that branch
        if len(MPO.parent_2_child[MP_focal])>0:
            MPO_focal = MPO.focus(MP_focal,verbose=False)
            focal_terms = MPO_focal.terms
        else:  # if the term has no children, just look at that term
            focal_terms=[MP_focal]
    else:
        MP_desc_focal = MPO.nodes[MP_focal]['name']
        if MPO.in_degree(MP_focal) > 0:
            focal_terms = focus_ontology(MPO, MP_focal, include_parents=include_parents, include_children=include_children)
        else:
                focal_terms=[MP_focal]
    return focal_terms, MP_desc_focal


def test_single_MPO_term(focal_terms, hier_df, hier_node, mgi_df, G_int_nodes, verbose=False,
                         min_genes=11, max_genes=2000, name=''):
    """
    Function to test for enrichment of genes resulting in selected phenotypes
    :param focal_terms: List of MPO phenotypes to check for enrichment against
    :param hier_df: Hierarchical DataFrame containing gene information
    :param hier_node: Node in the hierarchy to test
    :param mgi_df: DataFrame containing MGI gene information
    :param G_int_nodes: List of nodes in the interaction graph
    :param verbose: Whether to print detailed output
    :param min_genes: Minimum number of genes required for testing
    :param max_genes: Maximum number of genes allowed for testing
    :param name: Name of the test (for output purposes)
    :return: Tuple containing p-value, confidence interval, log odds ratio, and list of MGI genes
    """

    # check enrichment in root node
    focal_genes = hier_df['CD_MemberList'].loc[hier_node].split(' ')
    mgi_temp = mgi_df[mgi_df['MP'].isin(focal_terms)]
    mgi_temp = mgi_temp.dropna(subset=['human_ortholog'])
    mgi_genes = list(np.unique(mgi_temp['human_ortholog']))

    #new_index=[g.upper() for g in mgi_temp.index.tolist()]
    #mgi_temp.index=new_index

    #mgi_genes = mgi_temp.index.tolist()
    mgi_genes = list(np.intersect1d(mgi_genes,list(G_int_nodes)))

    # only test if there are at least min genes, and fewer than max genes
    if (len(mgi_genes)>=min_genes) and (len(mgi_genes) < max_genes):
        OR_p_temp, OR_CI_temp, log_OR_temp = perform_log_odds_z_test(mgi_genes, focal_genes, G_int_nodes, verbose=verbose, name=name)
    else:
        OR_p_temp = np.nan
        OR_CI_temp = (np.nan, np.nan)
        log_OR_temp = np.nan    
        mgi_genes = []
    return OR_p_temp, OR_CI_temp, log_OR_temp, mgi_genes


def perform_hypergeometric_test(hier_df, hier_node, G_int, mgi_genes):
    focal_genes = hier_df['CD_MemberList'].loc[hier_node].split(' ')
    n=len(focal_genes)
    M=len(list(G_int.nodes()))
    N=len(mgi_genes)
    x = len(np.intersect1d(focal_genes,mgi_genes))
    if x > 0:
        return hypergeom.sf(k=x-1,M=M,n=n,N=N), x, ' '.join(list(np.intersect1d(focal_genes,mgi_genes)))
    else:
        return np.nan, 0, ''

def perform_log_odds_z_test(mgi_genes, focal_genes, G_int_nodes, verbose=False, name=''):
    # is there anything to test?
    if (len(mgi_genes) == 0) or (len(focal_genes) == 0) or (len(G_int_nodes) == 0):
        return np.nan, (np.nan, np.nan), np.nan
    
    q00 = len(np.intersect1d(mgi_genes, focal_genes))
    q01 = len(mgi_genes)-q00
    q10 = len(focal_genes)-q00
    q11 = len(G_int_nodes)-q00-q01-q10
    table_temp = [[q00,q01],[q10,q11]]

    CT= contingency_tables.Table2x2(table_temp)
    OR_p_temp = CT.log_oddsratio_pvalue()
    OR_CI_temp = CT.log_oddsratio_confint()
    log_OR_temp = CT.log_oddsratio
    if verbose:
        print('\n'+name)
        print('number of genes in root node = '+str(len(focal_genes)))
        print('number of genes in focal MPO term = '+str(len(mgi_genes)))
        print('number overlapping genes = '+str(q00))
        print(OR_p_temp)
        print(OR_CI_temp)
        print(log_OR_temp)
    return OR_p_temp, OR_CI_temp, log_OR_temp

def MPO_enrichment_full(hier_df,MPO,mgi_df,MP_focal_list,G_int, use_ddot=False,
                        min_genes=10, max_genes=2000, verbose=False):
    """
    Function to test for enrichment of genes resulting in selected phenotypes
    when knocked out in every NetColoc system (not just root)

    The returned :py:class:`pandas.DataFrame` will have these columns:

    * **log(OR_p)** - -log10(Odds ratio p-vlaue)
    * **log_OR** - natural log odds ratio
    * **num_genes** - number of genes in MPO term overlapping with focal system
    * **gene_ids** - list of overlapping genes between MPO term and focal system

    :param hier_df: NetColoc systems map (processed output from cdaps_util)
    :type hier_df: :py:class:`pandas.DataFrame`
    :param MPO: DDOT ontology containing the parsed mammalian phenotype ontology
    :type MPO: :py:class:`ddot.Ontology`
    :param mgi_df: parsed MGI knockout dataframe
    :type mgi_df: :py:class:`pandas.DataFrame`
    :param MP_focal_list: List of MPO phenotypes to check for enrichment against
    :type MP_focal_list: list
    :param G_int: Background interactome
    :type G_int: :py:class:`networkx.Graph`
    :param use_ddot: Use the deprecated DDOT package to load the MGI ontology. Default False to use obonet
    :type use_ddot: boolean
    :return: Dataframe containing enrichment results
    :rtype: :py:class:`pandas.DataFrame`
    """
    MP_full_results_df=pd.DataFrame(index=hier_df.index.tolist())

    for MP_focal in MP_focal_list:
        focal_terms, MP_desc_focal = get_focal_terms(MPO, MP_focal, use_ddot=use_ddot, include_children=True, include_parents=False)

        hyper_p_list = []
        num_genes_list = []
        genes_id_list = []

        OR_p_list,OR_CI_list,log_OR_list=[],[],[]
        for focal_cluster in hier_df.index.tolist():
            
            OR_p_temp, OR_CI_temp, log_OR_temp, mgi_genes = test_single_MPO_term(focal_terms, hier_df, focal_cluster, mgi_df, list(G_int.nodes()), verbose=verbose,
                                                                             min_genes=min_genes, max_genes=max_genes, name=MP_desc_focal)
            
            if OR_p_temp is not None and not np.isnan(OR_p_temp):
                OR_p_list.append(OR_p_temp)
                OR_CI_list.append(OR_CI_temp)
                log_OR_list.append(log_OR_temp)
            else:
                OR_p_list.append(1)
                OR_CI_list.append(0)
                log_OR_list.append(0)
            
            hyper_p, x, genes_id = perform_hypergeometric_test(hier_df, focal_cluster, G_int, mgi_genes)
            
            if hyper_p is not None and not np.isnan(hyper_p):
                hyper_p_list.append(hyper_p)
                num_genes_list.append(x)
                genes_id_list.append(genes_id)
            else:
                hyper_p_list.append(1)
                num_genes_list.append(0)
                genes_id_list.append('')   
        MP_focal_df = pd.DataFrame({MP_desc_focal+':-log(OR_p)':-np.log10(OR_p_list),
                                    MP_desc_focal+':log_OR':log_OR_list,
                                    MP_desc_focal+':num_genes':num_genes_list,
                                    MP_desc_focal+':gene_ids':genes_id_list,
                                    MP_desc_focal+':hyper_p':hyper_p_list,
                                    },index=hier_df.index.tolist())

        # check that the column isn't already present in MP_full_results_df
        if MP_desc_focal+':-log(OR_p)' not in MP_full_results_df.columns.tolist():
            MP_full_results_df=MP_full_results_df.join(MP_focal_df)

    return MP_full_results_df


def find_related_terms(MPO, keywords, use_ddot=False):
    """Function to find terms in the MPO that match a list of keywords.
    
    :param MPO: DDOT ontology containing the parsed mammalian phenotype ontology
    :type MPO: :py:class:`networkx.MultiDiGraph`
    :param keywords: List of keywords to search for in the term descriptions
    :type keywords: list
    :param use_ddot: Use the deprecated DDOT package to load the MGI ontology. Default False to use obonet
    :type use_ddot: boolean
    :return: List of terms that match the keywords
    :rtype: list
    """
    focal_list = []
    if use_ddot:
        warnings.warn('Use of ddot will be deprecated from version 0.1.9', DeprecationWarning)
        assert DDOT_LOADED, 'DDOT not successfully loaded. Please check installation or set use_ddot=False'
        for t in MPO.node_attr.index.tolist():
            descr_temp = MPO.node_attr.loc[t]['description']
            if check_keywords(keywords, descr_temp):
                focal_list.append(t)
    else:
        for t in MPO.nodes:
            descr_temp = MPO.nodes[t]['name']
            if check_keywords(keywords, descr_temp):
                focal_list.append(t)
    return focal_list
            
            
def check_keywords(keywords, description):
    """Function to check if any of the keywords are present in the description.
    :param keywords: List of keywords to search for
    :type keywords: list
    :param description: Description to search in
    :type description: str
    :return: True if any keyword is found, False otherwise
    :rtype: bool
    """
    for keyword in keywords:
        if description.find(keyword) >-1:
            return True
    return False


def get_MP_description(term_id, ontology, use_ddot=False, include_definition=False):
    """Function to get the description of a given MPO term.
    :param term_id: ID of the MPO term
    :type term_id: str
    :param ontology: The ontology to get the term description from
    :type ontology: :py:class:`ddot.Ontology` or :py:class:`networkx.MultiDiGraph`
    :param use_ddot: Use the deprecated DDOT package to load the MGI ontology. Default False to use obonet
    :type use_ddot: boolean
    :param include_definition: If True, include the definition in the description
    :type include_definition: bool
    :return: Description of the MPO term
    :rtype: str
    """
    if use_ddot:
        warnings.warn('Use of ddot will be deprecated from version 0.1.9', DeprecationWarning)
        assert DDOT_LOADED, 'DDOT not successfully loaded. Please check installation or set use_ddot=False'
        return ontology.node_attr.loc[term_id, "description"]
    else:
        term_info = ontology.nodes[term_id]
        if include_definition:    
            return f'{term_info["name"]}: {term_info["def"]}'
        else:
            return term_info["name"]