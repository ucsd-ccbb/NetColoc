# -*- coding: utf-8 -*-

# External library imports
import ndex2
import networkx as nx
import numpy as np
import pandas as pd
import os

# Internal module convenience imports
from .netcoloc_utils import *
from .netprop import *

def __init__(self):
    pass

def netprop_zscore(seed_gene_file, seed_gene_file_delimiter=None, num_reps=10, alpha=0.5, interactome_file=None, interactome_uuid='f93f402c-86d4-11e7-a10d-0ac135e8bacf', ndex_server='public.ndexbio.org', ndex_user=None, ndex_password=None, out_name='out', save_final_heat=False, save_random_final_heats=False):
    '''Performs network heat propagation on the given interactome with the given
    seed genes, then returns the z-scores of the final heat values of each node
    in the interactome.

    The z-scores are calculated based on a null model, which is built by running
    the network propagation multiple times using randomly selected seed genes
    with similar degree distributions to the original seed gene set.

    Args:
        seed_gene_file (str): Location of file containing a delimited list of 
            seed genes.
        seed_gene_file_delimiter (str): Delimiter used to separate genes in seed
            gene file. (Default: any whitespace)
        num_reps (int): Number of times the network propagation algorithm should
            be run using random seed genes in order to build the null model.
            (Default: 10)
        alpha (float): Number between 0 and 1. Denotes the importance of the 
            propagation step in the network propagation, as opposed to the step
            where heat is added to seed genes only. Recommended to be 0.5 or
            greater. (Default: 0.5)
        interactome_file (str): Location of file containing the interactome in
            NetworkX gpickle format. Either the interactome_file argument or the
            interactome_uuid argument must be defined.
        interactome_uuid (str): UUID of the interactome on NDEx. Either the
            interactome_file argument or the interactome_uuid argument must be
            defined. (Default: The UUID of PCNet, the Parsimonious Composite
            Network: f93f402c-86d4-11e7-a10d-0ac135e8bacf)
        ndex_server (str): The NDEx server on which the interactome is stored.
            Only needs to be defined if interactome_uuid is defined. (Default: 
            ndexbio.org)
        ndex_user (str): The NDEx user that the interactome belongs to. Only
            needs to be defined if interactome_uuid is defined, and the 
            interactome is private.
        ndex_password (str): The password of the NDEx user's account. Only needs
            to be defined if interactome_uuid is defined, and the interactome is
            private.
        out_name (str): Prefix for saving output files. (Default: out)
        save_final_heat (bool): If this is set to true, then the raw network
            propagation heat scores for the original seed gene set will be saved
            in the form of a csv file in the current directory. (Default: False)
        save_random_final_heats (bool): If this is set to true, then the raw 
            network propagation heat scores for every repetition of the 
            algorithm using random seed genes will be saved in the form of a csv
            file in the current directory. (Beware: This can be a large file if 
            num_reps is large.) (Default: False)
    '''

    # TODO: INTEGRATE OUTPUT WITH network_localization.py, and network_colocalization.py
    # TODO: IMPROVE COMMENTS  

    # Process arguments

    # seed_gene_file
    seed_gene_file = os.path.abspath(seed_gene_file)
    
    #num_reps
    try:
        num_reps = int(num_reps)
    except:
        raise TypeError("The num_reps argument should be an integer")
    #int_file and int_uuid
    if interactome_file is None and interactome_uuid is None:
        raise TypeError("Either interactome_file or interactome_uuid argument must be provided")

    # Load interactome
    if interactome_file is not None:
        interactome_file = os.path.abspath(interactome_file)
        interactome = nx.Graph()
        interactome = nx.read_gpickle(interactome_file)
    else:
        interactome = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=interactome_uuid
        ).to_networkx()
    if 'None' in interactome.nodes():
        interactome.remove_node('None')
    nodes = list(interactome.nodes)
        
    # print out interactome num nodes and edges for diagnostic purposes
    print('number of nodes:')
    print(len(interactome.nodes))
    print('\nnumber of edges:')
    print(len(interactome.edges))

    # Load seed genes
    seed_file = open(seed_gene_file, 'r')
    seed_genes = list(np.intersect1d(nodes, seed_file.read().split(seed_gene_file_delimiter)))
    print('\nnumber of seed genes in interactome:')
    print(len(seed_genes))

    # Calculate individual_heats_matrix from interactome
    print('\ncalculating w_prime')
    w_prime = get_normalized_adjacency_matrix(interactome, conserve_heat=True)
    print('\ncalculating individual_heats_matrix')
    individual_heats_matrix = get_individual_heats_matrix(w_prime, alpha)

    # Calculate the z-score
    print('\nCalculating z-scores: ' + seed_gene_file)
    z_scores, final_heat, random_final_heats = calc_zscore_heat(individual_heats_matrix, nodes, dict(interactome.degree), seed_genes, num_reps=num_reps)

    # Save z-score results
    z_scores.to_csv('z_' + out_name + '_' + str(num_reps) + '_reps_.tsv', sep='\t')

    # If save_final_heat is true, save out the final heat vectore
    if save_final_heat == 'True':
        final_heat.to_csv('final_heat_' + out_name + '_' + str(num_reps) + '_reps_.tsv', sep='\t')

    # If save_random_final_heats is true, save out the vector of randoms (this can be a large file)
    if save_random_final_heats=='True': 
        pd.DataFrame(random_final_heats).to_csv('Fnew_'+out_name+'_rand'+str(num_reps)+'_reps_.tsv',sep='\t')

    return z_scores, random_final_heats
    
def calc_zscore_heat(individual_heats_matrix, nodes, degrees, seed_genes, num_reps=10, alpha=0.5,min_in_bin=10):
    '''
    Helper function to calculate the z-score of heat values from one input seet of genes
    '''

    final_heat = network_prop.network_propagation(individual_heats_matrix, nodes, seed_genes)   
    random_final_heats = np.zeros([num_reps, len(final_heat)])

    bins, actual_degree_to_bin_index = get_degree_binning(degrees, min_in_bin)
    
    for i in range(num_reps):
        if (i%250)==0:
            print(i)

        random_seed_genes = []
        for gene in seed_genes:
            degree = degrees[gene]

            # Find genes with similar degrees to focal gene degree
            genes_of_similar_degree = bins[actual_degree_to_bin_index[degree]]
            np.random.shuffle(genes_of_similar_degree) # shuffle them

            index = 0
            while genes_of_similar_degree[index] in random_seed_genes: # make sure the gene isn't already in the list
                index += 1
            
            random_seed_genes.append(genes_of_similar_degree[index]) # build the random_seed_genes list

        random_final_heat = network_propagation(individual_heats_matrix, nodes, random_seed_genes)
        random_final_heat.loc[random_seed_genes]=np.nan # set seeds to nan so they don't bias results
        random_final_heats[i] = random_final_heat

    z_scores = (np.log(final_heat) - np.nanmean(np.log(random_final_heats), axis=0)) / np.nanstd(np.log(random_final_heats), axis=0)
    
    return z_scores, final_heat, random_final_heats