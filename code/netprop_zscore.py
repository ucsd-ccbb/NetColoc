import numpy as np
import networkx as nx
import pandas as pd
import random
import string
import scipy.stats
import network_prop
import sys
import ndex2
import time

def main(seed_gene_file, num_reps=10, interactome_file=None, out_name='out', alpha=0.5, seed_gene_file_delimiter=None, interactome_uuid='f93f402c-86d4-11e7-a10d-0ac135e8bacf', ndex_server='public.ndexbio.org', ndex_user=None, ndex_password=None, save_final_heat=False, save_random_final_heats=False):
    '''
    Calculate z-scores for heat propagation
    
    num_reps: Number of random samplings to perform
    seed_gene_file: Location of seed gene file (delimited list of genes)
    int_file: Location of interactome (networkx gpickle format)
    out_name: Prefix for saving output
    rand_method: Type of randomization (degree_binning should be used most often, degree_ks_test deprecated)
    save_fnew_rand: Whether to save the full randomization output (beware can be a large file if large num_reps)
    seed_gene_file_delimiter: The delimiter to use when parsing the seed gene file (default: any whitespace)
    int_uuid: NDEx UUID of interactome (can be used instead of int_file)
    int_server: NDEx server of interactome (default: public.ndexbio.org)
    
    Example command:
    python netprop_zscore.py 10 HC_genes/ASD_HC_no_shared_200114.tsv ../interactomes/G_PCnet.gpickle ASD degree_binning single
    '''
    # TODO: INTEGRATE OUTPUT WITH network_localization.py, and network_colocalization.py
    # TODO: Improve efficiency (currently takes hours to run with num_reps=5000)
    # TODO: IMPROVE COMMENTS  

    # Process arguments
    #num_reps
    try:
        num_reps = int(num_reps)
    except:
        raise TypeError("num_reps argument should be an integer")
    #int_file and int_uuid
    if interactome_file is None and interactome_uuid is None:
        raise TypeError("Either int_file or int_uuid argument must be provided to the netprop_zscore function")

    # Load interactome
    if interactome_file is not None:
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

    # Calculate w_double_prime from interactome
    print('\ncalculating w_prime')
    w_prime = network_prop.normalized_adj_matrix(interactome, conserve_heat=True)
    print('\ncalculating w_double_prime')
    w_double_prime = network_prop.get_w_double_prime(w_prime, alpha)

    # Calculate the z-score
    print('\nCalculating z-scores: ' + seed_gene_file)
    z_scores, final_heat, random_final_heats = calc_zscore_heat(w_double_prime, nodes, dict(interactome.degree), seed_genes, num_reps=num_reps)

    # Save z-score results
    z_scores.to_csv('z_' + out_name + '_' + str(num_reps) + '_reps_.tsv', sep='\t')

    # If save_final_heat is true, save out the final heat vectore
    if save_final_heat == 'True':
        final_heat.to_csv('final_heat_' + out_name + '_' + str(num_reps) + '_reps_.tsv', sep='\t')

    # If save_random_final_heats is true, save out the vector of randoms (this can be a large file)
    if save_random_final_heats=='True': 
        pd.DataFrame(random_final_heats).to_csv('Fnew_'+out_name+'_rand'+str(num_reps)+'_reps_.tsv',sep='\t')
    
def calc_zscore_heat(w_double_prime, nodes, degrees, seed_genes, num_reps=10, alpha=0.5):
    '''
    Helper function to calculate the z-score of heat values from one input seet of genes
    rand_method = 'degree_ks_test', or 'degree_binning'.  select the type of randomization
    '''

    final_heat = network_prop.network_propagation(w_double_prime, nodes, seed_genes)   
    random_final_heats = np.zeros([num_reps, len(final_heat)])

    start = time.time()
    bins, actual_degree_to_bin_index = get_degree_binning(degrees, 10)

    t0=time.time()
    for i in range(num_reps):
        if (i%10)==0:
            print(i)
            print(time.time()-t0)
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

        random_final_heat = network_prop.network_propagation(w_double_prime, nodes, random_seed_genes)
        random_final_heat.loc[random_seed_genes]=np.nan # set seeds to nan so they don't bias results
        random_final_heats[i] = random_final_heat

    z_scores = (np.log(final_heat) - np.nanmean(np.log(random_final_heats), axis=0)) / np.nanstd(np.log(random_final_heats), axis=0)
    
    return z_scores, final_heat, random_final_heats

def calc_zscore_heat_continuous(w_double_prime, nodes, degrees, seed_genes, num_reps=10, alpha=0.5):
    '''
    
    Same function, but allows continuous input
    
    Helper function to calculate the z-score of heat values from one input seet of genes
    rand_method = 'degree_ks_test', or 'degree_binning'.  select the type of randomization
    '''

    final_heat = network_prop.network_propagation_continuous(w_double_prime, nodes, seed_genes)  
    random_final_heats = np.zeros([num_reps, len(final_heat)])

    start = time.time()
    bins, actual_degree_to_bin_index = get_degree_binning(degrees, 10)

    # loop over bins and shuffle seed heat values within bins
    t0 = time.time()
    for i in range(num_reps):
        if (i%10)==0:
            print(i)
            print(time.time()-t0)
        seed_random = pd.Series(np.zeros(len(nodes)),index=nodes)
        for b in np.arange(len(bins)):
            bin_genes = bins[b]

            vals_temp = seed_genes.loc[bin_genes].tolist()
            np.random.shuffle(vals_temp)
            seed_random.loc[bin_genes]=vals_temp

        #print(seed_random.head())
        random_final_heat = network_prop.network_propagation_continuous(w_double_prime, nodes, seed_random)
        #random_final_heat.loc[random_seed_genes]=np.nan # set seeds to nan so they don't bias results
        random_final_heats[i] = random_final_heat

    z_scores = (np.log(final_heat) - np.nanmean(np.log(random_final_heats), axis=0)) / np.nanstd(np.log(random_final_heats), axis=0)
    
    return z_scores, final_heat, random_final_heats

def get_degree_binning(node_to_degree_dict, min_bin_size, lengths=None):
    '''
    This function comes from network_utilities.py of emregtoobox.  
    '''
    degree_to_nodes = {}
    for node, degree in node_to_degree_dict.items():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    
    degrees = degree_to_nodes.keys()
    degrees = list(degrees)
    degrees.sort()

    bins = []
    bins_boundaries = []
    degree_to_bin_index = {}
    i = 0
    while i < len(degrees):
        low = degrees[i]
        nodes_of_certain_degree = degree_to_nodes[low]
        while len(nodes_of_certain_degree) < min_bin_size:
            i += 1
            if i == len(degrees):
                break
            nodes_of_certain_degree.extend(degree_to_nodes[degrees[i]])
        if i == len(degrees):
            i -= 1
        high = degrees[i]
        if len(nodes_of_certain_degree) < min_bin_size:
            nodes_ = bins[-1]
            low_, high_ = bins_boundaries[-1]
            bins[-1] = nodes_ + nodes_of_certain_degree
            bins_boundaries[-1] = (low_, high)
            for d in range(high_ + 1, high + 1):
                degree_to_bin_index[d] = len(bins) - 1
        else:
            for d in range(low, high + 1):
                degree_to_bin_index[d] = len(bins) 
            bins.append(nodes_of_certain_degree)
            bins_boundaries.append((low, high))
        i += 1 
    return bins, degree_to_bin_index

if __name__ == "__main__":
    main(*sys.argv[1:])
