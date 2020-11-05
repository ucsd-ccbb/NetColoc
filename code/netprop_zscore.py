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

# TODO: Change argument order to put seed_gene_file first, since it's required? (Will break things)
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

    # Load seed genes
    seed_file = open(seed_gene_file, 'r')
    seed_genes = list(np.intersect1d(nodes, seed_file.read().split(seed_gene_file_delimiter)))

    # Calculate w_double_prime from interactome
    w_prime = network_prop.normalized_adj_matrix(interactome, conserve_heat=True)
    w_double_prime = network_prop.get_w_double_prime(w_prime, alpha)

    # Calculate the z-score
    print('Calculating z-scores: ' + seed_gene_file)
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

    bins = get_degree_binning(degrees, 10)
    min_degree, max_degree, genes_binned = zip(*bins)
    bin_df = pd.DataFrame({'min_degree':min_degree,'max_degree':max_degree,'genes_binned':genes_binned})

    # Create a lookup table for degree and index
    actual_degree_to_bin_df_idx = {}
    for i in range(0, bin_df['max_degree'].max() + 1):
        idx_temp = bin_df[ (bin_df['min_degree'].lt(i + 1)) & (bin_df['max_degree'].gt(i - 1)) ].index.tolist()
        if len(idx_temp) > 0: # there are some degrees which aren't represented in the graph
            actual_degree_to_bin_df_idx[i] = idx_temp[0]

    for i in range(num_reps):
        random_seed_genes = []
        for gene in seed_genes:
            degree = degrees[gene]

            # Find genes with similar degrees to focal gene degree
            genes_of_similar_degree = bin_df.loc[actual_degree_to_bin_df_idx[degree]]['genes_binned']
            np.random.shuffle(genes_of_similar_degree) # shuffle them

            index = 0
            while genes_of_similar_degree[index] in random_seed_genes: # make sure the gene isn't already in the list
                index += 1
            
            random_seed_genes.append(genes_of_similar_degree[index]) # build the seed_D1_random list

        random_final_heat = network_prop.network_propagation(w_double_prime, nodes, random_seed_genes)
        random_final_heat.loc[random_seed_genes]=np.nan # set seeds to nan so they don't bias results
        random_final_heats[i] = random_final_heat

    z_scores = (np.log(final_heat) - np.nanmean(np.log(random_final_heats), axis=0)) / np.nanstd(np.log(random_final_heats), axis=0)
    
    return z_scores, final_heat, random_final_heats

def get_degree_binning(degrees, bin_size, lengths=None):
    '''
    This function comes from network_utilities.py of emregtoobox.  
    '''
    degree_to_nodes = {}
    for node, degree in degrees.items():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    
    values = degree_to_nodes.keys()
    values = list(values)
    values.sort()

    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins



if __name__ == "__main__":
    main(*sys.argv[1:])
