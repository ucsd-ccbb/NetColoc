import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#import seaborn as sns
import networkx as nx
import pandas as pd
import random
import string
import scipy.stats
import network_prop
import sys
#import visJS2jupyter.visJS_module
#import visJS2jupyter.visualizations

# for parallel processing
#from joblib import Parallel, delayed
#import multiprocessing

def main(num_reps=10, seed_gene_file='HC_genes/ASD_HC_no_shared_200114.tsv',int_file='../interactomes/G_PCnet.gpickle', out_name='ASD',rand_method = 'degree_binning',save_fnew_rand=False):
    '''
    
    Calculate z-scores for heat propagation
    
    num_reps: number of random samplings to perform
    seed_gene_file: location of seed gene file
    int_file: location of interactome (networkx gpickle format)
    out_name: prefix for saving output
    rand_method: type of randomization (degree_binning should be used most often, degree_ks_test deprecated)
    save_fnew_rand: whether to save the full randomization output (beware can be a large file if large num_reps)
    
    Example command:
    python netprop_zscore.py 10 HC_genes/ASD_HC_no_shared_200114.tsv ../interactomes/G_PCnet.gpickle ASD degree_binning single

    
    '''
    
    # TODO: IMPROVE GENE INPUT
    # TODO: INTEGRATE OUTPUT WITH network_localization.py, and network_colocalization.py
    # TODO: Improve efficiency (currently takes hours to run with num_reps=5000)
    # TODO: IMPROVE COMMENTS
    
    
    print('number of randomizations = '+str(num_reps))
    print('background interactome = ' + int_file)
    print('randomization method = ' + rand_method)
    print('save Fnew rand = '+save_fnew_rand)
    
    num_reps = int(num_reps)
    # load interactome and select focal interactome
    Gint = nx.Graph()
    Gint = nx.read_gpickle(int_file)
    if 'None' in Gint.nodes():
        Gint.remove_node('None')


    # load HC genes
    HC_genes_temp = pd.read_csv(seed_gene_file,sep='\t',index_col='Unnamed: 0')
    seed_HC = [str(g[1:-1]).strip("'") for g in HC_genes_temp['seed_genes'].tolist()[0][1:-1].split(', ')]
  
    print(seed_gene_file+':')
    print(len(seed_HC))
    seed_HC = list(np.intersect1d(Gint.nodes(),seed_HC))
    print(len(seed_HC))
    
    # calculate the z-score
    # calc Wprime from Gint
    Wprime = network_prop.normalized_adj_matrix(Gint,conserve_heat=True)


    print('calculating z-scores: '+seed_gene_file)
    z_seed,Fnew_rand_seed = calc_zscore_heat(Gint,Wprime,seed_HC,num_reps=num_reps,rand_method=rand_method)
    z_seed.to_csv('z_'+out_name+'_'+str(num_reps)+'_reps_'+rand_method+'.tsv',sep='\t')
    if save_fnew_rand=='True': # if true, save out the vector of randoms (this can be a large file)
        pd.DataFrame(Fnew_rand_seed).to_csv('Fnew_'+out_name+'_rand'+str(num_reps)+'_reps_'+rand_method+'.tsv',sep='\t')

    
    
def calc_zscore_heat(Gint,Wprime,genes_D1,num_reps=10,ks_sig = 0.3,rand_method = 'degree_binning'):
    '''
    Helper function to calculate the z-score of heat values from one input seet of genes
    
    rand_method = 'degree_ks_test', or 'degree_binning'.  select the type of randomization
    '''
    seed_D1 = list(np.intersect1d(list(genes_D1),Gint.nodes()))
    Fnew_D1 = network_prop.network_propagation(Gint,Wprime,seed_D1,alpha=.5,num_its=20)
    
    num_focal_edges=len(nx.subgraph(Gint,seed_D1).edges())
    
    Fnew_rand_D1 = np.zeros([num_reps,len(Fnew_D1)])
    if rand_method == 'degree_ks_test':
        for r in range(num_reps):
            if (r%50)==0:
                print(r)
            # UPDATE 8/23/17 -- replace with randomly selecting seed nodes, checking for degree distribution equivalence

            p=0
            # resample until degree distributions are not significantly different
            while p<ks_sig:
                seed_D1_random = Gint.nodes()
                np.random.shuffle(seed_D1_random)
                seed_D1_random = seed_D1_random[0:len(seed_D1)]
                ks_stat,p=scipy.stats.ks_2samp(pd.Series(Gint.degree(seed_D1)),pd.Series(Gint.degree(seed_D1_random)))


            Fnew_rand_tmp = network_prop.network_propagation(Gint,Wprime,seed_D1_random,alpha=.5,num_its=20)
            Fnew_rand_tmp.loc[seed_D1_random]=np.nan # set seeds to nan so they don't bias results
            Fnew_rand_D1[r] = Fnew_rand_tmp.loc[Fnew_D1.index.tolist()]

    elif rand_method == 'degree_binning':
        bins = get_degree_binning(Gint,10)
        min_degree, max_degree, genes_binned = zip(*bins)
        bin_df = pd.DataFrame({'min_degree':min_degree,'max_degree':max_degree,'genes_binned':genes_binned})

        # create a lookup table for degree and index
        actual_degree_to_bin_df_idx = {}
        for i in range(0, bin_df['max_degree'].max() + 1):
            idx_temp = bin_df[ (bin_df['min_degree'].lt(i + 1)) & (bin_df['max_degree'].gt(i - 1)) ].index.tolist()
            if len(idx_temp) > 0: # there are some degrees which aren't represented in the graph
                actual_degree_to_bin_df_idx[i] = idx_temp[0]
        for r in range(num_reps):
            if (r%50)==0:
                print(r)
            # UPDATE 1/30/18 -- sample from degree bins
            seed_D1_random = []
            for g in seed_D1:
                degree_temp = nx.degree(Gint,g)
                # find genes with similar degrees to focal gene degree
                genes_temp = bin_df.loc[actual_degree_to_bin_df_idx[degree_temp]]['genes_binned']
                
                np.random.shuffle(genes_temp) # shuffle them
                while genes_temp[0] in seed_D1_random: # make sure the gene isn't already in the list
                    np.random.shuffle(genes_temp) # shuffle them
                seed_D1_random.append(genes_temp[0]) # build the seed_D1_random list
 

            Fnew_rand_tmp = network_prop.network_propagation(Gint,Wprime,seed_D1_random,alpha=.5,num_its=20)
            Fnew_rand_tmp.loc[seed_D1_random]=np.nan # set seeds to nan so they don't bias results
            Fnew_rand_D1[r] = Fnew_rand_tmp.loc[Fnew_D1.index.tolist()]


    z_score_D1 = (np.log(Fnew_D1)-np.nanmean(np.log(Fnew_rand_D1),axis=0))/np.nanstd(np.log(Fnew_rand_D1),axis=0)
    
    return z_score_D1, Fnew_rand_D1


def get_degree_binning(g, bin_size, lengths=None):
    '''
    This function comes from network_utilities.py of emregtoobox.  
    '''
    degree_to_nodes = {}
    for node, degree in g.degree:
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    #Change later?
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
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
