"""
-------------------------------------------
Authors: Brin Rosenthal, Sophie Liu
Date: 10/07/20
-------------------------------------------
"""


# -------------------- NETWORK CO-LOCALIZATION ---------------------------------#


 def calculate_network_overlap(z1,z2,zthresh=3):
    '''
    Function to determine which genes overlap
    TODO: improve documentation
    '''
    z_merged = z1.join(z2['z'],lsuffix='_1',rsuffix='_2')
    z_combined = z_merged['z_1']*z_merged['z_2']*(z_merged['z_1']>0)*(z_merged['z_2']>0)
    high_z_genes = z_combined[z_combined>=zthresh].index.tolist()
    
    return(high_z_genes)

 def calculate_network_overlap_subgraph(Gint,z1,z2,zthresh=3):
    '''
    Function to return subgraph of network intersection
    TODO: improve documentation
    '''
    z_merged = z1.join(z2['z'],lsuffix='_1',rsuffix='_2')
    z_combined = z_merged['z_1']*z_merged['z_2']*(z_merged['z_1']>0)*(z_merged['z_2']>0)
    high_z_genes = z_combined[z_combined>=zthresh].index.tolist()
    
    G_comb_sub = nx.subgraph(Gint,high_z_genes)
    
    return(G_comb_sub)

def calculate_expected_overlap(d1,d2,z1,z2,plot=False,zthresh=4,numreps=1000):    
    '''
    Function to determine size of expected network overlap
    d1: name of geneset 1
    d2: name of geneset 2
    z1: zscore for geneset 1
    z2: zscore for geneset 2
    
    TODO: improve documentation
    
    '''

    z_d1d2_size=len(calculate_network_overlap(z1,z2,zthresh=zthresh))

    high_z_rand = []
    for r in np.arange(numreps):
        # use permutation shuffling method instead of Fnew comparison
        d1_shuf_genes = z1.index.tolist()
        np.random.shuffle(d1_shuf_genes)
        d1_shuf=z1[:]
        d1_shuf.index=d1_shuf_genes

        d2_shuf_genes = z2.index.tolist()
        np.random.shuffle(d2_shuf_genes)
        d2_shuf=z2[:]
        d2_shuf.index=d2_shuf_genes

        high_z_temp = len(calculate_network_overlap(d1_shuf,d2_shuf,zthresh=zthresh))
        high_z_rand.append(high_z_temp)

    if plot==True:
        sns.distplot(high_z_rand,label='expected network intersection size')
        plt.plot([z_d1d2_size,z_d1d2_size],[0,0.015],label='observed '+d1+'-'+d2+' network intersection size')
        plt.xlabel('size of proximal subgraph, z>='+str(zthresh),fontsize=16)
        plt.legend(fontsize=12)
    return z_d1d2_size,high_z_rand


# TODO: add plotting functionality