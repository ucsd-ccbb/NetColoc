# -*- coding: utf-8 -*-

'''Utility functions useful across multiple modules.
'''

def get_degree_binning(node_to_degree_dict, min_bin_size, lengths=None):
    '''
    This function comes from network_utilities.py of emregtoobox.  
    '''
    # Create dictionary mapping degrees to nodes
    degree_to_nodes = {}
    for node, degree in node_to_degree_dict.items():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    
    # Get sorted list of degrees
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