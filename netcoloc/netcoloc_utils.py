# -*- coding: utf-8 -*-

'''Utility functions useful across multiple modules.
'''

def get_degree_binning(node_to_degree_dict, min_bin_size, lengths=None):
    '''
    Groups nodes by degree into similarly sized bins. This function comes from
    network_utilities.py of emregtoobox.

    Args:
        node_to_degree_dict (dict): Dictionary mapping nodes to their degrees.
        min_bin_size (int): The minimum number of nodes each bin should contain.
        lengths (list): List of nodes to bin. If lengths is equal to None, then
            all nodes will be binned.

    Returns:
        bins (list of lists): A list of bins. Each bin contains a list of nodes
            of similar degree.
        degree_to_bin_index (dict): A dictionary mapping a degree to the index 
            of the bin in the bins list which contains nodes of that degree.
            For example, in order to retreive the bin containing nodes of degree
            5: bins[degree_to_bin_index[5]]
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