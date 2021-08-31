# -*- coding: utf-8 -*-

'''Utility functions useful across multiple modules.
'''


def get_degree_binning(node_to_degree_dict, min_bin_size, lengths=None):
    """
    Groups nodes by degree into similarly sized bins. This function
    comes from
    `network_utilities.py of emreg00/toolbox <https://github.com/emreg00/toolbox/blob/master/network_utilities.py>`__


    Returns a tuple with following two values:

    * **list of bins** where each bin contains a list of nodes of similar degree
    * **mapping of degree to index of bin** dict mapping a degree to the index
      of the bin in the bins list which contains nodes of that degree

    :param node_to_degree_dict: Map of nodes to their degrees
    :type node_to_degree_dict: dict
    :param min_bin_size: minimum number of nodes each bin should contain.
    :type min_bin_size: int
    :param lengths: List of nodes to bin. If lengths is equal to None, then
                    all nodes will be binned
    :type lengths: list
    :return: (list of bins, mapping of degree to index of bin)
    :rtype: tuple
    """
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

    degree_index = 0
    while degree_index < len(degrees):
        # Add nodes of each degree to bin until bin reaches minimum bin size
        low = degrees[degree_index]
        nodes_of_certain_degree = degree_to_nodes[low]
        while len(nodes_of_certain_degree) < min_bin_size:
            degree_index += 1
            if degree_index == len(degrees):
                degree_index -= 1
                break
            nodes_of_certain_degree.extend(degree_to_nodes[degrees[degree_index]])

        high = degrees[degree_index]
        if len(nodes_of_certain_degree) >= min_bin_size:
            # For each degree represented in bin, set degree to bin index
            for deg in range(low, high + 1):
                degree_to_bin_index[deg] = len(bins)
            bins.append(nodes_of_certain_degree)
            bins_boundaries.append((low, high))
        else:
            # Combine last bin with second last bin, if last bin is too small
            bins[-1].extend(nodes_of_certain_degree)
            low_of_previous_bin, high_of_previous_bin = bins_boundaries[-1]
            bins_boundaries[-1] = (low_of_previous_bin, high)
            for deg in range(high_of_previous_bin + 1, high + 1):
                degree_to_bin_index[deg] = len(bins) - 1
            
        degree_index += 1 

    return bins, degree_to_bin_index
