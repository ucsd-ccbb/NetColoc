# -*- coding: utf-8 -*-

'''Functions for performing network localization
'''

# External library imports
import copy
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tqdm.auto import tqdm
import warnings

# Internal module convenience imports
from .netcoloc_utils import *
#from netcoloc_utils import *


def __init__(self):
    pass

# -------------------- LOCALIZATION ---------------------------------#


def netprop_localization(z_scores, random_final_heats,
                         seed_genes, z_score_threshold=3, plot=True):
    """
    Calculates the size of a gene set's network proximity based on the
    z-score results of network propagation, and evaluates if it is larger than
    chance by repeating the algorithm using random degree-matched seed genes.
    Returns a z-score.

    :param z_scores: The output of the
                     :py:func:`~netcoloc.netprop_zscore.netprop_zscore`
                     method. A pandas Series where the index
                     contains gene names, and the data contains
                     the z-scores associated with each gene after
                     network propagation.
    :type z_scores: :py:class:`pandas.Series`
    :param random_final_heats: The output of the
                               :py:func:`~netcoloc.netprop_zscore.netprop_zscore`
                               method. A matrix containing ``z_scores``,
                               where each column corresponds to a gene, and
                               each row corresponds to a network propagation
                               performed using random seed genes.
    :type random_final_heats: :py:class:`numpy.ndarray`
    :param seed_genes: list of seed genes used to generate the z_scores
            and random_final_heats arguments.
    :type seed_genes: list
    :param z_score_threshold: The threshold after which genes are
                              considered significantly proximal to each
                              other.
    :type z_score_threshold: float
    :param plot: If True, then the distribution will be plotted
    :type plot: bool
    :return: The z-score representing how much larger the network
             proximity of this seed gene is than would be expected
             by chance.
    :rtype: float
    """
    nan_row = [np.nan] * len(random_final_heats[0])

    random_proximal_network_sizes = []
    for row in tqdm(range(len(random_final_heats))):
        # Calculate mean and standard deviation, excluding focal row
        focal_row = copy.deepcopy(random_final_heats[row])
        random_final_heats[row] = nan_row
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            log = np.log(random_final_heats)
            mean = np.nanmean(log, axis=0)
            standard_deviation = np.nanstd(log, axis=0)
            # Calculate z-scores for each gene in random network
            random_z_score = (np.log(focal_row) - mean) / standard_deviation
        # Count gene in proximal network, not including seed genes
        random_proximal_network_sizes.append(sum((random_z_score > z_score_threshold)) - len(seed_genes))
        # Replace focal row in random final heats matrix
        random_final_heats[row] = copy.deepcopy(focal_row)

    # Calculate the size of the true proximal subgraph, not including seed genes
    proximal_network_size = sum((z_scores > z_score_threshold)) - len(seed_genes)

    # Calculate z-score
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        z_score = (np.log(proximal_network_size) - np.nanmean(np.log(random_proximal_network_sizes))) / np.nanstd(np.log(random_proximal_network_sizes))

    # Plot figure
    if plot:
        plt.figure(figsize=(5, 4))
        dfig = sns.histplot(random_proximal_network_sizes, label='Random')
        plt.vlines(proximal_network_size, ymin=0, ymax=dfig.dataLim.bounds[3], color='r', label='True gene set')
        plt.xlabel('Size of proximal network, z > ' + str(z_score_threshold), fontsize=16)
        plt.ylabel('Count', fontsize=16) 
        plt.legend(loc='upper left')

    return z_score, proximal_network_size, random_proximal_network_sizes, 
