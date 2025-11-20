# -*- coding: utf-8 -*-

'''Utility functions useful across multiple modules.
'''

import warnings
import pandas as pd
import os
import numpy as np

class Seeds:
    def __init__(self, inputdata, option='score_file', agg_method='mean'):
        if isinstance(inputdata, str):
            # check that file exists
            if not os.path.exists(inputdata):
                raise FileNotFoundError(f'File {inputdata} not found')
            self.datafile = inputdata
            try:
                self.data = self._read_data(option='score_file')
                assert 'gene' in self.data.columns
            except AssertionError:
                self.data = self._read_data(option='list')
                assert 'gene' in self.data.columns
                
        elif isinstance(inputdata, pd.DataFrame):
            self.data = inputdata
        elif isinstance(inputdata, dict):
            self.data = pd.DataFrame({'gene': list(inputdata.keys()), 'score': list(inputdata.values())})
        elif isinstance(inputdata, list):
            self.data = pd.DataFrame({'gene': inputdata, 'score': [1]*len(inputdata)})
        
        # check if there are duplicate genes in self.data
        if (len(self.data['gene']) != len(set(self.data['gene']))):
            self._aggregate_scores(aggfunc=agg_method)
        
        self.genes = set(self.data['gene'])
        self.original_genes = self.genes.copy()
        self.scores = dict(zip(self.data['gene'], self.data['score']))    
        self.original_scores = self.scores.copy()
        self.describe_scores()
        
    def _read_data(self, gene_col='Entrez', score_col='P-value', option='score_file'):
        assert option in ['score_file', 'list'], f'Invalid option {option}. Must be `score_file` or `list`'
        if option == 'score_file':
            data = pd.read_csv(self.datafile, sep='\t')
            data.rename(columns={gene_col: 'gene', score_col: 'score'}, inplace=True)
            return data
        if option == 'list':
            data = pd.read_csv(self.datafile, header=None, names=['gene'])
            data['score'] = 1e-10
            return data
    
    def get_genes(self):
        return self.genes
    
    def get_score(self):
        return self.scores
    
    def describe_scores(self):
        return pd.Series(self.scores).describe()
    
    def get_seed_dict(self):
        return self.scores
    
    def transform_scores(self, method='log'):
        assert method in ['log', 'log10', 'log2', 'neglog10'], f'Invalid method {method}. Must be one of `log`, `log10`, `log2`, `neglog10`'
        assert all(score >= 0 for score in self.scores.values()), 'Scores must be positive for log transformation'
        
        if method == 'log10':
            self.scores = {gene: np.log10(score) for gene, score in self.scores.items()}
        elif method == 'log2':
            self.scores = {gene: np.log2(score) for gene, score in self.scores.items()}
        elif method == 'neglog10':
            self.scores = {gene: -1 * np.log10(max(score, 1e-250)) for gene, score in self.scores.items()}
        
    def normalize_scores(self, method='max', score_cap=None):

        assert method in ['max', 'minmax', 'zscore', 'sum', 'log'], f'Invalid method {method}. Must be one of `max`, `minmax`, `zscore`, `sum`, `log`'
        if score_cap is not None:
            self.scores = {gene: min(score, score_cap) for gene, score in self.scores.items()}
        if method == 'max':
            max_score = max(self.scores.values())
            self.scores = {gene: score/max_score for gene, score in self.scores.items()}
        elif method == 'minmax':
            min_score = min([x for x in self.scores.values() if x > 0])
            max_score = max(self.scores.values())
            self.scores = {gene: (score - min_score)/(max_score - min_score) for gene, score in self.scores.items()}
        elif method == 'sum':
            sum_score = sum(self.scores.values())
            self.scores = {gene: score/sum_score for gene, score in self.scores.items()}
        elif method == 'log':
            self.scores = {gene: max(np.log(score), 0) for gene, score in self.scores.items()}
            max_score = max(self.scores.values())
            self.scores = {gene: score/max_score for gene, score in self.scores.items()}
        elif method == 'zscore':
            mean_score = np.mean(list(self.scores.values()))
            std_score = np.std(list(self.scores.values()))
            self.scores = {gene: (score - mean_score)/std_score for gene, score in self.scores.items()}
            # remove negative scores (so only keep the top ~50% of genes)
            self.scores = {gene: score if score > 0  else 0 for gene, score in self.scores.items()}
        
    def reset_seeds(self):
        self.genes = self.original_genes.copy()
        self.scores = self.original_scores.copy()
    
    def filter_seeds_by_network(self, network_nodes):
        filtered_seeds = [seed for seed in list(self.genes) if seed in network_nodes]
        print(f'{len(self.genes) - len(filtered_seeds)}/{len(self.genes)} seeds not found in network')
        self.genes = set(filtered_seeds)
        self.scores = {gene: self.scores[gene] for gene in self.genes}
    
    def filter_seeds_by_score(self, min_score=None, max_score=None):
        if min_score is not None:
            filtered_seeds_min = [seed for seed in self.genes if self.scores[seed] >= min_score]
        else:
            filtered_seeds_min = self.genes.copy()
        if max_score is not None:
            filtered_seeds_max = [seed for seed in self.genes if self.scores[seed] <= max_score]
        else:
            filtered_seeds_max = self.genes.copy()
        filtered_seeds = set(filtered_seeds_min).intersection(filtered_seeds_max)
        
        print(f'{len(self.genes) - len(filtered_seeds)}/{len(self.genes)} seeds removed due to low/high score(s)')
        self.genes = set(filtered_seeds)
        self.scores = {gene: self.scores[gene] for gene in self.genes}
    
    def _aggregate_scores(self, aggfunc='mean'):
        # Check and aggregate multiple scores for the same gene?
        duplicate_genes = self.data[self.data.duplicated(subset='gene', keep=False)]['gene']
        if len(duplicate_genes) > 0:
            print(f'{len(duplicate_genes)} genes have multiple scores. Aggregating scores using `{aggfunc}`. To change aggregation method, use `agg_method` parameter')
            duplicated_data = self.data[self.data['gene'].isin(duplicate_genes)]
            unduplicated_data = self.data[~self.data['gene'].isin(duplicate_genes)]
            agg_data = duplicated_data.groupby('gene').agg(aggfunc).reset_index()
            self.data = pd.concat([unduplicated_data, agg_data], axis=0).reset_index(drop=True)
    
    def get_top_ranked_genes(self, k=500, ascending=False):
        return self.data.sort_values(by='score', ascending=ascending).gene[:k]


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
    assert min_bin_size <= len(node_to_degree_dict), f'Minimum bin size must be less than number of nodes {len(node_to_degree_dict)}'

    if lengths is not None:
        if len(lengths) == 0:
            warnings.warn("Lengths is empty. Returning empty bins and degree to bin index dictionary.")
            return [], {}
        missing_nodes = [x for x in lengths if x not in node_to_degree_dict]
        if len(missing_nodes) > 0:
            warnings.warn(f"The following nodes are not in the degree dictionary: {missing_nodes}")

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

from datetime import datetime

class Timer:
    """ This class is a simple timer that allows for the tracking of time elapsed between starting and ending tasks.
    It is useful for tracking time spent in different parts of a program or for timing different processes.
    """
    def __init__(self):
        self.start_times = {}
        self.finish_times = {}
        self.elapsed_times = {}
        self.tasks = []
        self.current_task_stack = []
        self.indents = {}
        
    def start(self, taskstr):
        if taskstr in self.start_times.keys():
            taskstr=taskstr + "1"
            i=1
            while taskstr in self.start_times.keys():
                i += 1
                taskstr = taskstr[0:-1] + str(i) 
        self.current_task_stack.append(taskstr)
        self.indents[taskstr] = len(self.current_task_stack) - 1
        self.tasks.append(taskstr)
        self.start_times[taskstr] = datetime.now()
        
    def end(self, taskstr):
        if taskstr in self.finish_times:
            matching_tasks = [task for task in self.start_times.keys() if taskstr in task]
            taskstr = matching_tasks[-1]
        self.current_task_stack.remove(taskstr)
        self.finish_times[taskstr] = datetime.now()
        self.elapsed_times[taskstr] = str(self.finish_times[taskstr] - self.start_times[taskstr])
        
    def print_all_times(self):
        try:
            for task in self.tasks:
                if task not in self.elapsed_times:
                    self.end(task)
                if self.indents[task] > 0:
                    print("".join(["|", "---"*self.indents[task], ">"]),self.elapsed_times[task], task)
                else:
                    print(self.elapsed_times[task], task)
        except:
            print(self.elapsed_times)
