=====
Usage
=====

To use netcoloc in a project

.. code-block:: python

    import netcoloc
    
Example Notebooks
-----------------

For example usage, please see `NetColoc/example_notebooks <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks>`__.  

For detailed commentary on the NetColoc pipeline see the associated protocal publicaiton at Rosenthal, S.B., Wright, S.N., Liu, S. et al. **Mapping the common gene networks that underlie related diseases.** Nat Protoc 18, 1745–1759 (2023). `https://doi.org/10.1038/s41596-022-00797-1 <https://doi.org/10.1038/s41596-022-00797-1>`__

**Background**

Here we introduce NetColoc, a tool which evaluates the extent to which two gene sets are related in network space, i.e. the extent to which they are colocalized in a molecular interaction network, and interrogates the underlying biological pathways and processes using multiscale community detection. This framework may be applied to any number of scenarios in which gene sets have been associated with a phenotype or condition, including rare and common variants within the same disease, genes associated with two comorbid diseases, genetically correlated GWAS phenotypes, GWAS across two different species, or gene expression changes after treatment with two different drugs, to name a few. NetColoc relies on a dual network propagation approach to identify the region of network space which is significantly proximal to both input gene sets, and as such is highly effective for small to medium input gene sets.


NetColoc analysis of rare variants in Autism spectrum disorder (ASD) and Congenital Heart Disease (CHD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of NetColoc workflow on genes associated with rare exome variants across two conditions (ASD and CHD)

See `example_notebooks/ASD_CHD_NetColoc_analysis.ipynb <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks/ASD_CHD_NetColoc_analysis.ipynb>`__.


NetColoc analysis of rare variants in Autism spectrum disorder (ASD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of NetColoc workflow on genes associated with rare exome variants in ASD only, to demonstrate workflow from single gene set

See `example_notebooks/ASD_NetColoc_analysis.ipynb <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks/ASD_NetColoc_analysis.ipynb>`__.


An example notebook to evaluate NetColoc scores for a range of thresholds on scored input gene lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See `example_notebooks/Evalute_scored_input_gene_lists.ipynb <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks/Evalute_scored_input_gene_lists.ipynb>`__.


Introducing Quantitative NetColoc (QNetColoc)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This notebook demonstrates the implementation of quantitative network colocalization which allows input seed genes to be assigned individual scores prior to propagation. 

See `example_notebooks/Quantitative_NetColoc_QNetColoc_example.ipynb <https://github.com/ucsd-ccbb/NetColoc/tree/v110/example_notebooks/Quantitative_NetColoc_QNetColoc_example.ipynb>`__.


