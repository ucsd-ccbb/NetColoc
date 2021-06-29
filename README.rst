===============================
NetColoc
===============================

.. image:: https://img.shields.io/pypi/v/netcoloc.svg
        :target: https://pypi.python.org/pypi/netcoloc

.. image:: https://img.shields.io/travis/ceofy/netcoloc.svg
        :target: https://travis-ci.org/ceofy/netcoloc

.. image:: https://readthedocs.org/projects/netcoloc/badge/?version=latest
        :target: https://netcoloc.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://requires.io/github/ceofy/netcoloc/requirements.svg?branch=master
        :target: https://requires.io/github/ceofy/netcoloc/requirements?branch=master
        :alt: Dependencies


Description
-----------

* Free software: MIT license
* Documentation: https://netcoloc.readthedocs.io.

Here we introduce NetColoc, a tool which evaluates the extent to which two gene sets are related in network space, i.e. the extent to which they are colocalized in a molecular interaction network, and interrogates the underlying biological pathways and processes using multiscale community detection. This framework may be applied to any number of scenarios in which gene sets have been associated with a phenotype or condition, including rare and common variants within the same disease, genes associated with two comorbid diseases, genetically correlated GWAS phenotypes, GWAS across two different species, or gene expression changes after treatment with two different drugs, to name a few. NetColoc relies on a dual network propagation approach to identify the region of network space which is significantly proximal to both input gene sets, and as such is highly effective for small to medium input gene sets.


Documentation
-------------

For a quick-start on NetColoc's functionality, please see the example notebooks (LINK). 

Installation
------------

NetColoc requires the following python packages:

* click
* matplotlib
* ndex2
* networkx
* numpy
* seaborn
* tqdm

Additional requirements for full functionality of example notebook (INCLUDE VERSIONS):

* getpass
* ndex2
* scipy
* json
* cdapsutil
* ddot
* requests
* mygene
* gprofiler

Citing NetColoc
---------------


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
