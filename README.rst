===============================
NetColoc
===============================

.. image:: https://img.shields.io/pypi/v/netcoloc.svg
        :target: https://pypi.python.org/pypi/netcoloc

.. image:: https://app.travis-ci.com/ucsd-ccbb/NetColoc.svg?token=nGm2MDyYBLqUCSB4exJZ&branch=main
    :target: https://app.travis-ci.com/ucsd-ccbb/NetColoc

.. image:: https://readthedocs.org/projects/netcoloc/badge/?version=latest
        :target: https://netcoloc.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6773330.svg
        :target: https://doi.org/10.5281/zenodo.6773330


Description
-----------

Here we introduce NetColoc, a tool which evaluates the extent to
which two gene sets are related in network space, i.e. the extent
to which they are colocalized in a molecular interaction network,
and interrogates the underlying biological pathways and processes
using multiscale community detection.

This framework may be applied to any number of scenarios in which
gene sets have been associated with a phenotype or condition,
including rare and common variants within the same disease,
genes associated with two comorbid diseases, genetically
correlated GWAS phenotypes, GWAS across two different species,
or gene expression changes after treatment with two different
drugs, to name a few.

NetColoc relies on a dual network propagation
approach to identify the region of network space which is
significantly proximal to both input gene sets, and as such is
highly effective for small to medium input gene sets.


Installation
--------------

NetColoc is available on `PyPI <https://pypi.org/>`__

::

     pip install netcoloc


Version Notes
~~~~~~~~~~~~~

Prior version `netcoloc v0.1.6 <https://doi.org/10.5281/zenodo.6654561>`__ was utilized in the *NetColoc* publication: Rosenthal, S. B. et al.,
"`Mapping the common gene networks that underlie related diseases <https://doi.org/10.1038/s41596-022-00797-1>`__." Nature Protocols (2023).
To install this version, please use the following command:

::

     pip install netcoloc==0.1.6

And follow the additional installation instructions at `<https://pypi.org/project/netcoloc/0.1.6/>`__.

The original source code and example notebooks can be acquired from Zenodo: `DOI:6654561 <https://doi.org/10.5281/zenodo.6654561>`__, or from GitHub:

::

        git clone git@github.com:ucsd-ccbb/NetColoc.git
        git checkout -b v0.1.6 tags/v0.1.6


Documentation
-------------

For a quick-start on NetColoc's functionality, please see the
`example notebooks <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks>`__.

**Usage Note:** Please follow steps in example notebooks for correct usage of NetColoc. At this time, individual functionalities have not been tested for independent use.

Dependencies
--------------

NetColoc requires the following python packages (automatically installed via ``pip install netcoloc``)


* `click <https://pypi.org/project/click>`__
* `matplotlib <https://pypi.org/project/matplotlib>`__
* `ndex2 <https://pypi.org/project/ndex2>`__
* `networkx <https://pypi.org/project/networkx>`__
* `numpy <https://pypi.org/project/numpy>`__
* `seaborn <https://pypi.org/project/seaborn>`__
* `tqdm <https://pypi.org/project/tqdm>`__
* `mygene >= 3.2.2 <https://pypi.org/project/mygene/>`__
* `scipy >= 1.5.3 <https://pypi.org/project/scipy/>`__
* `statsmodels <https://pypi.org/project/statsmodels/>`__
* `gprofiler-official >= 1.0.0 <https://pypi.org/project/gprofiler-official/>`__
* `ipywidgets <https://pypi.org/project/ipywidgets>`__
* `ipycytoscape <https://ipycytoscape.readthedocs.io/en/latest>`__
* `obonet <https://pypi.org/project/obonet/>`__
* `cdapsutil <https://pypi.org/project/cdapsutil/>`__




License
--------

* Free software: MIT license

Citing NetColoc
---------------

Rosenthal, Sara Brin, Sarah N. Wright, Sophie Liu, Christopher Churas, Daisy Chilin-Fuentes, Chi-Hua Chen, Kathleen M. Fisch, Dexter Pratt, Jason F. Kreisberg, and Trey Ideker.
"Mapping the common gene networks that underlie related diseases." Nature protocols (2023): 1-15. `<https://doi.org/10.1038/s41596-022-00797-1>`__

Other publications utilizing NetColoc:

* Rosenthal, S. B. et al. A convergent molecular network underlying autism and congenital heart disease. Cell Syst. 12, 1094-1107.e6 (2021). `10.1016/j.cels.2021.07.009 <http://dx.doi.org/10.1016/j.cels.2021.07.009>`__
* Wright, S. N. et al. Genome-wide association studies of human and rat BMI converge on synapse, epigenome, and hormone signaling networks. Cell Rep. 42, 112873 (2023). `10.1016/j.celrep.2023.112873 <http://dx.doi.org/10.1016/j.celrep.2023.112873>`__



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
