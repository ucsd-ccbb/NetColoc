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


Documentation
-------------

For a quick-start on NetColoc's functionality, please see the
`example notebooks <https://github.com/ucsd-ccbb/NetColoc/tree/main/example_notebooks>`__.

Dependencies
--------------

NetColoc requires the following python packages:

* click
* matplotlib
* ndex2
* networkx
* numpy
* seaborn
* tqdm
* `ddot <https://github.com/Ceofy/ddot>`__
* `mygene >= 3.2.2 <https://pypi.org/project/mygene/>`__
* `scipy >= 1.5.3 <https://pypi.org/project/scipy/>`__
* `statsmodels <https://pypi.org/project/statsmodels/>`__

Additional requirements for full functionality of example notebook:

* `cdapsutil >= 0.2.0a1 <https://pypi.org/project/cdapsutil/>`__
* `gprofiler-official >= 1.0.0 <https://pypi.org/project/gprofiler-official/>`__




Installation
--------------

NetColoc is available on `PyPI <https://pypi.org/>`__

::

     pip install netcoloc

License
--------

* Free software: MIT license

Citing NetColoc
---------------

Coming soon...

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
