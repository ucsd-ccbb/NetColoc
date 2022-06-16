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
        
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6654561.svg
   :target: https://doi.org/10.5281/zenodo.6654561


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

**Usage Note:** Please follow steps in example notebooks for correct usage of NetColoc. At this time, individual functionalities have not been tested for independent use. 

Dependencies
--------------

NetColoc requires the following python packages:

.. note:: All of the following packages minus `DDOT <https://github.com/idekerlab/ddot>`__ and `cdapsutil <https://pypi.org/project/cdapsutil>`__ will be automatically installed via ``pip install netcoloc``

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

* `python3 branch of DDOT <https://github.com/idekerlab/ddot/tree/python3>`__

  `DDOT <https://github.com/idekerlab/ddot>`__ can be installed one of
  two ways:

  1. To install `DDOT <https://github.com/idekerlab/ddot>`__ by downloading
     the zip file of the source tree:

     .. code-block::

        wget https://github.com/idekerlab/ddot/archive/refs/heads/python3.zip
        unzip python3.zip
        cd ddot-python3
        python setup.py bdist_wheel
        pip install dist/ddot*py3*whl

  2. To install `DDOT <https://github.com/idekerlab/ddot>`__ by cloning the repo:

     .. code-block::

        git clone --branch python3 https://github.com/idekerlab/ddot.git
        cd ddot
        python setup.py bdist_wheel
        pip install dist/ddot*py3*whl


Additional requirements for full functionality of example notebook:

* `cdapsutil >= 0.2.0a1 <https://pypi.org/project/cdapsutil/>`__


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
