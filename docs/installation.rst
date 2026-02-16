============
Installation
============

To install netcoloc, run this command in your terminal:

.. code-block:: console

    $ pip install netcoloc

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

Prior Versions
--------------

Prior version `netcoloc v0.1.6 <https://doi.org/10.5281/zenodo.6654561>`__ was utilized in the *NetColoc* publication: Rosenthal, S. B. et al.,
"`Mapping the common gene networks that underlie related diseases <https://doi.org/10.1038/s41596-022-00797-1>`__." Nature Protocols (2023).
To install this version, please use the following command:

.. code-block:: console

     pip install netcoloc==0.1.6

And follow the additional installation instructions for DDOT at `<https://pypi.org/project/netcoloc/0.1.6/>`__.

The original source code and example notebooks can be acquired from Zenodo: `DOI:6654561 <https://doi.org/10.5281/zenodo.6654561>`__, or from GitHub:

.. code-block:: console

        git clone git@github.com:ucsd-ccbb/NetColoc.git
        git checkout -b v0.1.6 tags/v0.1.6


From sources
------------

The sources for NetColoc can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/ucsd-ccbb/NetColoc

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/ucsd-ccbb/NetColoc/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/ucsd-ccbb/NetColoc
.. _tarball: https://github.com/ucsd-ccbb/NetColoc/tarball/master
.. _DDOT: https://github.com/idekerlab/ddot
