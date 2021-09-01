============
Installation
============

netcoloc depends on `DDOT`_ which has
to be manually installed.

`DDOT`_ can be installed one of two ways:

  1. To install `DDOT`_ by downloading
     the zip file of the source tree:

     .. code-block::

        wget https://github.com/idekerlab/ddot/archive/refs/heads/python3.zip
        unzip python3.zip
        cd ddot-python3
        python setup.py bdist_wheel
        pip install dist/ddot*py3*whl

  2. To install `DDOT`_ by cloning the repo:

     .. code-block::

        git clone --branch python3 https://github.com/idekerlab/ddot.git
        cd ddot
        python setup.py bdist_wheel
        pip install dist/ddot*py3*whl

Stable release
--------------

To install netcoloc, run this command in your terminal:

.. code-block:: console

    $ pip install netcoloc

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


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
