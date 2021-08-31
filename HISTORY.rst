=======
History
=======

0.1.4 (2021-08-31)
--------------------

* If import of DDOT package fails, only a warning message will be
  displayed unless user invokes ``netcoloc.validation.load_MPO()``
  in which case an ``ImportError`` is raised

* Fixed bug where ``z1_threshold`` parameter was being passed to ``z2_threshold`` parameter in
  ``netcoloc.network_cololcalization.calcualte_network_overlap`` method called by ``netcoloc.network_colocalization.calculate_network_overlap_subgraph`` method

0.1.3 (2021-08-18)
--------------------

* Added dependency `gprofiler-official <https://pypi.org/project/gprofiler-official>`__
  to ``setup.py`` and ``requirements.txt`` because this is used by
  ``network_colocalization.py``

0.1.2 (2021-08-17)
--------------------

* Added new `validation.py` module containing mouse knockout database
  functionality

0.1.1 (2021-08-06)
-------------------

* Fixed netcoloc imports in `netprop_zcore.py`


0.1.0 (2021-03-10)
------------------

* First release on PyPI.
