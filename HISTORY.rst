=======
History
=======

0.1.6 (2022-06-16)
--------------------

* `ipycytoscape <https://ipycytoscape.readthedocs.io/en/latest>`__ added as a dependency

* `ipywidgets <https://ipywidgets.readthedocs.io/en/latest>`__ added as a dependency

* Added `network_colocalization.sweep_input_pvals()` to sweep p-values and scores

* Added `network_colocalization.calculate_network_enrichment()` to sweep over z-score thresholds

* `netprop.get_individual_heats_matrix()` can take a networkx `Graph` object and internally call
  `netprop.get_normalized_adjancency_matrix()`. Documentation updated in both methods to note
  that the resulting matrices can be saved via `numpy.save()` and retrieved via `numpy.load()`

* `example_notebooks/ASD_CHD_NetColoc_analysis.ipynb` now visualizes hierarchy using
  `ipycytoscape <https://ipycytoscape.readthedocs.io/en/latest>`__

* `example_notebooks/ASD_CHD_NetColoc_analysis.ipynb` updated with a note about using `numpy.save()`
  and `numpy.load()` to save and retrieve result from `netprop.get_individual_heats_matrix()`
  




0.1.5 (2022-03-09)
--------------------

* Fixed divide by zero error seen when calculating cosine distance by updating `netprop.get_normalized_adjancency_matrix()`
  to properly normalize an adjacency matrix that is asymetric (UD-1863)

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
