#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Analyses -- :mod:`pbxplore.analysis`
====================================

Build occurence matrix
----------------------

.. autofunction:: count_matrix

.. autofunction:: read_occurence_file

.. autofunction:: plot_map

Compare protein block sequences
-------------------------------

.. autofunction:: compare

Visualize deformability
-----------------------

.. autofunction:: compute_neq

.. autofunction:: plot_neq

.. autofunction:: generate_weblogo

Cluster protein block sequences
-------------------------------

.. autofunction: distance_matrix

.. autofunction: hclust

Utils
-----

.. autofunction:: substitution_score

.. autofunction:: compute_freq_matrix

.. autofunction:: compute_score_by_position
"""

from .clustering import hclust, distance_matrix, RError
from .compare import compare
from .count import count_matrix, read_occurence_file
from .neq import compute_neq
from .utils import substitution_score, compute_freq_matrix, compute_score_by_position
from .visualization import *
