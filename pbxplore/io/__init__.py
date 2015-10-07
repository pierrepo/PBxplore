#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input/Output in files --- :mod:`pbxplore.io`
============================================

Deals with writing and reading of files in various formats.

Fasta
-----

.. autofunction:: read_fasta

.. autofunction:: read_several_fasta

.. autofunction:: write_fasta

.. autofunction:: write_fasta_entry

Results of an assignation
-------------------------

.. autofunction:: write_phipsi

.. autofunction:: write_phipsi_entry

.. autofunction:: write_flat

Results af analyses
-------------------

.. autofunction:: write_count_matrix

.. autofunction:: write_neq
"""

from .fasta import read_fasta, read_several_fasta, write_fasta, write_fasta_entry
from .write import write_phipsi, write_phipsi_entry, write_flat, write_count_matrix, write_neq
