#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Protein structure manipulation --- :mod:`pbxplore.structure`
============================================================

Constants
---------

* :data:`PDB_EXTENSIONS`
   list of file extensions corresponding to PDB files

* :data:`PDBx_EXTENSIONS`
   list of file extensions corresponding to PDBx/mmCIF files

Functions
---------

.. autofunction:: chains_from_files

.. autofunction:: chains_from_trajectory

Objects
-------

.. autoclass:: pbxplore.structure.structure.Chain

   .. automethod:: pbxplore.structure.structure.Chain.get_phi_psi_angles

.. autoclass:: pbxplore.structure.structure.Atom

.. autoclass:: pbxplore.structure.PDB.PDB

Exceptions
----------

.. autoexception:: pbxplore.structure.structure.ChainError

.. autoexception:: pbxplore.structure.structure.AtomError
"""

from .PDB import PDB_EXTENSIONS, PDBx_EXTENSIONS
from .loader import *
