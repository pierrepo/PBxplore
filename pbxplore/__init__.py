#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. function:: pbxplore.chains_from_files(path_list)

   See :func:`pbxplore.structure.chains_from_files`

.. function:: pbxplore.chains_from_trajectory(trajectory, topology)

   See :func:`pbxplore.structure.chains_from_trajectory` 

.. function:: pbxplore.assign(dihedrals)

   See :func:`pbxplore.assignment.assign`
"""

__version__ = "1.1"

from .structure.loader import *
from .assignment import assign
from . import PB
from . import io
from . import structure
from . import analysis
