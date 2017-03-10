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

__version__ = "1.3.5"

from .structure.loader import *
from .assignment import assign
from . import PB
from . import io
from . import structure
from . import analysis


def test():
    import os

    try:
        import nose
    except ImportError:
        raise ImportError("Nose have to be installed for tests")

    # find the directory where the test package lives
    from . import tests
    test_dir = os.path.dirname(tests.__file__)
    argv = [test_dir]

    #Get informations about system
    tests.system_info()
    print("nose version {}".format(nose.__version__))

    # run nose
    try:
        return nose.main(defaultTest=argv)
    except SystemExit as e:
        return e.code
