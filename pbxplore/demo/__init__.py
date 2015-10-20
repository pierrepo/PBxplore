"""
Demonstration files --- :mod:`pbxplore.demo`
============================================

PBxplore bundles a set of demonstration files. This module ease the access to
these files.

The path to the demonstration files is stored in :const:`DEMO_DATA_PATH`. This
constant can be accessed as :const:`pbxplore.demo.DEMO_DATA_PATH`, or as
:const:`pbxplore.DEMO_DATA_PATH`.

A list of the available demonstration files is available with the
:func:`list_demo_files` function.

.. autofunction:: pbxplore.demo.list_demo_files
"""

import os

DEMO_DATA_PATH=os.path.abspath(os.path.dirname(__file__))


def list_demo_files():
    """
    List the names of the bundled demo files
    """
    return os.listdir(DEMO_DATA_PATH)
