"""
Demonstration files --- :mod:`pbxplore.demo`
============================================

PBxplore bundles a set of demonstration files. This module ease the access to
these files.

The path to the demonstration files is stored in :const:`DEMO_DATA_PATH`. This
constant can be accessed as :const:`pbxplore.demo.DEMO_DATA_PATH`, or as
:const:`pbxplore.DEMO_DATA_PATH`.

A list of the available demonstration files is available with the
:func:`list_demo_files` function. The same list with absolute path instead of
file names is provided by :func:`list_demo_files_absolute`.

.. autofunction:: pbxplore.demo.list_demo_files

.. autofunction:: pbxplore.demo.list_demo_files_absolute
"""

import os

DEMO_DATA_PATH=os.path.abspath(os.path.dirname(__file__))


def list_demo_files():
    """
    List the names of the bundled demo files

    File names starting with _ or . are not listed. This allows to omit
    __init__.py, and hiden files.
    """
    return [demo_file for demo_file in os.listdir(DEMO_DATA_PATH)
            if not demo_file[0] in '_.']

def list_demo_files_absolute():
    """
    List the absolute path to the bundled demo files

    File names starting with _ or . are not listed. This allows to omit
    __init__.py, and hiden files.
    """
    return [os.path.join(DEMO_DATA_PATH, demo_file) for demo_file
            in list_demo_files()]
