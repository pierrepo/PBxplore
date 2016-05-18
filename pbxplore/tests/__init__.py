#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import importlib
from .. import __version__


def module_info(module_name):
    """
    Print information regarding a module, i.e version and install directory

    Parameters
    ----------
    module_name: string
        the name of the module
    """

    try:
        module = importlib.import_module(module_name)
        if hasattr(module, '__version__'):
            print("{0} version {1}".format(module_name, module.__version__))
        module_dir = os.path.dirname(module.__file__)
        print("{0} is installed in {1}".format(module_name, module_dir))
    except:
        print("{0} is not installed".format(module_name))


def system_info():
    """"
    Print information regarding the dependencies of PBxplore and Python version
    """

    print("PBxplore version {}".format(__version__))

    module_info("numpy")
    module_info("MDAnalysis")
    module_info("matplotlib")
    module_info("weblogolib")

    py_version = sys.version.replace('\n', '')
    print("Python version {0}".format(py_version))
