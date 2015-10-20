#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Print info on the demonstration files bundled with PBxplore

2013 - P. Poulain, A. G. de Brevern
"""


# Use print as a function for python 3 compatibility
from __future__ import print_function

# Standard modules
import os
import sys
import glob
import argparse

# Local modules
import pbxplore as pbx

def user_inputs():
    """
    Handle the user parameter from the command line
    """
    parser = argparse.ArgumentParser(
        description='Print info on the demonstration files bundled '
                    'with PBxplore')

    # arguments
    parser.add_argument('--list', '--ls', '-l',
                        action="store_true", default=False,
                        help="List the available files instead of "
                             "the directory path.")
    parser.add_argument('--list-abs', '-L',
                        action="store_true", default=False,
                        help="List the absolute path to the available files "
                             "instead of the directory path.")
    # get all arguments
    options = parser.parse_args()
    return options


def pbdata_cli():
    options = user_inputs()
    if options.list:
        for demo_file in pbx.demo.list_demo_files():
            print(demo_file)
    elif options.list_abs:
        for demo_file in pbx.demo.list_demo_files_absolute():
            print(demo_file)
    else:
        print(pbx.DEMO_DATA_PATH)

