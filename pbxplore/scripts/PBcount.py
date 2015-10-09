#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute PB frequency along protein sequence.

2013 - P. Poulain, A. G. de Brevern
"""


# Use print as a function for python 3 compatibility
from __future__ import print_function

# Standard modules
import os
import sys
import argparse

# Local module
import pbxplore as pbx

# Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


# MAIN - program starts here

def user_input():
    """
    Handle the user parameters for PBcount.py.

    Returns
    -------
    options : the parsed arguments as parsed by `argparse`.

    """
    parser = argparse.ArgumentParser(
        description='Compute PB frequency along protein sequence.')
    # mandatory arguments
    parser.add_argument("-f", action="append", required=True,
                        help="name(s) of the PBs file (in fasta format)")
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")
    # optional arguments
    parser.add_argument("--first-residue", action="store", type=int, default=1,
                        dest="first_residue",
                        help="define first residue number (1 by default)")

    # parse the arguments
    options = parser.parse_args()

    # check options
    if options.first_residue and options.first_residue < 0:
        print("Warning: first residue is < 1.")

    # check input files
    for name in options.f:
        if not os.path.isfile(name):
            sys.exit("{0}: not a valid file. Bye.".format(name))

    return options


def pbcount_cli():
    """
    PBcount command line.
    """
    options = user_input()

    # read PBs files
    pb_name, pb_seq = pbx.io.read_several_fasta(options.f)

    # count PBs at each position of the sequence
    try:
        pb_count = pbx.analysis.count_matrix(pb_seq)
    except pbx.PB.SizeError:
        sys.exit("cannot compute PB frequencies / different sequence lengths")
    except pbx.PB.InvalidBlockError as e:
        msg = "'{0}' is not a valid protein block (abcdefghijklmnop)".format
        sys.exit(msg(e.block))

    # write PBs count file
    count_file_name = options.o + ".PB.count"
    with open(count_file_name, 'w') as outfile:
        pbx.io.write_count_matrix(pb_count, outfile, options.first_residue)
    print("wrote {0}".format(count_file_name))


if __name__ == '__main__':
    pbcount_cli()
