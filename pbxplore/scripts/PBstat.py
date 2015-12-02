#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Statistical analysis and graphical representations of PBs.

Compute Neq, PBs distribution and draw logo representation of PBs.

2013 - P. Poulain, A. G. de Brevern
"""

# ============================================================================
# Modules
# ============================================================================
# Use print as a function for python 3 compatibility
from __future__ import print_function

# Standard modules
import os
import sys
import argparse

# Local module
import pbxplore as pbx

# Weblogolib is an optional requirement
try:
    import weblogolib
except:
    IS_WEBLOGO = False
else:
    IS_WEBLOGO = True

# matplotlib is an optional requirement
try:
    import matplotlib
except:
    IS_MATPLOTLIB = False
else:
    IS_MATPLOTLIB = True

# ============================================================================
# Python2/Python3 compatibility
# ============================================================================

# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


def user_inputs():
    """
    Handle the user parameters for PBstat.py.

    Returns
    -------
    options : the parsed arguments as parsed by `argparse`.
    """

    parser = argparse.ArgumentParser(
        description="Statistical analysis and graphical representations of PBs.")

    # mandatory arguments
    parser.add_argument("-f", action="store", required=True,
                        help="name of file that contains PBs frequency (count)")
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")

    # optional arguments
    parser.add_argument("--map", action="store_true", default=False, dest="mapdist",
                        help="generate map of the distribution of PBs along protein sequence")
    parser.add_argument("--neq", action="store_true", default=False, dest="neq",
                        help="compute Neq and generate Neq plot along protein sequence")
    parser.add_argument("--logo", action="store_true", default=False, dest="logo",
                        help="generate logo representation of PBs frequency along protein sequence")
    parser.add_argument("--image-format", action='store', type=str,
                        dest='image_format', default='png',
                        choices=['pdf', 'png', 'jpg'],
                        help='File format for all image output.')
    parser.add_argument("--residue-min", action="store", type=int,
                        dest="residue_min", help="defines lower bound of residue frame")
    parser.add_argument("--residue-max", action="store", type=int,
                        dest="residue_max", help="defines upper bound of residue frame")

    # get all parameters
    options = parser.parse_args()

    # check file
    if not os.path.isfile(options.f):
        parser.error("{0}: not a valid file".format(options.f))

    # Check residues min/max
    residues = [options.residue_min, options.residue_max]
    for residue in residues:
        if residue is not None and residue < 0:
            parser.error("residue argument must be >=0")
    if None not in residues and options.residue_min >= options.residue_max:
        parser.error("residue-min must be < residue-max.")

    # Check weblogo
    if options.logo:
        if not IS_WEBLOGO:
            parser.error("Weblogo is not installed; cannot generate the logo image.")

    # Check matplotlib
    if options.mapdist or options.neq:
        if not IS_MATPLOTLIB:
            parser.error("matplotlib is not installed; cannot generate plots.")

    return options


def check_residue_range(residues, residue_min, residue_max):
    """"
    Ensure that the lower bound and the upper bound parameters are in the range of
    the list `residues`.

    Parameters
    ----------
    residues : list
        the list of residues indexes
    residue_min:
        the lower bound of residue
    residue_max:
        the upper bound of residue

    Exceptions
    ----------
    IndexError : if `residue_min` or `residue_max` is not in the range
    """

    if residue_min is None:
        residue_min = residues[0]

    if residue_max is None:
        residue_max = residues[-1]

    if residue_min not in residues:
        raise IndexError("residue_min does not belong to the residue range")

    if residue_max not in residues:
        raise IndexError("residue_max does not belong to the residue range")

    if residue_min >= residue_max:
        raise IndexError("Lower bound > upper bound")

    return residue_min, residue_max


def pbstat_cli():
    """
    PBstat command line.
    """

    options = user_inputs()

    try:
        count, residues = pbx.analysis.read_occurence_file(options.f)
        residue_min, residue_max = check_residue_range(residues,
                                                       options.residue_min, options.residue_max)
    except Exception as e:
        sys.exit("ERROR: {0}".format(e))

    print("Index of first residue is: {0}".format(residue_min))

    # Handle output file name...
    output_file_name = options.o + ".PB.{0}"
    if options.residue_min or options.residue_max:
        output_file_name = "{0}.{1}-{2}".format(output_file_name, residue_min, residue_max)
    # ... and figure name
    output_fig_name = output_file_name + "." + options.image_format

    # -------------------------------------------------------------------------------
    # generates map of the distribution of PBs along protein sequence
    # -------------------------------------------------------------------------------
    if options.mapdist:
        file_fig_name = output_fig_name.format("map")
        pbx.analysis.plot_map(file_fig_name, count, residue_min, residue_max)
        print("wrote " + file_fig_name)

    # -------------------------------------------------------------------------------
    # computes Neq and generates neq plot along protein sequence
    # -------------------------------------------------------------------------------
    if options.neq:
        # compute Neq
        neq = pbx.analysis.compute_neq(count)

        # write Neq
        neq_file_name = output_file_name.format("Neq")
        with open(neq_file_name, "w") as outfile:
            pbx.io.write_neq(outfile, neq, residue_min, residue_max)
        print("wrote {0}".format(neq_file_name))

        # draw Neq
        file_fig_name = output_fig_name.format("Neq")
        pbx.analysis.plot_neq(file_fig_name, neq, residue_min, residue_max)
        print("wrote {}".format(file_fig_name))

    # -------------------------------------------------------------------------------
    # generates logo representation of PBs frequency along protein sequence
    # -------------------------------------------------------------------------------
    if options.logo:
        file_fig_name = output_fig_name.format("logo")
        title = options.f.replace(".PB.count", "")
        pbx.analysis.generate_weblogo(file_fig_name, count, residue_min, residue_max, title)
        print("wrote {}".format(file_fig_name))


if __name__ == '__main__':
    pbstat_cli()
