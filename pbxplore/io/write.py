#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
# import os

# Local modules
from .. import PB
from ..analysis import utils


def write_count_matrix(pb_count, outfile, first=1):
    """
    Write a PB occurence matrix in a file.

    Parameters
    ----------
    pb_count: an occurence matrix as a 2D numpy array.
    outfile: an open file where to write the matrix.
    first: the residue number of the first position.
    """
    # write the header (PB names)
    print("    " + "".join(["%6s" % name for name in PB.NAMES]), file=outfile)
    # write the data table
    for residue_idx, residue_pb in enumerate(pb_count):
        print("%-5d" % (residue_idx + first) +
              " ".join("%5d" % i for i in residue_pb), file=outfile)


def write_neq(outfile, neq_array, residue_min=1, residue_max=None):
    """
    Write the Neq matrix in an open file

    Parameters
    ----------
    outfile : file descriptor
        The file descriptor to write in. It must allow writing.
    neq_array : numpy array
        a 1D array containing the neq values.
    residue_min: int
        the lower bound of residue frame
    residue_max: int
        the upper bound of residue frame

    """

    # Slice
    neq = utils._slice_matrix(neq_array, residue_min, residue_max)

    print("%-6s %8s " % ("resid", "Neq"), file=outfile)
    for (res, neq) in enumerate(neq):
        print("%-6d %8.2f " % (res + residue_min, neq), file=outfile)
