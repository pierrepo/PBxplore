#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
import math

# Third-party module
import numpy

# Local module
from . import utils


def _neq_per_residue(row_freq):
    """
    Compute the Neq of a vector of frequencies coresponding to one residue.
    The vector is a row of the frequency matrix.
    This function is intended to be called from the `numpy.apply_along_axis` function.

    Parameters
    ----------
    row_freq : 1D numpy array
        vector of frequencies

    Return
    ------
        The value of the Neq
    """

    H = sum([math.log(freq)*freq for freq in row_freq if freq != 0])

    return math.exp(-H)


def compute_neq(count_mat):
    """
    Compute the Neq for each residue from an occurence matrix.

    Parameters
    ----------
    count_mat : numpy array
        an occurence matrix returned by `count_matrix`.

    Returns
    -------
    neq_array : numpy array
        a 1D array containing the neq values
    """

    # get the frequency matrix
    freq = utils.compute_freq_matrix(count_mat)

    # Compute neq
    neq_array = numpy.apply_along_axis(_neq_per_residue, 1, freq)

    return neq_array
