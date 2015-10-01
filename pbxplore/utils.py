#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import


# Local module
from .core import PB

# Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


def _slice_matrix(mat, residue_min=1, residue_max=None):
    """
    Slice a matrix given the lower and upper bond in parameters.
    The matrix has to be a numpy array with one row per residue.
    The slice will occurs on the rows and the sub-matrix is returned.

    Parameters
    ----------
    mat : numpy array
        the matrix to slice
    residue_min: int
        the lower bound of residue frame
    residue_max: int
        the upper bound of residue frame

    Returns
    -------
    sub_mat: numpy 2D array
        the matrix sliced

    Exceptions
    ----------
    IndexError : when something is wrong about the lower/upper bound
    """

    if residue_max is None:
        residue_max = mat.shape[0]

    if residue_min <= 0 or residue_max <= 0:
        raise IndexError("Index start at 1")

    # range of indexes of the matrix
    residues_idx = range(1, mat.shape[0] + 1)

    if residue_min not in residues_idx or residue_max not in residues_idx:
        raise IndexError("Index out of range")

    if residue_min >= residue_max:
        raise IndexError("Lower bound > upper bound")

    return mat[residue_min - 1: residue_max]


def compute_freq_matrix(count_mat):
    """
    Compute a PB frequency matrix from an occurence matrix.

    The frequency matrix has one row per sequence, and one column per block.
    The columns are ordered in as PB.NAMES.

    Parameters
    ----------
    count_mat : numpy array
        an occurence matrix returned by `count_matrix`.

    Returns
    -------
    freq : numpy array
        The frequency matrix
    """

    # Assert the occurence matrix is in the good format
    nb_pb = count_mat.shape[1]
    if nb_pb != len(PB.NAMES):
        raise ValueError("Wrong number of PB({} should be {}".format(nb_pb, len(PB.NAMES)))

    # Retrieve the number of sequences compiled
    # use the sum of all residue at position 3 since positions 1 and 2 have no PBs assignement
    nb_sequences = sum(count_mat[2, :])

    return count_mat / float(nb_sequences)
