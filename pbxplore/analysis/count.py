#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import


# Third-party module
import numpy

# Local module
from .. import PB


def _assert_same_size(sequences):
    """
    Raise an exception is all sequence are not the same length.

    Parameters
    ----------
    sequences
        a list of sequences

    Raises
    ------
    pbxplore.PB.SizeError
        not all the sequences are the same length.
    """
    seq_size = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_size:
            raise PB.SizeError


def count_matrix(pb_seq):
    """
    Count the occurences of each block at each position.

    The occurence matrix has one row per sequence, and one column per block.
    The columns are ordered in as :const:`pbxplore.PB.NAMES`.

    Parameters
    ----------
    pb_seq
        a list of PB sequences.

    Returns
    -------
    pb_count : numpy array
        The occurence matrix.

    Raises
    ------
    pbxplore.PB.InvalidBlockError
        encountered an unexpected PB
    """
    _assert_same_size(pb_seq)
    pb_count = numpy.zeros((len(pb_seq[0]),  len(PB.NAMES)))
    for seq in pb_seq:
        for idx, block in enumerate(seq):
            if block in PB.NAMES:
                pb_count[idx, PB.NAMES.index(block)] += 1.0
            elif block not in ["Z", "z"]:
                raise PB.InvalidBlockError(block=block)
    return pb_count


def read_occurence_file(name):
    """
    Read an occurence matrix from a file.
    It will return the matrix as a numpy array and the indexes of residues.

    Parameters
    ----------
    name : str
        Name of the file.

    Returns
    -------
    count_mat : numpy array
        the occurence matrix without the residue number
    residues: list
        the list of residues indexes

    Raises
    ------
    ValueError
        when something is wrong about the file
    """

    # load count file
    # skip first row that contains PBs labels
    try:
        count = numpy.loadtxt(name, dtype=int, skiprows=1)
    except:
        raise ValueError("ERROR: {0}: wrong data format".format(name))

    # determine number of sequences compiled
    # use the sum of all residue at position 3
    # since positions 1 and 2 have no PBs assignement
    # and begin at 1 to not sum the index of the line (here is 3)
    sequence_number = sum(count[2, 1:])
    if sequence_number == 0:
        raise ValueError("ERROR: counting 0 sequences!")

    # read residues number
    residues = count[:, 0]
    # remove residue numbers (first column)
    count = count[:, 1:]

    # get index of first residue
    try:
        int(residues[0])
    except:
        raise ValueError("""ERROR: cannot read index of first residue.
                         Wrong data format in {0}""".format(name))

    return count, residues
