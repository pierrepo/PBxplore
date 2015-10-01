#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Third-party module
import numpy

# Local module
from ..core import PB


def matrix_to_single_digit(matrix):
    """
    Convert a similarity substitution matrix to a single digit distance
    substitution matrix

    This function takes a substitution matrix expressed as similarity scores,
    and expresses it as integer distances in the range [0; 9]. Where 0 means
    similar PBs, and 9 means different PBs.

    Such single digit substitution matrix is convenient to display sequences
    of substitution scores.

    ..note:

        If a substitution matrix expressed with distances is provided
        rather than a matrix expressed with similarity scores, then the
        returned matrix is expressed with similarity scores.
    """
    mini = numpy.min(matrix)
    maxi = numpy.max(matrix)
    # Normalize between 0 and 1
    mat_modified = (matrix - mini)/(maxi - mini)
    # Convert similarity scores to distances and change the range to [0; 9]
    mat_modified = 9 * (1 - mat_modified)
    # Convert to integers
    mat_modified = mat_modified.astype(int)
    # Set diagonal to 0
    numpy.fill_diagonal(mat_modified, 0)
    return mat_modified


def compare_to_first_sequence(headers, sequences, substitution_mat):
    """
    Compare all sequence to the first one
    """
    iheaders = iter(headers)
    isequences = iter(sequences)
    ref_name = next(iheaders)
    ref_seq = next(isequences)
    for target_header, target_seq in zip(iheaders, isequences):
        header = "%s vs %s" % (ref_name, target_header)
        score_lst = PB.compute_score_by_position(substitution_mat, ref_seq, target_seq)
        yield header, score_lst
