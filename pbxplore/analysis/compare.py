#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Third-party module
import numpy

# Local module
from . import utils
from ..io import fasta


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
        score_lst = utils.compute_score_by_position(substitution_mat, ref_seq, target_seq)
        yield header, score_lst


def compare(header_lst, seq_lst, substitution_mat, fname):
    """
    Command line wrapper for the comparison of all sequences with the first one

    When the --compare option is given to the command line, the program
    compares all the sequences to the first one and writes these comparison as
    sequences of digits. These digits represent the distance between the PB
    in the target and the one in the reference at the same position. The digits
    are normalized in the [0; 9] range.

    This function run the comparison, write the result in a fasta file, and
    display on screen informations about the process.

    Parameters
    ----------
    header_lst: list of strings
        The list of sequence headers ordered as the sequences
    seq_lst: list of strings
        The list of sequences ordered as the headers
    substitution_mat: numpy.array
        A substitution matrix expressed as similarity scores
    fname: str
        The output file name
    """
    ref_name = header_lst[0]
    substitution_mat_modified = matrix_to_single_digit(substitution_mat)
    print("Normalized substitution matrix (between 0 and 9)")
    print(substitution_mat_modified)
    print("Compare first sequence ({0}) with others".format(ref_name))
    with open(fname, 'w') as outfile:
        for header, score_lst in compare_to_first_sequence(header_lst, seq_lst,
                                                           substitution_mat_modified):
            seq = "".join([str(s) for s in score_lst])
            fasta.write_fasta_entry(outfile, seq, header)
    print("wrote {0}".format(fname))
