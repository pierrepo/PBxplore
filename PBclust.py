#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cluster protein structures based on their PB sequences.

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# Modules
#===============================================================================
## Use print as a function for python 3 compatibility
from __future__ import print_function, division

## standard modules
import collections
import sys
import os
import subprocess
import argparse

## third-party module
import numpy

## local module
import PBlib as PB

## Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


def user_input():
    """
    Handle PBclust command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Cluster protein structures based on their PB sequences.")

    # mandatory arguments
    parser.add_argument("-f", action="append", required=True, 
        help="name(s) of the PBs file (in fasta format)")
    parser.add_argument("-o", action="store", required=True,
        help="name for results")
    parser.add_argument("-c", action="store", required=True, type=int,
        help="number of wanted clusters")  

    # optional arguments
    parser.add_argument("--compare", action="store_true", default=False,
        dest="compare", help="compare the first sequence versus all others")

    # get all parameters
    options = parser.parse_args()

    # test the validity of the arguments
    if options.c <= 0:
        parser.error("number of clusters must be > 0.")

    # check if the input files exist
    for name in options.f:
        if not os.path.isfile(name):
            sys.exit( "{0}: not a valid file. Bye".format(name) )

    return options


def display_clust_report(cluster_id):
    """
    Display a quick report on the clustering

    Display the number of structures in each cluster, and the fraction of the
    overall sequence set they represent.
    """
    nclusters = len(cluster_id)
    cluster_count = collections.Counter(cluster_id)
    for cluster, count in cluster_count.most_common():
        print("cluster {}: {} sequences ({:>2.0f}%)"
              .format(cluster, count, 100*count/nclusters))


def write_clusters(fname, cluster_id, medoid_id, seq_names):
    """
    Write the result of a clustering in a file

    The output file contains two types of lines:

    * first, lines that start with SEQ_CLU link each sequence header to a
      cluster ID; these lines are written n the same order as the input fasta
      file(s)
    * then, lines that start with MED_CLU link an input sequence to a cluster
      as its medoid; these lines are ordered as the cluster IDs so the first
      medoid is the medoid of the first cluster. The sequence index given in
      these lines start at 1.

    Parameters
    ----------
    fname : str
        The path to the file to write in
    cluster_id : list of int
        The cluster ID for each sequence ordered like the sequences
    medoid_id : list of int
        The index of the medoid for each group in the list of sequences
    seq_names: list of str
        The header for each sequence
    """
    with open(fname, "w") as outfile:
        for name, cluster in zip(seq_names, cluster_id):
            outfile.write('SEQ_CLU  "{}"  {} \n'.format(name, cluster))
        for idx, med in enumerate(medoid_id, start=1):
            outfile.write('MED_CLU  "{}"  {} \n'.format(seq_names[med], idx))


def write_distance_matrix(distance_matrix, fname):
    """
    Write a distance matrix in a file

    Parameters
    ----------
    distance_matrix : 2D numpy array
        The matrix to write
    fname : str
        The path to the file to write in
    """
    numpy.savetxt(fname, distance_matrix)


def compare(header_lst, seq_lst, substitution_mat, fname):
    """
    Wrap, for the command line, the comparison of all sequences with the first one

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
    substitution_mat_modified = PB.matrix_to_single_digit(substitution_mat)
    print("Normalized substitution matrix (between 0 and 9)")
    print(substitution_mat_modified)
    print("Compare first sequence ({0}) with others".format(ref_name))
    with open(fname, 'w') as outfile:
        for header, score_lst in PB.compare_to_first_sequence(header_lst, seq_lst,
                                                              substitution_mat_modified):
            seq = "".join([str(s) for s in score_lst])
            PB.write_fasta_entry(outfile, seq, header)
    print("wrote {0}".format(fname))


def pbclust_cli():
    """
    Run the PBclust command line
    """
    # Read user inputs
    options = user_input()
    header_lst, seq_lst = PB.read_several_fasta(options.f)

    # Load subtitution matrix
    substitution_mat = PB.load_substitution_matrix(PB.SUBSTITUTION_MATRIX_NAME)

    # --compare option
    # compare the first sequence (in the fasta file) versus all others
    if options.compare:
        compare_file_name = options.o + ".PB.compare.fasta"
        compare(header_lst, seq_lst, substitution_mat, compare_file_name)
        sys.exit(0)

    # Compute the distance matrix for the clustering
    distance_mat = PB.distance_matrix(seq_lst, substitution_mat)
    distance_fname = options.o + ".PB.dist"
    write_distance_matrix(distance_mat, distance_fname)
    print("wrote {0}".format(distance_fname))

    # Carry out the clustering
    try:
        cluster_id, medoid_id = PB.hclust(distance_mat, nclusters=options.c)
    except PB.RError as e:
        sys.exit('Error with R:\n' + str(e))
    display_clust_report(cluster_id)
    output_fname = options.o + ".PB.clust"
    write_clusters(output_fname, cluster_id, medoid_id, header_lst)
    print("wrote {0}".format(output_fname))


if __name__ == '__main__':
    pbclust_cli()
