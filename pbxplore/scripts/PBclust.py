#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cluster protein structures based on their PB sequences.

2013 - P. Poulain, A. G. de Brevern
"""


# Use print as a function for python 3 compatibility
from __future__ import print_function, division

# Standard modules
import collections
import sys
import os
import argparse

# Third-party module
import numpy

# Local module
import pbxplore as pbx
from pbxplore.PB import SUBSTITUTION_MATRIX_NAME

# Python2/Python3 compatibility
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

    # --clusters or --compare arguments
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--clusters", action="store", type=int,
                       help="number of wanted clusters")
    # optional arguments
    group.add_argument("--compare", action="store_true", default=False,
                       help="compare the first sequence versus all others")

    # get all arguments
    options = parser.parse_args()

    # test if the number of clusters is valid
    if options.clusters is not None and options.clusters <= 0:
            parser.error("Number of clusters must be > 0.")

    # check if input files exist
    for name in options.f:
        if not os.path.isfile(name):
            sys.exit("{0}: not a valid file. Bye".format(name))

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


def pbclust_cli():
    """
    Run the PBclust command line
    """
    # Read user inputs
    options = user_input()
    header_lst, seq_lst = pbx.io.read_several_fasta(options.f)

    # Load subtitution matrix
    try:
        substitution_mat = pbx.PB.load_substitution_matrix(SUBSTITUTION_MATRIX_NAME)
    except ValueError:
        sys.exit("Substitution matrix is not symetric.")
    except IOError:
        sys.exit("Error reading the substitution matrix.")

    # --compare option
    # compare the first sequence (in the fasta file) versus all others
    if options.compare:
        compare_file_name = options.o + ".PB.compare.fasta"
        pbx.analysis.compare(header_lst, seq_lst, substitution_mat, compare_file_name)
        sys.exit(0)

    # Compute the distance matrix for the clustering
    try:
        distance_mat = pbx.analysis.distance_matrix(seq_lst, substitution_mat)
    except pbx.PB.InvalidBlockError as e:
        sys.exit('Unexpected PB in the input ({})'.format(e.block))
    distance_fname = options.o + ".PB.dist"
    write_distance_matrix(distance_mat, distance_fname)
    print("wrote {0}".format(distance_fname))

    # Carry out the clustering
    try:
        cluster_id, medoid_id = pbx.analysis.hclust(distance_mat, nclusters=options.clusters)
    except pbx.analysis.RError as e:
        sys.exit('Error with R:\n' + str(e))
    display_clust_report(cluster_id)
    output_fname = options.o + ".PB.clust"
    write_clusters(output_fname, cluster_id, medoid_id, header_lst)
    print("wrote {0}".format(output_fname))


if __name__ == '__main__':
    pbclust_cli()
