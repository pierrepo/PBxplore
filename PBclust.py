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
from __future__ import print_function

## standard modules
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
    #-------------------------------------------------------------------------------
    # get arguments
    #-------------------------------------------------------------------------------
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

    #-------------------------------------------------------------------------------
    # check options
    #-------------------------------------------------------------------------------
    if options.c <= 0:
        parser.error("number of clusters must be > 0.")

    #-------------------------------------------------------------------------------
    # check input files
    #-------------------------------------------------------------------------------
    for name in options.f:
        if not os.path.isfile(name):
            sys.exit( "{0}: not a valid file. Bye".format(name) )

    return options


def display_clust_report(cluster_id):
    # count number of sequences in clusters
    cluster_count = {}
    for idx in cluster_id:
        cluster_count[idx] = cluster_count.get(idx, 0) + 1
    for idx in sorted(cluster_count):
        print("cluster {}: {} sequences ({:>2.0f}%)"
              .format(idx, cluster_count[idx],
                      1.0*cluster_count[idx]/len(cluster_id)*100))


def write_clusters(fname, seq_id, cluster_id, medoid_id, seq_names):
    with open(fname, "w") as outfile:
        for seq, cluster in zip(seq_id, cluster_id):
            outfile.write('SEQ_CLU  "{}"  {} \n'.format(seq_names[seq], cluster))
        for idx, med in enumerate(medoid_id):
            outfile.write('MED_CLU  "{}"  {} \n'.format(seq_names[med], idx+1))


def write_distance_matrix(distance_matrix, fname):
    numpy.savetxt(fname, distance_matrix)


def compare(header_lst, seq_lst, substitution_mat, fname):
    PB.clean_file(fname)
    ref_name = header_lst[0]
    ref_seq = seq_lst[0]
    substitution_mat_modified = PB.matrix_to_single_digit(substitution_mat)
    print("Normalized substitution matrix (between 0 and 9)")
    print(substitution_mat_modified)
    print("Compare first sequence ({0}) with others".format(ref_name))
    for header, score_lst in PB.compare_to_first_sequence(header_lst, seq_lst,
                                                          substitution_mat_modified):
        seq = "".join([str(s) for s in score_lst])
        PB.write_fasta(fname, seq, header)
    print("wrote {0}".format(fname))


def pbclust_cli():
    options = user_input()
    #-------------------------------------------------------------------------------
    # read PBs files
    #-------------------------------------------------------------------------------
    header_lst, seq_lst = PB.read_several_fasta(options.f)
    pb_seq = numpy.array(list(zip(header_lst, seq_lst)))

    #-------------------------------------------------------------------------------
    # load subtitution matrix
    #-------------------------------------------------------------------------------
    substitution_mat = PB.load_substitution_matrix(PB.SUBSTITUTION_MATRIX_NAME)

    #-------------------------------------------------------------------------------
    # --compare option
    # compare the first sequence (in the fasta file) versus all others
    #-------------------------------------------------------------------------------
    if options.compare:
        compare_file_name = options.o + ".PB.compare.fasta"
        compare(header_lst, seq_lst, substitution_mat, compare_file_name)
        sys.exit(0)

    similarity_mat = PB.distance_matrix(pb_seq, substitution_mat, PB.substitution_score)
    distance_mat = PB.similarity_to_distance(similarity_mat)
    seq_id, cluster_id, medoid_id = PB.hclust(distance_mat, nclusters=options.c)
    distance_fname = options.o + ".PB.dist"
    write_distance_matrix(distance_mat, distance_fname)
    print("wrote {0}".format(distance_fname))
    display_clust_report(cluster_id)

    output_fname = options.o + ".PB.clust"
    write_clusters(output_fname, seq_id, cluster_id, medoid_id, header_lst)
    print("wrote {0}".format(output_fname))


if __name__ == '__main__':
    pbclust_cli()
