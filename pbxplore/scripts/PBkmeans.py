#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read PDB structures and assign protein blocs (PBs).

2013 - P. Poulain, A. G. de Brevern
"""


# Use print as a function for python 3 compatibility
from __future__ import print_function

import argparse
import os
import sys

import pbxplore as pbx
from pbxplore.analysis import kmeans
from pbxplore.scripts import PBclust

try:
    import builtins
except ImportError:
    import itertools
    zip = itertools.izip
else:
    del builtins


def write_kmeans_clusters(fname, headers, groups):
    with open(fname, 'w') as infile:
        for header, group in zip(headers, groups):
            print('{0}\t{1}'.format(header, group), file=infile)


def user_input():
    """
    Handle PBkmeans command line arguments
    """
    parser = argparse.ArgumentParser(
        description=("Cluster protein structures based "
                     "on their PB sequences using the k-means "
                     "algorithm."))

    # mandatory arguments
    parser.add_argument("-f", action="append", required=True,
                        help="name(s) of the PBs file (in fasta format)")
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")
    parser.add_argument("--clusters", action="store", type=int, required=True,
                       help="number of wanted clusters")

    # options
    parser.add_argument("--max-iter", dest='max_iter', type=int, default=100,
                        help="maximum number of iterations")

    # get all arguments
    options = parser.parse_args()

    # test if the number of clusters is valid
    if options.clusters <= 0:
        parser.error("Number of clusters must be > 0.")

    # check if input files exist
    for name in options.f:
        if not os.path.isfile(name):
            parse.error("{0}: not a valid file. Bye".format(name))

    return options


def pbkmeans_cli():
    # Read user inputs
    options = user_input()
    header_lst, seq_lst = pbx.io.read_several_fasta(options.f)

    # Do the clustering
    groups, _, convergeance = kmeans.k_means(seq_lst,
                                             ngroups=options.clusters,
                                             max_iter=options.max_iter)

    # Write the output
    output_fname = options.o + ".PB.kmeans"
    write_kmeans_clusters(output_fname, header_lst, seq_lst)

    # Write the report
    PBclust.display_clust_report(groups)

    return 0

if __name__ =='__main__':
    sys.exit(pbkmeans_cli())
