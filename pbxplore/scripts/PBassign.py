#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read PDB structures and assign protein blocs (PBs).

2013 - P. Poulain, A. G. de Brevern
"""


# Use print as a function for python 3 compatibility
from __future__ import print_function

# Standard modules
import os
import sys
import glob
import argparse

# Local modules
import pbxplore as pbx

# MDAnalysis is an optional requirement
try:
    import MDAnalysis
except ImportError:
    IS_MDANALYSIS = False
else:
    IS_MDANALYSIS = True


# Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


def user_inputs():
    """
    Handle the user parameter from the command line

    Parse the command line parameter and build the list of input files.
    """
    parser = argparse.ArgumentParser(
        description='Read PDB structures and assign protein blocs (PBs).')

    # arguments
    parser.add_argument("-p", action="append",
                        help=("name of a pdb file "
                              "or name of a directory containing pdb files"))
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")
    # arguments for MDanalysis
    group = parser.add_argument_group(
        title='other options [if MDanalysis module is available]')
    group.add_argument("-x", action="store",
                       help="name of xtc file (Gromacs)")
    group.add_argument("-g", action="store",
                       help="name of gro file (Gromacs)")

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')
    # get all arguments
    options = parser.parse_args()

    # check options
    if not options.p:
        if not IS_MDANALYSIS:
            parser.error("MDAnalysis is not installed; options -x and -d are not available")
        if not options.x:
            parser.print_help()
            parser.error("use at least option -p or -x")
        elif not options.g:
            parser.print_help()
            parser.error("option -g is mandatory, with use of option -x")

    # check files
    pdb_name_lst = []
    if options.p:
        for name in options.p:
            # input is a file: store file name
            if os.path.isfile(name):
                pdb_name_lst.append(name)
            # input is a directory: list and store all PDB and PDBx/mmCIF files
            elif os.path.isdir(name):
                for extension in (pbx.structure.PDB_EXTENSIONS + pbx.structure.PDB.PDBx_EXTENSIONS):
                    pdb_name_lst += glob.glob(os.path.join(name, "*" + extension))
            # input is not a file neither a directory: say it
            elif (not os.path.isfile(name) or not os.path.isdir(name)):
                parser.error("{0}: not a valid file or directory".format(name))
    else:
        if not os.path.isfile(options.x):
            sys.exit("{0}: not a valid file".format(options.x))
        elif not os.path.isfile(options.g):
            sys.exit("{0}: not a valid file".format(options.g))

    return options, pdb_name_lst


def pbassign_cli():
    """
    PBassign command line.
    """
    options, pdb_name_lst = user_inputs()

    if options.p:
        if pdb_name_lst:
            print("{} PDB file(s) to process".format(len(pdb_name_lst)))
        else:
            print('Nothing to do. Good bye.')
            return
        # PB assignement of PDB structures
        chains = pbx.chains_from_files(pdb_name_lst)
    else:
        # PB assignement of a Gromacs trajectory
        chains = pbx.chains_from_trajectory(options.x, options.g)

    all_comments = []
    all_sequences = []
    for comment, chain in chains:
        dihedrals = chain.get_phi_psi_angles()
        sequence = pbx.assign(dihedrals)
        all_comments.append(comment)
        all_sequences.append(sequence)

    fasta_name = options.o + ".PB.fasta"
    with open(fasta_name, 'w') as outfile:
        pbx.io.write_fasta(outfile, all_sequences, all_comments)

    print("wrote {0}".format(fasta_name))


if __name__ == '__main__':
    pbassign_cli()
