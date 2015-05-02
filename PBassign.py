#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read PDB structures and assign protein blocs (PBs).

2013 - P. Poulain, A. G. de Brevern
"""


## Use print as a function for python 3 compatibility
from __future__ import print_function

## standard modules
import os
import sys
import glob
import argparse

## local module
import PBlib as PB
import PDBlib as PDB

# MDAnalysis is an optional requirement
try:
    import MDAnalysis
except ImportError:
    IS_MDANALYSIS = False
else:
    IS_MDANALYSIS = True


## Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


class _NullContext(object):
    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        pass


def PB_assign(pb_ref, structure, comment,
              fasta_file, flat_file=_NullContext(), phipsi_file=_NullContext()):
    """assign Protein Blocks (PB) from phi and psi angles
    """
    dihedrals = structure.get_phi_psi_angles()
    pb_seq = PB.assign(dihedrals, pb_ref)

    # write PBs in fasta file
    PB.write_fasta_entry(fasta_file, pb_seq, comment)
    # write phi and psi angles
    if not isinstance(phipsi_file, _NullContext):
        PB.write_phipsi_entry(phipsi_file, dihedrals, comment)
    # write PBs in flat file
    if not isinstance(flat_file, _NullContext):
        print(pb_seq, file=flat_file)
 
    print("PBs assigned for {0}".format(comment))


def user_inputs():
    parser = argparse.ArgumentParser(
        description='Read PDB structures and assign protein blocs (PBs).')

    # arguments
    parser.add_argument("-p", action="append",
                        help=("name of a pdb file "
                              "or name of a directory containing pdb files")
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")
    # arguments for MDanalysis
    group = parser.add_argument_group(
        title='other options [if MDanalysis module is available]')
    group.add_argument("-x", action="store",
                       help="name of xtc file (Gromacs)")
    group.add_argument("-g", action="store",
                       help="name of gro file (Gromacs)")
    # optional arguments
    parser.add_argument("--phipsi", action="store_true", default=False,
                        help="writes phi and psi angle")
    parser.add_argument("--flat", action="store_true", default=False,
                        help="writes one PBs sequence per line")
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
                for extension in (PDB.PDB_EXTENSIONS + PDB.PDBx_EXTENSIONS):
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
    options, pdb_name_lst = user_inputs()

    if options.p:
        if pdb_name_lst:
            print("{} PDB file(s) to process".format(len(pdb_name_lst)))
        else:
            print('Nothing to do. Good bye.')
            return

    # prepare fasta file for output
    fasta_name = options.o + ".PB.fasta"
    fasta_file = open(fasta_name, 'w')
    # prepare phi psi file for output
    phipsi_file = _NullContext()
    if options.phipsi:
        phipsi_name = options.o + ".PB.phipsi"
        phipsi_file = open(phipsi_name, 'w')
    # prepare flat file for output
    flat_file = _NullContext()
    if options.flat:
        flat_name = options.o + ".PB.flat"
        flat_file = open(flat_name, 'w')

    if options.p:
        # PB assignement of PDB structures
        chains = PDB.chains_from_files(pdb_name_lst)
    else:
        # PB assignement of a Gromacs trajectory
        chains = PDB.chains_from_trajectory(options.x, options.g)

    with fasta_file, phipsi_file, flat_file:
        for comment, chain in chains:
                PB_assign(PB.REFERENCES, chain, comment,
                          fasta_file=fasta_file, flat_file=flat_file,
                          phipsi_file=phipsi_file)

    print("wrote {0}".format(fasta_name))
    if options.flat:
        print("wrote {0}".format(flat_name))
    if options.phipsi:
        print("wrote {0}".format(phipsi_name))

if __name__ == '__main__':
    pbassign_cli()
