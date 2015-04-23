#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read PDB structures and assign protein blocs (PBs).

2013 - P. Poulain, A. G. de Brevern 
"""


#===============================================================================
# Modules
#===============================================================================
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


#===============================================================================
# Python2/Python3 compatibility
#===============================================================================

# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass


def PB_assign(pb_ref, structure, comment,
              fasta_name, flat_name=None, phipsi_name=None):
    """assign Protein Blocks (PB) from phi and psi angles
    """
    # get phi and psi angles from structure
    dihedrals = structure.get_phi_psi_angles()
    #print(dihedrals)
    # write phi and psi angles
    if not phipsi_name is None:
        PB.write_phipsi(phipsi_name, dihedrals, comment)

    pb_seq = PB.assign(dihedrals, pb_ref)
    
    # write PBs in fasta file
    PB.write_fasta(fasta_name, pb_seq, comment)
    
    # write PBs in flat file
    if not flat_name is None:
        PB.write_flat(flat_name, pb_seq)
 
    print("PBs assigned for {0}".format(comment))


#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# manage parameters
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description = 'Read PDB structures and assign protein blocs (PBs).')

# arguments
parser.add_argument("-p", action="append",
    help="name of a pdb file or name of a directory containing pdb files")
parser.add_argument("-o", action="store", required = True,
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
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

# get all arguments
options = parser.parse_args()

# check options
if not options.p:
    if not options.x:
        parser.print_help()
        parser.error("use at least option -p or -x")
    elif not options.g:
        parser.print_help()
        parser.error("option -g is mandatory, with use of option -x")


#-------------------------------------------------------------------------------
# check files
#-------------------------------------------------------------------------------
pdb_name_lst = []
if options.p:
    for name in options.p:
        # input is a file: store file name
        if os.path.isfile(name):
            pdb_name_lst.append(name)
        # input is a directory: list and store all PDB and PDBx/mmCIF files
        elif os.path.isdir(name):
            for extension in (PDB.PDB_EXTENSIONS + PDB.PDBx_EXTENSIONS):
                pdb_name_lst += glob.glob(os.path.join(name,  "*" + extension))
        # input is not a file neither a directory: say it
        elif (not os.path.isfile(name) or not os.path.isdir(name)):
            print("{0}: not a valid file or directory".format(name))
    
    if pdb_name_lst:
        print("{0} PDB file(s) to process".format( len(pdb_name_lst) ))
    else:
        sys.exit("Nothing to do. Bye.")
else:   
    if not os.path.isfile(options.x):
        sys.exit("{0}: not a valid file".format(options.x))
    elif not os.path.isfile(options.g):
        sys.exit("{0}: not a valid file".format(options.g))


#-------------------------------------------------------------------------------
# prepare fasta file for output
#-------------------------------------------------------------------------------
fasta_name = options.o + ".PB.fasta"
PB.clean_file(fasta_name)

#-------------------------------------------------------------------------------
# prepare phi psi file for output
#-------------------------------------------------------------------------------
phipsi_name = None
if options.phipsi:
    phipsi_name = options.o + ".PB.phipsi"
    PB.clean_file(phipsi_name)
 
#-------------------------------------------------------------------------------
# prepare flat file for output
#-------------------------------------------------------------------------------
flat_name = None
if options.flat:
    flat_name = options.o + ".PB.flat"
    PB.clean_file(flat_name)

#-------------------------------------------------------------------------------
# read PDB files
#-------------------------------------------------------------------------------


# PB assignement of PDB structures
if options.p:
    for pdb_name in pdb_name_lst:
        pdb = PDB.PDB(pdb_name)
        for chain in pdb.get_chains():
            # build comment 
            comment = pdb_name 
            if chain.model:
                comment += " | model %s" % (chain.model)
            if chain.name:
                comment += " | chain %s" % (chain.name)
            # assign PBs
            PB_assign(PB.REFERENCES, chain, comment,
                      fasta_name=fasta_name, flat_name=flat_name,
                      phipsi_name=phipsi_name)



# PB assignement of a Gromacs trajectory
if not options.p:
    try:
        import MDAnalysis
    except ImportError:
        sys.exit("Error: failed to import MDAnalysis")

    model = ""
    chain = ""
    comment = ""

    conf = options.g
    traj = options.x

    universe = MDAnalysis.Universe(conf, traj)

    for ts in universe.trajectory:
        structure = PDB.Chain()
        selection = universe.selectAtoms("backbone")
        for atm in selection:
            atom = PDB.Atom()
            atom.read_from_xtc(atm)
            # append structure with atom
            structure.add_atom(atom)
            # define structure comment
            # when the structure contains 1 atom
            if structure.size() == 1:
                comment = "%s | frame %s" % (options.x, ts.frame)
        # assign structure after end of frame
        if structure.size() != 0 :
            PB_assign(PB.REFERENCES, structure, comment,
                      fasta_name=fasta_name, flat_name=flat_name,
                      phipsi_name=phipsi_name)

print( "wrote {0}".format(fasta_name) )
if options.flat:
    print( "wrote {0}".format(flat_name) )
if options.phipsi:
    print( "wrote {0}".format(phipsi_name) )

