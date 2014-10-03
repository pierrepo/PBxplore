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
import math
import glob
import argparse 

## third-party module
import numpy 

## local module
import PBlib as PB

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

#===============================================================================
# Classes
#===============================================================================

class PdbAtom:
    """class for atoms in PDB format"""
    def __init__(self):
        """default constructor"""
        self.id = 0
        self.name = None
        self.resalt = ' '
        self.resname = None
        self.chain = None
        self.resid = 0
        self.resins = ''
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occupancy = 0.0
        self.tempfactor = 0.0
        self.element = '  '
        self.charge =  '  '
        
    def read(self, line):
        """read ATOM data from a PDB file line"""
        self.id = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resid = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip()) 

    def from_xtc(self, atm):
        """fill atom data from a .xtc mdanalysis selections"""
        self.id = atm.id
        self.name = atm.name
        self.resname = atm.resname
        self.chain = ""
        self.resid = atm.resid
        self.x = atm.pos[0]#.get_positions()[index][0]
        self.y = atm.pos[1]#.get_positions()[index][1]
        self.z = atm.pos[2]#.get_positions()[index][2]

    def __repr__(self):
        """representation for atom"""
        return 'atom %4d %4s in %4d %3s' \
        % (self.id, self.name, self.resid, self.resname)
        
    def format(self):
        """atom data formated as PDB"""
        return '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' \
        % ('ATOM  ', self.id, self.name, self.resalt,self.resname, self.chain, 
        self.resid, self.resins, self.x, self.y, self.z, 
        self.occupancy, self.tempfactor, self.element, self.charge)
        
    def coord(self):
        return [self.x, self.y, self.z]
    

class PdbStructure:
    """
    Class for a Protein Data Bank structure
    """
    def __init__(self):
        """default constructor for PDB structure"""
        self.chain = ""
        self.atoms = []
        self.comment = ""

    def add_atom(self, atom):
        """add atom to the structure"""
        # update chain name
        if not self.atoms:
            self.chain = atom.chain
        elif self.chain != atom.chain:
            print("WARNING: several chains in the same structure")
        # add atom to structure
        self.atoms.append(atom)

    def clean(self):
        """clean up chain, atoms and comment"""
        self.chain = ""
        self.atoms = []
        self.comment = ""
        
    def size(self):
        """get number of atoms from a structure"""
        return len(self.atoms)
               
    def get_all_dihedral(self):
        """compute phi and psi angles for a protein"""
        # extract backbone atoms
        backbone = {}
        for atom in self.atoms:
            if atom.name in ["CA", "C", "O", "N"]:
                resid = atom.resid
                if resid in backbone:
                    backbone[resid][atom.name] = atom
                else:
                    backbone[resid] = {atom.name: atom}
        
        # get dihedrals 
        phi_psi_angles = {}
        for res in sorted(backbone):
            # phi : C(i-1) - N(i) - CA(i) - C(i)
            try:
                phi = PB.get_dihedral(backbone[res-1]["C" ].coord(), 
                                   backbone[res  ]["N" ].coord(), 
                                   backbone[res  ]["CA"].coord(), 
                                   backbone[res  ]["C" ].coord())
            except:
                phi = None
            # psi : N(i) - CA(i) - C(i) - N(i+1)
            try:
                psi = PB.get_dihedral(backbone[res  ]["N" ].coord(), 
                                   backbone[res  ]["CA"].coord(), 
                                   backbone[res  ]["C" ].coord(), 
                                   backbone[res+1]["N" ].coord())
            except:
                psi = None
            #print res, phi, psi
            phi_psi_angles[res] = {"phi":phi, "psi":psi}
        return phi_psi_angles


#===============================================================================
# functions
#===============================================================================

def read_pb_definitions(pb_angles_string):
    """
    Read angle definitions of PBs
    """
    pb_angles = {}
    for line in pb_angles_string.split("\n"):
        if line and "#" not in line:
            items = line.split()
            pb_angles[items[0]] = numpy.array([float(items[i]) for i in range(1, len(items))])
    print("read PB definitions: %d PBs x %d angles " \
          % (len(pb_angles), len(pb_angles["a"])))
    return pb_angles

#-------------------------------------------------------------------------------
def angle_modulo_360(angle):
    """keep angle in the range -180 / +180 [degrees]
    """
    if angle > 180.0:
        return angle - 360.0
    elif angle < -180.0:
        return angle + 360.0
    else:
        return angle
    
#-------------------------------------------------------------------------------
def write_phipsi(name, torsion, com):
    """save phi and psi angles
    """
    f_out = open(name, "a")
    for res in sorted(torsion):
        try:
            phi_str = "%8.2f" % torsion[res]["phi"]
        except:
            phi_str = "    None"
        try:
            psi_str = "%8.2f" % torsion[res]["psi"]
        except:
            psi_str = "    None"
        f_out.write("%s %6d %s %s \n" % (com, res, phi_str, psi_str))
    f_out.close()

#-------------------------------------------------------------------------------
def write_flat(name, seq):
    """write flat sequence to file 
    """
    f_out = open(name, "a")
    f_out.write(seq + "\n")
    f_out.close()

#-------------------------------------------------------------------------------
def PB_assign(pb_ref, structure, comment):
    """assign Protein Blocks (PB) from phi and psi angles
    """
    # get phi and psi angles from structure
    dihedrals = structure.get_all_dihedral()
    # write phi and psi angles
    if options.phipsi:
        write_phipsi(phipsi_name, dihedrals, comment)

    pb_seq = ""
    # iterate over all residues
    for res in sorted(dihedrals):
        angles = []
        # try to get all eight angles required for PB assignement
        try:
            angles.append(dihedrals[res-2]["psi"])
            angles.append(dihedrals[res-1]["phi"])
            angles.append(dihedrals[res-1]["psi"])
            angles.append(dihedrals[res  ]["phi"])
            angles.append(dihedrals[res  ]["psi"])
            angles.append(dihedrals[res+1]["phi"])
            angles.append(dihedrals[res+1]["psi"])
            angles.append(dihedrals[res+2]["phi"])
            # check for bad angles 
            # (error while calculating torsion: missing atoms)
            if None in angles:
                pb_seq += "Z"
                continue 
           
        # cannot get required angles (Nter, Cter or missign residues)
        # -> cannot assign PB
        # jump to next residue
        except:
            pb_seq += "Z"
            continue
        
        # convert to array
        angles = numpy.array(angles)

        # compare to reference PB angles
        rmsda_lst = {}
        for block in pb_ref:
            diff = pb_ref[block] - angles
            diff2 = angle_modulo_360_vect(diff)
            rmsda = numpy.sum(diff2**2)
            rmsda_lst[rmsda] = block
        pb_seq += rmsda_lst[min(rmsda_lst)]

    # write PBs in fasta file
    PB.write_fasta(fasta_name, pb_seq, comment)
    
    # write PBs in flat file
    if options.flat:
        write_flat(flat_name, pb_seq)
 
    print("PBs assigned for", comment)
             
#-------------------------------------------------------------------------------
# vertorize function
#-------------------------------------------------------------------------------
angle_modulo_360_vect = numpy.vectorize(angle_modulo_360)

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
        if os.path.isfile(name):
            pdb_name_lst.append(name)
        elif os.path.isdir(name):
            pdb_name_lst += glob.glob(name + os.pathsep + "*.pdb")
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
# read PB definitions
#-------------------------------------------------------------------------------
pb_def = read_pb_definitions(PB.DEFINITIONS)

#-------------------------------------------------------------------------------
# prepare fasta file for output
#-------------------------------------------------------------------------------
fasta_name = options.o + ".PB.fasta"
PB.clean_file(fasta_name)

#-------------------------------------------------------------------------------
# prepare phi psi file for output
#-------------------------------------------------------------------------------
if options.phipsi:
    phipsi_name = options.o + ".PB.phipsi"
    PB.clean_file(phipsi_name)
 
#-------------------------------------------------------------------------------
# prepare flat file for output
#-------------------------------------------------------------------------------
if options.flat:
    flat_name = options.o + ".PB.flat"
    PB.clean_file(flat_name)

#-------------------------------------------------------------------------------
# read PDB files
#-------------------------------------------------------------------------------
structure = PdbStructure()

# PB assignement of PDB structures
if options.p:

    model = ""
    chain = " "
    comment = ""

    for pdb_name in pdb_name_lst:
        print(pdb_name)
        f_in = open(pdb_name, 'r')
        for line in f_in:
            flag = line[0:6].strip()
            if flag == "MODEL":
                model = line.split()[1]
            if flag == "ATOM":
                atom = PdbAtom()
                atom.read(line)
                # assign structure upon new chain
                if structure.size() != 0 and structure.chain != atom.chain:
                    PB_assign(pb_def, structure, comment)
                    structure.clean()
                # append structure with atom
                structure.add_atom(atom)
                # define structure comment
                # when the structure contains 1 atom
                if structure.size() == 1:
                    comment = pdb_name 
                    if model:
                        comment += " | model %s" % (model)
                        model = ""
                    if atom.chain:
                        comment += " | chain %s" % (atom.chain)
                        atom.chain = ""
            # assign structure after end of model (or chain)
            if structure.size() != 0 and flag in ["TER", "ENDMDL"]:
                PB_assign(pb_def, structure, comment)
                structure.clean()
        # assign last structure
        if structure.size() != 0:
            PB_seq = PB_assign(pb_def, structure, comment)

        f_in.close()   

# PB assignement of a Gromacs trajectory
if not options.p:
    try:
        import MDAnalysis
    except:
        sys.exit("Error: failed to import MDAnalysis")

    model = ""
    chain = ""
    comment = ""

    conf = options.g
    traj = options.x

    universe = MDAnalysis.Universe(conf, traj)

    for ts in universe.trajectory:
        selection = universe.selectAtoms("backbone")
        for atm in selection:
            atom = PdbAtom()        
            atom.from_xtc(atm)
            # append structure with atom
            structure.add_atom(atom)
            # define structure comment
            # when the structure contains 1 atom
            if structure.size() == 1:
                comment = "%s | frame %s" % (options.x, ts.frame)
        # assign structure after end of frame
        if structure.size() != 0 :
            PB_assign(pb_def, structure, comment)
            structure.clean()

print( "wrote {0}".format(fasta_name) )
if options.flat:
    print( "wrote {0}".format(flat_name) )
if options.phipsi:
    print( "wrote {0}".format(phipsi_name) )

