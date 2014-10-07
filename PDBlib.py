#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python library to handle Protein Data Bank files: PDB and PDBx/mmCIF

2014 - P. Poulain
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
import gzip

## third-party modules
import numpy

#===============================================================================
# Functions
#===============================================================================
def get_dihedral(atomA, atomB, atomC, atomD):
    """
    Compute dihedral angle between 4 atoms (A, B, C, D)
    each atom is represented as a list of three coordinates [x, y, z]
    output is in degree in the range -180 / +180
    
    Note: this function is not part of any class to ease its reusability
    """
    
    # convert lists to Numpy objects
    A = numpy.array(atomA)
    B = numpy.array(atomB)
    C = numpy.array(atomC)
    D = numpy.array(atomD)
 
    # vectors
    AB = B - A 
    BC = C - B 
    CD = D - C 

    # normal vectors
    n1 = numpy.cross(AB, BC)
    n2 = numpy.cross(BC, CD)

    # normalize normal vectors
    n1 /= numpy.linalg.norm(n1)
    n2 /= numpy.linalg.norm(n2)
    
    # angle between normals
    cosine = numpy.sum(n1*n2) / (numpy.linalg.norm(n1) * numpy.linalg.norm(n2))
    try :
        torsion = math.acos(cosine)
    except:
        cosine = int(cosine) #+0.0001
        torsion = math.acos(cosine)

    # convert radion to degree
    torsion = torsion * 180.0 / math.pi 

    # find if the torsion is clockwise or counterclockwise
    #if numpy.sum(n1 * CD) < 0.0:
    if numpy.dot(n1, CD) < 0.0:
        torsion = 360 - torsion
    if torsion == 360.0:
        torsion = 0.0
    
    # get range -180 / +180
    if torsion > 180.0:
        torsion = torsion - 360
    if torsion < -180.0:
        torsion = torsion + 360
   
    return torsion
    
#===============================================================================
# Classes
#===============================================================================

class Atom:
    """
    Class for atoms in PDB format
    """
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
        self.model = None
        
    def read_from_PDB(self, line):
        """read ATOM data from a PDB file line"""
        self.id = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resid = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip()) 

    def read_from_PDBx(self, line, fields):
        """
        Read ATOM data from a PDBx/mmCIF file line
        for more details regarding this format see:
        http://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html
        """
        try:
            dic = dict(zip( fields, line.split() ))
        except:
            print("Something went wrong in reading\n{0}".format(line))
            raise
        try:
            self.id = int( dic['id'] )
            self.name = dic['label_atom_id']
            self.resname = dic['label_comp_id']
            self.chain = dic['label_asym_id']
            self.resid = int( dic['label_seq_id'] )
            self.x = float( dic['Cartn_x'] )
            self.y = float( dic['Cartn_y'] )
            self.z = float( dic['Cartn_z'] ) 
            self.model = dic['pdbx_PDB_model_num']
        except:
            print("Something went wrong in data convertion\n{0}".format(dic))
            raise
        
    def read_from_xtc(self, atm):
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
    

class Chain:
    """
    Class for a Protein Data Bank structure./chain
    """
    def __init__(self):
        """default constructor for PDB structure"""
        self.name = ""
        self.model = ""
        self.atoms = []

    def __repr__(self):
        """Representation of chain"""
        return "Chain {0} / model {1}: {2} atoms".format(self.name, 
                                                         self.model, 
                                                         len(self.atoms))
        
    def add_atom(self, atom):
        """add atom to the structure"""
        # update chain name when first atom is stored
        if not self.atoms:
            self.name = atom.chain
        # check that chain name is always the same
        elif self.name != atom.chain:
            print("WARNING: several chains in the same structure")
        # add atom to structure
        self.atoms.append(atom)
    
    def set_model(self, model):
        """Set model number"""
        self.model = model
       
    def size(self):
        """get number of atoms from a structure"""
        return len(self.atoms)
               
    def get_phi_psi_angles(self):
        """
        Compute phi and psi angles for a protein
        """
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
            # phi: angle between C(i-1) - N(i) - CA(i) - C(i)
            try:
                phi = get_dihedral(backbone[res-1]["C" ].coord(), 
                                   backbone[res  ]["N" ].coord(), 
                                   backbone[res  ]["CA"].coord(), 
                                   backbone[res  ]["C" ].coord())
            except:
                phi = None
            # psi: angle between N(i) - CA(i) - C(i) - N(i+1)
            try:
                psi = get_dihedral(backbone[res  ]["N" ].coord(), 
                                   backbone[res  ]["CA"].coord(), 
                                   backbone[res  ]["C" ].coord(), 
                                   backbone[res+1]["N" ].coord())
            except:
                psi = None
            #print(res, phi, psi)
            phi_psi_angles[res] = {"phi":phi, "psi":psi}
        return phi_psi_angles


class PDBFile:
    """
    Class to read structures from PDB files
    """
    def __init__(self, name):
        """
        Default constructor for PDB file
        """
        self.filename = name
        self.chains = []

    def read_chains_from_PDB(self):
        """
        Read PDB file
        """
        # check that file exists
        if not os.path.isfile(self.filename):
            sys.exit("Cannot read {}: does not exist".format(self.filename))
        # create new chain
        chain = Chain()
        # get chains from file
        # A PDB file can have several models 
        # that can have several chains themselves.
        f_in = open(self.filename, 'r')
        for line in f_in:
            flag = line[0:6].strip()
            if flag == "MODEL":
                chain.set_model( line.split()[1] )
            if flag == "ATOM":
                atom = Atom()
                atom.read_from_PDB(line)
                # store current chain and clean object
                if chain.size() != 0 and chain.name != atom.chain:
                    self.chains.append( chain )
                    chain = Chain()
                # append structure with atom
                chain.add_atom(atom)
            # store chain after end of model or chain
            if chain.size() != 0 and flag in ["TER", "ENDMDL"]:
                self.chains.append( chain )
                chain = Chain()
        # store last chain
        if chain.size() != 0:
            self.chains.append( chain )
        f_in.close()
        print("Read a total of {0} chain(s) in {1}"
              .format(len(self.chains), self.filename))

    
    def read_chains_from_PDBx(self):
        """
        Read chains from PDBx/mmCIF files
        """
        # check that file exists
        if not os.path.isfile(self.filename):
            sys.exit("Cannot read {}: does not exist".format(self.filename))
        # create new chain
        chain = Chain()
        # get chains from file
        # A PDB file can have several models 
        # that can have several chains themselves.
        atom_fields = []
        atom_coordinates = []
        if self.filename.endswith( ('.gz', '.GZ') ):
            f_in = gzip.open(self.filename, 'r')
        else:
            f_in = open(self.filename, 'r')
        for line in f_in:
            item = line.strip()
            # then store atom field definitions
            if item.startswith("_atom_site."):
                atom_fields.append( item.replace("_atom_site.", "") )
            # then store atom coordinates
            if atom_fields and item.startswith('ATOM'):
                atom_coordinates.append( item )
        f_in.close()
        # separate all chains and store atoms
        chain = Chain()
        for atom_line in atom_coordinates:
            atom = Atom()
            atom.read_from_PDBx(atom_line, atom_fields)
            # define model at first atom
            if chain.size() == 1:
                chain.set_model(atom.model)
            # store current chain when chain name changed
            if chain.size() != 0 and chain.name != atom.chain:
                # store model number only if there is more than one model
                if chain.model == atom.model:
                    chain.set_model("")
                print(chain)
                self.chains.append( chain )
                chain = Chain()
            # store current chain when model number changed
            if chain.size() != 0 and chain.model != atom.model:
                self.chains.append( chain )
                chain = Chain()            
            # append structure with atom
            chain.add_atom(atom)
        # store last chain
        if chain.size() != 0:
            self.chains.append( chain )
            

    def get_chains(self):
        """
        Give chains, one at a time
        """
        for chain in self.chains:
            yield chain
            
