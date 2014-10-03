#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python library to handle Protein Data Bank files: PDB and PDBx/mmCIF

2014 - P. Poulain
"""
#===============================================================================
# Modules
#===============================================================================
## standard modules
import os
import sys
import math

## third-party modules
import numpy

#===============================================================================
# Classes
#===============================================================================

class Atom:
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
    

class Structure:
    """
    Class for a Protein Data Bank structure (chain)
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
            print "WARNING: several chains in the same structure"
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
        for res in sorted(backbone.iterkeys()):
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


class PDBFile:
    def __init__(self):
        """
        Default constructor for PDB file
        """
        self.filename = ""

class PDBxFile:
    def __init__(self):
        """
        Default constructor for PDBx/mmCIF file
        """
        self.filename = ""

