#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
import os
import gzip

# Local module
from .structure import Atom, Chain

# =============================================================================
# Data
# =============================================================================
# file extensions for PDB and PDBx/mmCIF files
PDB_EXTENSIONS = ('.pdb', '.PDB', '.ent', '.ENT4')
PDBx_EXTENSIONS = ('.cif', '.CIF', '.cif.gz', '.CIF.GZ')


class PDB:
    """
    Class to read PDB files.

    """
    def __init__(self, name):
        """
        Default constructor for PDB file.
        """
        self.filename = name
        self.chains = []
        # check that file exists
        if not os.path.isfile(self.filename):
            raise IOError("Cannot read {}: does not exist or is not a file."
                          .format(self.filename))
        if self.filename.endswith(PDB_EXTENSIONS):
            self.__read_PDB()
        if self.filename.endswith(PDBx_EXTENSIONS):
            self.__read_PDBx()

    def __read_PDB(self):
        """
        Read PDB file.
        """
        # create new chain
        chain = Chain()
        # get chains from file
        # A PDB file can have several models
        # that can have several chains themselves.
        if self.filename.endswith(('.gz', '.GZ')):
            # for compressed file
            f_in = gzip.open(self.filename, 'rt')
        else:
            f_in = open(self.filename, 'rt')
        for line in f_in:
            flag = line[0:6].strip()
            if flag == "MODEL":
                chain.set_model(line.split()[1])
            if flag == "ATOM":
                atom = Atom.read_from_PDB(line)
                # store current chain and clean object
                if chain.size() != 0 and chain.name != atom.chain:
                    self.chains.append(chain)
                    chain = Chain()
                # append structure with atom
                chain.add_atom(atom)
            # store chain after end of model or chain
            if chain.size() != 0 and flag in ["TER", "ENDMDL"]:
                self.chains.append(chain)
                chain = Chain()
        # store last chain
        if chain.size() != 0:
            self.chains.append(chain)
        f_in.close()
        print("Read {0} chain(s) in {1}"
              .format(len(self.chains), self.filename))

    def __read_PDBx(self):
        """
        Read PDBx/mmCIF file
        """
        # create new chain
        chain = Chain()
        # get chains from file
        # A PDBx file can have several models
        # that can have several chains themselves.
        atom_fields = []
        atom_coordinates = []
        if self.filename.endswith(('.gz', '.GZ')):
            # for compressed file
            f_in = gzip.open(self.filename, 'rt')
        else:
            f_in = open(self.filename, 'rt')
        for line in f_in:
            item = line.strip()
            # then store atom field definitions
            if item.startswith("_atom_site."):
                atom_fields.append(item.replace("_atom_site.", ""))
            # then store atom coordinates
            if atom_fields and item.startswith('ATOM'):
                atom_coordinates.append(item)
        f_in.close()
        # separate all chains and store atoms
        chain = Chain()
        for atom_line in atom_coordinates:
            atom = Atom.read_from_PDBx(atom_line, atom_fields)
            # define model at first atom
            if chain.size() == 1:
                chain.set_model(atom.model)
            # store current chain when chain name changed
            if chain.size() != 0 and chain.name != atom.chain:
                # store model number only if there is more than one model
                if chain.model == atom.model:
                    chain.set_model("")
                print(chain)
                self.chains.append(chain)
                chain = Chain()
            # store current chain when model number changed
            if chain.size() != 0 and chain.model != atom.model:
                self.chains.append(chain)
                chain = Chain()
            # append structure with atom
            chain.add_atom(atom)
        # store last chain
        if chain.size() != 0:
            self.chains.append(chain)

    def get_chains(self):
        """
        Give chains, one at a time.

        Returns
        -------
        generator
            Chains in PDB structure.
        """
        for chain in self.chains:
            yield chain
