#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

# Local module
from .structure import Chain, Atom
from .PDB import PDB

# Conditional import
try:
    import MDAnalysis
except ImportError:
    IS_MDANALYSIS = False
else:
    IS_MDANALYSIS = True


# Create the __all__ keyword according to the conditional import
__all__ = ['chains_from_files']
if IS_MDANALYSIS:
    __all__ += ['chains_from_trajectory']


def chains_from_files(path_list):
    for pdb_name in path_list:
        pdb = PDB(pdb_name)
        for chain in pdb.get_chains():
            # build comment
            comment = pdb_name
            if chain.model:
                comment += " | model %s" % (chain.model)
            if chain.name:
                comment += " | chain %s" % (chain.name)
            yield comment, chain


def chains_from_trajectory(trajectory, topology):
    universe = MDAnalysis.Universe(topology, trajectory)
    selection = universe.select_atoms("backbone")

    #Initialize structure with the selection
    structure = Chain()
    for atm in selection:
        atom = Atom.read_from_xtc(atm)
        # append structure with atom
        structure.add_atom(atom)

    for ts in universe.trajectory:
        #Update only with new coordinates
        structure.set_coordinates(selection.positions)

        # define structure comment
        comment = "%s | frame %s" % (trajectory, ts.frame)
        yield comment, structure
