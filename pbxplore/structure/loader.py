#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function
import sys

# Local module
from .structure import Chain, Atom
from .PDB import PDB


# load MDAnalysis with limited support for Python 3
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis


# Create the __all__ keyword according to the conditional import
__all__ = ['chains_from_files', 'chains_from_trajectory']


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

        print("Read {0} chain(s) in {1}".format(pdb.nb_chains, pdb_name), file=sys.stderr)


def chains_from_trajectory(trajectory, topology):
    universe = MDAnalysis.Universe(topology, trajectory)
    selection = universe.select_atoms("backbone")

    #Initialize structure with the selection
    structure = Chain()
    for atm in selection:
        atom = Atom.read_from_xtc(atm)
        # append structure with atom
        structure.add_atom(atom)

    nb_frames = len(universe.trajectory)

    # Print the first frame
    print("Frame {}/{}.".format(1, nb_frames), file=sys.stderr)

    for ts in universe.trajectory:
        #Update only with new coordinates
        structure.set_coordinates(selection.positions)

        # define structure comment
        comment = "%s | frame %s" % (trajectory, ts.frame)
        yield comment, structure

        # Progress bar
        # Print one frame every 100.
        if ((ts.frame + 1) % 100 == 0):
            print("Frame {}/{}.".format(ts.frame + 1, nb_frames), file=sys.stderr)

    # Print the last frame
    print("Frame {}/{}.".format(nb_frames, nb_frames), file=sys.stderr)
