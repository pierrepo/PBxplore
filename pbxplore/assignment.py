#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Protrein Block assignation --- :mod:`pbxplore.assignment`
=========================================================

.. autofunction:: assign
"""

from __future__ import print_function, absolute_import

# Third-party module
import numpy

# Local module
from . import PB


def assign(dihedrals, pb_ref=PB.REFERENCES):
    """
    Assign Protein Blocks.

    Dihedral angles are provided as a dictionnary with one item per residue.
    The key is the residue number, and the value is a dictionnary with phi and
    psi as keys, and the dihedral angles as values.

    The protein block definitions are provided as a dictionnary. Each key is a
    block name, the values are lists of dihedral angles.

    Parameters
    ----------
    dihedrals : dict
        Phi and psi dihedral angles for each residue.
    pb_ref : dict
        The definition of the protein blocks.
    """
    pb_seq = ""

    # Transform the dict into a numpy array with the right order
    ref = numpy.array([pb_ref[key] for key in sorted(pb_ref)])

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
        except KeyError:
            pb_seq += "Z"
            continue

        # Compute the RMSD between all reference angles and angles of the current residue
        # (16*8)*(1,8) vectorization
        # Note (ref - ang + 180) % 360 - 180) ensure the real difference between 2 angles
        rmsda = numpy.sum(((ref - angles + 180) % 360 - 180)**2, axis=1)

        # Find the PB with the lowest RMSD
        pb_seq += PB.NAMES[numpy.argmin(rmsda)]

    return pb_seq
