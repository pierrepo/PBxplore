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


def _angle_modulo_360(angle):
    """
    Keep angle in the range -180 / +180 [degrees]
    """
    if angle > 180.0:
        return angle - 360.0
    elif angle < -180.0:
        return angle + 360.0
    else:
        return angle


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

        # convert to array
        angles = numpy.array(angles)

        # compare to reference PB angles
        rmsda_lst = {}
        for block in pb_ref:
            diff = pb_ref[block] - angles
            diff2 = _angle_modulo_360_vect(diff)
            rmsda = numpy.sum(diff2**2)
            rmsda_lst[rmsda] = block
        pb_seq += rmsda_lst[min(rmsda_lst)]
    return pb_seq


# vertorize function
_angle_modulo_360_vect = numpy.vectorize(_angle_modulo_360)
_angle_modulo_360_vect.__doc__ = _angle_modulo_360.__doc__
