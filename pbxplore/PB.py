#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Protein Blocks related data --- :mod:`pbxplore.PB`
==================================================

Protein blocks definition
-------------------------

* :data:`REFERENCES`
  The definition of each block as a dictionary; each key is a block name (as
  a lower case letter), and each value is a list of the dihedral angle values
  that define the block.
* :data:`NAMES`
  The names of all the protein blocks.

Substitution matrix
-------------------

* :data:`SUBSTITUTION_MATRIX_NAME`
  The absolute path to the file that contains the subtitution matrix

.. autofunction:: load_substitution_matrix

Exceptions
----------

.. autoexception:: SizeError

.. autoexception:: InvalidBlockError
"""

from __future__ import print_function, absolute_import

# Standard modules
import os

# Third-party modules
import numpy

# Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass

# Data
# Protein Blocks reference angles
# taken from A. G. de Brevern, C. Etchebest and S. Hazout.
# "Bayesian probabilistic approach for predicting backbone structures
# in terms of protein blocks"
# Proteins, 41: 271-288 (2000)
REFERENCES = {
    'a': [ 41.14,   75.53,  13.92,  -99.80,  131.88,  -96.27, 122.08,  -99.68],
    'b': [108.24,  -90.12, 119.54,  -92.21,  -18.06, -128.93, 147.04,  -99.90],
    'c': [-11.61, -105.66,  94.81, -106.09,  133.56, -106.93, 135.97, -100.63],
    'd': [141.98, -112.79, 132.20, -114.79,  140.11, -111.05, 139.54, -103.16],
    'e': [133.25, -112.37, 137.64, -108.13,  133.00,  -87.30, 120.54,   77.40],
    'f': [116.40, -105.53, 129.32,  -96.68,  140.72,  -74.19, -26.65,  -94.51],
    'g': [  0.40,  -81.83,   4.91, -100.59,   85.50,  -71.65, 130.78,   84.98],
    'h': [119.14, -102.58, 130.83,  -67.91,  121.55,   76.25,  -2.95,  -90.88],
    'i': [130.68,  -56.92, 119.26,   77.85,   10.42,  -99.43, 141.40,  -98.01],
    'j': [114.32, -121.47, 118.14,   82.88, -150.05,  -83.81,  23.35,  -85.82],
    'k': [117.16,  -95.41, 140.40,  -59.35,  -29.23,  -72.39, -25.08,  -76.16],
    'l': [139.20,  -55.96, -32.70,  -68.51,  -26.09,  -74.44, -22.60,  -71.74],
    'm': [-39.62,  -64.73, -39.52,  -65.54,  -38.88,  -66.89, -37.76,  -70.19],
    'n': [-35.34,  -65.03, -38.12,  -66.34,  -29.51,  -89.10,  -2.91,   77.90],
    'o': [-45.29,  -67.44, -27.72,  -87.27,    5.13,   77.49,  30.71,  -93.23],
    'p': [-27.09,  -86.14,   0.30,   59.85,   21.51,  -96.30, 132.67,  -92.91],
}
#   PB  psi(n-2) phi(n-1) psi(n-1)   phi(n)  psi(n)  phi(n+1) psi(n+1) phi(n+2)

NAMES = 'abcdefghijklmnop'  # name of the 16 PBs
SUBSTITUTION_MATRIX_NAME = os.path.join(os.path.dirname(__file__),
                                        "PBs_substitution_matrix.dat")


class InvalidBlockError(ValueError):
    """
    Exception raised when encounter an invalid protein block.
    """
    def __init__(self, block=None):
        super(InvalidBlockError, self).__init__(self)
        self.block = block

    def __repr__(self):
        if self.block is None:
            return "Invald block"
        else:
            return "Ivalid block '{}'".format(self.block)


class SizeError(AssertionError):
    """
    Exception raised when a sequence does not have the expected length.
    """
    pass


def load_substitution_matrix(name):
    """
    Load PB substitution matrix.

    The matrix must be 16x16.

    Parameters
    ----------
    name : str
        Name of the file containing the PBs susbtitution matrix.

    Returns
    -------
    mat : numpy array
        Array of floats.

    Raises
    ------
    InvalidBlockError
        encountered an unexpected PB
    """
    mat = numpy.loadtxt(name, dtype=float, skiprows=2)
    assert mat.shape == (16, 16), 'wrong substitution matrix size'
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] != mat[j][i]:
                raise ValueError("Matrix is not symetric - idx {} and {}".format(i, j))
    return mat
