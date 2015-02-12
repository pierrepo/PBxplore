#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python library to handle Protein Blocks

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

## third-party modules
import numpy

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
# Data
#===============================================================================
# real module directory
__location__ = os.path.realpath(
               os.path.join(os.getcwd(), os.path.dirname(sys.argv[0])) )

# Protein Blocks angles definitions
# taken from A. G. de Brevern, C. Etchebest and S. Hazout. 
# "Bayesian probabilistic approach for predicting backbone structures 
# in terms of protein blocks"
# Proteins, 41: 271-288 (2000)
DEFINITIONS = """
#PB psi(n-2) phi(n-1)  psi(n-1)   phi(n)   psi(n)   phi(n+1)  psi(n+1)  phi(n+2) 
a    41.14      75.53     13.92   -99.80   131.88     -96.27    122.08    -99.68  
b   108.24     -90.12    119.54   -92.21   -18.06    -128.93    147.04    -99.90  
c   -11.61    -105.66     94.81  -106.09   133.56    -106.93    135.97   -100.63 
d   141.98    -112.79    132.20  -114.79   140.11    -111.05    139.54   -103.16 
e   133.25    -112.37    137.64  -108.13   133.00     -87.30    120.54     77.40   
f   116.40    -105.53    129.32   -96.68   140.72     -74.19    -26.65    -94.51  
g     0.40     -81.83      4.91  -100.59    85.50     -71.65    130.78     84.98   
h   119.14    -102.58    130.83   -67.91   121.55      76.25     -2.95    -90.88  
i   130.68     -56.92    119.26    77.85    10.42     -99.43    141.40    -98.01  
j   114.32    -121.47    118.14    82.88  -150.05     -83.81     23.35    -85.82  
k   117.16     -95.41    140.40   -59.35   -29.23     -72.39    -25.08    -76.16  
l   139.20     -55.96    -32.70   -68.51   -26.09     -74.44    -22.60    -71.74  
m   -39.62     -64.73    -39.52   -65.54   -38.88     -66.89    -37.76    -70.19  
n   -35.34     -65.03    -38.12   -66.34   -29.51     -89.10     -2.91     77.90   
o   -45.29     -67.44    -27.72   -87.27     5.13      77.49     30.71    -93.23  
p   -27.09     -86.14      0.30    59.85    21.51     -96.30    132.67    -92.91 
"""

# names of the 16 PBs
NAMES = ["a", "b", "c", "d", "e", "f", "g", "h",
           "i", "j", "k", "l", "m", "n", "o", "p"]
NUMBER = len(NAMES)

print(__location__)
SUBSTITUTION_MATRIX_NAME = os.path.join(__location__, "PBs_substitution_matrix.dat")


# line width for fasta format
FASTA_WIDTH = 60


#===============================================================================
# Functions
#===============================================================================
def get_dihedral(atomA, atomB, atomC, atomD):
    """
    Compute dihedral angle between 4 atoms (A, B, C, D)
    each atom is represented as a list of three coordinates [x, y, z]
    output is in degree in the range -180 / +180
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

#-------------------------------------------------------------------------------    
def read_fasta(name):
    """
    Read fasta file and output sequences in a list
    """
    assert os.path.exists(name), name + ' does not exist'
    sequence_lst = []
    header_lst = []
    header = ""
    sequence = ""
    f_in = open(name, "r")
    for line in f_in:
        data = line.strip()
        # jump empty lines
        if not data:
            continue
        # store header and sequence when a new header (i.e. sequence) is found
        if sequence and header and data.startswith(">"):
            header_lst.append(header)
            sequence_lst.append(sequence)
            # reset header and sequence
            header = ""
            sequence = ""
        # save header of sequence
        if data.startswith(">"):
            header = data[1:]
        # save sequence
        if ">" not in data:
            sequence += data
    f_in.close()
    # save last sequence
    if header and sequence:
        header_lst.append(header)
        sequence_lst.append(sequence)
    # outputs
    assert len(header_lst) == len(sequence_lst), \
           "cannot read same number of headers and sequences"
    print("read %d sequences in %s" % (len(sequence_lst), name))
    if len(sequence_lst) == 0:
        print("WARNING: %s seems empty of sequence" %(name))
    return header_lst, sequence_lst

#-------------------------------------------------------------------------------
def load_substitution_matrix(name):
    """
    Load PB substitution matrix
    """
    try:
        # mat = numpy.loadtxt(name, dtype=int, skiprows=2)
        mat = numpy.loadtxt(name, dtype=float, skiprows=2)
    except:
        sys.exit("ERROR: cannot read %s" % name)
    assert len(mat) == 16, 'wrong substitution matrix size'
    assert len(mat[0]) == 16, 'wrong substitution matrix size'
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] != mat[j][i]:
                print(i, j)
                print(mat[i][j], mat[j][i])
                sys.exit("ERROR: matrix is not symetric - idx %i and %i" % (i, j))
    print("read substitution matrix")
    return mat

#-------------------------------------------------------------------------------
def clean_file(name):
    """
    Clean existing file
    """
    if os.path.exists(name):
        os.remove(name)

#-------------------------------------------------------------------------------
def write_fasta(name, seq, comment):
    """
    Format seq and comment to fasta format
    and write file
    """
    fasta_content  = ">"+comment+"\n"
    fasta_content += "\n".join( [seq[i:i+FASTA_WIDTH] for i in range(0, len(seq), FASTA_WIDTH)] )
    fasta_content += "\n"
    f_out = open(name, "a")
    f_out.write(fasta_content)
    f_out.close()
