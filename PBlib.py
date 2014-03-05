#! /usr/bin/env python

"""
PBlib.py 

Python library to handle Protein Blocks

2013 - P. Poulain, A. G. de Brevern 
"""
#===============================================================================
# modules import
#===============================================================================
import os

#===============================================================================
# data
#===============================================================================

# Protein Blocks angle definitions
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

SUBSTITUTION_MATRIX_NAME = "PBs_substitution_matrix.dat"


# line width for fasta format
FASTA_WIDTH = 60



#===============================================================================
# functions
#===============================================================================
def read_fasta(name):
    """read fasta file and output sequences in a list
    """
    assert os.path.exists(name), name + ' does not exist'
    sequence_lst = []
    header_lst = []
    header = ""
    sequence = ""
    f_in = open(name, "r")
    for line in f_in:
        data = line.strip()
        if data and ">" == data[0]:
            header = data[1:]
        if data and ">" not in data:
            sequence += data
        if sequence and data and ">" == data[0]:
            header_lst.append(header)
            sequence_lst.append(sequence)
            header = ""
            sequence = ""
    f_in.close()
    # save last sequence
    header_lst.append(header)
    sequence_lst.append(sequence)
    # outputs
    assert len(header_lst) == len(sequence_lst), \
           "cannot read same number of headers and sequences"
    print "read %d sequences in %s" % (len(sequence_lst), name)
    return header_lst, sequence_lst


def clean_file(name):
    """clean existing file
    """
    if os.path.exists(name):
        os.remove(name)

#-------------------------------------------------------------------------------
def write_fasta(name, seq, comment):
    """format seq and comment to fasta format
       and write file
    """
    fasta_content  = ">"+comment+"\n"
    fasta_content += "\n".join( [seq[i:i+FASTA_WIDTH] for i in xrange(0, len(seq), FASTA_WIDTH)] )
    fasta_content += "\n"
    f_out = open(name, "a")
    f_out.write(fasta_content)
    f_out.close()