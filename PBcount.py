#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute PB frequency along protein sequence.

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
import argparse

## third-party module
import numpy 

## local module
import PBlib as PB

#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# get options
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description = 'Compute PB frequency along protein sequence.')

# mandatory arguments
parser.add_argument("-f", action="append", required = True,
    help="name(s) of the PBs file (in fasta format)")
parser.add_argument("-o", action="store", required = True,
    help="name for results")


# optional arguments
parser.add_argument("--residue-shift", action="store", type=int,
    dest = "residue_shift", help="shift to adjust residue number")
parser.add_argument("--first-frame", action="store", type=int,
    dest = "first_frame", help="lower index of trajectory frame (default = 1)")

# get all arguments
options = parser.parse_args()

#-------------------------------------------------------------------------------
# check options
#-------------------------------------------------------------------------------
if options.residue_shift and options.residue_shift < 0:
	parser.error("residue shift must be positive")

if not options.first_frame :
    index_first_frame = 0
elif options.first_frame <= 0:
    parser.error("lower index of trajectory frame must be positive")
else:
    index_first_frame = options.first_frame-1

#-------------------------------------------------------------------------------
# check input files
#-------------------------------------------------------------------------------
for name in options.f:
    if not os.path.isfile(name):
        sys.exit( "{0}: not a valid file. Bye.".format(name) )
    
#-------------------------------------------------------------------------------
# read PBs files
#-------------------------------------------------------------------------------
pb_seq = []
pb_name = []
for name in options.f:
    header, seq = PB.read_fasta(name)
    pb_name += header
    pb_seq += seq

#-------------------------------------------------------------------------------
# check all sequences have the same size
#-------------------------------------------------------------------------------
pb_seq_size = len(pb_seq[0])
for seq in pb_seq:
    if len(seq) != pb_seq_size:
        sys.exit("cannot compute PB frequencies / different sequence lengths")

#-------------------------------------------------------------------------------
# count PBs at each position of the sequence
#-------------------------------------------------------------------------------
pb_count = numpy.zeros((pb_seq_size, len(PB.NAMES)))

for seq in pb_seq:
    for idx, block in enumerate(seq):
        if block in PB.NAMES:
            pb_count[idx, PB.NAMES.index(block)] += 1.0
        elif block not in ["Z", "z"]:
            msg = "{0} is not a valid protein block (abcdefghijklmnop)".format
            sys.exit( msg(block) )

#-------------------------------------------------------------------------------
# write PBs count file
#-------------------------------------------------------------------------------
shift = 0
if options.residue_shift:
	shift = options.residue_shift
	print( "first residue will be numbered {0}".format(shift + 1) )

count_file_name = options.o + ".PB.count"
content = "    "
# build header (PB names)
content += "".join(["%6s" % name for name in PB.NAMES]) + "\n"
# build data table
for residue_idx, residue_pb in enumerate(pb_count):
    content += "%-5d" % (residue_idx + 1 + shift) + " ".join("%5d" % i for i in residue_pb) + "\n"
# write data
count_file = open(count_file_name, "w")
count_file.write(content)
count_file.close()
print( "wrote {0}".format(count_file_name) )

