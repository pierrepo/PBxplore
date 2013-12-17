#! /usr/bin/env python

"""
PBcount.py reads protein blocks (PBs) sequence files in fasta format
and computes PBs frequency along protein sequence

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# load modules
#===============================================================================
import PBlib as PB
import optparse 
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python installations are steal using Python < 2.7
import os
import sys
import numpy 


#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# get options
#-------------------------------------------------------------------------------
parser = optparse.OptionParser(
    usage="%prog -f file_1.PB.fasta [options] -o output_root_name",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-f", action="append", type="string", 
help="name(s) of the PBs file (in fasta format)")
mandatory_opts.add_option("-o", action="store", type="string", 
help="root name for results")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--residue-shift", action="store", type="int",
    dest = "residue_shift", help="shift to adjust residue number")
optional_opts.add_option("--first-frame", action="store", type="int",
    dest = "first_frame", help="lower index of trajectory frame (default = 1)")
parser.add_option_group(optional_opts)
# get all parameters
(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
# check options
#-------------------------------------------------------------------------------
if not options.f:
    parser.print_help()
    parser.error("option -f is mandatory")

if not options.o:
    parser.print_help()
    parser.error("option -o is mandatory")

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
        sys.exit("%s does not appear to be a valid file.\nBye" % name)
    
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
            sys.exit("%s is not a valid protein block (abcdefghijklmnop)" 
            % block)

#-------------------------------------------------------------------------------
# write PBs count file
#-------------------------------------------------------------------------------
shift = 0
if options.residue_shift:
	shift = options.residue_shift
	print "first residue will be numbered %d" % (shift + 1)

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
print "wrote %s" % count_file_name

