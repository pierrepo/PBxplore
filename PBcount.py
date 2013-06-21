#! /usr/bin/env python

"""
PBcount.py reads protein blocks (PBs) sequence files in fasta format
and computes PBs frequency along protein sequence

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# load modules
#===============================================================================
import optparse 
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python installations are steal using Python < 2.7
import os
import sys
import numpy 

#===============================================================================
# data
#===============================================================================
# 16 PBs
PB_DIC = {"a":0, "b":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8,
"j":9, "k":10, "l":11, "m":12, "n":13, "o":14, "p":15}


#===============================================================================
# functions
#===============================================================================
def read_fasta(name):
    """read fasta file and output sequences in a list"""
    seqLst = []
    seq = ""
    f_in = open(name, "r")
    for line in f_in:
        data = line.strip()
        if data and ">" not in data:
            seq += data
        if seq and data and ">" == data[0]:
            seqLst.append(seq)
            seq = ""
    f_in.close()
    # save last sequence
    if seq:
        seqLst.append(seq)
    # outputs
    print "read %d sequences in %s" % (len(seqLst), name)
    return seqLst

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
	parser.error("residue shift must positive")

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
for name in options.f:
    pb_seq += read_fasta(name)

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
pb_count = numpy.zeros((pb_seq_size, len(PB_DIC)))

for seq in pb_seq:
    for idx, block in enumerate(seq):
        if block in PB_DIC:
            pb_count[idx, PB_DIC[block]] += 1.0
        elif block not in ["Z", "z"]:
            sys.exit("%s is not a valid protein block (abcdefghijklmnop)" 
            % block)

#-------------------------------------------------------------------------------
# write PBs count file
#-------------------------------------------------------------------------------
shift = 0
if options.residue_shift:
	shift = options.residue_shift
	print "first residue will be numbered %d" % shift

count_file_name = options.o + ".PB.count"
content = "    "
# build header (PB names)
for PB_name in sorted(PB_DIC):
    content += "%6s" % PB_name
content += "\n"
# build data table
for residue_idx, residue_pb in enumerate(pb_count):
    content += "%-5d" % (residue_idx + 1 + shift) + " ".join("%5d" % i for i in residue_pb) + "\n"
# write data
count_file = open(count_file_name, "w")
count_file.write(content)
count_file.close()
print "wrote %s" % count_file_name


