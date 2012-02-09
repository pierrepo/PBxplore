#! /usr/bin/env python

"""
PBstat.py reads protein blocks (PBs) sequence files in fasta format
and compute statisiques (count, frequency and Neq).

2012 - P. Poulain, A. G. de Brevern 
"""


#===============================================================================
# load modules
#===============================================================================
from optparse import OptionParser
import os
import sys
import numpy 
import math

#===============================================================================
# data
#===============================================================================
# 16 PBs + block "Z"
PBdic = {"a":0, "b":1, "c":2,"d":3, "e":4, "f":5, "g":6, "h":7, "i":8, \
"j":9, "k":10, "l":11, "m":12, "n":13, "o":14, "p":15}
PBdicSize = len(PBdic)


#===============================================================================
# functions
#===============================================================================
def readFasta(name):
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
    print "read %d sequences in %s" %(len(seqLst), name)
    return seqLst
#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# get options
#-------------------------------------------------------------------------------
parser = OptionParser(usage="%prog -f PB_1.PB [-f PB_2.PB] -o output_root_name")
parser.add_option("-f", action="append", type="string", dest="PBNameLst", help="name(s) of the PB file (in fasta format)")
parser.add_option("-o", action="store", type="string", dest="outName", help="root name for results")

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
# check options
#-------------------------------------------------------------------------------
if not options.PBNameLst:
    parser.print_help()
    parser.error("option -f is mandatory")

if not options.outName:
    parser.print_help()
    parser.error("option -o is mandatory")


PBNameLst = options.PBNameLst
outName = options.outName

#-------------------------------------------------------------------------------
# check files
#-------------------------------------------------------------------------------
for name in PBNameLst:
    if not os.path.isfile(name):
        sys.exit("%s does not appear to be a valid file.\nBye" %(name))
    
#-------------------------------------------------------------------------------
# read PB files
#-------------------------------------------------------------------------------
PBseq = []
for name in PBNameLst:
    PBseq += readFasta(name)

#-------------------------------------------------------------------------------
# check all sequences have the same size
#-------------------------------------------------------------------------------
PBseqSize = len(PBseq[0])
for seq in PBseq:
    if len(seq) != PBseqSize:
        sys.exit("cannot compute PB frequencies / different sequence lengths")

#-------------------------------------------------------------------------------
# count PB at each position of the sequence
#-------------------------------------------------------------------------------
PBseqNb = len(PBseq)
PBcount = numpy.zeros((PBseqSize, PBdicSize))

for seq in PBseq:
    for idx, block in enumerate(seq):
        if block in PBdic:
            PBcount[idx, PBdic[block]] += 1.0

#        else:
#            print "%s is not a valid protein block" %(block)
#            print "skipping this sequence"
#            break

#-------------------------------------------------------------------------------
# write PB counts
#-------------------------------------------------------------------------------
countName = outName + ".PB.count"
content = "    "
# build header (PB names)
for PB_name in sorted(PBdic):
    content += "%6s" %(PB_name)
content += "\n"
# build data table
for residue_idx, residue_PB in enumerate(PBcount):
    content += "%-5d" % (residue_idx + 1) + " ".join("%5d" % i for i in residue_PB) + "\n"
# write data
f_in = open(countName, "w")
f_in.write(content)
f_in.close()
print "wrote %s" %(countName)

#-------------------------------------------------------------------------------
# compute PB frequencies
#-------------------------------------------------------------------------------
PBfreq = PBcount / PBseqNb

#-------------------------------------------------------------------------------
# compute Neq
#-------------------------------------------------------------------------------
Neq = []
for pos in xrange(PBseqSize):
    H = 0.0
    for b in xrange(PBdicSize):
        f = PBfreq[pos, b] 
        if f != 0:
            H += f * math.log(f)
    Neq.append( math.exp(-H))

#-------------------------------------------------------------------------------
# write Neq
#-------------------------------------------------------------------------------
NeqName = outName + ".PB.Neq"
content = "%6s %8s \n" %("#resid", "Neq")
for idx in xrange(len(Neq)):
    content += "%-6d %8.2f \n" %(idx+1, Neq[idx])
f_in = open(NeqName, "w")
f_in.write(content)
f_in.close()
print "wrote %s" %(NeqName)

