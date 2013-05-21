#! /usr/bin/env python

"""
PBstat.py reads protein blocks (PBs) sequence files in fasta format,
computes statistics (count and Neq) and optionally displays Neq.

2012 - P. Poulain, A. G. de Brevern 
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
import math
import subprocess

#===============================================================================
# data
#===============================================================================
# 16 PBs
PB_DIC = {"a":0, "b":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8,
"j":9, "k":10, "l":11, "m":12, "n":13, "o":14, "p":15}
PB_DIC_SIZE = len(PB_DIC)


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
    usage="%prog -f file_1.PB.fasta [-f file_2.PB.fasta] -o output_root_name",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-f", action="append", type="string", 
help="name(s) of the PB file (in fasta format)")
mandatory_opts.add_option("-o", action="store", type="string", 
help="root name for results")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--neq-lower", action="store", type="int",
    dest = "neq_lower", help="lower bound for Neq display")
optional_opts.add_option("--neq-upper", action="store", type="int",
    dest = "neq_upper", help="upper bound for Neq display")
optional_opts.add_option("--neq-shift", action="store", type="int",
    dest = "neq_shift", help="shift to adjust residue number")
optional_opts.add_option("--no-neq", action="store_true",
    dest = "no_neq", help="disables Neq display")
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

#-------------------------------------------------------------------------------
# check input files
#-------------------------------------------------------------------------------
for name in options.f:
    if not os.path.isfile(name):
        sys.exit("%s does not appear to be a valid file.\nBye" % name)
    
#-------------------------------------------------------------------------------
# read PB files
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
# count PB at each position of the sequence
#-------------------------------------------------------------------------------
pb_seq_nb = len(pb_seq)
pb_count = numpy.zeros((pb_seq_size, PB_DIC_SIZE))

for seq in pb_seq:
    for idx, block in enumerate(seq):
        if block in PB_DIC:
            pb_count[idx, PB_DIC[block]] += 1.0
        elif block not in ["Z", "z"]:
            sys.exit("%s is not a valid protein block (abcdefghijklmnop)" 
            % block)

#-------------------------------------------------------------------------------
# write PB counts
#-------------------------------------------------------------------------------
count_file_name = options.o + ".PB.count"
content = "    "
# build header (PB names)
for PB_name in sorted(PB_DIC):
    content += "%6s" % PB_name
content += "\n"
# build data table
for residue_idx, residue_pb in enumerate(pb_count):
    content += "%-5d" % (residue_idx + 1) + " ".join("%5d" % i for i in residue_pb) + "\n"
# write data
count_file = open(count_file_name, "w")
count_file.write(content)
count_file.close()
print "wrote %s" % count_file_name

#-------------------------------------------------------------------------------
# compute PB frequencies
#-------------------------------------------------------------------------------
pb_freq = pb_count / pb_seq_nb

#-------------------------------------------------------------------------------
# compute Neq
#-------------------------------------------------------------------------------
neq = []
for pos in xrange(pb_seq_size):
    H = 0.0
    for b in xrange(PB_DIC_SIZE):
        f = pb_freq[pos, b] 
        if f != 0:
            H += f * math.log(f)
    neq.append( math.exp(-H))

#-------------------------------------------------------------------------------
# write Neq
#-------------------------------------------------------------------------------
neq_file_name = options.o + ".PB.Neq"
content = "%-6s %8s \n" % ("resid", "Neq")
for idx in xrange(len(neq)):
    content += "%-6d %8.2f \n" % (idx+1, neq[idx])
neq_file = open(neq_file_name, "w")
neq_file.write(content)
neq_file.close()
print "wrote %s" % neq_file_name

#-------------------------------------------------------------------------------
# display Neq (generates a png picture)
#-------------------------------------------------------------------------------
if options.no_neq:
    print "no Neq display generated"
    sys.exit(0)


#-------------------------------------------------------------------------------
# prepare data
#-------------------------------------------------------------------------------
res = range(1, len(neq)+1)

if options.neq_lower:
    lower = options.neq_lower
else:
    lower = min(res)
    
if options.neq_upper:
    upper = options.neq_upper
else:
    upper = max(res)

if options.neq_shift:
    shift = options.neq_shift
else:
    shift = 0
    
# create numpy array for export
# and select data range
data = numpy.transpose([res, neq])
data = data[lower-1:upper, :]

#  convert array to string for further import in R
numpy.set_printoptions(threshold=numpy.inf)
str_data = numpy.array_str(data, max_line_width = 100000, precision = 4).translate(None, '[]')
#  max_line_width : big enough to avoid unneeded newline
#  precision : float with 4 digits

#-------------------------------------------------------------------------------
# build R script
#-------------------------------------------------------------------------------

R_script="""
connector = textConnection("%s")
data = read.table(connector, header = FALSE)
lower = %d
upper = %d
shift = %d
""" % (str_data, lower, upper, shift)

#  graphical parameters
R_script += """
png(filename='%s', width = 1600, height = 1200)
par(
    # default margins are: 5.1 4.1 4.1 2.1
    # extend bottom margin for text (+5 line)
    mar = c(5.1, 5.1, 4.1, 2.1),
    oma = c(2,0,0,0), # 2 lines for comments: 0 to 4
    lwd=3,            # line width
    bty = 'o',        # type of box around graphic
    font.lab = 2,     # axis label font (bold)
    font.axis = 2,    # axis font (bold)
    cex.lab=2.5,      # axis label width
    cex.axis=2.0      # axis width
)
""" % (options.o + ".PB.Neq.png")

# plot data
R_script += """
plot(data, type= 'l', 
    xlab = 'Residue number', ylab = 'Neq', 
    xaxt="n", ylim=c(1,max(round(data[,2]))+2))
axis(1, lower:upper, (lower:upper)+shift)
"""

#-------------------------------------------------------------------------------
# execute R script
#-------------------------------------------------------------------------------
command="R --vanilla --slave"
proc = subprocess.Popen(command, shell = True, 
stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
(out, err) = proc.communicate(R_script)
if err:
    print "ERROR:", err
code = proc.wait()
if code:
    print "ERROR: exit code != 0"
    print "exit code:", code
else:
    print "wrote %s" % options.o + ".PB.Neq.png"


print out
