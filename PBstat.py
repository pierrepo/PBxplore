#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Statistical analysis and graphical representations of PBs.

Compute Neq, PBs distribution and draw logo representation of PBs.

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
import subprocess
import argparse

## third-party module
import numpy 

## local module
import PBlib as PB

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
# Functions
#===============================================================================
def array_to_string(ar):
    """
    Convert numpy array to string for further import in R
    """
    numpy.set_printoptions(threshold=numpy.inf)
    #  max_line_width : big enough to avoid unneeded newline
    #  precision : float with 4 digits
    return numpy.array_str(ar, max_line_width = 100000, precision = 4).translate(None, '[]')

#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# manage parameters
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description = "Statistical analysis and graphical representations of PBs.")

# mandatory arguments
parser.add_argument("-f", action="store", required=True,
    help="name of file that contains PBs frequency (count)")
parser.add_argument("-o", action="store", required=True,
    help="name for results")

# optional arguments
parser.add_argument("--map", action="store_true", default=False, dest="mapdist",
    help="generate map of the distribution of PBs along protein sequence")
parser.add_argument("--neq", action="store_true", default=False, dest="neq", 
    help="compute Neq and generate Neq plot along protein sequence")
parser.add_argument("--logo", action="store_true", default=False, dest="logo", 
    help="generate logo representation of PBs frequency along protein sequence")
parser.add_argument("--residue-min", action="store", type=int,
    dest="residue_min", help="defines lower bound of residue frame")
parser.add_argument("--residue-max", action="store", type=int,
    dest="residue_max", help="defines upper bound of residue frame")

# get all parameters
options = parser.parse_args()

#-------------------------------------------------------------------------------
# load and check data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  check count file is readable
if not os.path.isfile(options.f):
    sys.exit( "ERROR: {0}: not a valid file".format(options.f) )

# load count file 
# skip first row that contains PBs labels
try:
    freq = numpy.loadtxt(options.f, dtype=int, skiprows=1)
except:
    sys.exit( "ERROR: {0}: wrong data format".format(options.f) )

# check format
# 17 columns (residue number + 16 PBs) should be present
if len(freq[0,:]) != (PB.NUMBER + 1):
    sys.exit( "ERROR: {0}: wrong data format".format(options.f) )

# read residue numbers
residue_lst = list(freq[:, 0])

#-------------------------------------------------------------------------------
# check / define residue min / max
if options.residue_min:
    residue_min = options.residue_min
else:
    residue_min = min(freq[:,0])

if options.residue_max:
    residue_max = options.residue_max
else:
    residue_max = max(freq[:,0])

if residue_min < 0:
    sys.exit("ERROR: residue_min must be >= 0")

if residue_max < 0:
    sys.exit("ERROR: residue_max must be >= 0")

if residue_min >= residue_max:
    sys.exit("ERROR: residue_min must be > residue_max")

if residue_min not in residue_lst:
    sys.exit( "ERROR: residue_min does not belong to the residue range \
    in {0}".format(options.f) )

if residue_max not in residue_lst:
    sys.exit( "ERROR: residue_max does not belong to the residue range \
    in {0}".format(options.f) )    


# get index of first residue
try:
    first_residue_index = int(freq[0, 0])
except:
    sys.exit( "ERROR: cannot read index of first residue. \
    Wrong data format in {0}".format(options.f) )
print( "Index of first residue is: {0}".format(first_residue_index) )

#  slice data to the required frame
freq = freq[residue_min - first_residue_index : residue_max - first_residue_index + 1, : ]

# determine number of sequences compiled
# use the sum of all residue at position 3
# since positions 1 and 2 have no PBs assignement
# and begin at 1 to not sum the index of the line (here is 3)
sequence_number = sum(freq[2, 1:])
if sequence_number == 0:
    sys.exit("ERROR: counting 0 sequences!")

# update residue numbers
residue_lst = list(freq[:, 0])

# remove residue number    
# extract and normalize PBs frequencies
freq = freq[:, 1:] / float(sequence_number)

#-------------------------------------------------------------------------------
# generates map of the distribution of PBs along protein sequence
#-------------------------------------------------------------------------------
if options.mapdist:
    #  convert array to string for further import in R
    freq_string = array_to_string(freq)

    # build R script
    #-------------------------------------------------------------------------------
    # data
    R_script="""
connector = textConnection("%s")
freq = read.table(connector, header = FALSE)
xnames = %d:%d
ynames = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')
    """ % (freq_string, residue_min, residue_max)

    # define output file name
    map_file_name = options.o + ".PB.map.png"
    if options.residue_min or options.residue_max:
        map_file_name = options.o + ".PB.map.%i-%i.png" % (residue_min, residue_max)

    #  graphical parameters
    R_script += """
png(filename='%s', width = log(length(xnames))*500, height = 1000)
par(
    # default margins are: 5.1 4.1 4.1 2.1
    # extend bottom margin for text (+5 line)
    mar = c(5.1, 5.1, 4.1, 2.1),
    oma = c(2,0,0,0), # 2 lines for comments: 0 to 4
    lwd=3,            # line width
    bty = 'o',        # type of box around graphic
    font.lab = 2,     # axis label font (bold)
    font.axis = 2,    # axis font (bold)
    cex.lab=1.9,      # axis label width
    cex.axis=1.5      # axis width
)
    """ % (map_file_name)

    # color gradient goes 
    # from dark blue (freq = 0) to green/yellow (freq ~ 0.5) to red (freq = 1)
    R_script += """
grad = matrix(nrow=848, ncol=3)
grad[1,1] = 20
grad[1,2] = 20
grad[1,3] = 232
for(i in 2:212){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]+1
grad[i,3] = grad[i-1,3]
}
for(i in 213:424){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]
grad[i,3] = grad[i-1,3]-1
}
for(i in 425:636){
grad[i,1] = grad[i-1,1]+1
grad[i,2] = grad[i-1,2]
grad[i,3] = grad[i-1,3]
}
for(i in 637:848){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]-1
grad[i,3] = grad[i-1,3]
}
colorpal = rgb(grad[,1]/255,grad[,2]/255,grad[,3]/255)
    """

    # plot data map
    R_script += """
layout(matrix(1:2, 1, 2), width=c(log(length(xnames))*3, 1))

par(mar = c(5.1, 7.1, 4.1, 1.1))
image(as.matrix(freq), axes=FALSE, xlab="Residue number", col=colorpal, zlim = c(0, 1))
box()
axis(1, seq(0, 1, 1/(length(xnames)-1)), xnames)
axis(2, seq(0, 1, 1/(length(ynames)-1)), ynames, font = 4)
mtext('PBs', side = 2, line = 5, cex=1.9, font=2)
mtext(bquote(beta~'strand'), side = 2, line = 3, at = 3*1/15, cex=1.5)
mtext('coil', side = 2, line = 3, at = 7*1/15, cex=1.5)
mtext(bquote(alpha~'helix'), side = 2, line = 3, at = 12*1/15, cex=1.5)

par(mar = c(5.1, 1.1, 4.1, 5.1))
image(t(seq(1, 848)), col=colorpal, axes=FALSE)
axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2))
mtext("intensity", side = 4, line = 3, cex = 1.7, font = 2)
box()
    """

    # execute R script
    #-------------------------------------------------------------------------------
    command="R --vanilla --slave"
    proc = subprocess.Popen(command, shell = True, 
    stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    (out, err) = proc.communicate(R_script)
    if err:
        print( "ERROR: {0}".format(err) )
    code = proc.wait()
    if code:
        print( "ERROR: exit code != 0" )
        print( "exit code: {0}".format(code) )
    else:
        print( "wrote {0}".format(map_file_name) )
    print(out)

#-------------------------------------------------------------------------------
# computes Neq and generates neq plot along protein sequence
#-------------------------------------------------------------------------------
if options.neq:
    # compute Neq
    #-------------------------------------------------------------------------------
    neq_array = numpy.zeros((len(residue_lst), 2))
    neq_array[:, 0] = residue_lst
    for idx in range(len(residue_lst)):
        H = 0.0
        for b in range(PB.NUMBER):
            f = freq[idx, b] 
            if f != 0:
                H += f * math.log(f)
        neq_array[idx, 1] = math.exp(-H)
    
    # define output file name
    #-------------------------------------------------------------------------------
    neq_file_name = options.o + ".PB.Neq"
    if options.residue_min or options.residue_max:
        neq_file_name = options.o + ".PB.Neq.%i-%i" % (residue_min, residue_max)
    
    # write Neq
    #-------------------------------------------------------------------------------
    content = "%-6s %8s \n" % ("resid", "Neq")
    for (res, neq) in neq_array:
        content += "%-6d %8.2f \n" % (res, neq)
    neq_file = open(neq_file_name, "w")
    neq_file.write(content)
    neq_file.close()
    print( "wrote {0}".format(neq_file_name) )
    
    #  convert array to string for further import in R
    neq_array_string = array_to_string(neq_array)

    # build R script
    #-------------------------------------------------------------------------------

    R_script="""
connector = textConnection("%s")
neq = read.table(connector, header = FALSE)
    """ % (neq_array_string)

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
    """ % (neq_file_name + ".png")

    # plot data
    R_script += """
plot(neq, type= 'l', 
xlab = 'Residue number', ylab = 'Neq', 
ylim=c(1,max(round(neq[,2]))+2))
    """

    # execute R script
    #-------------------------------------------------------------------------------
    command="R --vanilla --slave"
    proc = subprocess.Popen(command, shell = True, 
    stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    (out, err) = proc.communicate(R_script)
    if err:
        print( "ERROR: {0}".format(err) )
    code = proc.wait()
    if code:
        print( "ERROR: exit code != 0" )
        print( "exit code: {0}".format(code) )
    else:
        print( "wrote {0}.png".format(neq_file_name) )
    
    print( out )

#-------------------------------------------------------------------------------
# generates logo representation of PBs frequency along protein sequence
#
# G. E. Crooks, G. Hon, J.-M. Chandonia, and S. E. Brenner. 
# 'WebLogo: A Sequence Logo Generator.'
# Genome Research 14:1188–90 (2004)
# doi:10.1101/gr.849004.
# http://weblogo.threeplusone.com/
#-------------------------------------------------------------------------------
if options.logo:
    # read count file
    #-------------------------------------------------------------------------------
    f_in = open(options.f, 'r')
    count_content = f_in.readlines()
    f_in.close()

    # convert a table of PB frequencies into transfac format as required by weblogo
    # http://meme.sdsc.edu/meme/doc/transfac-format.html
    #-------------------------------------------------------------------------------
    residue_lst = []
    transfac_content  = "ID %s\n" % options.f
    transfac_content += "BF unknown\n"
    transfac_content += "P0" + count_content[0][2:]
    for line in count_content[1:]:
        item = line.split()
        residue = int(item[0])
        residue_lst.append(residue)
        transfac_content += "%05d" % residue + line[5:-1] +  "    X" + "\n"
    transfac_content += "XX\n"
    transfac_content += "//"

    # write transfac file (debug only)
    #-------------------------------------------------------------------------------
    debug = False
    if debug:
        transfac_name = options.o + ".PB.transfac"
        f_out = open(transfac_name, 'w')
        f_out.write(transfac_content)
        f_out.close()
        print( "wrote {0}".format(transfac_name) )

    # define output file name
    #-------------------------------------------------------------------------------
    logo_file_name = options.o + ".PB.logo.pdf"
    if options.residue_min or options.residue_max:
        logo_file_name = options.o + ".PB.logo.%i-%i.pdf" % (residue_min, residue_max)

    # call weblogo
    #-------------------------------------------------------------------------------
    command = """weblogo \
--format pdf \
--errorbars NO \
--fineprint "PBlogo" \
--title %s \
--color "#1240AB" d      "strand main" \
--color "#1240AB" abcdef "strand others" \
--color "#0BD500" ghij "coil" \
--color "#FD0006" m     "helix main" \
--color "#FD0006" klnop "helix others" \
--composition none \
--datatype transfac \
-s large \
-o %s \
--lower %i \
--upper %i 
    """ % (options.f.replace(".PB.count", ""), logo_file_name, 
    residue_min, residue_max)

    proc = subprocess.Popen(command, shell = True,
    stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    (out, err) = proc.communicate(transfac_content)
    if err:
        print( "ERROR: {0}".format(err) )
    code = proc.wait()
    if code:
        print( "ERROR: exit code != 0" )
        print( "exit code: {0}".format(code) )
    else:
        print( "wrote {0}".format(logo_file_name) )
    print( out )

