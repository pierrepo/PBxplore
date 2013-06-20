#! /usr/bin/env python

"""
PBmap.py reads protein blocks (PBs) counts
and outputs frequency map for all PBs.

2013 - P. Poulain, A. G. de Brevern 
"""
#===============================================================================
# load modules
#===============================================================================
import optparse
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python setups are still using Python < 2.7
import subprocess
import numpy
import sys 
import os

#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# manage parameters
#-------------------------------------------------------------------------------
parser = optparse.OptionParser(
    usage="%prog -f file.PB.count [--min residue --max residue] -o output_root_name",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-f", action="store", type="string",
    help="name of file that contains PB frequency")
mandatory_opts.add_option("-o", action="store", type="string",
    help="root name for results")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--min", action="store", type="int",
    help="lower residue number to display PB map")
optional_opts.add_option("--max", action="store", type="int",
    help="upper residue number to display PB map")
optional_opts.add_option("--shift", action="store", type="int",
    help="shift to adjust residue number")
parser.add_option_group(optional_opts)
# get all parameters
(options, args) = parser.parse_args()

#  check options
if not options.f:
    parser.print_help()
    parser.error("option -f in mandatory")

name = options.f

if not options.o:
    parser.print_help()
    parser.error("option -o is mandatory")

#-------------------------------------------------------------------------------
# load and check data
#-------------------------------------------------------------------------------

#  check file
if not os.path.isfile(options.f):
    sys.exit("ERROR: %s is not a valid file" % name)

try:
    data = numpy.loadtxt(name, dtype=int, skiprows=1)
except:
    sys.exit("ERROR: wrong data format in %s" % name)

# 17 rows (residue number + 16 PBs) should be present
if len(data[0,:]) != 17:
    sys.exit("ERROR: wrong data format in %s" % name)

if options.min:
    lower = options.min
else:
    lower = min(data[:,0])

if options.max:
    upper = options.max
else:
    upper = max(data[:,0])

if options.shift:
    shift = options.shift
else:
    shift = 0


#-------------------------------------------------------------------------------
# prepare data
#-------------------------------------------------------------------------------
#  slice data
data = data[lower-1:upper, :]

# remove first column (contains residue number)
data = data[:, 1:]

# normalize data
total_count = sum(data[2,:])
data = data / float(total_count)

#  convert array to string for further import in R
#str_data = " ".join(str(d)+"\n" for d in data).translate(None, '[]')
numpy.set_printoptions(threshold=numpy.inf)
str_data = numpy.array_str(data, max_line_width = 1000, precision = 4).translate(None, '[]')
#  max_line_width : big enough to avoid unneeded newline
#  precision : float with 4 digits

#-------------------------------------------------------------------------------
# build R script
#-------------------------------------------------------------------------------

# data
R_script="""
connector = textConnection("%s")
freq = read.table(connector, header = FALSE)
lower = %d
upper = %d
xnames = lower:upper
ynames = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')
""" % (str_data, lower+shift, upper+shift)

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
""" % (options.o + ".PB.map.png")

# color gradient goes 
# from dark blue (freq = 0) to green/yellow to red (freq = 1)
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
image(as.matrix(freq), axes=FALSE, xlab="Residue number", col=colorpal)
box()
axis(1, seq(0, 1, 1/(length(xnames)-1)), xnames)
axis(2, seq(0, 1, 1/(length(ynames)-1)), ynames, font = 4)
mtext('PB', side = 2, line = 5, cex=1.9, font=2)
mtext(bquote(beta~'strand'), side = 2, line = 3, at = 3*1/15, cex=1.5)
mtext('coil', side = 2, line = 3, at = 7*1/15, cex=1.5)
mtext(bquote(alpha~'helix'), side = 2, line = 3, at = 12*1/15, cex=1.5)

par(mar = c(5.1, 1.1, 4.1, 5.1))
image(t(seq(1, 848)), col=colorpal, axes=FALSE)
axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2))
mtext("intensity", side = 4, line = 3, cex = 1.7, font = 2)
box()

void = dev.off()
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
    print "wrote %s" % options.o + ".PB.map.png"

print out

