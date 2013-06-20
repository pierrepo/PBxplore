#! /usr/bin/env python

"""
PBlogo.py 
calls 'weblogo' to build logo representations of PBs
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
import subprocess

#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# get options
#-------------------------------------------------------------------------------
parser = optparse.OptionParser(
    usage="%prog -f file.PB.count [options]",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-f", action="store", type="string",
help="name of the file that contains PB counts")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--write-transfac", action="store_true", 
    dest = "write_transfac", help = "write intermediate transfac data to file")
optional_opts.add_option("--residue-lower", action="store", type="int",
    dest = "residue_lower", help="lower bound residue frame")
optional_opts.add_option("--residue-upper", action="store", type="int",
    dest = "residue_upper", help="upper bound for residue frame")
parser.add_option_group(optional_opts)
# get all parameters
(options, args) = parser.parse_args()

# check options
if not options.f:
    parser.print_help()
    parser.error("option -f is mandatory")

#-------------------------------------------------------------------------------
#
# convert a table of PB counts into transfac format that is required by weblogo
# http://meme.sdsc.edu/meme/doc/transfac-format.html
#
#-------------------------------------------------------------------------------

# read count file
#-------------------------------------------------------------------------------
if options.f and not os.path.isfile(options.f):
    sys.exit("%s does not appear to be a valid file" % (options.count_name))
f_in = open(options.f, 'r')
count_content = f_in.readlines()
f_in.close()

# convert to transfac format
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

# write transfac file (optional)
#-------------------------------------------------------------------------------
if options.write_transfac:
    transfac_name = options.f.replace('.count', '.transfac')
    f_out = open(transfac_name, 'w')
    f_out.write(transfac_content)
    f_out.close()
    print "wrote %s" % (transfac_name)

#-------------------------------------------------------------------------------
#
# call weblogo
#
#-------------------------------------------------------------------------------

lower = residue_lst[0]
if options.residue_lower:
    lower = options.residue_lower

upper = residue_lst[-1]
if options.residue_upper:
    upper = options.residue_upper

logo_name = options.f.replace(".count", ".logo.pdf")
if options.residue_lower or options.residue_upper:
    logo_name = options.f.replace(".count", ".logo.%i-%i.pdf" % (lower, upper))


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
""" % (options.f.replace(".PB.count", ""), logo_name, lower, upper)

proc = subprocess.Popen(command, shell = True,
stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
(out, err) = proc.communicate(transfac_content)
if err:
    print "ERROR:", err
code = proc.wait()
if code:
    print "ERROR: exit code != 0"
    print "exit code:", code
else:
    print "wrote %s" % logo_name

print out


