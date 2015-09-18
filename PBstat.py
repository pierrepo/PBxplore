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
import matplotlib
import matplotlib.pyplot as plt

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


def cmd_exists(cmd):
    """
    Check whether a command/program exists and can be executed.

    It uses the subprocess library to call the command.
    The return boolean is based on the ENOENT (No such file or directory) symbol.


    Parameters
    ----------
    cmd : list of str
        The list of the command and its arguments

    Returns
    -------
        bool
            True if the command can be executed, False otherwise.
    """
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc.communicate()
        proc.wait()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


class CommandAction(argparse.Action):
    """Handle command line argument which will call an external program.

    It means to check wether the external program is available and
    can be executed before running the main program. The check is
    performed through the `cmd_exists` function.

    It behaves the same way as the "store_true" action but it:
        - doesn't require the `default` keyword (always at False)
        - add a new keyword `command`
        - add a new optional keyword `cmd_help`

    Keyword argument
    ----------------
    command: str
        The name of the external program to call.
        It can contains arguments of the external program (like --help).
        The string will be split to satisfy the `cmd_exists` function.

    cmd_help: str (optional)
        A useful information to help the user to install the external program.
        It could be the right command line or a link to a documentation.
        This string will be part of the argparse.ArgumentError message.

    Raise
    -----
    ValueError: when the command value is not str.

    argparse.ArgumentError: when the external program can not be executed.
    This exception will be handle automatically by argparse during the
    parsing step.

    Examples
    --------
    parser = argparse.ArgumentParser()

    # External program with its own argument
    parser.add_argument("--logo", action=CommandAction,
                        command="weblogo --help", dest="logo")

    # External program without argument
    parser.add_argument("--ls", action=CommandAction, command="ls", dest="ls")

    # External program with a help message
    parser.add_argument("--logo", dest="logo", action=CommandAction,
                        command="weblogo --help",
                        cmd_help="Use 'pip install weblogo' to install")
    """

    def __init__(self, option_strings, command, dest, cmd_help=None, help=None):
        if not isinstance(command, str):
            raise ValueError("command must be a string")

        if cmd_help is not None and not isinstance(cmd_help, str):
            raise ValueError("cmd_help must be a string")

        self.command = command.split()
        self.cmd_help = cmd_help
        super(CommandAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            const=True,
            default=False,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        if not cmd_exists(self.command):
            msg="{0} was not found.".format(self.command[0])
            if self.cmd_help is not None:
                msg += " {0}".format(self.cmd_help)
            # Special exception so argparse can handle it himself
            # during the parsing step
            raise argparse.ArgumentError(self,msg)

        # Set with the self.const (which is for argument without value)
        # and not with values     (which is for argument with value)
        setattr(namespace, self.dest, self.const)


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
parser.add_argument("--logo", action=CommandAction, command="weblogo --help", dest="logo",
    cmd_help="(See https://github.com/pierrepo/PBxplore/blob/master/doc/installation.md)",
    help="generate logo representation of PBs frequency along protein sequence")
parser.add_argument("--image-format", action='store', type=str,
                    dest='image_format', default='png',
                    choices=['pdf', 'png', 'jpg'],
                    help='File format for all image output.')
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
if len(freq[0,:]) != (len(PB.NAMES) + 1):
    sys.exit( "ERROR: {0}: wrong data format".format(options.f) )

# read residues
residues = freq[:, 0]

#-------------------------------------------------------------------------------
# check / define residue min / max
if options.residue_min:
    residue_min = options.residue_min
else:
    residue_min = min(residues)

if options.residue_max:
    residue_max = options.residue_max
else:
    residue_max = max(residues)

if residue_min < 0:
    sys.exit("ERROR: residue_min must be >= 0")

if residue_max < 0:
    sys.exit("ERROR: residue_max must be >= 0")

if residue_min >= residue_max:
    sys.exit("ERROR: residue_min must be > residue_max")

if residue_min not in residues:
    sys.exit( "ERROR: residue_min does not belong to the residue range \
    in {0}".format(options.f) )

if residue_max not in residues:
    sys.exit( "ERROR: residue_max does not belong to the residue range \
    in {0}".format(options.f) )    


# get index of first residue
try:
    first_residue_index = int(residues[0])
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

# update residues
residues = freq[:, 0]

# remove residue numbers (first column)
freq = freq[:, 1:]
# normalize PBs frequencies
freq = freq / float(sequence_number)

#-------------------------------------------------------------------------------
# generates map of the distribution of PBs along protein sequence
#-------------------------------------------------------------------------------
if options.mapdist:

    # define output file name
    map_file_name = options.o + ".PB.map." + options.image_format
    if options.residue_min or options.residue_max:
        map_file_name = "{}.PB.map.{}-{}.{}".format(
                        options.o, residue_min, residue_max,
                        options.image_format)

    # define ticks for x-axis
    x_step = 5
    xticks = residues[::x_step]
    # trying to round ticks: 5, 10, 15 instead of 6, 11, 16...
    if xticks[0] == 1:
        xticks = xticks-1
        xticks[0] += 1

    yticks = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 
              'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

    fig = plt.figure(figsize=(2.0*math.log(len(residues)), 4))
    ax = fig.add_axes([0.1, 0.1, 0.75, 0.8])
    
    # Color scheme inspired from ColorBrewer 
    # http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=5
    # This color scheme is colorblind safe
    colors = [
    (44.0 / 255, 123.0 / 255, 182.0 / 255),
    (171.0 / 255, 217.0 / 255, 233.0 / 255),
    (255.0 / 255, 255.0 / 255, 191.0 / 255),
    (253.0 / 255, 174.0 / 255, 97.0 / 255), 
    (215.0 / 255, 25.0 / 255, 28.0 / 255)
    ]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('ColBrewerRdYlBu', colors)
    
    img = ax.imshow(numpy.transpose(freq[: , :]), interpolation='none', vmin=0, vmax=1, origin='lower', aspect='auto', cmap=cmap)
    
    ax.set_xticks( xticks - numpy.min(xticks) )
    ax.set_xticklabels(xticks)
    ax.set_yticks(range(len(yticks)))
    ax.set_yticklabels(yticks, style='italic', weight='bold')
    colorbar_ax = fig.add_axes([0.87, 0.1, 0.03, 0.8])
    fig.colorbar(img, cax=colorbar_ax)
    # print "beta-strand", "coil" and "alpha-helix" text 
    # only if there is more than 20 residues
    if len(residues) >= 20:
        # center alpha-helix: PB m (13th PB out of 16 PBs)
        # center coil: PB h and i (8th and 9th PBs out of 16 PBs)
        # center beta-sheet: PB d (4th PB out of 16 PBs) 
        fig.text(0.05, 4.0/16*0.8+0.075, r"$\beta$-strand", rotation = 90,  va = 'center', transform=ax.transAxes)
        fig.text(0.05, 8.5/16*0.8+0.075, r"coil", rotation = 90, va = 'center')
        fig.text(0.05, 13.0/16*0.8+0.075, r"$\alpha$-helix", rotation = 90, va = 'center', transform=ax.transAxes)
    fig.text(0.01, 0.5, "PBs", rotation = 90, weight='bold', size = 'larger', transform=ax.transAxes)
    fig.text(0.4, 0.01, "Residue number", weight="bold")
    fig.text(0.96, 0.6, "Intensity", rotation = 90, weight="bold")
    fig.savefig(map_file_name, dpi=300)
    print("wrote " + map_file_name)

#-------------------------------------------------------------------------------
# computes Neq and generates neq plot along protein sequence
#-------------------------------------------------------------------------------
if options.neq:
    # compute Neq
    #-------------------------------------------------------------------------------
    neq_array = numpy.zeros((len(residues), 2))
    neq_array[:, 0] = residues
    for idx in range(len(residues)):
        H = 0.0
        for b in range(len(PB.NAMES)):
            f = freq[idx, b] 
            if f != 0:
                H += f * math.log(f)
        neq_array[idx, 1] = math.exp(-H)
    
    # define output file name
    #-------------------------------------------------------------------------------
    neq_file_name = options.o + ".PB.Neq"
    if options.residue_min or options.residue_max:
        neq_file_name = "{}.PB.Neq.{}-{}".format(
                        options.o, residue_min, residue_max)
    
    # write Neq
    #-------------------------------------------------------------------------------
    content = "%-6s %8s \n" % ("resid", "Neq")
    for (res, neq) in neq_array:
        content += "%-6d %8.2f \n" % (res, neq)
    neq_file = open(neq_file_name, "w")
    neq_file.write(content)
    neq_file.close()
    print( "wrote {0}".format(neq_file_name) )
    
    # draw Neq
    #-------------------------------------------------------------------------------
    neq_fig_name = neq_file_name + '.' + options.image_format
    fig = plt.figure(figsize=(2.0*math.log(len(residues)), 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim([0, round(max(neq_array[:, 1]), 0)+1])
    ax.plot(neq_array[:, 0], neq_array[:, 1])
    ax.set_xlabel('Residue number', fontsize=18)
    ax.set_ylabel('Neq', fontsize=18, style='italic')
    fig.savefig(neq_fig_name)
    print("wrote {}".format(neq_fig_name))


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
    transfac_content = PB.count_to_transfac(options.f, count_content)

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
    logo_file_name = options.o + ".PB.logo." + options.image_format
    if options.residue_min or options.residue_max:
        logo_file_name = "{}.PB.logo.{}-{}.{}".format(
            options.o, residue_min, residue_max, options.image_format)

    # call weblogo
    #-------------------------------------------------------------------------------
    # If the output file format is 'jpg', then the --format argument of
    # weblogo should be 'jpeg'.
    logo_format = options.image_format
    if logo_format == 'jpg':
        logo_format = 'jpeg'
    command = """weblogo \
--format %s \
--errorbars NO \
--fineprint "" \
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
    """ % (logo_format, options.f.replace(".PB.count", ""),
           logo_file_name, residue_min, residue_max)

    proc = subprocess.Popen(command, shell = True,
    stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    (out, err) = proc.communicate(transfac_content.encode('UTF-8'))
    if err:
        print( "ERROR: {0}".format(err) )
    code = proc.wait()
    if code:
        print( "ERROR: exit code != 0" )
        print( "exit code: {0}".format(code) )
    else:
        print( "wrote {0}".format(logo_file_name) )
    print( out )

