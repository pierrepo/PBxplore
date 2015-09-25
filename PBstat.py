#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Statistical analysis and graphical representations of PBs.

Compute Neq, PBs distribution and draw logo representation of PBs.

2013 - P. Poulain, A. G. de Brevern
"""

# ============================================================================
# Modules
# ============================================================================
# Use print as a function for python 3 compatibility
from __future__ import print_function

# Standard modules
import os
import sys
import math
import subprocess
import argparse

# Third-party module
import numpy
import matplotlib
import matplotlib.pyplot as plt

# Local module
import PBlib as PB

# ============================================================================
# Python2/Python3 compatibility
# ============================================================================

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
            msg = "{0} was not found.".format(self.command[0])
            if self.cmd_help is not None:
                msg += " {0}".format(self.cmd_help)
            # Special exception so argparse can handle it himself
            # during the parsing step
            raise argparse.ArgumentError(self, msg)

        # Set with the self.const (which is for argument without value)
        # and not with values     (which is for argument with value)
        setattr(namespace, self.dest, self.const)


def user_inputs():
    """
    Handle the user parameters for PBstat.py.

    Returns
    -------
    options : the parsed arguments as parsed by `argparse`.
    """

    parser = argparse.ArgumentParser(
        description="Statistical analysis and graphical representations of PBs.")

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

    # check file
    if not os.path.isfile(options.f):
        parser.error("{0}: not a valid file".format(options.f))

    # Check residues min/max
    residues = [options.residue_min, options.residue_max]
    for residue in residues:
        if residue is not None and residue < 0:
            parser.error("residue argument must be >=0")
    if None not in residues and options.residue_min >= options.residue_max:
        parser.error("residue-min must be < residue-max.")

    return options


def generate_weblogo(count_file, residue_min, residue_max, title,
                     logo_format, outfile, debug=False):
    """
    Generates logo representation of PBs frequency along protein sequence through
    the weblogo binary

    Parameters
    ----------
    count_file: str
        The file containing the PB frequency
    residue_min: int
        the lower bound of residue frame
    residue_max: int
        the upper bound of residue frame
    title: str
        the title of the weblogo
    logo_format: str
        the format of the image output
    outfile : str
        The path to the file to write in
    debug: bool
        If True, write the transfac object into a file
    """

    with open(count_file, 'r') as f_in:
        count_content = f_in.readlines()

    # convert a table of PB frequencies into transfac format as required by weblogo
    # http://meme.sdsc.edu/meme/doc/transfac-format.html
    transfac_content = PB.count_to_transfac(count_file, count_content)

    # write transfac file (debug only)
    if debug:
        # change extension
        transfac_name = os.path.splitext(outfile)[0] + ".transfac"
        with open(transfac_name, 'w') as f_out:
            f_out.write(transfac_content)
        print("wrote {0}".format(transfac_name))

    # If the output file format is 'jpg', then the --format argument of
    # weblogo should be 'jpeg'.
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
              """ % (logo_format, title, outfile, residue_min, residue_max)

    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    (out, err) = proc.communicate(transfac_content.encode('UTF-8'))
    if err:
        print("ERROR: {0}".format(err))
    code = proc.wait()
    if code:
        print("ERROR: exit code != 0")
        print("exit code: {0}".format(code))
    else:
        print("wrote {0}".format(outfile))
    print(out)


def read_occurence_file(name):
    """
    Read an occurence matrix from a file.
    It will return the matrix as a numpy array and the indexes of residues.

    Parameters
    ----------
    name : str
        Name of the file.

    Returns
    -------
    count_mat : numpy array
        the occurence matrix without the residue number
    residues: list
        the list of residues indexes

    Exceptions
    ----------
    ValueError : when something is wrong about the file
    """

    # load count file
    # skip first row that contains PBs labels
    try:
        count = numpy.loadtxt(name, dtype=int, skiprows=1)
    except:
        raise ValueError("ERROR: {0}: wrong data format".format(name))

    # determine number of sequences compiled
    # use the sum of all residue at position 3
    # since positions 1 and 2 have no PBs assignement
    # and begin at 1 to not sum the index of the line (here is 3)
    sequence_number = sum(count[2, 1:])
    if sequence_number == 0:
        raise ValueError("ERROR: counting 0 sequences!")

    # read residues number
    residues = count[:, 0]
    # remove residue numbers (first column)
    count = count[:, 1:]

    # get index of first residue
    try:
        int(residues[0])
    except:
        raise ValueError("""ERROR: cannot read index of first residue.
                         Wrong data format in {0}""".format(name))

    return count, residues


def check_residue_range(residues, residue_min, residue_max):
    """"
    Ensure that the lower bound and the upper bound parameters are in the range of
    the list `residues`.

    Parameters
    ----------
    residues : list
        the list of residues indexes
    residue_min:
        the lower bound of residue
    residue_max:
        the upper bound of residue

    Exceptions
    ----------
    IndexError : if `residue_min` or `residue_max` is not in the range
    """

    if residue_min is None:
        residue_min = residues[0]

    if residue_max is None:
        residue_max = residues[-1]

    if residue_min not in residues:
        raise IndexError("residue_min does not belong to the residue range")

    if residue_max not in residues:
        raise IndexError("residue_max does not belong to the residue range")

    if residue_min >= residue_max:
        raise IndexError("Lower bound > upper bound")

    return residue_min, residue_max


def pbstat_cli():
    """
    PBstat command line.
    """

    options = user_inputs()

    try:
        count, residues = read_occurence_file(options.f)
        residue_min, residue_max = check_residue_range(residues,
                                                       options.residue_min, options.residue_max)
    except Exception as e:
        sys.exit("ERROR: {0}".format(e))

    print("Index of first residue is: {0}".format(residue_min))

    # Handle output file name...
    output_file_name = options.o + ".PB.{0}"
    if options.residue_min or options.residue_max:
        output_file_name = "{0}.{1}-{2}".format(output_file_name, residue_min, residue_max)
    # ... and figure name
    output_fig_name = output_file_name + "." + options.image_format

    # -------------------------------------------------------------------------------
    # generates map of the distribution of PBs along protein sequence
    # -------------------------------------------------------------------------------
    if options.mapdist:
        file_fig_name = output_fig_name.format("map")
        PB.plot_map(file_fig_name, count, residue_min, residue_max)
        print("wrote " + file_fig_name)

    # -------------------------------------------------------------------------------
    # computes Neq and generates neq plot along protein sequence
    # -------------------------------------------------------------------------------
    if options.neq:
        # compute Neq
        neq = PB.compute_neq(count)

        # write Neq
        neq_file_name = output_file_name.format("Neq")
        with open(neq_file_name, "w") as outfile:
            PB.write_neq(outfile, neq, residue_min, residue_max)
        print("wrote {0}".format(neq_file_name))

        # draw Neq
        file_fig_name = output_fig_name.format("Neq")
        PB.plot_neq(file_fig_name, neq, residue_min, residue_max)
        print("wrote {}".format(file_fig_name))

    # -------------------------------------------------------------------------------
    # generates logo representation of PBs frequency along protein sequence
    #
    # G. E. Crooks, G. Hon, J.-M. Chandonia, and S. E. Brenner.
    # 'WebLogo: A Sequence Logo Generator.'
    # Genome Research 14:1188–90 (2004)
    # doi:10.1101/gr.849004.
    # http://weblogo.threeplusone.com/
    #  -------------------------------------------------------------------------------
    if options.logo:
        title = options.f.replace(".PB.count", "")
        generate_weblogo(options.f, residue_min, residue_max, title,
                         options.image_format, output_fig_name.format("logo"))


if __name__ == '__main__':
    pbstat_cli()
