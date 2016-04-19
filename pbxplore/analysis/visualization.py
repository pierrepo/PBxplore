#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
import os
import math

# Third-party modules
import numpy

# Conditional imports
try:
    import matplotlib
    import matplotlib.pyplot as plt
except ImportError:
    IS_MATPLOTLIB = False
else:
    IS_MATPLOTLIB = True

try:
    import weblogolib
except ImportError:
    IS_WEBLOGO = False
else:
    IS_WEBLOGO = True

# Local modules
from .. import PB
from . import utils


# Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass

# Create the __all__ keyword according to the conditional imports
__all__ = []
if IS_MATPLOTLIB:
    __all__ += ['plot_neq', 'plot_map']
if IS_WEBLOGO:
    __all__ += ['generate_weblogo']


def plot_neq(fname, neq_array, residue_min=1, residue_max=None):
    """
    Generate the Neq plot along the protein sequence

    Parameters
    ----------
    neq_array : numpy array
        an array containing the neq value associated to the residue number
    fname : str
        The path to the file to write in
    residue_min: int
        the lower bound of the protein sequence
    residue_max: int
        the upper bound of the protein sequence
    """

    neq = utils._slice_matrix(neq_array, residue_min, residue_max)
    nb_residues = neq.shape[0]

    # Residue number with good offset given the slice
    x = numpy.arange(residue_min, residue_min + nb_residues)

    fig = plt.figure(figsize=(2.0*math.log(nb_residues), 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim([0, round(max(neq), 0)+1])
    ax.plot(x, neq)
    ax.set_xlabel('Residue number', fontsize=18)
    ax.set_ylabel('Neq', fontsize=18, style='italic')
    fig.savefig(fname)


def plot_map(fname, count_mat, residue_min=1, residue_max=None):
    """
    Generate a map of the distribution of PBs along protein sequence from
    an occurence matrix.

    Parameters
    ----------
    fname : str
        The path to the file to write in
    count_mat : numpy array
        an occurence matrix returned by `count_matrix`.
    residue_min: int
        the lower bound of the protein sequence
    residue_max: int
        the upper bound of the protein sequence
    """
    # Get the frequency matrix
    freq_mat = utils.compute_freq_matrix(count_mat)
    # Slice it
    freq = utils._slice_matrix(freq_mat, residue_min, residue_max)
    nb_residues = freq.shape[0]
    # Residue number with good offset given the slice
    x = numpy.arange(residue_min, residue_min + nb_residues)

    # Define a scaling factor to handle nice rendering for small and large proteins
    # This is empirical!
    scaling_factor = math.log(nb_residues)

    # define ticks for x-axis
    x_step = 5
    # space ticks for large proteins
    if nb_residues > 100:
        x_step = 10
    if nb_residues > 200:
        x_step = int( scaling_factor) * 5
    xticks = x[::x_step]
    # trying to round ticks: 5, 10, 15 instead of 6, 11, 16...
    if xticks[0] == 1:
        xticks = xticks-1
        xticks[0] += 1

    # define ticks for y-axis
    yticks = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
              'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

    fig = plt.figure(figsize=(2.0 * scaling_factor, 4))
    ax = fig.add_axes([0.1, 0.1, 0.75, 0.8])

    # Color scheme inspired from ColorBrewer
    # http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=5
    # This color scheme is colorblind safe
    colors = [(44.0 / 255, 123.0 / 255, 182.0 / 255),
              (171.0 / 255, 217.0 / 255, 233.0 / 255),
              (255.0 / 255, 255.0 / 255, 191.0 / 255),
              (253.0 / 255, 174.0 / 255, 97.0 / 255),
              (215.0 / 255, 25.0 / 255, 28.0 / 255)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('ColBrewerRdYlBu', colors)

    img = ax.imshow(numpy.transpose(freq[:, :]), interpolation='none', vmin=0, vmax=1,
                    origin='lower', aspect='auto', cmap=cmap)

    ax.set_xticks(xticks - numpy.min(xticks))
    ax.set_xticklabels(xticks)
    ax.set_yticks(range(len(yticks)))
    ax.set_yticklabels(yticks, style='italic', weight='bold')
    colorbar_ax = fig.add_axes([0.87, 0.1, 0.03, 0.8])
    fig.colorbar(img, cax=colorbar_ax)
    # print "beta-strand", "coil" and "alpha-helix" text
    # only if there is more than 20 residues
    if nb_residues >= 20:
        text_shift = 0.0
        if nb_residues >= 100:
            text_shift = scaling_factor / 500
        # center alpha-helix: PB m (13th PB out of 16 PBs)
        # center coil: PB h and i (8th and 9th PBs out of 16 PBs)
        # center beta-sheet: PB d (4th PB out of 16 PBs)
        fig.text(0.05 + text_shift, 4.0/16*0.8+0.075, r"$\beta$-strand", rotation=90,
                 va='center', transform=ax.transAxes)
        fig.text(0.05 + text_shift, 8.5/16*0.8+0.075, r"coil", rotation=90,
                 va='center')
        fig.text(0.05 + text_shift, 13.0/16*0.8+0.075, r"$\alpha$-helix", rotation=90,
                 va='center', transform=ax.transAxes)

    fig.text(0.01 + text_shift * 2, 0.5, "PBs", rotation=90, weight="bold",
             size='larger', transform=ax.transAxes)
    fig.text(0.4, 0.01, "Residue number", weight="bold")
    fig.text(0.96 - text_shift, 0.6, "Intensity", rotation=90, weight="bold")
    fig.savefig(fname, dpi=300)


def generate_weblogo(fname, count_mat, residue_min=1, residue_max=None, title=""):
    """
    Generates logo representation of PBs frequency along protein sequence through
    the weblogo library.

    The weblogo reference:
    G. E. Crooks, G. Hon, J.-M. Chandonia, and S. E. Brenner.
    'WebLogo: A Sequence Logo Generator.'
    Genome Research 14:1188–90 (2004)
    doi:10.1101/gr.849004.
    http://weblogo.threeplusone.com/

    Parameters
    ----------
    fname : str
        The path to the file to write in
    count_mat : numpy array
        an occurence matrix returned by `count_matrix`.
    residue_min: int
        the lower bound of residue frame
    residue_max: int
        the upper bound of residue frame
    title: str
        the title of the weblogo. Default is empty.
    """

    # Slice the matrix
    count = utils._slice_matrix(count_mat, residue_min, residue_max)

    # Create a custom color scheme for PB
    colors = weblogolib.ColorScheme([weblogolib.ColorGroup("d", "#1240AB", "strand main"),
                                     weblogolib.ColorGroup("abcdef", "#1240AB", "strand others"),
                                     weblogolib.ColorGroup("ghij", "#0BD500", "coil"),
                                     weblogolib.ColorGroup("m",  "#FD0006", "helix main"),
                                     weblogolib.ColorGroup("klnop",  "#FD0006", "helix others")])

    # Load data from an occurence matrix
    data = weblogolib.LogoData.from_counts(PB.NAMES, count)

    # Create options
    options = weblogolib.LogoOptions(fineprint=False, logo_title=title, color_scheme=colors,
                                     stack_width=weblogolib.std_sizes["large"],
                                     first_residue=residue_min)

    # Generate weblogo
    logo = weblogolib.LogoFormat(data, options)

    # Retrieve image format
    image_format = os.path.splitext(fname)[1][1:]

    # Retrieve the right function given the image format
    try:
        if image_format == 'jpg':
            image_format = 'jpeg'
        formatter = weblogolib.formatters[image_format]
    except KeyError:
        raise ValueError("Invalid format image '{0}'."
                         " Valid ones are : eps, png, pdf, jpg/jpeg, svg".format(image_format))
    # Format the logo
    image = formatter(data, logo)

    # Write it
    with open(fname, "w") as f:
        print(image, file=f)
