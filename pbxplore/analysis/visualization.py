#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
import os
import math

# Third-party modules
import numpy

import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

try:
    import weblogolib
    # Weblogo compatibility
    # With version < 3.5, the color class is 'ColorGroup'. In version >= 3.5,
    # it is 'SymbolColor'. Here, we change to always have 'ColorGroup'.
    try:
        ColorGroup = weblogolib.SymbolColor
    except NameError:
        ColorGroup = weblogolib.ColorGroup
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
__all__ = ['plot_neq', 'plot_map']
if IS_WEBLOGO:
    __all__ += ['generate_weblogo']


def plot_neq(fname, neq_array, idx_first_residue=1, residue_min=1, residue_max=None):
    """
    Generate the Neq plot along the protein sequence

    Parameters
    ----------
    fname : str
        The path to the file to write in
    neq_array : numpy array
        an array containing the neq value associated to the residue number
    idx_first_residue: int
        the index of the first residue in the array
    residue_min: int
        the lower bound of the protein sequence
    residue_max: int
        the upper bound of the protein sequence
    """

    neq = utils._slice_matrix(neq_array, idx_first_residue, residue_min, residue_max)
    nb_residues = neq.shape[0]

    # Residue number with good offset given the slice
    x = numpy.arange(residue_min, residue_min + nb_residues)

    fig = plt.figure(figsize=(2.0 * math.log(nb_residues), 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim([0, round(max(neq), 0) + 1])
    ax.plot(x, neq)
    ax.set_xlabel('Residue number', fontsize=18)
    ax.set_ylabel('Neq', fontsize=18, style='italic')
    fig.savefig(fname)


def plot_map(fname, count_mat, idx_first_residue=1, residue_min=1, residue_max=None):
    """
    Generate a map of the distribution of PBs along protein sequence from
    an occurence matrix.

    Parameters
    ----------
    fname : str
        The path to the file to write in
    count_mat : numpy array
        an occurence matrix returned by `count_matrix`.
    idx_first_residue: int
        the index of the first residue in the matrix
    residue_min: int
        the lower bound of the protein sequence
    residue_max: int
        the upper bound of the protein sequence
    """
    # Get the frequency matrix
    freq_mat = utils.compute_freq_matrix(count_mat)
    # take a slice with min/max residues
    freq = utils._slice_matrix(freq_mat, idx_first_residue, residue_min, residue_max)
    nb_residues = freq.shape[0]
    # Residue number with good offset given the slice
    x = numpy.arange(residue_min, residue_min + nb_residues)

    # Define a scaling factor to handle nice rendering for small and large proteins
    # This is empirical!
    scaling_factor = math.log(nb_residues)

    # define ticks for x-axis
    x_step = 1
    # space ticks for large proteins
    if nb_residues > 12:
        x_step = 2
    if nb_residues > 25:
        x_step = 5
    if nb_residues > 100:
        x_step = 10
    if nb_residues > 200:
        x_step = int(scaling_factor) * 5
    xticks = x[::x_step]
    # trying to round ticks: 5, 10, 15 instead of 6, 11, 16...
    if xticks[0] == 1:
        xticks = xticks - 1
        xticks[0] += 1

    # define ticks for y-axis
    yticks = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
              'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

    fig = plt.figure(figsize=(2.0 * scaling_factor, 4))

    gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1, 2.0 * scaling_factor, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    # Color scheme inspired from ColorBrewer
    # http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=5
    # This color scheme is colorblind safe
    colors = [(44.0 / 255, 123.0 / 255, 182.0 / 255),
              (171.0 / 255, 217.0 / 255, 233.0 / 255),
              (255.0 / 255, 255.0 / 255, 191.0 / 255),
              (253.0 / 255, 174.0 / 255, 97.0 / 255),
              (215.0 / 255, 25.0 / 255, 28.0 / 255)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('ColBrewerRdYlBu', colors)

    img = ax2.imshow(numpy.transpose(freq[:, :]), interpolation='none', vmin=0, vmax=1,
                     origin='lower', aspect='auto', cmap=cmap)

    # add colorbar
    divider2 = make_axes_locatable(ax2)
    # cax2 = divider2.append_axes("right", size="5%", pad=0.08)
    cax2 = divider2.append_axes("right", size=0.15, pad=0.08)
    plt.colorbar(img, cax=cax2)

    # add ticks and labels
    ax2.set_xticks(xticks - numpy.min(xticks))
    ax2.set_xticklabels(xticks)
    ax2.set_yticks(range(len(yticks)))
    ax2.set_yticklabels(yticks, style='italic', weight='bold')
    ax2.set_xlabel("Residue number", weight="bold")

    # add secondary structures
    ax1.set_axis_off()
    ax1.text(0.2, 0.5, "PBs", rotation=90, weight="bold",
                       size='larger', ha='center', va='center')
    # center alpha-helix: PB m (13th PB out of 16 PBs)
    # center coil: PB h and i (8th and 9th PBs out of 16 PBs)
    # center beta-sheet: PB d (4th PB out of 16 PBs)
    ax1.text(0.5, 4.0 / 16, r"$\beta$-strand", rotation=90,
             ha='center', va='center')
    ax1.text(0.5, 8.5 / 16, r"coil", rotation=90,
             ha='center', va='center')
    ax1.text(0.5, 13.0 / 16, r"$\alpha$-helix", rotation=90,
             ha='center', va='center')

    # add "intensity"
    ax3.set_axis_off()
    ax3.text(0.7, 0.5, "Intensity", rotation=90, weight="bold", ha='center', va='center')

    # adjust margins and save
    fig.subplots_adjust(left=0.01, bottom=0.12, top=0.96, right=0.99, wspace=0)
    fig.savefig(fname, dpi=300)


def generate_weblogo(fname, count_mat, idx_first_residue=1, residue_min=1, residue_max=None, title=""):
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
    idx_first_residue: int
        the index of the first residue in the matrix
    residue_min: int
        the lower bound of residue frame
    residue_max: int
        the upper bound of residue frame
    title: str
        the title of the weblogo. Default is empty.
    """

    # Slice the matrix
    count = utils._slice_matrix(count_mat, idx_first_residue, residue_min, residue_max)

    # Create a custom color scheme for PB
    colors = weblogolib.ColorScheme([ColorGroup("d", "#1240AB", "strand main"),
                                     ColorGroup("abcdef", "#1240AB", "strand others"),
                                     ColorGroup("ghij", "#0BD500", "coil"),
                                     ColorGroup("m", "#FD0006", "helix main"),
                                     ColorGroup("klnop", "#FD0006", "helix others")])

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
    with open(fname, "wb") as f:
        f.write(image)
