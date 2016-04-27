#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Imports for the clustering
import random
import sys
import numpy
import collections
import copy

from ..PB import NAMES

def count_per_position_partial(sequences, seq_indices):
    """
    Count the population of each block at each position.

    Parameters
    ----------
    sequences
        a list of sequences
    seq_indices
        a list of indices; only the sequences from `sequences`
        which have their index in `seq_indices` will be taken into
        account.

    Returns
    -------
    counts: collections.defaultdict
        a defaultdict with the blocks as keys and their population as
        values; the default value is a 0 interger.
    """
    counts = [collections.defaultdict(int)
              for _ in range(len(sequences[0]))]
    for seq_idx in seq_indices:
        seq = sequences[seq_idx]
        for block, count in zip(seq, counts):
            count[block] += 1
    return counts


def make_profile_partial(sequences, seq_indices):
    """
    Calculate a frequency profile from a list of sequences.

    A profile is the frequency for each block to be a each position.

    Parameters
    ----------
    sequences
        a list of sequences
    seq_indices
        a list of indices; only the sequences from `sequences`
        which have their index in `seq_indices` will be taken into
        account.

    Returns
    -------
    profile: numpy array
        the frequency profile as a numpy array with a row for each
        block and a column for each position.
    """
    counts = count_per_position_partial(sequences, seq_indices)
    profile = numpy.zeros((16, len(counts)))
    for pos_idx, position in enumerate(counts):
        for block_idx, block in enumerate(NAMES):
            profile[block_idx, pos_idx] = position[block]
    profile /= len(seq_indices)
    return profile


def compatibility(profile, sequence):
    """
    Compute the compatibility of a sequence with a profile.

    Parameters
    ----------
    profile
        a frequency profile as a numpy array with a row for each
        block and a column for each position
    sequence: the block sequence to test as a string

    Returns
    -------
    probability: float
        cummulative probability of the given sequence given the profile
    """
    probabilities = numpy.zeros((profile.shape[1], ))
    for pos_idx, block in enumerate(sequence):
        block_idx = NAMES.find(block)
        probabilities[pos_idx] = profile[block_idx, pos_idx]
    return numpy.sum(probabilities)


def _argmax(values):
    """
    Return the index of the maximum value of a sequence
    """
    iter_values = iter(values)
    argmax = 0
    valmax = next(iter_values)
    for arg, val in enumerate(iter_values, start=1):
        if val > valmax:
            valmax = val
            argmax = arg
    return argmax


def attribute_center(centers, sequences):
    """
    Assign the more compatible center to each sequence.

    For a list of frequency profiles and a list of sequences, assign to
    each sequence the more compatible frequency profile.

    Parameters
    ----------
    centers
        a list of frequency profiles, each profile is a numpy array
        with rows for the blocks and collumns for the positions
    sequences: a list of block sequences

    Returns
    -------
    groups: list
        a list of the profile index for each sequence
    """
    groups = [0 for _ in range(len(sequences))]
    for seq_idx, sequence in enumerate(sequences):
        compatibilities = [compatibility(center, sequence)
                           for center in centers]
        groups[seq_idx] = _argmax(compatibilities)
    return groups


def update_centers(sequences, groups, ngroups):
    """
    Calculate the frequency profile for each group

    Parameters
    ----------
    sequences
        a list of block sequences
    groups
        a list with a group identifier for each sequence
    ngroups
        the number of groups; this number cannot be obtained with
        enough confidence from the `groups` list, indeed some group can
        be empty and therefore not appearing in the `group` list

    Returns
    -------
    centers: numpy array
        a frequency profile for each group
    """
    grouped_indices = [[] for _ in range(ngroups)]
    for seq_idx, group_idx in enumerate(groups):
        grouped_indices[group_idx].append(seq_idx)
    centers = []
    for group in grouped_indices:
        centers.append(make_profile_partial(sequences, group))
    return centers


def get_medoids(groups):
    raise NotImplemented


def initial_centers(sequences, ngroups):
    """
    Create single sequence frequency profiles from randomly chosen sequences

    Choose `ngroups` sequences and build a single sequence frequency
    profile for each of them. These profiles aims at being initial centers
    for the clustering.

    Parameters
    ----------
    sequences
        a list of block sequences
    ngroups
        the number of profiles to build

    Returns
    -------
    centers: list
        a list of frequency profiles
    """
    centers = []
    sample = random.sample(range(len(sequences)), ngroups)
    assert len(sample) == len(set(sample)), \
        'Redundances in the initial sampling'
    assert len([sequences[i] for i in sample]) \
        == len(set([sequences[i] for i in sample])), \
        'Redundances in the initial sampling'
    for sequence_idx in sample:
        centers.append(make_profile_partial(sequences, [sequence_idx, ]))
    return centers


def _is_converged(old_centers, new_centers):
    """
    Returns True is all `new_centers` are the same as the `old_centers`
    """
    for old, new in zip(old_centers, new_centers):
        if not numpy.all(old == new):
            return False
    return True

def k_means(sequences, ngroups, max_iter, logfile=sys.stdout):
    """
    Carry out a K-means clustering on block sequences

    Cluster the sequences into `ngroups` groups. The centers of each
    group are computed as a frequency profile for the sequences in the
    group. The distance metrix used to assign a center to the sequences
    is the probability to obtain a sequence from a given frequency profile.

    Parameters
    ----------
    sequences
        a list of block sequences to cluster
    ngroups
        the number of cluster to build
    max_iter
        the maximum number of iterations to run
    logfile
        a file descriptor where to write logs, stdout by default

    Returns
    -------
    groups: list
        a list of cluster identifier for each sequence
    convergeance: bool
        True if the clustering converged, else False
    """
    convergeance = False
    centers = initial_centers(sequences, ngroups)
    for iteration in range(1, max_iter + 1):
        groups = attribute_center(centers, sequences)
        print(iteration, collections.Counter(groups), len(centers),
              file=logfile)
        new_centers = update_centers(sequences, groups, ngroups)
        if _is_converged(centers, new_centers):
            print('Convergence reached in {} iterations'
                  .format(iteration), file=logfile)
            convergeance = True
            break
        centers = copy.copy(new_centers)
    else:
        print(('K-means reached {} iterations before '
               'reaching convergence.').format(max_iter),
              file=logfile)
    return groups, centers, convergeance
