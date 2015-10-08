#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard modules
import os
import textwrap

FASTA_WIDTH = 60  # line width for fasta format


def read_fasta(name):
    """
    Read fasta file and output sequences in a list.

    Parameters
    ----------
    name : str
        Name of file containing sequences in fasta format.

    Returns
    -------
    header_lst : list
        List of headers (str)
    sequence_lst : list
        List of sequences (str)
    """
    assert os.path.exists(name), name + ' does not exist'
    sequence_lst = []
    header_lst = []
    header = ""
    sequence = ""
    with open(name, "rt") as f_in:
        for line in f_in:
            data = line.strip()
            # jump empty lines
            if not data:
                continue
            # store header and sequence when a new header
            # (i.e. sequence) is found
            if sequence and header and data.startswith(">"):
                header_lst.append(header)
                sequence_lst.append(sequence)
                # reset header and sequence
                header = ""
                sequence = ""
            # save header of sequence
            if data.startswith(">"):
                header = data[1:]
            # save sequence
            if ">" not in data:
                sequence += data
    # save last sequence
    if header and sequence:
        header_lst.append(header)
        sequence_lst.append(sequence)
    # outputs
    assert len(header_lst) == len(sequence_lst), \
        "cannot read same number of headers and sequences"
    print("read %d sequences in %s" % (len(sequence_lst), name))
    if len(sequence_lst) == 0:
        print("WARNING: {} seems empty of sequence".format(name))
    return header_lst, sequence_lst


def read_several_fasta(input_files):
    """
    Read several fasta files

    Note that each fasta file may contain several sequences.

    Parameters
    ----------
    input_files: a list of fasta file paths.

    Returns
    -------
    pb_name: a list of the headers
    pb_seq: a list of the sequences
    """
    pb_seq = []
    pb_name = []
    for name in input_files:
        header, seq = read_fasta(name)
        pb_name += header
        pb_seq += seq
    return pb_name, pb_seq


def write_fasta_entry(outfile, sequence, comment, width=FASTA_WIDTH):
    """
    Write a fasta entry (header + sequence) in an open file

    Parameters
    ----------
    name : file descriptor
        The file descriptor to write in. It must allow writing.
    sequence : str
        Sequence to format.
    comment : str
        Comment to make header of sequence.
    width : int
        The width of a line. `FASTA_WIDTH` by default.
    """
    print('>' + comment, file=outfile)
    print(textwrap.fill(sequence, width=width), file=outfile)


def write_fasta(outfile, sequences, comments):
    """
    Write fasta entries (header + sequence) in an open file

    Parameters
    ----------
    name : file descriptor
        The file descriptor to write in. It must allow writing.
    header_lst : list
        List of headers (str)
    sequence_lst : list
        List of sequences (str)
    """
    for sequence, comment in zip(sequences, comments):
        write_fasta_entry(outfile, sequence, comment)
