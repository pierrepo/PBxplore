#!/usr/bin/env python

"""
Regression tests for PBxplore.

This script run the various PBxplore programs with various argument, and makes
sure the output is the expected one. The aim is to check that the programs are
not broken during development.

Be careful this script does not test that the output is right. It just test that
the output is the expected one based on previous a previous version.

To run this test suite, you can either run this script without arguments or use
the nose [1]_ package. The latter option gives a more readable output, as stdout
is captured, and is displayed only if a test fails.

.. [1] https://nose.readthedocs.org
"""

# Use print as a function like in python 3
from __future__ import print_function

from unittest import TestCase, main
from os import path
from uuid import uuid1
import os
import subprocess
import shutil

# Resources for the tests are stored in the following directory
REFDIR = "test-resources/"

# Temporary outputs will be written in this directory
OUTDIR = "test-outputs/"


class TestPBAssign(TestCase):
    """
    Regression tests for PBAssign.py
    """
    def setUp(self):
        """
        Run before each test.

        Make sure that the output directory exists.
        """
        if not path.isdir(OUTDIR):
            os.mkdir(OUTDIR)

    def test_fasta(self):
        """
        Run PBAssign on PDB files, and check the fasta output.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        _test_PBassign_options(references, ['{0}.PB.fasta'], [])


def _same_file_content(file_a, file_b):
    """
    Return True if two files are identical. Take file path as arguments.
    """
    with open(file_a) as f1, open(file_b) as f2:
        # Compare content line by line
        for f1_line, f2_line in zip(f1, f2):
            if f1_line != f2_line:
                return False
        # Check if one file is longer than the other; it would result as one
        # file iterator not completely consumed
        for infile in (f1, f2):
            try:
                f1.readline()
            except StopIteration:
                # The iterator is consumed
                pass
            else:
                # The iterator is not consumed, it means that this file is
                # longer than the other
                return False
    # If we reach this line, it means that we did not find any difference
    return True


def _assert_identical_files(file_a, file_b, message=''):
    """
    Raise an Assert exception if the two files are not identical.

    Take file path as arguments.
    """
    assert not _same_file_content(file_a, file_b), message


def _run_prog(program, pdbid, options, indir=REFDIR, outdir=OUTDIR):
    """
    Run a PBxplore program on a PDBID with the given options.

    options is expected to be a list that will be directly pass to subprocess,
    it must not contain the input or output options.
    """
    if program not in ('./PBassign.py',):
        raise NotImplementedError('_run_prog does not know how to run {0}'\
                                  .format(program))
    out_run_dir = path.join(OUTDIR, str(uuid1()))
    test_input = path.join(REFDIR, pdbid + '.pdb')
    out_basename = path.join(out_run_dir, pdbid)

    os.mkdir(out_run_dir)

    run_list = [program, '-p', test_input, '-o', out_basename] + options
    print(' '.join(run_list))
    status = subprocess.call(run_list, stdout=subprocess.PIPE)

    return status, out_run_dir


def _test_PBassign_options(basenames, outfiles, options, expected_exit = 0):
    for basename in basenames:
        status, out_run_dir = _run_prog('./PBassign.py', basename, options)
        assert status == expected_exit, \
               'PBassign stoped with a {0} exit code'.format(status)
        assert len(os.listdir(out_run_dir)) == len(outfiles),\
                ('PBassign did not produced the right number of files: '
                 '{0} files produced instead of {1}').format(
                    len(os.listdir(out_run_dir)), len(outfiles))
        for outfile in (template.format(basename) for template in outfiles):
            test_file = path.join(out_run_dir, outfile)
            ref_file = path.join(REFDIR, outfile)
            _assert_identical_files(test_file, ref_file)
        shutil.rmtree(out_run_dir)


if __name__ == '__main__':
    main()
