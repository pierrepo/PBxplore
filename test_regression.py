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
import os
import subprocess

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
        for pdbid in references:
            test_input = path.join(REFDIR, pdbid + '.pdb')
            out_basename = path.join(OUTDIR, pdbid)
            test_output = path.join(OUTDIR, pdbid + '.PB.fasta')
            reference_output = path.join(REFDIR, pdbid + '.PB.fasta')
            # We want to be sure the output does not already exists
            if path.exists(test_output):
                raise RuntimeError('{0} should not exists before the test '
                                   'actually run'.format(test_output))
            # Run the program
            run_list = ['./PBassign.py', '-p', test_input, '-o', out_basename]
            print(' '.join(run_list))
            subprocess.call(run_list, stdout=subprocess.PIPE)
            
            # Test that the output is identical to the expected one
            self.assertFalse(_same_file_content(test_output, reference_output),
                             "{0} and {1} should be identical".format(
                                 test_output, reference_output))

            # Clean the output
            os.remove(test_output)



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


if __name__ == '__main__':
    main()
