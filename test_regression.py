#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Regression tests for PBxplore.

This script run the various PBxplore programs with various argument, and makes
sure the output is the expected one. The aim is to check that the programs are
not broken during development.

Be careful this script does not test that the output is right. It just test
that the output is the expected one based on a previous version.

To run this test suite, you can either run this script without arguments or use
the nose [1]_ package. The latter option gives a more readable output, as
stdout is captured, and is displayed only if a test fails.

.. [1] https://nose.readthedocs.org
"""

# Use print as a function like in python 3
from __future__ import print_function

import unittest
from os import path
from uuid import uuid1
from functools import wraps
import os
import subprocess
import shutil

# Resources for the tests are stored in the following directory
REFDIR = "test_data/"

# Temporary outputs will be written in this directory
OUTDIR = "test_tmp/"


def _failure_test(method):
    """
    Decorate tests that are supposed to fail

    Some tests asses that things that should go wrong actually go wrong. These
    tests should fail. The _failure_test decorator reverse the assesment test
    so if the decorated test sucess it is reported as failed.

    This decorator differs to unittest.expectedFailure. The
    unittest.expectedFailure aims at decorate tests that are supposed to sucess
    but are known to failed.
    """
    @wraps(method)
    def wrapped(*args, **kwargs):
        try:
            method(*args, **kwargs)
        except AssertionError:
            pass
        else:
            raise AssertionError('Test should have failed.')
    return wrapped


class TestPBAssign(unittest.TestCase):
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

    def test_flat(self):
        """
        Run PBassign with the --flat option.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        _test_PBassign_options(references, ['{0}.PB.fasta', '{0}.PB.flat'],
                               ['--flat'])

    def test_phipsi(self):
        """
        Run PBassign with the --phipsi option.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        _test_PBassign_options(references, ['{0}.PB.fasta', '{0}.PB.phipsi'],
                               ['--phipsi'])

    def test_flat_phipsi(self):
        """
        Run PBassign with the both the --flat and --phipsi options.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        _test_PBassign_options(references,
                               ['{0}.PB.fasta', '{0}.PB.flat', '{0}.PB.phipsi'],
                               ['--flat', '--phipsi'])

    def test_multiple_inputs(self):
        """
        Run PBassign with multiple inputs.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        _test_PBassign_options(references,
                               ['{0}.PB.fasta', '{0}.PB.flat', '{0}.PB.phipsi'],
                               ['--flat', '--phipsi'], multiple='all')

    @_failure_test
    def test_missing_output(self):
        """
        Test if the tests properly fail if an output is missing.
        """
        references = ["1BTA"]
        _test_PBassign_options(references,
                               ['{0}.PB.fasta', '{0}.PB.flat', '{0}.PB.phipsi',
                                '{0}.missing'],
                               ['--flat', '--phipsi'])

    @_failure_test
    def test_too_many_outputs(self):
        """
        Test if the tests properly fail if an output is not expected.
        """
        references = ["1BTA"]
        _test_PBassign_options(references,
                               ['{0}.PB.fasta', '{0}.PB.flat'],
                               ['--flat', '--phipsi'])

    @_failure_test
    def test_different_outputs(self):
        """
        Test if the tests properly fail if an output content is different from
        expected.
        """
        references = ["test_fail"]
        _test_PBassign_options(references, ['{0}.PB.fasta'], [])


def _same_file_content(file_a, file_b):
    """
    Return True if two files are identical. Take file path as arguments.
    """
    with open(file_a) as f1, open(file_b) as f2:
        # Compare content line by line
        for f1_line, f2_line in zip(f1, f2):
            if f1_line != f2_line:
                print(file_a, file_b)
                print(f1_line, f2_line, sep='//')
                return False
        # Check if one file is longer than the other; it would result as one
        # file iterator not completely consumed
        for infile in (f1, f2):
            if infile.readline() != '':
                # The iterator is not consumed, it means that this file is
                # longer than the other
                print('File too long')
                return False
    # If we reach this line, it means that we did not find any difference
    return True


def _assert_identical_files(file_a, file_b):
    """
    Raise an Assert exception if the two files are not identical.

    Take file path as arguments.
    """
    assert _same_file_content(file_a, file_b), '{0} and {1} are not identical'\
                                               .format(file_a, file_b)


def _run_prog(program, pdbid, options,
              multiple=None, indir=REFDIR, outdir=OUTDIR):
    """
    Run a PBxplore program on a PDBID with the given options.

    options is expected to be a list that will be directly pass to subprocess,
    it must not contain the input or output options.
    """
    if program not in ('./PBassign.py',):
        raise NotImplementedError('_run_prog does not know how to run {0}'
                                  .format(program))
    out_run_dir = path.join(OUTDIR, str(uuid1()))
    if multiple is None:
        test_input = path.join(REFDIR, pdbid + '.pdb')
        out_basename = path.join(out_run_dir, pdbid)
        input_args = ['-p', test_input]
    else:
        input_args = []
        for basename in pdbid:
            input_args += ['-p', path.join(REFDIR, basename + '.pdb')]
        out_basename = path.join(out_run_dir, multiple)

    os.mkdir(out_run_dir)

    run_list = [program] + input_args + ['-o', out_basename] + options
    print(' '.join(run_list))
    status = subprocess.call(run_list, stdout=subprocess.PIPE)

    return status, out_run_dir


def _test_PBassign_options(basenames, outfiles, options,
                           multiple=None, expected_exit=0):
    if not multiple is None:
        basenames = [basenames]
        out_name = multiple
    for basename in basenames:
        status, out_run_dir = _run_prog('./PBassign.py', basename, options,
                                        multiple)
        assert status == expected_exit, \
               'PBassign stoped with a {0} exit code'.format(status)
        assert len(os.listdir(out_run_dir)) == len(outfiles),\
                ('PBassign did not produced the right number of files: '
                 '{0} files produced instead of {1}').format(
                    len(os.listdir(out_run_dir)), len(outfiles))
        out_name = basename if multiple is None else multiple
        for outfile in (template.format(out_name) for template in outfiles):
            test_file = path.join(out_run_dir, outfile)
            ref_file = path.join(REFDIR, outfile)
            _assert_identical_files(test_file, ref_file)
        shutil.rmtree(out_run_dir)


if __name__ == '__main__':
    unittest.main()
