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
import sys
import warnings
import contextlib

try:
    import MDAnalysis
except ImportError:
    IS_MDANALYSIS = False
    warnings.warn('MDanalysis is not available, tests will be run accordingly')
else:
    IS_MDANALYSIS = True

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


class TemplateTestCase(unittest.TestCase):
    """
    Template TestCase class for the other TestCase class to inherit from.

    Children class must overload the `_build_command_line` and the
    `_validate_output` methods.
    """
    def setUp(self):
        """
        Run before each test.

        Make sure that the output directory exists. And create a temporary
        directory to work in.
        """
        if not path.isdir(OUTDIR):
            os.mkdir(OUTDIR)
        self._temp_directory = os.path.join(OUTDIR, str(uuid1()))
        os.mkdir(self._temp_directory)

    def tearDown(self):
        """
        Run after each test.

        Delete the temporary directory except if the test failed. Note that, on
        python 3, the temporary directory is always deleted.
        """
        if ((sys.version_info[0] == 2
                 and sys.exc_info() == (None, None, None))
                or sys.version_info[0] == 3):
            # On python 2, sys.exc_info() is (None, None, None) only when a
            # test pass. On python 3, however, there is no difference in 
            # sys.exc_info() between a passing and a failing test. Here, on
            # python 2, we delete the temporary directory only is the test
            # passes; on python 3 we always delete the temporary directory.
            shutil.rmtree(self._temp_directory)

    def _run_program_and_validate(self, reference, **kwargs):
        """
        Run the program to test and validate its outputs.
        """
        # Build the command line to run. This relies on the _build_command_line
        # method that is a virtual method, which must be overloaded by the
        # child class.
        command = self._build_command_line(**kwargs)
        print(command)

        # Run the command.
        exe = subprocess.Popen(command,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = exe.communicate()
        print(out.decode('utf-8'))
        print(err.decode('utf-8'))

        # The return code should be 0.
        assert exe.returncode == 0, 'Program exited with a {} code.'.format(
            exe.returncode)

        # Validate the output files. This relies on the _validate_output
        # virtual method.
        self._validate_output(reference, **kwargs)

    def _build_command_line(self, **kwargs):
        """
        Build the command line to run.

        This is a virtual method. It must be overloaded by the child class.
        """
        raise NotImplementedError

    def _validate_output(self, reference, **kwargs):
        """
        Validate the output files.

        This is a virtual method. It must be overloaded by the child class.
        """
        raise NotImplementedError


class TestPBAssign(TemplateTestCase):
    """
    Regression tests for PBAssign.py
    """
    def _run_PBassign(self, pdbid, extension, options,
                      multiple=None, indir=REFDIR, outdir=OUTDIR):
        """
        Run a PBxplore program on a PDBID with the given options.

        `options` is expected to be a list that will be directly passed to
        subprocess, it must not contain the input or output options.
        """
        out_run_dir = os.path.join(self._temp_directory, str(uuid1()))
        os.mkdir(out_run_dir)
        if multiple is None:
            test_input = path.join(REFDIR, pdbid + extension)
            out_basename = path.join(out_run_dir, pdbid)
            input_args = ['-p', test_input]
        else:
            input_args = []
            for basename in pdbid:
                input_args += ['-p', path.join(REFDIR, basename + extension)]
            out_basename = path.join(out_run_dir, multiple)

        run_list = (['./PBassign.py'] + input_args + 
                    ['-o', out_basename + extension] + options)
        print(' '.join(run_list))
        exe = subprocess.Popen(run_list,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = exe.communicate()
        print(out.decode('utf-8'))
        print(err.decode('utf-8'))

        return exe.returncode, out_run_dir

    def _test_PBassign_options(self, basenames, extensions, outfiles, options,
                               multiple=None, expected_exit=0):
        if not multiple is None:
            basenames = [basenames]
            out_name = multiple
        for basename in basenames:
            for extension in extensions:
                status, out_run_dir = self._run_PBassign(basename, extension,
                                                         options, multiple)
                assert status == expected_exit, \
                       'PBassign stoped with a {0} exit code'.format(status)
                assert len(os.listdir(out_run_dir)) == len(outfiles),\
                        ('PBassign did not produced the right number of files: '
                         '{0} files produced instead of {1}').format(
                            len(os.listdir(out_run_dir)), len(outfiles))
                out_name = basename if multiple is None else multiple
                for outfile in (template.format(out_name + extension)
                                for template in outfiles):
                    test_file = path.join(out_run_dir, outfile)
                    ref_file = path.join(REFDIR, outfile)
                    _assert_identical_files(test_file, ref_file)

    def test_fasta(self):
        """
        Run PBAssign on PDB files, and check the fasta output.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        extensions = [".pdb", ".cif.gz"]
        self._test_PBassign_options(references, extensions,
                                    ['{0}.PB.fasta'], [])

    def test_flat(self):
        """
        Run PBassign with the --flat option.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        extensions = [".pdb", ".cif.gz"]
        self._test_PBassign_options(references, extensions,
                                    ['{0}.PB.fasta', '{0}.PB.flat'],
                                    ['--flat'])

    def test_phipsi(self):
        """
        Run PBassign with the --phipsi option.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        extensions = [".pdb", ".cif.gz"]
        self._test_PBassign_options(references, extensions,
                                    ['{0}.PB.fasta', '{0}.PB.phipsi'],
                                    ['--phipsi'])

    def test_flat_phipsi(self):
        """
        Run PBassign with the both the --flat and --phipsi options.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        extensions = [".pdb", ".cif.gz"]
        self._test_PBassign_options(references, extensions,
                                    ['{0}.PB.fasta', '{0}.PB.flat', '{0}.PB.phipsi'],
                                    ['--flat', '--phipsi'])

    def test_multiple_inputs(self):
        """
        Run PBassign with multiple inputs.
        """
        references = ["1BTA", "1AY7", "2LFU", "3ICH"]
        extensions = [".pdb", ".cif.gz"]
        self._test_PBassign_options(references, extensions,
                                    ['{0}.PB.fasta', '{0}.PB.flat', '{0}.PB.phipsi'],
                                    ['--flat', '--phipsi'], multiple='all')

    def test_xtc_input(self):
        """
        Run PBassign on a trajectory in the XTC format.

        This test should produce the righ output with python 2. With python 3,
        PBassign should fail as MDanalysis is not available.
        """
        name = 'barstar_md_traj'
        out_run_dir = self._temp_directory
        output_fname = name + '.PB.fasta'
        call_list = ['./PBassign.py',
                     '-x', os.path.join(REFDIR, name + '.xtc'),
                     '-g', os.path.join(REFDIR, name + '.gro'),
                     '-o', os.path.join(out_run_dir, name)]
        exe = subprocess.Popen(call_list,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        out, err = exe.communicate()
        status = exe.wait()
        print(out.decode('utf-8'))
        print(err.decode('utf-8'))
        if IS_MDANALYSIS:
            # MDanalysis is available, PBassign should run and produce the
            # correct output
            assert status == 0, 'PBassign exited with an error'
            _assert_identical_files(os.path.join(REFDIR, output_fname),
                                    os.path.join(out_run_dir, output_fname))
        else:
            # MDanalysis is not available, PBassign should fail
            assert status != 0, 'PBassign shoud not have exited with a 0 code'


    @_failure_test
    def test_missing_output(self):
        """
        Test if the tests properly fail if an output is missing.
        """
        references = ["1BTA"]
        extensions = [".pdb"]
        self._test_PBassign_options(references, extensions, 
                                    ['{0}.PB.fasta', '{0}.PB.flat',
                                     '{0}.PB.phipsi', '{0}.missing'],
                                    ['--flat', '--phipsi'])

    @_failure_test
    def test_too_many_outputs(self):
        """
        Test if the tests properly fail if an output is not expected.
        """
        references = ["1BTA"]
        extensions = [".pdb"]
        self._test_PBassign_options(references, extensions, 
                                    ['{0}.PB.fasta', '{0}.PB.flat'],
                                    ['--flat', '--phipsi'])

    @_failure_test
    def test_different_outputs(self):
        """
        Test if the tests properly fail if an output content is different from
        expected.
        """
        references = ["test_fail"]
        extensions = [".pdb"]
        self._test_PBassign_options(references, extensions, ['{0}.PB.fasta'], [])


class TestPBcount(TemplateTestCase):
    """
    Test running PBcount.
    """
    def _build_command_line(self, input_files, output, first_residue=None):
        output_full_path = os.path.join(self._temp_directory, output)
        command = ['./PBcount.py', '-o', output_full_path]
        for input_file in input_files:
            command += ['-f', os.path.join(REFDIR, input_file)]
        if not first_residue is None:
            command += ['--first-residue', str(first_residue)]
        return command

    def _validate_output(self, reference, output, **kwargs):
        reference_full_path = os.path.join(REFDIR, reference)
        output_full_path = os.path.join(self._temp_directory,
                                        output + '.PB.count')
        _assert_identical_files(output_full_path, reference_full_path)

    def test_single_file_single_model(self):
        """
        Run PBcount with a single input file that contains a single model.
        """
        input_files = ['count_single1.PB.fasta',]
        output = 'output'
        reference = 'count_single1.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output)

    def test_single_file_multiple_models(self):
        """
        Run PBcount with a single input file that contains multiple models.
        """
        input_files = ['count_multi1.PB.fasta',]
        output = 'output'
        reference = 'count_multi1.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output)

    def test_multiple_files_single_model(self):
        """
        Run PBcount with multiple input files that contain a single model.
        """
        input_files = ['count_single1.PB.fasta',
                       'count_single2.PB.fasta',
                       'count_single3.PB.fasta']
        output = 'output'
        reference = 'count_single123.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output)

    def test_multiple_files_multiple_models(self):
        """
        Run PBcount with multiple input files that contain multiple models each.
        """
        input_files = ['count_multi1.PB.fasta',
                       'count_multi2.PB.fasta',
                       'count_multi3.PB.fasta']
        output = 'output'
        reference = 'count_multi123.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output)

    def test_first_residue_positive(self):
        """
        Test PBcount on with the --first-residue option and a positive value.
        """
        input_files = ['count_multi1.PB.fasta',]
        output = 'output'
        reference = 'count_multi1_first20.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output,
                                       first_residue=20)

    def test_first_residue_negative(self):
        """
        Test PBcount on with the --first-residue option and a negative value.
        """
        input_files = ['count_multi1.PB.fasta',]
        output = 'output'
        reference = 'count_multi1_first-20.PB.count'
        self._run_program_and_validate(reference,
                                       input_files=input_files, output=output,
                                       first_residue=-20)


class TestPBclust(TemplateTestCase):
    def _build_command_line(self, input_files, output,
                            clusters=None, compare=False): 
        output_full_path = os.path.join(self._temp_directory, output)
        command = ['./PBclust.py', '-o', output_full_path]
        for input_file in input_files:
            command += ['-f', os.path.join(REFDIR, input_file)]
        if not clusters is None:
            command += ['-c', str(clusters)]
        if compare:
            command += ['--compare']
        return command

    def _validate_output(self, reference, input_files, output,
                         clusters=None, compare=False, **kwargs):
        output = os.path.join(self._temp_directory, output)
        if compare:
            # Asses the validity of the distance file
            reference_full_path = os.path.join(REFDIR,
                                               reference + '.PB.compare.fasta')
            output_full_path = output + '.PB.compare.fasta'
            _assert_identical_files(output_full_path, reference_full_path)
        else:
            # Asses the validity of the main output od PBclust (the clust file)
            reference_full_path = os.path.join(REFDIR, reference + '.PB.clust')
            output_full_path = output + '.PB.clust'
            _assert_identical_files(output_full_path, reference_full_path)

    def test_default_single_input(self):
        self._run_program_and_validate(reference='psi_md_traj_all',
                                       input_files=['psi_md_traj_all.PB.fasta',],
                                       output='output', clusters=3)

    def test_default_multi_input(self):
        self._run_program_and_validate(reference='psi_md_traj_all',
                                       input_files=['psi_md_traj_1.PB.fasta',
                                                    'psi_md_traj_2.PB.fasta',
                                                    'psi_md_traj_3.PB.fasta'],
                                       output='output', clusters=3)

    def test_nclusters(self):
        self._run_program_and_validate(reference='psi_md_traj_all_c5',
                                       input_files=['psi_md_traj_all.PB.fasta',],
                                       output='output',
                                       clusters=5)

    def test_compare(self):
        self._run_program_and_validate(reference='psi_md_traj_1',
                                       input_files=['psi_md_traj_1.PB.fasta',],
                                       output='output',
                                       compare=True, clusters=3)


def _same_file_content(file_a, file_b):
    """
    Return True if two files are identical. Take file path as arguments.
    """
    with open(file_a) as f1, open(file_b) as f2:
        # Compare content line by line
        for f1_line, f2_line in zip(f1, f2):
            if (f1_line != f2_line):
                print(file_a, file_b)
                print(f1_line, f2_line, sep='//')
                return False
        # Check if one file is longer than the other; it would result as one
        # file iterator not completely consumed
        for infile in (f1, f2):
            try:
                next(infile) 
            except StopIteration:
                pass
            else:
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


if __name__ == '__main__':
    unittest.main()
