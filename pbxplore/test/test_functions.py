#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for PBxplore.

Tests functions from different programs.

2014 - P. Poulain
"""

# =============================================================================
# Modules
# =============================================================================
import unittest
import collections
import os

import numpy

import pbxplore as pbx
from pbxplore.structure import structure
from pbxplore.analysis import kmeans

here = os.path.abspath(os.path.dirname(__file__))


# =============================================================================
# Classes for tests
# =============================================================================
class TestPDBlib(unittest.TestCase):
    """
    Tests for PDBlib
    """

    def test_get_dihedral(self):
        """
        Test for get_dihedral()
        """
        Result = collections.namedtuple('Result', ['A', 'B', 'C', 'D', 'torsion'])
        results = (Result((-7.28, -9.262, 5.077),
                          (-7.526, -10.643, 5.529),
                          (-6.221, -11.438, 5.555),
                          (-6.289, -12.685, 5.931),
                          -179.663656153),
                   Result((-1.373, -8.817, -4.389),
                          (-1.203, -8.335, -5.792),
                          (-1.891, -6.977, -5.927),
                          (-1.918, -6.429, -7.107),
                          -176.048770127),
                   Result((-0.533, -8.42, -3.47  ),
                          (-1.373, -8.817, -4.389),
                          (-1.203, -8.335, -5.792),
                          (-1.891, -6.977, -5.927),
                          -84.8356057692),
                   Result((-1.918, -6.429, -7.107),
                          (-2.609, -5.125, -7.305),
                          (-4.108, -5.392, -7.331),
                          (-4.469, -6.494, -7.911),
                          -36.8942888266),
                   Result((-11.285, 6.472, -7.44 ),
                          (-12.62, 5.829, -7.425 ),
                          (-13.585, 6.626, -6.544),
                          (-13.098, 7.621, -5.858),
                          -6.58786169376),
                   Result((-11.284, -0.971, -2.679),
                          (-12.65, -0.794, -3.226),
                          (-13.665, -1.664, -2.479),
                          (-13.262, -2.363, -1.452),
                          3.91626706556),
                   Result((-2.004, -10.892, -2.611),
                          (-1.87, -9.835, -1.853),
                          (-0.726, -8.877, -2.011),
                          (-0.533, -8.42, -3.47),
                          50.065196067),
                   Result((11.174, -6.725, 0.458),
                          (10.732, -7.258, -0.86),
                          (9.27, -6.869, -1.096),
                          (8.741, -7.185, -2.245),
                          175.872397707))

        for res in results:
            torsion = structure.get_dihedral(res.A, res.B, res.C, res.D)
            self.assertAlmostEqual(torsion, res.torsion)


class TestAtomClass(unittest.TestCase):
    """
    Tests for the Atom class in PDBlib
    """

    def test_read_from_PDB(self):
        """
        Tests for read_from_PDB()
        """
        a = structure.Atom.read_from_PDB("ATOM    512  N   GLU A  32      -1.870  -9.835  -1.853  1.00  0.56           N  ")
        self.assertAlmostEqual(a.coords, [-1.87, -9.835, -1.853])
        a = structure.Atom.read_from_PDB("ATOM   1424  CA  SER A  89       7.604  11.308   1.435  1.00  0.62           C  ")
        self.assertAlmostEqual(a.coords, [7.604, 11.308, 1.435])
        a = structure.Atom.read_from_PDB("ATOM   1167  CG2 VAL B  50       9.294  44.541  -4.830  1.00 27.62           C  ")
        self.assertAlmostEqual(a.coords, [9.294, 44.541, -4.83])

    def test_read_from_PDBx(self):
        """
        Tests for read_from_PDBx()
        """
        fields = ['group_PDB', 'id', 'type_symbol', 'label_atom_id',
                  'label_alt_id', 'label_comp_id', 'label_asym_id', 'label_entity_id',
                  'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
                  'occupancy', 'B_iso_or_equiv', 'Cartn_x_esd', 'Cartn_y_esd',
                  'Cartn_z_esd', 'occupancy_esd', 'B_iso_or_equiv_esd',
                  'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id',
                  'auth_atom_id', 'pdbx_PDB_model_num']
        line = "ATOM 4769  H HB   . ILE A 1 35  ? -20.422 5.104   -0.207  1.00 0.00 ? ? ? ? ? ? 277 ILE A HB   3"
        a = structure.Atom.read_from_PDBx(line, fields)
        self.assertAlmostEqual(a.coords, [-20.422, 5.104, -0.207])
        line = "ATOM 18201 H HG21 . THR A 1 140 ? 11.080  -12.466 -8.977  1.00 0.00 ? ? ? ? ? ? 382 THR A HG21 8"
        a = structure.Atom.read_from_PDBx(line, fields)
        self.assertAlmostEqual(a.coords, [11.08, -12.466, -8.977])
        line = "ATOM 23720 H HE2  . HIS A 1 193 ? 13.974  24.297  0.352   1.00 0.00 ? ? ? ? ? ? 435 HIS A HE2  10"
        a = structure.Atom.read_from_PDBx(line, fields)
        self.assertAlmostEqual(a.coords, [13.974, 24.297, 0.352])


class TestChainClass(unittest.TestCase):
    """
    Tests for Chain class in PDBlib
    """

    def test_size(self):
        """
        Tests for size()
        """
        lines = ("ATOM    840  C   ARG B  11      22.955  23.561  -4.012  1.00 28.07           C  ",
                 "ATOM    849  N   SER B  12      22.623  24.218  -2.883  1.00 24.77           N  ",
                 "ATOM    850  CA  SER B  12      22.385  23.396  -1.637  1.00 21.99           C  ",
                 "ATOM    851  C   SER B  12      21.150  24.066  -0.947  1.00 32.67           C  ",
                 "ATOM    855  N   ILE B  13      20.421  23.341  -0.088  1.00 30.25           N  ")
        ch = structure.Chain()
        for line in lines:
            at = structure.Atom.read_from_PDB(line)
            ch.add_atom(at)
        self.assertEqual(ch.size(), 5)

    def test_get_phi_psi_angles(self):
        """
        Tests for get_phi_psi_angles()
        """
        results = {11: {'phi': None, 'psi': None},
                   12: {'phi': -139.77684605036447, 'psi': 157.94348570201197},
                   13: {'phi': None, 'psi': None}}
        lines = ("ATOM    840  C   ARG B  11      22.955  23.561  -4.012  1.00 28.07           C  ",
                 "ATOM    849  N   SER B  12      22.623  24.218  -2.883  1.00 24.77           N  ",
                 "ATOM    850  CA  SER B  12      22.385  23.396  -1.637  1.00 21.99           C  ",
                 "ATOM    851  C   SER B  12      21.150  24.066  -0.947  1.00 32.67           C  ",
                 "ATOM    855  N   ILE B  13      20.421  23.341  -0.088  1.00 30.25           N  ")
        ch = structure.Chain()
        for line in lines:
            at = structure.Atom.read_from_PDB(line)
            ch.add_atom(at)

        phi_psi = ch.get_phi_psi_angles()
        for resid, angles in results.items():
            self.assertAlmostEqual(angles["phi"], phi_psi[resid]["phi"])
            self.assertAlmostEqual(angles["psi"], phi_psi[resid]["psi"])


class TestPBlib(unittest.TestCase):
    """
    Tests for PBlib
    """

    def test_read_fasta(self):
        headers, sequences = pbx.io.read_fasta(os.path.join(here, "test_data/1BTA.pdb.PB.fasta"))
        self.assertEqual(headers, ['test_data/1BTA.pdb | chain A'])
        self.assertEqual(sequences, ['ZZdddfklonbfklmmmmmmmmnopafklnoiakl'
                                     'mmmmmnoopacddddddehkllmmmmngoilmmmm'
                                     'mmmmmmmmnopacdcddZZ'])


class TestKMeans(unittest.TestCase):
    """
    Test individual functions for K-Means clustering
    """

    def test_count_per_position_partial(self):
        sequences = ['abcdef', 'bcdefg', 'cdefgh',
                     'defghi', 'efghij', 'fghijk',  # ignore in the test
                     'ghijkl', 'hijklm', 'ijklmn',
                     'ghijkl', 'hijklm', 'ijklmn',
                     'jklmno', 'klmnop', ]          # ignore in the test
        indices = [0, 1, 2, 6, 7, 8, 9, 10, 11]

        # Shorten defaultdict call to have a nicer table
        dd = collections.defaultdict
        # Create a list with the right number of elements. The type of each
        # element does not matter as the elements will be overwritten
        ref_count = ['' for _ in range(6)]
        # Fill the reference. Each element of the list corresponds to a
        # position in the sequences. Only the sequences that have their index
        # in the ``indices`` list will be used in the count.
        ref_count[0] = dd(int, {'a': 1, 'b': 1, 'c': 1, 'g': 2, 'h': 2, 'i': 2})
        ref_count[1] = dd(int, {'b': 1, 'c': 1, 'd': 1, 'h': 2, 'i': 2, 'j': 2})
        ref_count[2] = dd(int, {'c': 1, 'd': 1, 'e': 1, 'i': 2, 'j': 2, 'k': 2})
        ref_count[3] = dd(int, {'d': 1, 'e': 1, 'f': 1, 'j': 2, 'k': 2, 'l': 2})
        ref_count[4] = dd(int, {'e': 1, 'f': 1, 'g': 1, 'k': 2, 'l': 2, 'm': 2})
        ref_count[5] = dd(int, {'f': 1, 'g': 1, 'h': 1, 'l': 2, 'm': 2, 'n': 2})

        count = kmeans.count_per_position_partial(sequences, indices)
        self.assertEqual(count, ref_count)

    def test_make_profile_partial(self):
        sequences = ['abcdef', 'bcdefg', 'cdefgh',
                     'defghi', 'efghij', 'fghijk',  # ignore in the test
                     'ghijkl', 'hijklm', 'ijklmn',
                     'ghijkl', 'hijklm', 'ijklmn',
                     'jklmno', 'klmnop',            # ignore in the test
                     'ijklmn']
        # Using 10 sequences makes things easier
        indices = [0, 1, 2, 6, 7, 8, 9, 10, 11, 14]
        ref_profile = numpy.array([[0.1, 0.0, 0.0, 0.0, 0.0, 0.0],   # a
                                   [0.1, 0.1, 0.0, 0.0, 0.0, 0.0],   # b
                                   [0.1, 0.1, 0.1, 0.0, 0.0, 0.0],   # c
                                   [0.0, 0.1, 0.1, 0.1, 0.0, 0.0],   # d
                                   [0.0, 0.0, 0.1, 0.1, 0.1, 0.0],   # e
                                   [0.0, 0.0, 0.0, 0.1, 0.1, 0.1],   # f
                                   [0.2, 0.0, 0.0, 0.0, 0.1, 0.1],   # g
                                   [0.2, 0.2, 0.0, 0.0, 0.0, 0.1],   # h
                                   [0.3, 0.2, 0.2, 0.0, 0.0, 0.0],   # i
                                   [0.0, 0.3, 0.2, 0.2, 0.0, 0.0],   # j
                                   [0.0, 0.0, 0.3, 0.2, 0.2, 0.0],   # k
                                   [0.0, 0.0, 0.0, 0.3, 0.2, 0.2],   # l
                                   [0.0, 0.0, 0.0, 0.0, 0.3, 0.2],   # m
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.3],   # n
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # o
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])  # p
        profile = kmeans.make_profile_partial(sequences, indices)
        assert(numpy.allclose(ref_profile, profile))

    def test_compatibility_identity(self):
        profile = numpy.array([[1.0, 0.0, 0.0, 0.0],   # a
                               [0.0, 1.0, 0.0, 0.0],   # b
                               [0.0, 0.0, 1.0, 0.0],   # c
                               [0.0, 0.0, 0.0, 1.0],   # d
                               [0.0, 0.0, 0.0, 0.0],   # e
                               [0.0, 0.0, 0.0, 0.0],   # f
                               [0.0, 0.0, 0.0, 0.0],   # g
                               [0.0, 0.0, 0.0, 0.0],   # h
                               [0.0, 0.0, 0.0, 0.0],   # i
                               [0.0, 0.0, 0.0, 0.0],   # j
                               [0.0, 0.0, 0.0, 0.0],   # k
                               [0.0, 0.0, 0.0, 0.0],   # l
                               [0.0, 0.0, 0.0, 0.0],   # m
                               [0.0, 0.0, 0.0, 0.0],   # n
                               [0.0, 0.0, 0.0, 0.0],   # o
                               [0.0, 0.0, 0.0, 0.0]])  # p
        sequence = 'abcd'
        reference_compatibility = len(sequence)
        compatibility = kmeans.compatibility(profile, sequence)
        self.assertAlmostEqual(compatibility, reference_compatibility)

    def test_compatibility(self):
        # The profile correspond to the following sequences:
        # - abcd      - ioph    - alep
        # - jiop      - fgad    - hobe
        # - abjp      - piep
        # - hbno      - jojo
        profile = numpy.array([[0.3, 0.0, 0.1, 0.0],   # a
                               [0.0, 0.3, 0.1, 0.0],   # b
                               [0.0, 0.0, 0.1, 0.0],   # c
                               [0.0, 0.0, 0.0, 0.2],   # d
                               [0.0, 0.1, 0.2, 0.1],   # e
                               [0.1, 0.0, 0.0, 0.0],   # f
                               [0.0, 0.1, 0.0, 0.0],   # g
                               [0.2, 0.0, 0.0, 0.1],   # h
                               [0.1, 0.2, 0.0, 0.0],   # i
                               [0.2, 0.0, 0.2, 0.0],   # j
                               [0.0, 0.0, 0.0, 0.0],   # k
                               [0.0, 0.1, 0.0, 0.0],   # l
                               [0.0, 0.0, 0.0, 0.0],   # m
                               [0.0, 0.0, 0.1, 0.0],   # n
                               [0.0, 0.2, 0.1, 0.2],   # o
                               [0.1, 0.0, 0.1, 0.4]])  # p
        print(profile.sum(axis=0))
        assert(numpy.allclose(profile.sum(axis=0),
                              numpy.ones((profile.shape[1],),
                                         dtype=profile.dtype)))
        reference_sequences = (('abcd', 0.9), ('hipo', 0.7), ('bcda', 0.0),
                               ('aaaa', 0.4), ('kkkk', 0.0), ('oooo', 0.5),)
        for sequence, reference_compatibility in reference_sequences:
            compatibility = kmeans.compatibility(profile, sequence)
            self.assertAlmostEqual(compatibility, reference_compatibility)

    def test_argmax(self):
        reference = (([0, 1, 2, 3, 4], 4),  # Ordered
                     ([4, 3, 2, 1, 0], 0),  # Reverse ordered
                     ([1, 3, 2, 4, 0], 3),  # Random order
                     ([0, 0, 3, 4, 4], 3),  # Duplicates
                     ([4, 4, 4, 4, 4], 0),  # All the same
                     )
        for test_case, expectation in reference:
            self.assertEqual(kmeans._argmax(test_case), expectation)


if __name__ == '__main__':
    unittest.main()
