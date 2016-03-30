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

import pbxplore as pbx
from pbxplore.structure import structure

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

if __name__ == '__main__':
    unittest.main()
