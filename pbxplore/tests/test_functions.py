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
import collections
import os
import numpy

import pytest

import pbxplore as pbx
from pbxplore.structure import structure

import MDAnalysis

here = os.path.abspath(os.path.dirname(__file__))


# =============================================================================
# Classes for tests
# =============================================================================


Result = collections.namedtuple('Result', ['A', 'B', 'C', 'D', 'torsion'])


class TestStructurelib(object):
    """
    Tests for Structurelib
    """

    @pytest.mark.parametrize('result',
        (Result((-7.28, -9.262, 5.077),
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
     )
    def test_get_dihedral(self, result):
        """
        Test for get_dihedral()
        """
        torsion = structure.get_dihedral(result.A, result.B,
                                         result.C, result.D)
        assert torsion == pytest.approx(result.torsion)

    def test_loader_PDB(self):
        """
        Test for API loader function on PDBs
        """
        filename = os.path.join(here, "test_data/2LFU.pdb")
        comment, chain = list(pbx.chains_from_files([filename]))[0]

        ref_comment = "{0} | model 1 | chain A".format(filename)
        ref_chain = "Chain A / model 1: 2372 atoms"
        assert ref_comment == comment
        assert ref_chain == format(chain)

    def test_loader_xtc(self):
        """
        Test for API load function on xtc files
        """
        topol = os.path.join(here, "test_data/barstar_md_traj.gro")
        traj = os.path.join(here, "test_data/barstar_md_traj.xtc")
        chains = list(pbx.chains_from_trajectory(traj, topol))

        comment, chain = chains[0]
        ref_comment = "{0} | frame 0".format(traj)
        ref_chain = "Chain  / model : 355 atoms"
        assert ref_comment == comment
        assert ref_chain == format(chain)

        comment, chain = chains[-1]
        ref_comment = "{0} | frame 9".format(traj)
        assert ref_comment == comment
        assert ref_chain == format(chain)


class TestAtomClass(object):
    """
    Tests for the Atom class in PDBlib
    """

    PDBx_fields = ['group_PDB', 'id', 'type_symbol', 'label_atom_id',
                   'label_alt_id', 'label_comp_id', 'label_asym_id',
                   'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code',
                   'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
                   'B_iso_or_equiv', 'Cartn_x_esd', 'Cartn_y_esd',
                   'Cartn_z_esd', 'occupancy_esd', 'B_iso_or_equiv_esd',
                   'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id',
                   'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num']

    @pytest.mark.parametrize(
        'line, expected',
        (("ATOM    512  N   GLU A  32      -1.870  -9.835  -1.853  1.00  0.56           N  ", [-1.87, -9.835, -1.853]),
         ("ATOM   1424  CA  SER A  89       7.604  11.308   1.435  1.00  0.62           C  ", [7.604, 11.308, 1.435]),
         ("ATOM   1167  CG2 VAL B  50       9.294  44.541  -4.830  1.00 27.62           C  ", [9.294, 44.541, -4.83]))
    )
    def test_read_from_PDB(self, line, expected):
        """
        Tests for read_from_PDB()
        """
        atom = structure.Atom.read_from_PDB(line)
        assert atom.coords == pytest.approx(expected)

    def test_read_PDB_line_short(self):
        """
        Test when PDB line is too short
        """
        with pytest.raises(structure.AtomError):
            structure.Atom.read_from_PDB("ATOM    512  N   GLU A  32      -1.870  -9.835")

    @pytest.mark.parametrize(
        'line,coordinates',
        (("ATOM 4769  H HB   . ILE A 1 35  ? -20.422 5.104   -0.207  1.00 0.00 ? ? ? ? ? ? 277 ILE A HB   3", [-20.422, 5.104, -0.207]),
         ("ATOM 18201 H HG21 . THR A 1 140 ? 11.080  -12.466 -8.977  1.00 0.00 ? ? ? ? ? ? 382 THR A HG21 8", [11.08, -12.466, -8.977]),
         ("ATOM 23720 H HE2  . HIS A 1 193 ? 13.974  24.297  0.352   1.00 0.00 ? ? ? ? ? ? 435 HIS A HE2  10", [13.974, 24.297, 0.352]))
    )
    def test_read_from_PDBx(self, line, coordinates):
        """
        Tests for read_from_PDBx()
        """
        atom = structure.Atom.read_from_PDBx(line, self.PDBx_fields)
        assert atom.coords == pytest.approx(coordinates)

    @pytest.mark.parametrize(
        'line',
        ("ATOM 4769  H HB   . ILE A 1 35  ? -20.422 5.104",
         "ATOM 4769  H HB   . ILE A 1 XXX  ? -20.422 5.104   -0.207  1.00 0.00 ? ? ? ? ? ? 277 ILE A HB   3",)
    )
    def test_read_PDBx_failed_line(self, line):
        """
        Test when PDBx line is not correctly formated
        """
        with pytest.raises(structure.AtomError):
            structure.Atom.read_from_PDBx(line, self.PDBx_fields)

    def test_read_from_xtc(self):
        """
        Tests for read_from_xtc()
        """
        topology = os.path.join(here, "test_data/barstar_md_traj.gro")
        traj = os.path.join(here, "test_data/barstar_md_traj.xtc")
        universe = MDAnalysis.Universe(topology, traj)
        selection = universe.select_atoms("backbone")

        # First timeframe
        atom = structure.Atom.read_from_xtc(selection[0])
        assert atom.resid == 1
        assert atom.name == "N"
        for a, b in zip(atom.coords, [21.68, 33.87, 36.18]):
            assert a == pytest.approx(b, abs=1e-3)
        atom = structure.Atom.read_from_xtc(selection[-1])
        assert atom.resid == 89
        assert atom.name == "C"
        for a, b in zip(atom.coords, [40.14, 38.75, 28.42]):
            assert a == pytest.approx(b, abs=1e-3)

        #Last one
        ts = universe.trajectory[-1]
        atom = structure.Atom.read_from_xtc(selection[0])
        for a, b in zip(atom.coords, [20.63, 38.43, 32.09]):
            assert a == pytest.approx(b, abs=1e-3)
        atom = structure.Atom.read_from_xtc(selection[-1])
        for a, b in zip(atom.coords, [39.14, 39.77, 25.60]):
            assert a == pytest.approx(b, abs=1e-3)


class TestChainClass(object):
    """
    Tests for Chain class in PDBlib
    """

    @staticmethod
    @pytest.fixture
    def chain():
        """
        Run before each test.

        Create a chain object
        """
        lines = ("ATOM    840  C   ARG B  11      22.955  23.561  -4.012  1.00 28.07           C  ",
                 "ATOM    849  N   SER B  12      22.623  24.218  -2.883  1.00 24.77           N  ",
                 "ATOM    850  CA  SER B  12      22.385  23.396  -1.637  1.00 21.99           C  ",
                 "ATOM    851  C   SER B  12      21.150  24.066  -0.947  1.00 32.67           C  ",
                 "ATOM    855  N   ILE B  13      20.421  23.341  -0.088  1.00 30.25           N  ")
        chain = structure.Chain()
        for line in lines:
            atom = structure.Atom.read_from_PDB(line)
            chain.add_atom(atom)
        return chain

    def test_size(self, chain):
        """
        Tests for size()
        """
        assert chain.size() ==  5

    @pytest.mark.parametrize(
        'resid,angles',
        ((11, {'phi': None, 'psi': None}),
         (12, {'phi': -139.77684605036447, 'psi': 157.94348570201197}),
         (13, {'phi': None, 'psi': None}))
    )
    def test_get_phi_psi_angles(self, chain, resid, angles):
        """
        Tests for get_phi_psi_angles()
        """
        phi_psi = chain.get_phi_psi_angles()
        assert angles["phi"] == pytest.approx(phi_psi[resid]["phi"])

    def test_set_coordinates(self, chain):
        """
        Tests for coordinates update
        """
        new_coords = numpy.array([[1.00, 1.00, 1.00],
                                  [2.00, 2.00, 2.00],
                                  [3.00, 3.00, 3.00],
                                  [4.00, 4.00, 4.00],
                                  [5.00, 5.00, 5.00]])
        chain.set_coordinates(new_coords)

        for atom, ref_coords in zip(chain, new_coords):
            numpy.testing.assert_array_almost_equal(atom.coords, ref_coords)

        # Wrong shape
        wrong_coords = new_coords[:-1]
        with pytest.raises(ValueError):
            chain.set_coordinates(wrong_coords)


class TestPDBClass(object):
    """
    Tests for PDB class in Structurelib
    """

    @staticmethod
    @pytest.fixture
    def chains_1BTA():
        filename = os.path.join(here, "test_data/1BTA.pdb")
        pdb = pbx.structure.PDB.PDB(filename)
        return list(pdb.get_chains())

    @staticmethod
    @pytest.fixture
    def chains_1AY7_pdb():
        filename = os.path.join(here, "test_data/1AY7.pdb")
        pdb = pbx.structure.PDB.PDB(filename)
        return list(pdb.get_chains())

    @staticmethod
    @pytest.fixture
    def chains_1AY7_pdbx():
        filename = os.path.join(here, "test_data/1AY7.cif.gz")
        pdb = pbx.structure.PDB.PDB(filename)
        return list(pdb.get_chains())

    @staticmethod
    @pytest.fixture
    def chains_2LFU():
        filename = os.path.join(here, "test_data/2LFU.pdb")
        pdb = pbx.structure.PDB.PDB(filename)
        return list(pdb.get_chains())

    @pytest.mark.parametrize(
        'index,ref',
        ((0, "ATOM      1    N LYS A   1      -8.655   5.770   8.371  0.00  0.00              "),
         (-1, "ATOM   1434   HG SER A  89       6.663  12.440   4.229  0.00  0.00              "))
    )
    def test_read_single_PDB(self, chains_1BTA, index, ref):
        """
        Tests for single chain in one PDB
        """
        chain = chains_1BTA[0]
        assert chain[index].format() == ref

    @pytest.mark.parametrize(
        'chain_idx,ref_first,ref_last',
        ((0,
          "ATOM      1    N ASP A   1      11.860  13.207  12.724  0.00  0.00              ",
          "ATOM    751  OXT CYS A  96       9.922  16.291  36.110  0.00  0.00              "),
         (1,
          "ATOM    753    N LYS B   1      11.318  46.585   0.493  0.00  0.00              ",
          "ATOM   1489  OXT SER B  89      13.857  33.192 -16.133  0.00  0.00              "))
    )
    def test_read_multiple_PDB(self, chains_1AY7_pdb, chain_idx,
                               ref_first, ref_last):
        """
        Tests for multiple chains in one file.

        This test is called for both PDB abd PDBx format because the
        `chains_1AY7` fixture is parametrixed for both extensions.
        """
        chain = chains_1AY7_pdb[chain_idx]
        assert chain[0].format() == ref_first
        assert chain[-1].format() == ref_last

    def test_read_models_PDB(self, chains_2LFU):
        """
        Tests for multiple models in one PDB
        """
        #3 models of one chain
        assert len(chains_2LFU) == 3

        model_3 = chains_2LFU[2]
        ref = "ATOM      1    N ASN A 276     -21.874   9.349   4.010  0.00  0.00              "
        assert model_3[0].format() == ref
        assert model_3.model == "3"

    # This test could be factorized with `test_read_multiple_PDB` by
    # parametrizing the `chains_1AY7_*` fixture. Yet, the atom number of the
    # second chain are shifted by one between the PDB and the PDBx file. This
    # is due to the TER record counting in the sequence of atom numbers in the
    # PDB file.
    @pytest.mark.parametrize(
        'chain_idx,ref_first,ref_last',
        ((0,
          "ATOM      1    N ASP A   1      11.860  13.207  12.724  0.00  0.00              ",
          "ATOM    751  OXT CYS A  96       9.922  16.291  36.110  0.00  0.00              "),
         (1,
          "ATOM    752    N LYS B   1      11.318  46.585   0.493  0.00  0.00              ",
          "ATOM   1488  OXT SER B  89      13.857  33.192 -16.133  0.00  0.00              "))
    )
    def test_read_multiple_PDBx(self, chains_1AY7_pdbx, chain_idx,
                                ref_first, ref_last):
        """
        Tests for multiple chains in one file.

        This test is called for both PDB abd PDBx format because the
        `chains_1AY7` fixture is parametrixed for both extensions.
        """
        chain = chains_1AY7_pdbx[chain_idx]
        assert chain[0].format() == ref_first
        assert chain[-1].format() == ref_last


class TestIolib(object):
    """
    Tests for Iolib
    """

    def test_read_fasta(self):
        """
        Test for parsing mulitple fastas
        """
        filename = os.path.join(here, "test_data/1AY7.pdb.PB.fasta") 
        headers, sequences = pbx.io.read_fasta(filename)
        assert headers == ['test_data/1AY7.pdb | chain A',
                           'test_data/1AY7.pdb | chain B']
        assert sequences == ['ZZbjadfklmcfklmmmmmmmmnnpaafbfkgo'
                             'pacehlnomaccddehjaccdddddehklpnbja'
                             'dcdddfbehiacddfegolaccdddfkZZ',
                             'ZZcddfklpcbfklmmmmmmmmnopafklgoiakl'
                             'mmmmmmmmpacddddddehkllmmmmnnommmmmm'
                             'mmmmmmmmnopacddddZZ']
