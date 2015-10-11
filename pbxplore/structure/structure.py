#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard module
import math

# Third-party modules
import numpy


# =============================================================================
# Classes
# =============================================================================
class AtomError(Exception):
    """
    Exeption class for the Atom class.
    """
    pass


class ChainError(Exception):
    """
    Exeption class for the Chain class
    """
    pass


class Atom:
    """
    Class for atoms in PDB or PDBx/mmCIF format.
    """
    def __init__(self):
        """default constructor"""
        self.id = 0
        self.name = None
        self.resalt = ' '
        self.resname = None
        self.chain = None
        self.resid = 0
        self.resins = ''
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occupancy = 0.0
        self.tempfactor = 0.0
        self.element = '  '
        self.charge =  '  '
        self.model = None

    def read_from_PDB(self, line):
        """
        Read ATOM data from a PDB file line.

        Parameters
        ----------
        line : str
            Line from a PDB file starting with 'ATOM' or 'HETATM'.

        Raises
        ------
        AtomError
            If line is too short.

        Notes
        -----
        PDB format documentation:
        http://www.wwpdb.org/documentation/format33/v3.3.html
        """
        if len(line) < 55:
            raise AtomError("ATOM line too short:\n{0}".format(line))
        self.id = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resid = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())

    def read_from_PDBx(self, line, fields):
        """
        Read ATOM data from a PDBx/mmCIF file line

        Parameters
        ----------
        line : str
            Line from a PDBx/mmCIF file starting with 'ATOM' or 'HETATM'.
        fields : list
            List of str containing fields of data for PDBx/mmCIF format.

        Notes
        -----
        Format documentation:
        http://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html
        """
        try:
            dic = dict(zip(fields, line.split()))
        except:
            raise AtomError("Something went wrong in reading\n{0}".format(line))
        try:
            self.id = int(dic['id'])
            self.name = dic['label_atom_id']
            self.resname = dic['label_comp_id']
            self.chain = dic['label_asym_id']
            self.resid = int(dic['label_seq_id'])
            self.x = float(dic['Cartn_x'])
            self.y = float(dic['Cartn_y'])
            self.z = float(dic['Cartn_z'])
            self.model = dic['pdbx_PDB_model_num']
        except:
            raise AtomError("Something went wrong in data convertion\n{0}"
                            .format(dic))

    def read_from_xtc(self, atm):
        """
        Read ATOM date from a .xtc mdanalysis selection.

        Parameters
        ----------
        atm : atom object
        """
        self.id = atm.id
        self.name = atm.name
        self.resname = atm.resname
        self.chain = ""
        self.resid = atm.resid
        self.x = atm.pos[0]
        self.y = atm.pos[1]
        self.z = atm.pos[2]

    def __repr__(self):
        """
        Atom representation.
        """
        return 'atom {:4d} {:4s} in {:4d} {:3s}' \
               .format(self.id, self.name, self.resid, self.resname)

    def format(self):
        """
        Atom displayed in PDB format.
        """
        return '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' \
               % ('ATOM  ', self.id, self.name, self.resalt, self.resname, self.chain,
                  self.resid, self.resins, self.x, self.y, self.z,
                  self.occupancy, self.tempfactor, self.element, self.charge)

    def coords(self):
        """
        Return atom coordinates.

        """
        return [self.x, self.y, self.z]


class Chain:
    """
    Class to handle PDB chain
    """
    def __init__(self):
        """
        Constructor
        """
        self.name = ""
        self.model = ""
        self.atoms = []

    def __repr__(self):
        """
        Representation
        """
        return "Chain {0} / model {1}: {2} atoms".format(self.name,
                                                         self.model,
                                                         len(self.atoms))

    def add_atom(self, atom):
        """
        Add atom.

        Parameters
        ----------
        atom : object from Atom class
            Atom to be added to chain.

        Raises
        ------
        ChainError
            If the chain has several names.

        """
        # set chain name when first atom is stored
        if not self.atoms:
            self.name = atom.chain
        # check that chain name is always the same
        elif self.name != atom.chain:
            raise ChainError("Several chains are in the same structure")
        # add atom to structure
        self.atoms.append(atom)

    def set_model(self, model):
        """
        Set model number.

        Parameters
        ----------
        model : str
            Model identifier.
        """
        self.model = model

    def size(self):
        """
        Get number of atoms.
        """
        return len(self.atoms)

    def get_phi_psi_angles(self):
        """
        Compute phi and psi angles.

        Returns
        -------
        phi_psi_angles : dict
            Dict with residue number (int) as keys
            and a ``{'phi' : (float), 'psi' : (float)}`` dictionnary as values.

        Examples
        --------
        >>> lines = ("ATOM    840  C   ARG B  11      22.955  23.561  -4.012  1.00 28.07           C  ",
        ...          "ATOM    849  N   SER B  12      22.623  24.218  -2.883  1.00 24.77           N  ",
        ...          "ATOM    850  CA  SER B  12      22.385  23.396  -1.637  1.00 21.99           C  ",
        ...          "ATOM    851  C   SER B  12      21.150  24.066  -0.947  1.00 32.67           C  ",
        ...          "ATOM    855  N   ILE B  13      20.421  23.341  -0.088  1.00 30.25           N  ")
        >>>
        >>> import pbxplore as pbx
        >>> ch = pbx.structure.structure.Chain()
        >>> for line in lines:
        ...     at = pbx.structure.structure.Atom()
        ...     at.read_from_PDB(line)
        ...     ch.add_atom(at)
        ...
        >>> print(ch.get_phi_psi_angles())
        {11: {'phi': None, 'psi': None}, 12: {'phi': -139.77684605036447, 'psi': 157.94348570201197}, 13: {'phi': None, 'psi': None}}

        """
        # extract backbone atoms
        backbone = {}
        for atom in self.atoms:
            if atom.name in ["CA", "C", "O", "N"]:
                resid = atom.resid
                if resid in backbone:
                    backbone[resid][atom.name] = atom
                else:
                    backbone[resid] = {atom.name: atom}

        # get dihedrals
        phi_psi_angles = {}
        for res in sorted(backbone):
            # phi: angle between C(i-1) - N(i) - CA(i) - C(i)
            try:
                phi = get_dihedral(backbone[res-1]["C" ].coords(),
                                   backbone[res  ]["N" ].coords(),
                                   backbone[res  ]["CA"].coords(),
                                   backbone[res  ]["C" ].coords())
            except:
                phi = None
            # psi: angle between N(i) - CA(i) - C(i) - N(i+1)
            try:
                psi = get_dihedral(backbone[res  ]["N" ].coords(),
                                   backbone[res  ]["CA"].coords(),
                                   backbone[res  ]["C" ].coords(),
                                   backbone[res+1]["N" ].coords())
            except:
                psi = None
            # print(res, phi, psi)
            phi_psi_angles[res] = {"phi": phi, "psi": psi}
        return phi_psi_angles


# =============================================================================
# Functions
# =============================================================================
def get_dihedral(atomA, atomB, atomC, atomD):
    """
    Compute dihedral angle between 4 atoms (A, B, C, D).

    Parameters
    ----------
    atomA : list
        Coordinates of atom A as a list or tuple of floats [x, y, z].
    atomB : list
        Coordinates of atom B as a list or tuple of floats [x, y, z].
    atomC : list
        Coordinates of atom C as a list or tuple of floats [x, y, z].
    atomD : list
        Coordinates of atom D as a list or tuple of floats [x, y, z].

    Returns
    -------
    torsion : float
        Torsion angle defined by the atoms A, B, C and D. Angle is defined
        in degrees in the range -180, +180.

    Notes
    -----
    This function is on purpose not part of any class to ease its reusability.

    Examples
    --------
    >>> atom1 = (-1.918, -6.429, -7.107)
    >>> atom2 = (-2.609, -5.125, -7.305)
    >>> atom3 = (-4.108, -5.392, -7.331)
    >>> atom4 = (-4.469, -6.494, -7.911)
    >>> get_dihedral(atom1, atom2, atom3, atom4)
    -36.8942888266
    """

    # convert lists to Numpy objects
    A = numpy.array(atomA)
    B = numpy.array(atomB)
    C = numpy.array(atomC)
    D = numpy.array(atomD)

    # vectors
    AB = B - A
    BC = C - B
    CD = D - C

    # normal vectors
    n1 = numpy.cross(AB, BC)
    n2 = numpy.cross(BC, CD)

    # normalize normal vectors
    n1 /= numpy.linalg.norm(n1)
    n2 /= numpy.linalg.norm(n2)

    # angle between normals
    cosine = numpy.sum(n1*n2) / (numpy.linalg.norm(n1) * numpy.linalg.norm(n2))
    try:
        torsion = math.acos(cosine)
    except:
        cosine = int(cosine)  # +0.0001
        torsion = math.acos(cosine)

    # convert radion to degree
    torsion = torsion * 180.0 / math.pi

    # find if the torsion is clockwise or counterclockwise
    # if numpy.sum(n1 * CD) < 0.0:
    if numpy.dot(n1, CD) < 0.0:
        torsion = 360 - torsion
    if torsion == 360.0:
        torsion = 0.0

    # get range -180 / +180
    if torsion > 180.0:
        torsion = torsion - 360
    if torsion < -180.0:
        torsion = torsion + 360

    return torsion
