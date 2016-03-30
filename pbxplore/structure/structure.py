#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard module
import math
import operator

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
    def __init__(self, ident=0, name=None, resname=None, chain=None, resid=0,
                 x=0.0, y=0.0, z=0.0, model=None):
        """default constructor"""
        self.id = ident
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resid = resid
        self.x = x
        self.y = y
        self.z = z
        self.model = model

    @classmethod
    def read_from_PDB(cls, line):
        """
        Constructor from a PDB file line.

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
        ident = int(line[6:11].strip())
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resid = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        return cls(ident, name, resname, chain, resid, x, y, z)

    @classmethod
    def read_from_PDBx(cls, line, fields):
        """
        Constructor from a PDBx/mmCIF file line

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
            ident = int(dic['id'])
            name = dic['label_atom_id']
            resname = dic['label_comp_id']
            chain = dic['label_asym_id']
            resid = int(dic['label_seq_id'])
            x = float(dic['Cartn_x'])
            y = float(dic['Cartn_y'])
            z = float(dic['Cartn_z'])
            model = dic['pdbx_PDB_model_num']
        except:
            raise AtomError("Something went wrong in data convertion\n{0}"
                            .format(dic))

        return cls(ident, name, resname, chain, resid, x, y, z, model)

    @classmethod
    def read_from_xtc(cls, atm):
        """
        Constructor from a .xtc mdanalysis selection.

        Parameters
        ----------
        atm : atom object of MDAnlysis
        """
        x, y, z = atm.position
        return cls(atm.id, atm.name, atm.resname, "", atm.resid, x, y, z)

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
               % ('ATOM  ', self.id, self.name, ' ', self.resname, self.chain,
                  self.resid, '', self.x, self.y, self.z,
                  0.0, 0.0, '  ', '  ')

    @property
    def coords(self):
        """
        Return atom coordinates.
        """
        return [self.x, self.y, self.z]

    @coords.setter
    def coords(self, pos):
        """
        Set the cartesian coordinates of the atom.

        Parameters
        ----------
        pos: a list or numpy array of 3 elements
        """
        self.x, self.y, self.z = pos


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

    def __getitem__(self, i):
        return self.atoms[i]

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

    def set_coordinates(self, positions):
        """
        Update the coordinates of all atoms in a chain.

        Parameters
        ----------
        positions : a 2D numpy array with a shape of (number of atoms * 3)

        Raises
        ------
        TypeError
            If positions doesn't have the right shape

        """
        if numpy.shape(positions) != (self.size(), 3):
            raise ValueError("Coordinates array doesn't have the good shape.")

        for atm, coords in zip(self.atoms, positions):
            atm.coords = coords

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
                phi = get_dihedral(backbone[res-1]["C" ].coords,
                                   backbone[res  ]["N" ].coords,
                                   backbone[res  ]["CA"].coords,
                                   backbone[res  ]["C" ].coords)
            except:
                phi = None
            # psi: angle between N(i) - CA(i) - C(i) - N(i+1)
            try:
                psi = get_dihedral(backbone[res  ]["N" ].coords,
                                   backbone[res  ]["CA"].coords,
                                   backbone[res  ]["C" ].coords,
                                   backbone[res+1]["N" ].coords)
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

    # vectors
    AB = list(map(operator.sub, atomB, atomA))
    BC = list(map(operator.sub, atomC, atomB))
    CD = list(map(operator.sub, atomD, atomC))

    # normal vectors
    n1 = []
    n1.append(((AB[1] * BC[2]) - (AB[2] * BC[1])))
    n1.append(((AB[2] * BC[0]) - (AB[0] * BC[2])))
    n1.append(((AB[0] * BC[1]) - (AB[1] * BC[0])))
    n2 = []
    n2.append(((BC[1] * CD[2]) - (BC[2] * CD[1])))
    n2.append(((BC[2] * CD[0]) - (BC[0] * CD[2])))
    n2.append(((BC[0] * CD[1]) - (BC[1] * CD[0])))

    n1 = numpy.array(n1)
    n2 = numpy.array(n2)


    # normalize normal vectors
    n1 /= numpy.sqrt(n1.dot(n1))
    n2 /= numpy.sqrt(n2.dot(n2))

    # angle between normals
    cosine = n1.dot(n2)
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
