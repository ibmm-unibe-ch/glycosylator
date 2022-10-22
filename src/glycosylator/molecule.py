import math

import networkx as nx
import numpy as np
import prody
from prody import Atom, AtomGroup

from .glycan_representations import AtomGraph, ResidueGraph


class TooManyChains(Exception):
    """Raise when Molecule class has an AtomGroup with more than one chain/segment"""


class Molecule:
    def __init__(
        self,
        atom_group: AtomGroup,
        root_atom: Atom,
    ):
        _validate_atom_group(atom_group)
        self.atom_group: prody.AtomGroup = atom_group

        self.root_atom: prody.Atom = root_atom
        self.root_residue: prody.Residue = self.atom_group[
            self.root_atom.getChid(), self.root_atom.getResnum()
        ]

        self.bond_graph: AtomGraph = AtomGraph(self.atom_group, self.root_atom)
        self.residue_graph: ResidueGraph = ResidueGraph(self.bond_graph)

        self.guess_angles()
        self.guess_dihedrals()
        self.guess_torsionals()

    # @classmethod
    # def from_PDB(cls, pdb: str, root_atom: int = 1, **kwargs):
    #     atom_group = prody.parsePDB(pdb, **kwargs)
    #     return Molecule(atom_group, root_atom)

    def write_PDB(self, filename: str, selection: str = "all"):
        prody.writePDB(filename, self.atom_group.select(selection))

    def guess_angles(self):
        """Searches for all angles in a molecule based on the connectivity"""
        self.angles = []
        for node in self.bond_graph.nodes():
            self.angles.extend(_find_paths(self.bond_graph, node, 2))

    def guess_dihedrals(self):
        """Searches for all dihedrals in a molecule based on the connectivity"""
        self.dihedrals = []
        for node in self.bond_graph.nodes():
            self.dihedrals.extend(_find_paths(self.bond_graph, node, 3))

    def guess_torsionals(self, hydrogens=True):
        """Builds a list with all the torsional angles that can rotate
        Parameters:
            hydrogens: include torsional involving terminal hydrogens
        Initializes:
            torsionals: a list of serial number of atom defining a torsional angle (quadruplet)
        """
        cycles = nx.cycle_basis(self.bond_graph.to_undirected(as_view=True))
        # TODO: can cycle_id just be a flat set of all atoms in any cycles?
        cycle_id = {atom: i for i, cycle in enumerate(cycles) for atom in cycle}

        elements = self.atom_group.getElements()
        torsionals = []
        for dihedral in self.dihedrals:
            if dihedral[1] in cycle_id and dihedral[2] in cycle_id:
                # skip dihedral if both middle atoms in a cycle
                continue
            if not hydrogens:
                if elements[dihedral[0]] == "H" or elements[dihedral[-1]] == "H":
                    # skip dihedral if either of outer atoms are hydrogen
                    continue
            # use the direction of bond_graph to orient dihedral
            # to then check for uniqueness
            # i.e. ABCD and DCBA are same and should not be repeated
            if self.bond_graph.has_edge(dihedral[2], dihedral[1]):
                dihedral.reverse()
            if dihedral not in torsionals:
                torsionals.append(dihedral)

        self.torsionals = torsionals

    def measure_dihedral_angle(self, dihedral_atoms):
        """Calculates dihedral angle for 4 atoms.
        Parameters:
            torsional: list of atom serial numbers
        Returns:
            angle: dihedral angle in degrees
        """

        idx = np.argsort(dihedral_atoms)
        vec_sel = self.atom_group.select(f"index {' '.join(map(str, dihedral_atoms))}")
        c0, c1, c2, c3 = vec_sel.getCoords()[idx, :]

        q1 = c1 - c0
        q2 = c2 - c1
        q3 = c3 - c2

        q1xq2 = np.cross(q1, q2)
        q2xq3 = np.cross(q2, q3)

        n1 = q1xq2 / np.sqrt(np.dot(q1xq2, q1xq2))
        n2 = q2xq3 / np.sqrt(np.dot(q2xq3, q2xq3))

        u1 = n2
        u3 = q2 / (np.sqrt(np.dot(q2, q2)))
        u2 = np.cross(u3, u1)

        cos_theta = np.dot(n1, u1)
        sin_theta = np.dot(n1, u2)
        angle = np.degrees(-np.arctan2(sin_theta, cos_theta))

        return angle

    def get_all_torsional_angles(self):
        """Computes all the torsional angles of the molecule
        Return:
            angles: angles in degrees
        """
        return [self.measure_dihedral_angle(torsional) for torsional in self.torsionals]

    def rotate_bond(self, torsional, theta, absolute=False):
        """Rotate the molecule around a torsional angle. Atom affected are in direct graph.
        Parameters:
            torsional: index of torsional angle (in torsionals) or list of indices of atoms defining the torsional angle.
            theta: amount (degrees)
            absolute: defines if theta is the increment or the absolute value of the angle
        Returns:
            c_angle: angle before the rotation in degrees
        """
        if type(torsional) == int:
            torsional = self.torsionals[torsional]
        elif torsional not in self.torsionals:
            raise ValueError("Invalid Torsional")

        atoms = []
        a1 = torsional[-2]

        atoms = nx.descendants(self.bond_graph, a1)
        sel = self.atom_group.select(f"index {' '.join(map(str, atoms))}")
        t = torsional[1:-1]
        v1, v2 = self.atom_group.select(f"index {' '.join(map(str, t))}").getCoords()[
            np.argsort(t), :
        ]
        axis = v2 - v1
        c_angle = 0.0
        if absolute:
            c_angle = self.measure_dihedral_angle(torsional)
            theta = theta - c_angle

        coords = sel.getCoords()
        M = _rotation_matrix(axis, np.radians(theta))
        coords = M.dot(coords.transpose())
        sel.setCoords(coords.transpose() + v2 - np.dot(M, v2))


def _validate_atom_group(atom_group: prody.AtomGroup):
    """validates that the AtomGroup consists of a single physical molecule"""
    pass
    # if atom_group.numChains() != 1 or atom_group.numSegments() != 1:
    #     raise TooManyChains()


def _find_paths(G, node, length, excludeSet=None):
    """Finds all paths of a given length
    Parameters:
        G: graph (netwrokx)
        node: starting node
        length: length of path
        excludedSet: set
    Returns:
        paths: list of all paths of a length starting from node
    """
    if excludeSet == None:
        excludeSet = {node}
    else:
        excludeSet.add(node)

    if length == 0:
        return [[node]]
    paths = [
        [node] + path
        for neighbor in G.neighbors(node)
        if neighbor not in excludeSet
        for path in _find_paths(G, neighbor, length - 1, excludeSet)
    ]
    excludeSet.remove(node)
    return paths


def _rotation_matrix(axis, theta):
    """Computes the rotation matrix about an arbitrary axis in 3D
    Code from: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    Parameters:
        axis: axis
        theta: rotation angle
    Return:
        rotation matrix

    """
    # TODO: replace math with np
    # TODO: is there a single np function to do this?
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )
