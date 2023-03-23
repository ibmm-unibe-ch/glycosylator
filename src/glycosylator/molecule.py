import warnings
import networkx as nx
import numpy as np

from typing import Union

import Bio.PDB as bio

from glycosylator.graphs import AtomGraph, ResidueGraph
import glycosylator.utils as utils
import glycosylator.utils.structural as structural


class TooManyChains(Exception):
    """Raise when Molecule class has an AtomGroup with more than one chain/segment"""


class Molecule:
    def __init__(
        self,
        structure,
        root_atom: Union[str, int, bio.Atom.Atom] = None,
    ):
        self._struct = structure
        self.root_atom = root_atom

        self._AtomGraph = None
        self._ResidueGraph = None

        self._bonds = []

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        root_atom: Union[str, int] = None,
        id: str = None,
    ):
        """
        Read a Molecule from a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        return cls(struct, root_atom)

    @property
    def root_atom(self):
        return self._root_atom

    @root_atom.setter
    def root_atom(self, value):
        if value is None:
            self._root_atom = None
        self._root_atom = self._get_atom(value)

    @property
    def root_residue(self):
        if self._root_atom:
            return self.root_atom.get_parent()

    @property
    def structure(self):
        return self._struct

    @property
    def chains(self):
        return list(self._struct.get_chains())

    @property
    def residues(self):
        return list(self._struct.get_residues())

    @property
    def atoms(self):
        return list(self._struct.get_atoms())

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        raise NotImplementedError

    def get_root(self):
        return self.root_atom

    def set_root(self, atom):
        self.root_atom = atom

    def get_atom_graph(self):
        """
        Generate an AtomGraph for the Molecule
        """
        if self._AtomGraph is None:
            self._AtomGraph = AtomGraph(self._struct.id, self._bonds)
        return self._AtomGraph

    def get_residue_graph(self):
        """
        Generate a ResidueGraph for the Molecule
        """
        if self._ResidueGraph is None:
            self._ResidueGraph = ResidueGraph.from_AtomGraph(self.get_atom_graph())
        return self._ResidueGraph

    def infer_missing_atoms(self, _topology=None):
        """
        Infer missing atoms in the structure

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        """
        structural.fill_missing_atoms(self._struct, _topology)

    def get_standard_bonds(self, _topology=None) -> list:
        """
        Get the standard bonds for the structure

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.apply_standard_bonds(self._struct, _topology)
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

    def infer_bonds(self, max_bond_length: float = None, restrict_residues: bool = True) -> list:
        """
        Infer bonds between atoms in the structure

        Parameters
        ----------
        max_bond_length
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.
        restrict_residues
            Whether to restrict bonds to only those in the same residue.
            If False, bonds between atoms in different residues are also inferred.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.infer_bonds(self._struct, max_bond_length, restrict_residues)
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

    def infer_residue_connections(self, max_bond_length: float = None):
        """
        Infer bonds between atoms that connect different residues in the structure

        Parameters
        ----------
        max_bond_length
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.
        """
        bonds = structural.infer_residue_connections(self._struct, max_bond_length)
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

    def compute_angle(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        atom3: Union[str, int, bio.Atom.Atom],
    ):
        """
        Compute the angle between three atoms

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        atom3
            The third atom

        Returns
        -------
        float
            The angle in degrees
        """
        atom1 = self._get_atom(atom1)
        atom2 = self._get_atom(atom2)
        atom3 = self._get_atom(atom3)
        return structural.compute_angle(atom1, atom2, atom3)

    def _get_atom(self, _atom):
        """
        Get an atom based on its id, serial number or check if it is available...
        """
        if isinstance(_atom, str):
            atoms = {i.id: i for i in self._struct.get_atoms()}
            if not _atom in atoms:
                raise ValueError(f"Atom {_atom} not found in structure")
            return atoms[_atom]
        elif isinstance(_atom, int):
            atoms = {i.serial_number: i for i in self._struct.get_atoms()}
            if not _atom in atoms:
                raise ValueError(f"Atom {_atom} not found in structure")
            return atoms[_atom]
        elif isinstance(_atom, bio.Atom.Atom):
            if _atom in self._struct.get_atoms():
                return _atom
            else:
                warnings.warn(f"Atom {_atom} not found in structure")
                return _atom
        else:
            raise TypeError(f"Invalid type for atom: {type(_atom)}")


#         _validate_atom_group(atom_group)
#         self.atom_group: prody.AtomGroup = atom_group

#         self.root_atom: prody.Atom = root_atom
#         self.root_residue: prody.Residue = self.atom_group[self.root_atom.getChid(), self.root_atom.getResnum()]

#         self.atom_graph: AtomGraph = AtomGraph(self.atom_group, self.root_atom)
#         self.residue_graph: ResidueGraph = ResidueGraph.from_AtomGraph(self.atom_graph)

#         self.guess_angles()
#         self.guess_dihedrals()
#         self.guess_torsionals()

#     @classmethod
#     def from_PDB(cls, pdb: str, root_atom_serial: int, **kwargs):
#         atom_group = prody.parsePDB(pdb, **kwargs)
#         serials = list(atom_group.getSerials())
#         root_atom_index = serials.index(root_atom_serial)
#         return Molecule(atom_group, atom_group[root_atom_index])

#     @property
#     def glycan_topology(self):
#         return self.residue_graph.to_glycan_topology()

#     def write_PDB(self, filename: str, selection: str = "all"):
#         prody.writePDB(filename, self.atom_group.select(selection).toAtomGroup())

#     def guess_angles(self):
#         """Searches for all angles in a molecule based on the connectivity"""
#         self.angles = []
#         for node in self.atom_graph.nodes():
#             self.angles.extend(_find_paths(self.atom_graph, node, 2))

#     def guess_dihedrals(self):
#         """Searches for all dihedrals in a molecule based on the connectivity"""
#         self.dihedrals = []
#         for node in self.atom_graph.nodes():
#             self.dihedrals.extend(_find_paths(self.atom_graph, node, 3))

#     def guess_torsionals(self, hydrogens=True):
#         """Builds a list with all the torsional angles that can rotate
#         Parameters:
#             hydrogens: include torsional involving terminal hydrogens
#         Initializes:
#             torsionals: a list of serial number of atom defining a torsional angle (quadruplet)
#         """
#         cycles = nx.cycle_basis(self.atom_graph.to_undirected(as_view=True))
#         # TODO: can cycle_id just be a flat set of all atoms in any cycles?
#         cycle_id = {atom: i for i, cycle in enumerate(cycles) for atom in cycle}

#         elements = self.atom_group.getElements()
#         torsionals = []
#         for dihedral in self.dihedrals:
#             if dihedral[1] in cycle_id and dihedral[2] in cycle_id:
#                 # skip dihedral if both middle atoms in a cycle
#                 continue
#             if not hydrogens:
#                 if elements[dihedral[0]] == "H" or elements[dihedral[-1]] == "H":
#                     # skip dihedral if either of outer atoms are hydrogen
#                     continue
#             # use the direction of bond_graph to orient dihedral
#             # to then check for uniqueness
#             # i.e. ABCD and DCBA are same and should not be repeated
#             if self.atom_graph.has_edge(dihedral[2], dihedral[1]):
#                 dihedral.reverse()
#             if dihedral not in torsionals:
#                 torsionals.append(dihedral)

#         self.torsionals = torsionals

#     def measure_dihedral_angle(self, dihedral_atoms):
#         """Calculates dihedral angle for 4 atoms.
#         Parameters:
#             torsional: list of atom serial numbers
#         Returns:
#             angle: dihedral angle in degrees
#         """

#         idx = np.argsort(dihedral_atoms)
#         vec_sel = self.atom_group.select(f"index {' '.join(map(str, dihedral_atoms))}")
#         c0, c1, c2, c3 = vec_sel.getCoords()[idx, :]

#         q1 = c1 - c0
#         q2 = c2 - c1
#         q3 = c3 - c2

#         q1xq2 = np.cross(q1, q2)
#         q2xq3 = np.cross(q2, q3)

#         n1 = q1xq2 / np.sqrt(np.dot(q1xq2, q1xq2))
#         n2 = q2xq3 / np.sqrt(np.dot(q2xq3, q2xq3))

#         u1 = n2
#         u3 = q2 / (np.sqrt(np.dot(q2, q2)))
#         u2 = np.cross(u3, u1)

#         cos_theta = np.dot(n1, u1)
#         sin_theta = np.dot(n1, u2)
#         angle = np.degrees(-np.arctan2(sin_theta, cos_theta))

#         return angle

#     def get_all_torsional_angles(self):
#         """Computes all the torsional angles of the molecule
#         Return:
#             angles: angles in degrees
#         """
#         return [self.measure_dihedral_angle(torsional) for torsional in self.torsionals]

#     def rotate_bond(self, torsional, theta, absolute=False):
#         """Rotate the molecule around a torsional angle. Atom affected are in direct graph.
#         Parameters:
#             torsional: index of torsional angle (in torsionals) or list of indices of atoms defining the torsional angle.
#             theta: amount (degrees)
#             absolute: defines if theta is the increment or the absolute value of the angle
#         Returns:
#             c_angle: angle before the rotation in degrees
#         """
#         if type(torsional) == int:
#             torsional = self.torsionals[torsional]
#         elif torsional not in self.torsionals:
#             raise ValueError("Invalid Torsional")

#         atoms = []
#         a1 = torsional[-2]

#         atoms = nx.descendants(self.atom_graph, a1)
#         sel = self.atom_group.select(f"index {' '.join(map(str, atoms))}")
#         t = torsional[1:-1]
#         v1, v2 = self.atom_group.select(f"index {' '.join(map(str, t))}").getCoords()[np.argsort(t), :]
#         axis = v2 - v1
#         c_angle = 0.0
#         if absolute:
#             c_angle = self.measure_dihedral_angle(torsional)
#             theta = theta - c_angle

#         coords = sel.getCoords()
#         M = _rotation_matrix(axis, np.radians(theta))
#         coords = M.dot(coords.transpose())
#         sel.setCoords(coords.transpose() + v2 - np.dot(M, v2))


# def _validate_atom_group(atom_group: prody.AtomGroup):
#     """validates that the AtomGroup consists of a single physical molecule"""
#     pass
#     # if atom_group.numChains() != 1 or atom_group.numSegments() != 1:
#     #     raise TooManyChains()


# def _find_paths(G, node, length, excludeSet=None):
#     """Finds all paths of a given length
#     Parameters:
#         G: graph (netwrokx)
#         node: starting node
#         length: length of path
#         excludedSet: set
#     Returns:
#         paths: list of all paths of a length starting from node
#     """
#     if excludeSet == None:
#         excludeSet = {node}
#     else:
#         excludeSet.add(node)

#     if length == 0:
#         return [[node]]
#     paths = [[node] + path for neighbor in G.neighbors(node) if neighbor not in excludeSet for path in _find_paths(G, neighbor, length - 1, excludeSet)]
#     excludeSet.remove(node)
#     return paths
