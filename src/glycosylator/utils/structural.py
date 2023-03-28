"""
Auxiliary functions related to structure handling
"""

from collections import defaultdict
import os
import pickle
from typing import NamedTuple, Union
import warnings

import numpy as np
import Bio.PDB as bio
from openbabel import pybel
import numba

import glycosylator.utils.defaults as defaults
import glycosylator.utils.constants as constants
import glycosylator.utils.convert as convert

# =================================================================
# Structure related classes
# =================================================================


class Neighborhood:
    """
    This class handles graph connectivity of connected nodes and is the basis
    for the AtomNeighborhood and ResidueNeighborhood classes.

    Parameters
    ----------
    graph
        A networkx graph
    """

    # the method of node objects to use for mapping to an index-based dictionary
    # this assumes that the index is a **unique** identifier for the nodes
    __index_method__ = AttributeError

    # the method of node objects to use for mapping to an id-based dictionary
    # this allows for multiple nodes with the same id (e.g. multiple 'C1' atoms...)
    __id_method__ = AttributeError

    def __init__(self, graph):
        self._src = graph

        # each implementation may call this method before the super init
        # self._validate_getters()

        self._make_node_dicts()

    def get_neighbors(self, node, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of an atom

        Parameters
        ----------
        node : object
            The target node object.
        n : int
            The (maximal) number of edges that should
            separate the target from neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n edges away from the target or closer.
            If "at", this will get all neighbors that are exactly n edges
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the target node
        """
        if isinstance(node, list):
            return [self.get_neighbors(a, n, mode) for a in node]

        self._seen = set()
        if mode == "upto":
            return self._get_neighbors_upto(node, n) - {node}
        elif mode == "at":
            return self._get_neighbors_at(node, n) - {node}
        else:
            raise ValueError(f"Invalid mode: {mode}")

    def _make_node_dicts(self):
        """
        Make a dictionary of nodes, which are biopython objects (or anything that has an id and sequence identifier attribute) mapped to their serial number and ids
        """
        self._idx_nodes = {}
        self._ids_nodes = defaultdict(list)

        for node1, node2 in self._src.edges:

            self._idx_nodes.setdefault(self.__index_method__(node1), node1)
            self._idx_nodes.setdefault(self.__index_method__(node2), node2)

            id1 = self.__id_method__(node1)
            id2 = self.__id_method__(node2)

            if node1 not in self._ids_nodes[id1]:
                self._ids_nodes[id1].append(node1)
            if node2 not in self._ids_nodes[id2]:
                self._ids_nodes[id2].append(node2)

    def _get_node(self, node: Union[int, str]):
        """
        Get a node from its id or sequence index
        """
        if isinstance(node, int):
            node = self._idx_nodes.get(node)
        elif isinstance(node, str):
            node = self._ids_nodes.get(node)
            if node and len(node) > 1:
                warnings.warn(f"Found multiple nodes with the given id. Returning the full list!")
            elif node:
                node = node[0]
            else:
                warnings.warn(f"No nodes found with the given id!")
        else:
            raise TypeError(f"Invalid node type: {type(node)}")
        return node

    def _get_neighbors_upto(self, node, n: int):
        """
        Get all neighbors of a node that are n edges away from the target or closer
        """
        if n == 0:
            return {node}
        else:
            neighbors = set()
            for neighbor in self._src.adj[node]:
                if neighbor in self._seen:
                    continue
                if n >= 1:
                    neighbors.update(self._get_neighbors_upto(neighbor, n - 1))
                neighbors.add(neighbor)
                self._seen.add(neighbor)
            return neighbors

    def _get_neighbors_at(self, node, n: int):
        """
        Get all neighbors of a node that are exactly n edges away from the target
        """
        if n == 0:
            return {node}
        else:
            neighbors = set()
            for neighbor in self._src.adj[node]:
                if neighbor in self._seen:
                    continue
                if n >= 1:
                    neighbors.update(self._get_neighbors_upto(neighbor, n - 1))
                self._seen.add(neighbor)
            return neighbors


class AtomNeighborhood(Neighborhood):
    """
    This class handles the bond connectivity neighborhood of atoms in a structure
    and can be used to obtain bonded atom triplets.

    Parameters
    ----------
    graph
        An AtomGraph
    """

    __index_method__ = __index_method__ = lambda _, node: getattr(node, "serial_number")
    __id_method__ = lambda _, node: getattr(node, "id")

    @property
    def atoms(self):
        """
        Returns the atoms in the structure
        """
        return self._src.atoms

    @property
    def bonds(self):
        """
        Returns the bonds in the structure
        """
        return self._src.bonds

    def get_neighbors(self, atom, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of an atom

        Parameters
        ----------
        atom : int or str or Bio.PDB.Atom
            The serial number or the id of the atom, or the atom itself.
        n : int
            The (maximal) number of bonds that should
            separate the target from the neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n bonds away from the target or closer.
            If "at", this will get all neighbors that are exactly n bonds
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the atom
        """
        if not isinstance(atom, bio.Atom.Atom):
            atom = self._get_node(atom)
        if atom is None:
            return set()
        return super().get_neighbors(atom, n, mode)

    def get_atom(self, atom: Union[int, str]):
        """
        Get an atom from its id or serial number
        """
        return self._get_node(atom)


class ResidueNeighborhood(Neighborhood):
    """
    This class handles the residue connectivity neighborhood of residues in a structure
    and can be used to obtain residue triplets.

    Parameters
    ----------
    graph
        A ResidueGraph
    """

    __index_method__ = __index_method__ = lambda _, node: node.id[1]
    __id_method__ = lambda _, node: getattr(node, "resname")

    @property
    def residues(self):
        """
        Returns the residues in the structure
        """
        return self._src.residues

    @property
    def bonds(self):
        """
        Returns the bonds in the structure
        """
        return self._src.bonds

    def get_neighbors(self, residue, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of a residue

        Parameters
        ----------
        residue : int or str or Bio.PDB.Residue
            The index or the id (resname) of the residue, or the residue itself.
        n : int
            The (maximal) number of bonds that should
            separate the target from the neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n bonds away from the target or closer.
            If "at", this will get all neighbors that are exactly n bonds
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the residue
        """
        if not isinstance(residue, bio.Residue.Residue):
            residue = self._get_node(residue)
        if residue is None:
            return set()
        return super().get_neighbors(residue, n, mode)

    def get_residue(self, residue: Union[int, str]):
        """
        Get a residue from its id or index
        """
        return self._get_node(residue)


class Quartet:
    """
    An atom quartet that can be used to compute internal coordinates
    """

    def __init__(self, atom1, atom2, atom3, atom4, improper: bool = False) -> None:
        self._atoms = (atom1, atom2, atom3, atom4)
        self._improper = improper

    @property
    def atoms(self):
        return self._atoms

    @property
    def atom1(self):
        return self._atoms[0]

    @property
    def atom2(self):
        return self._atoms[1]

    @property
    def atom3(self):
        return self._atoms[2]

    @property
    def atom4(self):
        return self._atoms[3]

    @property
    def improper(self):
        return self._improper

    @property
    def center_atom(self):
        if self._improper:
            return self._atoms[2]
        else:
            raise TypeError("This quartet is not an improper")

    @property
    def bond_length_12(self):
        return np.linalg.norm(self._atoms[1].coord - self._atoms[0].coord)

    @property
    def bond_length_23(self):
        return np.linalg.norm(self._atoms[2].coord - self._atoms[1].coord)

    @property
    def bond_length_34(self):
        return np.linalg.norm(self._atoms[3].coord - self._atoms[2].coord)

    @property
    def bond_angle_123(self):
        return compute_angle(self._atoms[0], self._atoms[1], self._atoms[2])

    @property
    def bond_angle_234(self):
        return compute_angle(self._atoms[1], self._atoms[2], self._atoms[3])

    @property
    def dihedral(self):
        return compute_dihedral(self._atoms[0], self._atoms[1], self._atoms[2], self._atoms[3])

    def __hash__(self) -> int:
        return hash(tuple(sorted(self._atoms))) + hash(self._improper)

    def __eq__(self, other) -> bool:
        if isinstance(other, Quartet):
            return set(self._atoms) == set(other._atoms) and self._improper == other._improper
        elif isinstance(other, (tuple, list)):
            if len(other) == 4:
                return set(self._atoms) == set(other)
            elif len(other) == 5:
                return set(self._atoms) == set(other[:4]) and self._improper == other[4]
        return False

    def __repr__(self) -> str:
        if hasattr(self._atoms[0], "id"):
            return f"Quartet({self._atoms[0].id}, {self._atoms[1].id}, {self._atoms[2].id}, {self._atoms[3].id}, improper={self._improper})"
        else:
            return f"Quartet({self._atoms[0]}, {self._atoms[1]}, {self._atoms[2]}, {self._atoms[3]}, improper={self._improper})"

    def __iter__(self):
        return iter(self._atoms)

    def __getitem__(self, item):
        return self._atoms[item]


# =================================================================
# Structure related functions
# =================================================================


def read_smiles(smiles: str, add_hydrogens: bool = True):
    """
    Read a SMILES string into a Bio.PDB.Structure

    Parameters
    ----------
    smiles : str
        The SMILES string to read
    add_hydrogens : bool
        Whether to add hydrogens to the structure

    Returns
    -------
    structure : Bio.PDB.Structure
        The structure
    """
    mol = pybel.readstring("smi", smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES string {smiles}")
    mol.make3D()

    if not add_hydrogens:
        mol.removeh()

    converter = convert.PybelBioPythonConverter()
    mol = converter.pybel_molecule_to_biopython(mol)

    return mol


def infer_residue_connections(structure, bond_length: float = None):
    """
    Infer the connectivity graph of residues from the distances between atoms of residue pairs.
    This will establish only bonds between close-by atoms from non-identical residues.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host at least two Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
    bond_length : float
        The maximum distance between two atoms to be considered a bond.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms from different Residues that are bonded.
    """
    if bond_length is None:
        bond_length = defaults.DEFAULT_BOND_LENGTH

    bonds = []
    _seen_residues = set()
    for residue1 in structure.get_residues():
        for residue2 in structure.get_residues():
            if residue1 == residue2:
                continue
            elif residue2 in _seen_residues:
                continue

            atoms1 = list(residue1.get_atoms())
            atoms2 = list(residue2.get_atoms())

            _neighbors = bio.NeighborSearch(atoms1 + atoms2)
            _neighbors = _neighbors.search_all(radius=bond_length)

            bonds.extend(
                [n for n in _neighbors if n[0].get_parent() != n[1].get_parent()],
            )

        _seen_residues.add(residue1)

    # bonds = [Bond(i) for i in bonds]

    return bonds


def infer_bonds(structure, bond_length: float = None, restrict_residues: bool = True):
    """
    Generate a connectivity graph by inferring bonds from the distance between atoms.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
    bond_length : float
        The maximum distance between two atoms to be considered a bond.
    restrict_residues : bool
        If set to `True`, only bonds between atoms of the same residue will be considered.

    Returns
    -------
    bonds : list
        The connectivity graph of the molecule, storing tuples of `Bio.PDB.Atom` objects.
    """
    if bond_length is None:
        bond_length = defaults.DEFAULT_BOND_LENGTH

    if restrict_residues:
        bonds = []
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            _neighbors = bio.NeighborSearch(atoms)
            bonds.extend(
                _neighbors.search_all(radius=bond_length),
            )

    else:
        atoms = list(structure.get_atoms())
        _neighbors = bio.NeighborSearch(atoms)
        bonds = _neighbors.search_all(radius=bond_length)

    # # convert to Bond-tuples
    # bonds = [Bond(i) for i in bonds]

    return bonds


def apply_standard_bonds(structure, _topology=None):
    """
    Apply bonds according to loaded standard topologies. This will compute a list of tuples with bonded
    atoms from the same residue. It will not infer residue-residue connections! It is possible to provide a single
    atom as input structure, in which case all bonds that involve said atom are returned.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which hosts Atoms:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
        - `Bio.PDB.Atom`

    _topology : Topology
        The topology to use for applying the bonds.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms that are bonded.
    """
    if structure.level == "R":
        residue = structure

        if _topology is None:
            _topology = defaults.get_default_topology()

        if not _topology.has_residue(residue.resname):
            warnings.warn(f"[ignoring] No abstract residue found in Topology for {residue.resname}!")
            return

        atoms = residue.child_dict
        abstract_residue = _topology.get_residue(residue.resname)
        bonds = [(atoms.get(bond.atom1.id), atoms.get(bond.atom2.id)) for bond in abstract_residue.bonds]
        bonds = [bond for bond in bonds if bond[0] and bond[1]]  # make sure to have no None entries...
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_standard_bonds(residue, _topology)
        bonds = [bond for bond in bonds if bond[0].id == structure.id or bond[1].id == structure.id]
        return bonds

    else:
        bonds = []
        for residue in structure.get_residues():
            bonds.extend(apply_standard_bonds(residue, _topology))
        return bonds


def fill_missing_atoms(structure, _topology=None):
    """
    Fill missing atoms in residues.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`

    _topology
        A specific topology to use for references.
    """
    if structure.level in set(("S", "M", "C")):
        for residue in structure.get_residues():
            fill_missing_atoms(residue, _topology=_topology)
        return
    elif structure.level == "A":
        raise ValueError("Cannot fill missing atoms in an Atom object!")
    elif structure.level != "R":
        raise ValueError(f"Invalid level for structure: {structure.level}")

    if not _topology:
        _topology = defaults.get_default_topology()

    if not _topology.has_residue(structure.resname):
        warnings.warn(f"[ignoring] No abstract residue found in Topology for {structure.resname}!")
        return
    _abstract = _topology.get_residue(structure.resname)

    # first get all atoms that need to be added
    _atoms_to_add = _abstract.get_missing_atoms(structure)

    # now iterate through the atoms and impute them one by one
    # try to keep atoms that could not be imputed in a queue
    # in hopes of imputing them later, but abort if it's no good...
    _fail_timer = 0
    while len(_atoms_to_add) != 0:
        _atom = _atoms_to_add.pop(0)
        _success = compute_atom_coordinate(_atom, structure, _abstract)
        if not _success:
            _atoms_to_add.append(_atom)
            _fail_timer += 1
            if _fail_timer == len(_atoms_to_add) + 1:
                warnings.warn(f"Could not impute atoms: {[i.id for i in _atoms_to_add]}! Aborting...", RuntimeWarning)
                return
            continue
        _fail_timer = 0


def compute_atom_coordinate(atom, residue, abstract):
    """
    Compute the coordinates of an atom from internal coordinates from an AbstractResidue into a structure.

    Parameters
    ----------
    atom : AbstractAtom
        The atom to compute the coordinates for.
    residue : Bio.PDB.Residue
        The target residue wherein to add the atom.
    abstract : AbstractResidue
        The abstract residue to reference data from.

    Returns
    -------
    success: bool
        Whether the atom could successfully be integrated.
        It may be that a specific atom could not be computed due
        to lacking reference atoms. In this case it may be necessary
        to first impute other atoms.
    """

    # get the internal coordinates to impute the atom positions from
    # since we can either compute atoms at position 1 or 4, we will specifically
    # search for ICs with these properties...
    _coord_func = compute_atom1_from_others
    _ics = abstract.get_internal_coordinates(atom, None, None, None, mode="partial")
    if len(_ics) == 0:
        _coord_func = compute_atom4_from_others
        _ics = abstract.get_internal_coordinates(None, None, None, atom, mode="partial")
        if len(_ics) == 0:
            # warnings.warn(f"Cannot find suitable IC for {atom}!")
            return False

    # try finding an IC that has available reference atoms
    for ic in _ics:
        _ref = ic.get_reference_atoms(residue)
        if len(_ref) == 3:
            break
    else:
        # warnings.warn(f"Cannot find suitable IC for {atom}!")
        return False

    # get the coordinates of the reference atoms
    _ref = [i.coord for i in _ref]

    # now use the internal coordinates to compute the new atom coordinate
    _coord = _coord_func(*_ref, ic)

    # create the new atom
    _new_atom = atom.to_biopython()
    _new_atom.set_coord(_coord)

    # add the new atom to the residue
    residue.add(_new_atom)
    return True


def compute_atom1_from_others(coords2, coords3, coords4, ic):
    """
    Compute the coordinates of the first atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : AbstractInternalCoordinates
        The internal coordinates

    Returns
    -------
    coords1 : array-like
        The coordinates of the first atom
    """
    if ic.is_proper:
        return _IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=coords2,
            r=-ic.bond_length_12,
            theta=-np.radians(ic.bond_angle_123),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        _vec = _IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=np.full(3, 0),
            r=1,
            theta=np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
        if ic.bond_length_13:
            _vec *= ic.bond_length_13
        else:
            BC = np.linalg.norm(coords2 - coords3)
            AC = ic.bond_length_12
            AB = np.sqrt(BC**2 + AC**2 - 2 * BC * AC * np.cos(np.radians(ic.bond_angle_123)))
            _vec *= AB
        final = coords3 + _vec
        return final


def compute_atom4_from_others(coords1, coords2, coords3, ic):
    """
    Compute the coordinates of the fourth atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : AbstractInternalCoordinates
        The internal coordinates

    Returns
    -------
    coords4 : array-like
        The coordinates of the fourth atom
    """
    if ic.is_proper:
        return _IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        return _IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )


def center_of_gravity(masses, coords):
    """
    Compute the center of gravity of a molecule.

    Parameters
    ----------
    masses : array-like
        The masses of the atoms as an nx1 vector
    coords : array-like
        The coordinates of the atoms as an nx3 array

    Returns
    -------
    cog : array-like
        The center of gravity
    """
    return np.average(coords, axis=0, weights=masses)


def compute_angle(atom1, atom2, atom3):
    """
    Compute the angle between three atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom

    Returns
    -------
    angle : float
        The angle between the three atoms in degrees
    """
    a = atom1.coord - atom2.coord
    b = atom3.coord - atom2.coord
    return np.degrees(np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))))


def compute_dihedral(atom1, atom2, atom3, atom4):
    """
    Compute the dihedral angle between four atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom
    atom4 : Bio.PDB.Atom
        The fourth atom

    Returns
    -------
    dihedral : float
        The dihedral angle between the four atoms in degrees
    """
    ab = atom1.coord - atom2.coord
    bc = atom3.coord - atom2.coord
    cd = atom4.coord - atom3.coord

    # normalize bc so that it does not influence magnitude of vector
    # rejections that come next
    bc /= np.linalg.norm(bc)

    # vector rejections
    v = ab - np.dot(ab, bc) * bc
    w = cd - np.dot(cd, bc) * bc

    # angle between v and w in radians
    x = np.dot(v, w)
    y = np.dot(np.cross(bc, v), w)
    return np.degrees(np.arctan2(y, x))


def compute_torsional(atom1, atom2, atom3, atom4):
    """
    Compute the torsional angle between four atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom
    atom4 : Bio.PDB.Atom
        The fourth atom

    Returns
    -------
    torsional : float
        The torsional angle between the four atoms in degrees
    """
    ab = atom1.coord - atom2.coord
    bc = atom3.coord - atom2.coord
    cd = atom4.coord - atom3.coord

    # normalize b so that it does not influence magnitude of vector
    # rejections that come next
    bc /= np.linalg.norm(bc)

    # vector rejections
    v = ab - np.dot(ab, bc) * bc
    w = cd - np.dot(cd, bc) * bc

    # angle between v and w in radians
    x = np.dot(v, w)
    y = np.dot(np.cross(bc, v), w)
    return np.degrees(np.arctan2(y, x))


def compute_triplets(bonds: list):
    """
    Compute all possible triplets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds

    Returns
    -------
    triplets : list
        A list of triplets

    Examples
    --------
    ```
    ( 1 )---( 2 )---( 4 )
       \\
       ( 3 )
         |
       ( 5 )
    ```
    >>> bonds = [(1, 2), (1, 3), (2, 4), (3, 5)]
    >>> compute_triplets(bonds)
    [(2, 1, 3), (1, 2, 4), (1, 3, 5)]
    """
    triplets = []
    for i, bond1 in enumerate(bonds):
        atom_11, atom_12 = bond1
        for j, bond2 in enumerate(bonds[i + 1 :]):
            atom_21, atom_22 = bond2
            if atom_11 == atom_21:
                triplets.append((atom_12, atom_11, atom_22))
            elif atom_11 == atom_22:
                triplets.append((atom_12, atom_11, atom_21))
            elif atom_12 == atom_21:
                triplets.append((atom_11, atom_12, atom_22))
            elif atom_12 == atom_22:
                triplets.append((atom_11, atom_12, atom_21))
    return triplets


def compute_quartets(bonds: list):
    """
    Compute all possible quartets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds

    Returns
    -------
    quartets : list
        A list of quartets

    Examples
    --------
    ```
    ( 1 )---( 2 )---( 4 )
               \\
               ( 3 )
                 |
               ( 5 )
    ```
    >>> bonds = [(1, 2), (2, 3), (2, 4), (3, 5)]
    >>> compute_quartets(bonds)
    {Quartet(1, 2, 3, 5, improper=False), Quartet(5, 3, 2, 4, improper=False), Quartet(1, 3, 2, 4, improper=True)}
    """

    triplets = compute_triplets(bonds)
    quartets = set()
    for idx, triplet1 in enumerate(triplets):
        atom_1, atom_2, atom_3 = triplet1

        for triplet2 in triplets[idx + 1 :]:
            atom_4, atom_5, atom_6 = triplet2

            # decision tree to map atoms into quartets

            quartet = None
            if atom_2 == atom_4:

                if atom_1 == atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)

                elif atom_3 == atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)

            elif atom_2 == atom_5:

                if atom_1 == atom_4 or atom_3 == atom_4:
                    quartet = Quartet(atom_1, atom_3, atom_2, atom_6, True)

            elif atom_2 == atom_6:

                if atom_1 == atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)

                elif atom_3 == atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)

            if quartet:
                quartets.add(quartet)

    return quartets


# @numba.njit
def _rotation_matrix(axis, angle):
    """
    Compute the rotation matrix about an arbitrary axis in 3D

    Source
    ------
    Stackoverflow thread: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    Parameters
    ----------
    axis : array-like
        The axis to rotate around
    angle : float
        The angle to rotate by (in radians)

    Returns
    -------
    rotation_matrix : array-like
        The rotation matrix
    """
    # axis = np.asarray(axis)
    # angle = np.asarray(angle)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def _IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    compute the coordinates of a fourth atom from a proper internal coordinate
    system and the other three atom coordinates.

    Parameters
    ----------
    a, b, c : np.ndarray
        coordinates of the other three atoms
    anchor : np.ndarray
        coordinates of the anchor atom relative to which the new coordinate should be calculated
    r : float
        bond length of the new atom relative to the anchor
    theta : float
        bond angle between the new atom and its plane partners
    dihedral : float
        dihedral angle of the internal coordinates
    """
    ab = b - a
    bc = c - b

    # compute normalized bond vectors for available atoms
    ab /= np.linalg.norm(ab)
    bc /= np.linalg.norm(bc)

    # compute plane vector for atoms 1-2-3
    plane_abc = np.cross(ab, bc)
    plane_abc /= np.linalg.norm(plane_abc)

    # rotate the plane vector around the middle bond (2-3) to get the plane 2-3-4
    _rot = _rotation_matrix(bc, dihedral)
    plane_bcd = np.dot(_rot, plane_abc)
    plane_bcd /= np.linalg.norm(plane_bcd)

    # rotate the middle bond around the new plane
    _rot = _rotation_matrix(plane_bcd, theta)
    cd = np.dot(_rot, bc)
    cd /= np.linalg.norm(cd)

    # compute the coordinates of the fourth atom
    d = anchor + r * cd
    return d


if __name__ == '__main__':

    # MANNOSE = bio.PDBParser().get_structure("MAN", "./support/examples/MAN.pdb")

    # from glycosylator.graphs import AtomGraph

    # nei = AtomNeighborhood(AtomGraph.from_biopython(MANNOSE))

    # n2 = nei.get_neighbors("C1", 2)
    # print(n2)

    edges = [(2, 3), (3, 4), (4, 5), (4, 10), (5, 6), (5, 8), (10, 11), (6, 7), (8, 9), (11, 12)]
    edges = [(1, 2), (2, 3), (2, 4), (3, 5)]
    q = compute_quartets(edges)
    print(q)
