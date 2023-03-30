"""
Classes to handle the atom and residue neighborhoods of a structure
"""

import numpy as np
from collections import defaultdict
from ctypes import Union
import warnings

import Bio.PDB as bio
import glycosylator.utils.structural.base as base

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
        return base.compute_angle(self._atoms[0], self._atoms[1], self._atoms[2])

    @property
    def bond_angle_234(self):
        return base.compute_angle(self._atoms[1], self._atoms[2], self._atoms[3])

    @property
    def dihedral(self):
        return base.compute_dihedral(self._atoms[0], self._atoms[1], self._atoms[2], self._atoms[3])

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