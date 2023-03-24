"""
The basic Class for Molecular Graphs
"""

from abc import abstractmethod

import Bio.PDB as bio
import networkx as nx


class BaseGraph(nx.Graph):
    def __init__(self, id, bonds: list):
        if len(bonds) == 0:
            raise ValueError("Cannot create a graph with no bonds")
        super().__init__(bonds)
        self.id = id
        self._structure = None
        self._neighborhood = None

    @property
    def structure(self):
        """
        Returns the underlying `bio.PDB.Structure` object
        """
        if not self._structure:
            self._structure = self._get_structure()
        return self._structure

    @property
    def residues(self):
        """
        Returns the residues in the molecule
        """
        return list(self.structure.get_residues())

    @property
    def atoms(self):
        """
        Returns the atoms in the molecule
        """
        return list(self.structure.get_atoms())

    @property
    def bonds(self):
        """
        Returns the bonds in the molecule
        """
        return list(self.edges)

    def to_pdb(self, filename: str):
        """
        Save to a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        structure = bio.Structure.Structure(self.id)
        model = bio.Model.Model(0)
        chain = bio.Chain.Chain("A")
        structure.add(model)
        model.add(chain)

        for atom in self._structure:
            chain.add(atom)

        io = bio.PDBIO()
        io.set_structure(structure)
        io.save(filename)

    @abstractmethod
    def get_neighbors(self, node, n: int = 1, mode="upto"):
        """
        Get the neighbors of a node

        Parameters
        ----------
        node
            The target node

        n : int, optional
            The number of edges to separate the node from its neighbors.

        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` edges
            - "exact": get all neighbors exactly `n` edges away

        Returns
        -------
        set
            The neighbors of the node
        """
        raise NotImplementedError

    def _get_structure(self):
        """
        Get the underlying `bio.PDB.Structure` object
        """
        if not hasattr(list(self.nodes)[0], "get_parent"):
            raise ValueError("Nodes are not Biopython entities with linked parents!")
        structure = list(self.nodes)[0].get_parent()
        while not isinstance(structure, bio.Structure.Structure):
            structure = structure.get_parent()
        return structure
