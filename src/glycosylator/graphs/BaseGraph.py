"""
The basic Class for Molecular Graphs
"""

import Bio.PDB as bio
import networkx as nx


class BaseGraph(nx.Graph):
    def __init__(self, id, bonds: list):
        super().__init__(bonds)
        self.id = id
        self._structure = None

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

    def _get_structure(self):
        """
        Get the underlying `bio.PDB.Structure` object
        """
        structure = list(self.nodes)[0].get_parent()
        while not isinstance(structure, bio.Structure.Structure):
            structure = structure.get_parent()
        return structure
