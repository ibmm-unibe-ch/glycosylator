import networkx as nx
import numpy as np
import mdtraj as md
import Bio.PDB as bio

import glycosylator.utils.structural as struct
from glycosylator.graphs.BaseGraph import BaseGraph
import glycosylator.graphs.AtomGraph as AtomGraph


class ResidueGraph(BaseGraph):
    """
    A graph representation of residues bonded together as an abstraction of a large contiguous molecule.
    """

    def __init__(self, id, bonds: list):
        super().__init__(id, bonds)
        self._AtomGraph = None
        self._atomic_bonds = {}
        self._residues = {i.id: i for i in self.nodes}

    @classmethod
    def from_AtomGraph(cls, atom_graph):
        """
        Create a ResidueGraph from an AtomGraph.

        Parameters
        ----------
        atom_graph : AtomGraph
            The AtomGraph representation of the molecule

        Returns
        -------
        ResidueGraph
            The ResidueGraph representation of the molecule
        """
        id = atom_graph.id
        # we first get all the atom-resolution bonds between different residues,
        # then we get the participating parents (e.g. residues) and make the edges from them.
        bonds = struct.infer_residue_connections(atom_graph.structure)
        _bonds = [(p1.get_parent(), p2.get_parent()) for p1, p2 in bonds]
        new = cls(id, _bonds)
        new._AtomGraph = atom_graph
        new._atomic_bonds.update({_bond: bond for _bond, bond in zip(_bonds, bonds)})
        return new

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        id=None,
        apply_standard_bonds: bool = True,
        infer_bonds: bool = False,
        max_bond_length: float = None,
        restrict_residues: bool = True,
    ):
        """
        Create a ResidueGraph from a PDB of a single molecule

        Parameters
        ----------
        filename : str
            Path to the PDB file
        id : str
            The ID of the molecule. By default the filename is used.
        apply_standard_bonds : bool
            Whether to apply standard bonds from known molecule connectivities.
        infer_bonds : bool
            Whether to infer bonds from the distance between atoms.
        max_bond_length : float
            The maximum distance between atoms to infer a bond.
        restrict_residues : bool
            Whether to restrict to atoms of the same residue when inferring bonds.

        Returns
        -------
        ResidueGraph
            The ResidueGraph representation of the molecule
        """
        _atom_graph = AtomGraph.from_pdb(
            filename,
            id,
            apply_standard_bonds,
            infer_bonds,
            max_bond_length,
            restrict_residues,
        )

        new = cls.from_AtomGraph(_atom_graph)
        return new

    @property
    def residues(self):
        """
        Get the residues in the molecule.

        Returns
        -------
        list
            The residues in the molecule
        """
        return list(self._residues.values())

    @property
    def atomic_bonds(self):
        """
        Get the atomic-level bonds in the molecule.

        Returns
        -------
        dict
            The atomic-level bonds in the molecule
        """
        return self._atomic_bonds

    def to_AtomGraph(self):
        """
        Convert the ResidueGraph to an AtomGraph.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        return self._AtomGraph

    def get_atomic_bond(self, residue1, residue2) -> tuple:
        """
        Get the atomic-level bond between two residues.

        Parameters
        ----------
        residue1 : Residue or str
            The first residue or it's id
        residue2 : Residue or str
            The second residue or it's id

        Returns
        -------
        tuple
            The atomic bond between the two residues
        """
        if isinstance(residue1, str):
            residue1 = self._residues.get(residue1)
        if isinstance(residue2, str):
            residue2 = self._residues.get(residue2)

        bond = self._atomic_bonds.get((residue1, residue2))
        if bond is None:
            bond = self._atomic_bonds.get((residue2, residue1))
        return bond


if __name__ == '__main__':

    _man = "support/examples/MAN9.pdb"
    _man = AtomGraph.from_pdb(_man)
    man = ResidueGraph.from_AtomGraph(_man)

    import matplotlib.pyplot as plt
    import networkx as nx

    print(len(list(man.bonds)))
    nx.draw(man, with_labels=True, font_weight='bold')
    plt.show()
