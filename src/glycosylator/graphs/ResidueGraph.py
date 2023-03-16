import networkx as nx
import numpy as np
import mdtraj as md
import Bio.PDB as bio

import glycosylator.utils as utils
from glycosylator.graphs.BaseGraph import BaseGraph
import glycosylator.graphs.AtomGraph as AtomGraph


class ResidueGraph(BaseGraph):
    """
    A graph representation of residues bonded together as an abstraction of a large contiguous molecule.
    """

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
        bonds = utils.infer_residue_connections(atom_graph.structure)
        bonds = [(p1.get_parent(), p2.get_parent()) for p1, p2 in bonds]
        new = cls(id, bonds)
        new.__AtomGraph = atom_graph
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
        new.__AtomGraph = _atom_graph
        return new

    def to_AtomGraph(self):
        """
        Convert the ResidueGraph to an AtomGraph.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        if hasattr(self, "__AtomGraph"):
            return self.__AtomGraph
        else:
            raise AttributeError("This ResidueGraph does not have an associated AtomGraph.")


if __name__ == '__main__':

    _man = "support/examples/MAN9.pdb"
    man = ResidueGraph.from_pdb(_man)

    import matplotlib.pyplot as plt
    import networkx as nx

    print(len(list(man.bonds)))
    nx.draw(man)
    plt.show()
