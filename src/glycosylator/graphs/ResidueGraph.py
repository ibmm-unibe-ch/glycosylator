from typing import Union
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
        for r in self.nodes:
            r.coord = r.center_of_mass()
        nx.set_node_attributes(
            self, {i: i.center_of_mass() for i in self.nodes}, "coord"
        )

    @classmethod
    def from_molecule(cls, mol, detailed: bool = True, locked: bool = True):
        """
        Create a ResidueGraph from a molecule object.

        Parameters
        ----------
        mol : Molecule
            The molecule object
        detailed : bool
            Whether to make a "detailed" residue graph representation
            including the atomic-scale bonds between residues. If True,
            locked bonds can be directly migrated from the molecule.
        locked : bool
            Whether to migrate locked bonds from the molecule. This
            is only possible if detailed is True.

        Returns
        -------
        ResidueGraph
            The ResidueGraph representation of the molecule
        """
        atomgraph = mol.make_atom_graph()
        new = cls.from_AtomGraph(atomgraph)
        if detailed:
            new.make_detailed()
            if locked:
                new._locked_edges.update(
                    (i for i in mol._locked_bonds if i in new.edges)
                )
        return new

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
        main_connections = struct.infer_residue_connections(
            atom_graph.structure, triplet=False
        )
        bonds = struct.infer_residue_connections(atom_graph.structure, triplet=True)
        _bonds = [(p1.get_parent(), p2.get_parent()) for p1, p2 in main_connections]
        new = cls(id, _bonds)
        new._AtomGraph = atom_graph

        for bond in bonds:
            parent_1 = bond[0].get_parent()
            parent_2 = bond[1].get_parent()
            new._atomic_bonds.setdefault((parent_1, parent_2), []).append(bond)

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

    def make_detailed(self):
        """
        Use a detailed representation of the residues in the molecule by adding the specific atoms
        that connect the residues together. This is useful for visualization and analysis.

        Note
        ----
        This function is not reversible.
        """

        for residue_pair, bonds in self._atomic_bonds.items():
            if residue_pair in self.edges:
                self.remove_edge(*residue_pair)

            for atom1, atom2 in bonds:

                # add edges between the atoms
                self.add_edge(atom1, atom2)

                p1 = atom1.get_parent()
                p2 = atom2.get_parent()

                # and selectively connect either atom to their center of mass residue
                if p1 != p2:
                    self.add_edge(atom2, p2)
                else:
                    self.add_edge(atom1, p1)

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

    def get_neighbors(
        self, residue: Union[int, str, bio.Residue.Residue], n: int = 1, mode="upto"
    ):
        """
        Get the neighbors of a residue

        Parameters
        ----------
        residue : int, str, bio.Residue.Residue
            The target residue
        n : int, optional
            The number of connections to separate the residue from its neighbors.
        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` bonds
            - "exact": get all neighbors exactly `n` bonds away

        Returns
        -------
        set
            The neighbors of the residue
        """
        if not self._neighborhood:
            self._neighborhood = struct.ResidueNeighborhood(self)
        return self._neighborhood.get_neighbors(residue, n, mode)

    def centers_of_mass(self):
        """
        Get the centers of mass of the residues in the molecule.

        Returns
        -------
        dict
            The centers of mass of the residues in the molecule
        """
        return {residue.id: residue.center_of_mass() for residue in self.residues}


if __name__ == "__main__":

    _man = "support/examples/MAN9.pdb"
    _man = AtomGraph.from_pdb(_man)
    man = ResidueGraph.from_AtomGraph(_man)

    import matplotlib.pyplot as plt
    import networkx as nx

    man.make_detailed()

    import glycosylator.utils.visual as vis

    v = vis.MoleculeViewer3D(man)
    v.show()

    print(len(list(man.bonds)))
    nx.draw(man, with_labels=True, font_weight="bold", pos=nx.spectral_layout(man))
    plt.show()
