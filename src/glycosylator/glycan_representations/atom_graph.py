from abc import abstractmethod

import networkx as nx
import prody


class AtomGraph(nx.DiGraph):
    """A directed graph representation of atoms and bonds in a contiguous molecule.

    Each node is an atom; each edge is a bond between atoms.
    The root of the directed graph is the atom in the molecule (glycan) that would link to a protein.

    This is a subclass of networkx.DiGraph

    :param atom_group: a prody AtomGroup or Selection containing a single contiguous molecule
    :type atom_group: prody.AtomGroup | prody.Selection
    :param root_atom: an atom within atom_group that will be the root of the directed graph
    :type root_atom: prody.Atom
    :param attrs: graph attributes as allowed by networkx.DiGraph
    """

    def __init__(
        self,
        atom_group: prody.AtomGroup | prody.Selection,
        root_atom: prody.Atom,
        **attrs
    ):
        if type(atom_group) == prody.AtomGroup:
            atom_selection = atom_group.select("all")
        elif type(atom_group) == prody.Selection:
            atom_selection = atom_group
            atom_group = atom_selection.getAtomGroup()
        else:
            raise ValueError("atom_group of wrong type")

        root_atom_index = root_atom.getIndex()

        if atom_group.getBonds() is None:
            atom_group.inferBonds()

        # create bond graph to pass to initialiser
        directed_bond_graph = AtomGraph._create_directed_bond_graph(
            atom_selection, root_atom_index
        )
        super().__init__(
            incoming_graph_data=directed_bond_graph,
            atom_group=atom_group,
            atom_selection=atom_selection,
            root_atom=root_atom,
            root_atom_index=root_atom_index,
            **attrs
        )

        self._set_all_node_attributes()

    @classmethod
    def from_PDB(cls, file_path: str, root_atom_serial: int):
        """Create an AtomGraph from a PDB of a single molecule

        :param file_path: a PDB file path
        :type file_path: str
        :param root_atom_serial: an atom's serial number that will become the root of the directed graph
        :type root_atom_serial: int
        :return: AtomGraph of the molecule
        :rtype: AtomGraph
        """
        atom_group = prody.parsePDB(file_path)
        serials = list(atom_group.getSerials())
        root_atom_index = serials.index(root_atom_serial)
        return cls(atom_group, atom_group[root_atom_index])

    def _set_all_node_attributes(self):
        atom_group: prody.AtomGroup = self.graph["atom_group"]
        # need to run getResindices once to initialise these values within the atomgroup
        atom_group.getResindices()
        for atom_index, data in self.nodes(data=True):
            atom = atom_group[atom_index]
            data.update({label: atom.getData(label) for label in atom.getDataLabels()})
            data.update({"atom": atom, "frozen": False})

    @staticmethod
    def _create_directed_bond_graph(
        atom_selection: prody.Selection, root_atom_index: int
    ) -> nx.DiGraph:
        bonds = [bond.getIndices() for bond in atom_selection.iterBonds()]
        bond_graph = nx.Graph()
        bond_graph.add_nodes_from(
            [atom.getIndex() for atom in atom_selection.iterAtoms()]
        )
        bond_graph.add_edges_from(bonds)
        bond_graph = nx.bfs_tree(bond_graph, root_atom_index)
        return bond_graph

    def freeze_bonds(self, bonds: str | list[tuple[int, int]] = "all"):
        """Set the edge attribute 'frozen' to True for bonds

        If bonds is 'all', all edges get 'frozen' set to True
        Otherwise, bonds is a list of the graph's edges to be frozen

        :param bonds: list of edges or 'all', defaults to "all"
        :type bonds: str | list[tuple[int, int]], optional
        """
        if bonds == "all":
            for node1, node2, data in self.edges(data=True):
                data["frozen"] = True
        else:
            for node1, node2 in bonds:
                self[(node1, node2)]["frozen"] = True


if __name__ == "__main__":
    ag = prody.parsePDB("support/examples/man9.pdb")
    x = AtomGraph(ag, ag[0])
    x = AtomGraph.from_PDB("support/examples/man9.pdb", 1)
    print("done")
