from abc import abstractmethod

import networkx as nx
import prody


class AtomGraph(nx.DiGraph):
    def __init__(
        self,
        atom_group: prody.AtomGroup,
        root_atom: prody.Atom = None,
        root_atom_index: int = None,
        **attr
    ):
        if root_atom is None:
            root_atom = atom_group[root_atom_index]
        if root_atom_index is None:
            root_atom_index = root_atom.getIndex()

        # numFragments needs bond information;
        # .inferBonds() will set bond information if none is present
        if atom_group.getBonds() is None:
            atom_group.inferBonds()
        if atom_group.numFragments() != 1:
            raise ValueError("atom_group contains more than one fragment")

        # create bond graph to pass to initialiser
        directed_bond_graph = AtomGraph._create_directed_graph(
            atom_group, root_atom_index
        )

        # initialise self
        # this will populate all the edge, node info
        super().__init__(
            incoming_graph_data=directed_bond_graph,
            atom_group=atom_group,
            root_atom=root_atom,
            root_atom_index=root_atom_index,
            **attr
        )

    # def _edge_attr_dict_factory(self):
    #     ...
    # edge_attr_dict_factory=_edge_attr_dict_factory

    @staticmethod
    def _create_directed_graph(
        atom_group: prody.AtomGroup, root_atom_index: int
    ) -> nx.DiGraph:
        bonds = [bond.getIndices() for bond in atom_group.getBonds()]
        bond_graph = nx.Graph()
        bond_graph.add_edges_from(bonds)
        bond_graph = nx.bfs_tree(bond_graph, root_atom_index)
        return bond_graph
