import networkx as nx
import prody
from glycosylator.glycan_representations.atom_graph import AtomGraph
from prody import Atom, AtomGroup, Residue


class ResidueGraph(nx.DiGraph):
    def __init__(self, atom_graph: AtomGraph, **attr):
        directed_residue_graph = ResidueGraph._create_directed_residue_graph(atom_graph)
        super().__init__(
            incoming_graph_data=directed_residue_graph,
            atom_group=atom_graph.graph["atom_group"],
            atom_graph=atom_graph,
            **attr
        )

        self._set_all_node_attributes()

        root_residue_index = atom_graph.graph["root_atom"].getResindex()
        self.graph.update(
            {
                "root_residue_index": root_residue_index,
                "root_residue": self.nodes(data="residue")[root_residue_index],
            }
        )

    def _set_all_node_attributes(self):
        atom_group: AtomGroup = self.graph["atom_group"]
        atom_graph: AtomGraph = self.graph["atom_graph"]
        residues: list[Residue] = [residue for residue in atom_group.iterResidues()]
        for residue_index, data in self.nodes(data=True):
            residue = residues[residue_index]
            data.update(
                {label: residue.getData(label) for label in residue.getDataLabels()}
            )
            data.update({"residue": residue}),

    @classmethod
    def from_AtomGroup(
        cls, atom_group: AtomGroup, root_atom: Atom = None, root_atom_index: int = None
    ):
        return cls(AtomGraph(atom_group, root_atom, root_atom_index))

    @staticmethod
    def _create_directed_residue_graph(atom_graph: AtomGraph) -> nx.DiGraph:
        res_indices = atom_graph.graph["atom_group"].getResindices()
        inter_res_bonds = [
            (res_indices[atom_i], res_indices[atom_j])
            for atom_i, atom_j in atom_graph.edges
            if res_indices[atom_i]
            != res_indices[atom_j]  # atoms that are part of different residues
        ]
        residue_graph = nx.DiGraph()
        residue_graph.add_edges_from(inter_res_bonds)
        return residue_graph


if __name__ == "__main__":
    ag = prody.parsePDB("support/examples/man9.pdb")
    x = ResidueGraph.from_AtomGroup(ag, ag[0])
    print("done")
