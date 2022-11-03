import networkx as nx
import prody
from glycosylator.file_parsers.glycan_topology import GlycanTopology
from glycosylator.glycan_representations.atom_graph import AtomGraph


class ResidueGraph(nx.DiGraph):
    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)

    @staticmethod
    def _set_all_node_attributes(res_digraph: nx.DiGraph):
        atom_group: prody.AtomGroup = res_digraph.graph["atom_group"]
        atom_graph: AtomGraph = res_digraph.graph["atom_graph"]
        residues: list[prody.Residue] = [
            residue for residue in atom_group.iterResidues()
        ]
        for residue_index, data in res_digraph.nodes(data=True):
            residue = residues[residue_index]
            data.update(
                {label: residue.getData(label) for label in residue.getDataLabels()}
            )
            data.update({"residue": residue}),

    @classmethod
    def from_AtomGroup(
        cls, atom_group: prody.AtomGroup | prody.Selection, root_atom: prody.Atom
    ):
        atom_graph = AtomGraph(atom_group, root_atom)
        return cls.from_AtomGraph(atom_graph)

    @classmethod
    def from_AtomGraph(cls, atom_graph: AtomGraph):
        directed_residue_graph = ResidueGraph._create_directed_residue_graph(atom_graph)
        root_residue_index = atom_graph.graph["root_atom"].getResindex()
        root_residue = directed_residue_graph.nodes(data="residue")[root_residue_index]
        return cls(
            incoming_graph_data=directed_residue_graph,
            root_residue_index=root_residue_index,
            root_residue=root_residue,
        )

    @classmethod
    def from_PDB(cls, file_path: str, root_atom_serial: int):
        atom_graph = AtomGraph.from_PDB(file_path, root_atom_serial)
        return cls.from_AtomGraph(atom_graph)

    @classmethod
    def from_GlycanTopology(cls, glycan_topology: GlycanTopology):
        name = glycan_topology.name
        paths = sorted(glycan_topology.paths, key=len)
        patches_to_index = {tuple(path.patches): i for i, path in enumerate(paths)}
        nodes = []
        edges = []
        for i, path in enumerate(paths):
            node_attrs = {
                "res_name": path.residue_name,
                "linking_atom": path.linking_atom,
            }
            node = (i, node_attrs)
            nodes.append(node)

            patches = path.patches
            if len(patches) == 0:
                # there are no edges for root residue
                continue
            *prev_patches, last_patch = patches
            previous_i = patches_to_index[tuple(prev_patches)]

            edge_attrs = {"patch": last_patch}
            edge = (previous_i, i, edge_attrs)
            edges.append(edge)

        g = nx.DiGraph(incoming_graph_data=None, glycan_name=name, root_residue_index=0)
        g.add_nodes_from(nodes)
        g.add_edges_from(edges)
        return cls(g)

    @staticmethod
    def _create_directed_residue_graph(atom_graph: AtomGraph) -> nx.DiGraph:
        atoms = {atom.getIndex(): atom for atom in atom_graph.graph["atom_group"]}
        res_indices = {
            atom.getIndex(): atom.getResindex()
            for atom in atom_graph.graph["atom_group"]
        }
        inter_res_bonds = [
            (
                res_indices[atom_i],
                res_indices[atom_j],
                {"atoms": (atoms[atom_i], atoms[atom_j])},
            )
            for atom_i, atom_j in atom_graph.edges
            if res_indices[atom_i]
            != res_indices[atom_j]  # atoms that are part of different residues
        ]
        residue_graph = nx.DiGraph(
            incoming_graph_data=None,
            atom_group=atom_graph.graph["atom_group"],
            atom_selection=atom_graph.graph["atom_selection"],
            atom_graph=atom_graph,
        )
        residue_graph.add_nodes_from(res_indices.values())
        residue_graph.add_edges_from(inter_res_bonds)
        return residue_graph

    # TODO: rework maybe?
    @property
    def linkage_paths(self):
        paths = nx.shortest_path(self, source=self.graph["root_residue_index"])
        # TODO: exclude the root_residue path as it's the protein backbone residue
        linkage_paths = {self.graph["root_residue_index"]: " "}
        for path_target, path in paths.items():
            edges = zip(path[:-1], path[1:])
            value = [self[node1][node2]["patch"] for node1, node2 in edges]
            linkage_paths[path_target] = " ".join(value)
        return linkage_paths


if __name__ == "__main__":
    ag = prody.parsePDB("support/examples/man9.pdb")
    x = ResidueGraph.from_AtomGroup(ag, ag[0])
    x = ResidueGraph.from_PDB("support/examples/man9.pdb", 1)
    print("done")
