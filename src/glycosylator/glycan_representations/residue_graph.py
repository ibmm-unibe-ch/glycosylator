import networkx as nx
import prody

from glycosylator.glycan_representations import AtomGraph
from glycosylator.glycan_representations.glycan_topology import GlycanTopology


class ResidueGraph(nx.DiGraph):
    """A directed graph representation a molecule's residues and their connections

    Each node is a residue; each edge is bond between residues.
    The root of the graph is the residue that would link to a protein or is itself a protein's amino acid.

    This is a subclass of networkx.DiGraph
    """

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
        """Creates a ResidueGraph directly from a prody.AtomGroup or prody.Selection containing a single molecule

        :param atom_group: a prody AtomGroup or Selection containing a single contiguous molecule
        :type atom_group: prody.AtomGroup | prody.Selection
        :param root_atom: the atom within atom_group defining the residue that will be the root of the directed graph
        :type root_atom: prody.Atom
        :return: ResidueGraph instance of the molecule in atom_group
        :rtype: ResidueGraph
        """
        atom_graph = AtomGraph(atom_group, root_atom)
        return cls.from_AtomGraph(atom_graph)

    @classmethod
    def from_AtomGraph(cls, atom_graph: AtomGraph):
        """Create a ResidueGraph directly from an AtomGraph

        :param atom_graph: AtomGraph respresentation of a molecule
        :type atom_graph: AtomGraph
        :return: ResidueGraph instance of the molecule in atom_graph
        :rtype: ResidueGraph
        """
        directed_residue_graph = ResidueGraph._create_directed_residue_graph(atom_graph)
        cls._set_all_node_attributes(directed_residue_graph)
        root_residue_index = atom_graph.graph["root_atom"].getResindex()
        root_residue = directed_residue_graph.nodes[root_residue_index]["residue"]
        return cls(
            incoming_graph_data=directed_residue_graph,
            root_residue_index=root_residue_index,
            root_residue=root_residue,
        )

    @classmethod
    def from_PDB(cls, file_path: str, root_atom_serial: int):
        """Create a ResidueGraph instance directly from a PDB file

        :param file_path: PDB file path containing a single molecule
        :type file_path: str
        :param root_atom_serial: an atom's serial number defining the residue that will become the root of the directed graph
        :type root_atom_serial: int
        :return: ResidueGraph instance of the molecule
        :rtype: ResidueGraph
        """
        atom_graph = AtomGraph.from_PDB(file_path, root_atom_serial)
        return cls.from_AtomGraph(atom_graph)

    @classmethod
    def from_GlycanTopology(cls, glycan_topology: GlycanTopology):
        """Create a ResidueGraph instance from a GlycanTopology

        Use the information within GlycanTopology to create a ResidueGraph. The resultant graph will not contain any associated atomic coordinates.

        :param glycan_topology: GlycanTopology instance
        :type glycan_topology: GlycanTopology
        :return: ResidueGraph instance of the glycan described by the topology
        :rtype: ResidueGraph
        """
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

    def identify_patches(self, glycosylator):
        for res1, res2, data in self.edges(data=True):
            atom_names = {atom.getName() for atom in data["atoms"]}
            patch = glycosylator.find_patch(atom_names)
            data["patch"] = patch

    def path_to_index(self, path: list[str]):
        res_i = self.graph["root_residue_index"]
        for reference_patch in path:
            for neighbour_i, data in self.succ[res_i].items():
                if data["patch"] == reference_patch:
                    res_i = neighbour_i
                    break
            # runs only if loop doesn't break
            else:
                raise ValueError("Patch not found")
        return res_i

    def to_glycan_topology(self, exclude_protein_root: bool = True):
        """Create a GlycanTopology instance

        Option to exclude the amino acid root if present.

        :param exclude_protein_root: Boolean to exclude or include the amino acid root, defaults to True
        :type exclude_protein_root: bool, optional
        :return: GlycanTopology instance
        :rtype: GlycanTopology
        """
        # TODO: exclude the root_residue path as it's the protein backbone residue
        root_res_name = self.graph["root_residue"].getResname()
        # default value UnnamedGlycan is name not present
        glycan_name = self.graph.get("name", "UnnamedGlycan")

        l = []
        paths = nx.shortest_path(self, source=self.graph["root_residue_index"])
        for path_target, path in paths.items():
            res_name = self.nodes[path_target]["residue"].getResname()
            edges = zip(path[:-1], path[1:])
            patches = [self.edges[node1, node2]["patch"] for node1, node2 in edges]
            if len(path) < 2:
                linking_atom = ""
            else:
                linking_atom = self.edges[path[-2], path[-1]]["atoms"][1].getName()
            l.append((res_name, linking_atom, patches))

        # if root is an amino acid e.g. ASN, remove it from dict, and shorten paths by 1
        if exclude_protein_root and self.graph["root_residue"].getSegname() == "PROT":
            l.pop(0)
            for _, _, patches in l:
                # 0th patch is patch linking protein to glycan
                patches.pop[0]

        return GlycanTopology.from_tuples(glycan_name, l)


if __name__ == "__main__":
    ag = prody.parsePDB("support/examples/man9.pdb")
    x = ResidueGraph.from_AtomGroup(ag, ag[0])
    x = ResidueGraph.from_PDB("support/examples/man9.pdb", 1)
    print("done")
