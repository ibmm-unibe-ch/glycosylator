import copy
import re
import sqlite3
from itertools import zip_longest

import networkx as nx
import prody
from prody import AtomGroup, Residue, Selection
from typing_extensions import Self

from glycosylator.file_parsers import glycan_topology
from glycosylator.file_parsers.glycan_topology import GlycanTopology
from glycosylator.glycan_representations import ResidueGraph
from glycosylator.molecule import Molecule
from glycosylator.molecule_builder import MoleculeBuilder


class Glycosylator:
    def __init__(
        self,
        charmm_topology_path: str,
        charm_parameter_path: str,
    ):
        self.charmm_topology_path = charmm_topology_path
        self.charmm_parameter_path = charm_parameter_path
        self.builder = MoleculeBuilder(charmm_topology_path, charm_parameter_path)

    def load_glycoprotein_from_PDB(self, protein_file_path: str):
        glycoprotein_atom_group = prody.parsePDB(protein_file_path)
        self.load_glycoprotein_from_AtomGroup(glycoprotein_atom_group)

    def load_glycoprotein_from_AtomGroup(
        self, glycoprotein_atom_group: prody.AtomGroup
    ):
        self.glycoprotein = glycoprotein_atom_group
        self.residues = [residue for residue in self.glycoprotein.iterResidues()]
        self.linking_residues = self.find_linking_residues()
        if self.glycoprotein.getBonds() is None:
            self._infer_glycan_bonds()

    def find_linking_residues(self) -> list[Residue]:
        sequon_pattern = re.compile("(N[^P][S|T])")
        linking_residues = []
        for chain in self.glycoprotein.iterChains():
            residues = [res for res in chain.iterResidues()]
            matches = sequon_pattern.finditer(chain.getSequence())
            linking_residues.extend([residues[match.start()] for match in matches])
        return linking_residues

    def find_existing_glycans(self, freeze_bonds=True):
        # this graph will contain many disjoint subgraphs,
        # each containing a continuous glycan+linking_residue
        glycan_graphs = nx.Graph()
        glycan_graphs.add_nodes_from(
            ((atom.getIndex()) for atom in self._glycans_and_links_sel.iterAtoms())
        )
        glycan_graphs.add_edges_from(
            (bond.getIndices() for bond in self._glycans_and_links_sel.iterBonds())
        )
        # separate the contiguous components into separate graphs
        glycan_graphs = [
            glycan_graphs.subgraph(c).copy()
            for c in nx.connected_components(glycan_graphs)
        ]
        glycans = []
        for graph in glycan_graphs:
            atoms = [self.glycoprotein[index] for index in graph]
            sel_string = (
                f"(serial {' '.join([str(atom.getSerial()) for atom in atoms])})"
            )
            selection = self.glycoprotein.select(sel_string)

            # the selection contains only one amino acid, so this returns just one atom
            root_atom = next(selection.select("protein and name CA").iterAtoms())
            glycan = Molecule(selection, root_atom=root_atom)
            glycan.atom_graph.freeze_bonds("all")
            glycans.append(glycan)

        return glycans

    def identify_glycan(self, glycan: Molecule):
        residue_graph = glycan.residue_graph
        linkage_paths = glycan.residue_graph.linkage_paths

        target = {
            f"{residue_graph.nodes[residue_index]['resname']} {path}"
            for residue_index, path in linkage_paths.items()
        }

        for glycan_name, glycan_linkage_paths in self.glycan_keys.items():
            if target == set(glycan_linkage_paths):
                glycan.name = glycan_name
                return glycan_name
        return ""

    def check_protein_glycan_linkage(self, glycan: Molecule, patches: list[str]):
        residue_graph = glycan.residue_graph
        root_residue_index = residue_graph.graph["root_residue_index"]
        try:
            next_residue_index = next(residue_graph.successors(root_residue_index))
        # if no successors, return False
        # happens when the glycan is just the protein linking residue e.g. ASN
        except StopIteration:
            return False

        linking_atoms = residue_graph[root_residue_index][next_residue_index]["atoms"]
        atom_names = {atom.getName() for atom in linking_atoms}
        for patch in patches:
            # patch_atoms example: ['2C1', '1ND2']
            # str starting with "1" is atom from a protein residue
            # str starting with "2" is atom from a glycan residue
            patch_atoms = self.builder.Topology.patches[patch]["BOND"][0:2]
            patch_atoms = set([atom[1:] for atom in patch_atoms])
            if atom_names == patch_atoms:
                return True
        return False

    def assign_patches(self, glycan: Molecule, remove_unknown=True):
        residue_graph = glycan.residue_graph
        for res1, res2, data in residue_graph.edges(data=True):
            atom_names = {atom.getName() for atom in data["atoms"]}
            patch = self.find_patch(atom_names)
            data["patch"] = patch
            if remove_unknown:
                residue_graph.remove_edge(res1, res2)
                # glycan.del_atoms

    def find_patch(self, atom_names: set[str]) -> str | None:
        # TODO: add anomer support
        key = "-".join(sorted(atom_names))
        if key in self.builder.Topology.atomnames_to_patch:
            return self.builder.Topology.atomnames_to_patch[key][0]
        else:
            return None

    def _infer_glycan_bonds(self):
        glycans_sel = "(not protein and not water)"
        linking_res_sel = f"(resnum {' '.join([str(residue.getResnum()) for residue in self.linking_residues])})"
        self._glycans_and_links_sel = self.glycoprotein.select(
            f"{glycans_sel} or {linking_res_sel}"
        )
        # self.glycoprotein.setBonds(
        #     [
        #         bond.getIndices()
        #         for bond in self._glycans_and_links_sel.toAtomGroup().inferBonds()
        #     ]
        # )
        # kdtree = prody.KDTree(self._glycans_and_links_sel.getCoords())
        # kdtree.search(1.7)
        # self.glycoprotein.setBonds(kdtree.getIndices())

        # can't get bond search for selection to work, so have to brute force whole glycoprotein
        # which is slow
        self.glycoprotein.inferBonds()

    def glycosylate(
        self,
        protein_residue: Residue = None,
        protein_link_patch: str = None,
        glycan: Molecule = None,
        template_glycan: str | GlycanTopology = None,
    ) -> Molecule:
        """Build a glycan

        Args:
            protein_residue (Residue, optional): protein residue to be glycosylated. If None, just the glycan is returned.
            protein_link_patch (str, optional): patch to connect protein_residue to glycan. Necessary if ab initio glycosylating a protein residue.
            glycan (Molecule, optional): glycan Molecule instance to be modified. If none, a glycan is built ab initio.
            template_glycan (_type_): target glycan to be created. Can be by name or by linkage topology.

        Side effects:
            if protein_residue is provided, modifications will be done in place. Provide a deep copy of protein_residue if this is undesired
            if glycan is provided, unknown glycan patches will be removed
        Returns:
            Molecule: the built glycan
        """
        # Step 1: get the data for template_glycan
        # Case 1: template provided by name
        if type(template_glycan) == str:
            template_topology: GlycanTopology = self.glycan_topologies[template_glycan]
        # Case 2: template provided as linkage topology
        elif type(template_glycan) == GlycanTopology:
            template_topology: GlycanTopology = template_glycan

        template_patches_to_i = {
            tuple(path.patches): i for i, path in enumerate(template_topology.paths)
        }

        # Step 2: prepare root
        # Case 1: root is an amino acid
        if protein_residue is not None:
            root = protein_residue
            ...
        # Case 2: no root provided -> DUMMY root residue
        elif protein_residue is None:
            root = self.builder.build_from_DUMMY(...)
            ...

        # Step 3: modify glycan
        # Case 1: glycan was provided
        # TODO: compare root of template and glycan -> should be same
        if glycan is not None:
            atom_group = glycan.atom_group
            glycan_graph: ResidueGraph = glycan.residue_graph
            glycan_topology = ...
            glycan_patches_to_i = ...
        # Case 2: no glycan provided -> ab initio building
        elif glycan is None:
            atom_group = prody.AtomGroup()
            glycan_graph = ResidueGraph()

        # template_edges = nx.edge_bfs(
        #     template_glycan_graph,
        #     source=template_glycan_graph.graph["root_residue_index"],
        # )
        # glycan_edges = nx.edge_bfs(
        #     glycan_graph, source=glycan_graph.graph["root_residue_index"]
        # )
        # TODO: guaranteed bug - if bfs traverses a short branch first in one graph, and a long branch first in the other graph,
        # then template residues go out of sync with glycan residues
        # solution: use .linkage_paths property, and iterate over template paths, then check if it exists in glycan paths

        for path in sorted(template_topology.paths, key=len):
            # TODO: handle root residue case where there are no patches
            *prev_patches, patch = path.patches
            prev_template_i, template_i = template_patches_to_i.get(
                tuple(prev_patches)
            ), template_patches_to_i.get(tuple(path.patches))
            prev_res_i, res_i = (
                glycan_patches_to_i.get(tuple(prev_patches)),
                glycan_patches_to_i.get(tuple(path.patches)),
            )

            # Case 1: residue exists and may need fixing
            if res_i is not None:
                residue = glycan_graph[res_i]["residue"]
                resnum, chain, segname = (
                    residue.getResnum(),
                    residue.getChid(),
                    residue.getSegname(),
                )
                # build new complete residue
                newres, missing_atoms, bonds = self.builder.find_missing_atoms(
                    residue, resnum, chain, segname
                )
                # set the coordinates for missing atoms using template ICs
                ics = self.builder.Topology[path.residue_name]["IC"]
                self.builder.build_missing_atom_coord(newres, missing_atoms, ics)
                # apply_patch ?

            # Case 2: second residue needs to be created ab initio
            else:
                # previous residue info
                residue1 = glycan_graph[prev_res_i]["residue"]
                chain, segname = (
                    residue1.getChid(),
                    residue1.getSegname(),
                )
                # new res name
                resname = path.residue_name
                # new res num
                resnum = max(atom_group.getResnums()) + 1
                res_i = max(atom_group.getResindices()) + 1
                # make new residue ab initio
                newres, new_dele_atoms, bonds = self.builder.build_from_patch(
                    residue1, resnum, resname, chain, segname, patch
                )
                # need to update graph...

            # delete atoms
            # tinker with bonds?

        return glycan

    def load_glycan_topologies(self, file_path):
        try:
            glycan_topologies = glycan_topology.read_glycan_topology(file_path)
            self.glycan_topologies = glycan_topologies
            return
        except UnicodeDecodeError:
            pass

        try:
            glycan_topologies = glycan_topology.import_connectivity_topology(file_path)
            self.glycan_topologies = glycan_topologies
            return
        except sqlite3.DatabaseError:
            pass

        raise ValueError("File format of file not recognised.")
