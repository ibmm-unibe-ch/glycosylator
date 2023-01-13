import copy
import re
import sqlite3

import networkx as nx
import prody

from glycosylator.file_parsers import glycan_topology
from glycosylator.file_parsers.glycan_topology import GlycanTopology
from glycosylator.glycan_representations import ResidueGraph
from glycosylator.molecule import Molecule
from glycosylator.molecule_builder import HighLevelBuilder


class Glycosylator:
    def __init__(
        self,
        charmm_topology_path: str,
        charm_parameter_path: str,
    ):
        self.charmm_topology_path = charmm_topology_path
        self.charmm_parameter_path = charm_parameter_path
        self.builder = HighLevelBuilder(charmm_topology_path, charm_parameter_path)
        self.glycan_topologies: dict[str, GlycanTopology] = dict()

    def load_glycoprotein_from_PDB(self, protein_file_path: str):
        glycoprotein_atom_group = prody.parsePDB(protein_file_path)
        self.load_glycoprotein_from_AtomGroup(glycoprotein_atom_group)

    def load_glycoprotein_from_AtomGroup(
        self, glycoprotein_atom_group: prody.AtomGroup
    ):
        self.glycoprotein = glycoprotein_atom_group
        self.residues = [residue for residue in self.glycoprotein.iterResidues()]
        # self.linking_residues = self.find_linking_residues()
        # if self.glycoprotein.getBonds() is None:
        #     self._infer_glycan_bonds()

    def write_glycoprotein_to_PDB(self, file_path: str):
        prody.writePDB(file_path, self.glycoprotein)

    def find_linking_residues(
        self, sequon_pattern: str = "(N[^P][S|T])"
    ) -> list[prody.Residue]:
        sequon_pattern = re.compile(sequon_pattern)
        linking_residues = []
        for chain in self.glycoprotein.iterChains():
            residues = [res for res in chain.iterResidues()]
            matches = sequon_pattern.finditer(chain.getSequence())
            linking_residues.extend([residues[match.start()] for match in matches])
        return linking_residues

    def find_existing_glycans(self, freeze_bonds=True) -> list[Molecule]:
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
            if freeze_bonds:
                glycan.atom_graph.freeze_bonds("all")
            glycans.append(glycan)

        self.glycans = glycans
        return glycans

    def identify_glycan(self, glycan: Molecule):
        glycan_paths = set(glycan.glycan_topology.paths)

        for glycan_name, template_topology in self.glycan_topologies.items():
            template_paths = set(template_topology.paths)
            if glycan_paths == template_paths:
                return glycan_name

        return "UnknownGlycan"

    def find_glycan(self, amino_acid_root: prody.Residue | str):
        try:
            glycans = self.glycans
        except AttributeError:
            glycans = self.find_existing_glycans()

        for glycan in glycans:
            if glycan.root_residue == amino_acid_root:
                return glycan

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

    def assign_patches(self, glycan: Molecule):
        residue_graph = glycan.residue_graph
        for res1, res2, data in residue_graph.edges(data=True):
            atom_names = {atom.getName() for atom in data["atoms"]}
            patch = self.find_patch(atom_names)
            data["patch"] = patch

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
        glycan: Molecule = None,
        template_glycan: str | GlycanTopology = None,
        protein_residue: prody.Residue = None,
        protein_link_patch: str = None,
        inplace: bool = False,
    ) -> Molecule:
        """Build a glycan

        Args:
            protein_residue (Residue, optional): protein residue to be glycosylated. If None, just the glycan is returned.
            protein_link_patch (str, optional): patch to connect protein_residue to glycan. Necessary if ab initio glycosylating a protein residue.
            glycan (Molecule, optional): glycan Molecule instance to be modified. If none, a glycan is built ab initio.
            template_glycan (str | GlycanTopology): target glycan to be created. Can be by name or by linkage topology.

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
        # Case 2: template provided as GlycanTopology instance
        elif type(template_glycan) == GlycanTopology:
            template_topology: GlycanTopology = template_glycan
        else:
            raise TypeError(
                "Argument template_glycan should either be a valid key for self.glycan_topologies or a GlycanTopology instance"
            )

        if glycan is not None:
            old_glycan = glycan
            old_glycan_patches_to_i = {
                path.patches: i
                for i, path in enumerate(old_glycan.residue_graph.to_glycan_topology())
            }
        else:
            old_glycan_patches_to_i = dict()

        new_glycan, new_root_atom = self._initialise_new_glycan(
            glycan, template_topology, protein_residue, protein_link_patch
        )
        new_glycan_graph = ResidueGraph.from_AtomGroup(new_glycan, new_root_atom)
        new_glycan_patches_to_i = {
            path.patches: i
            for i, path in enumerate(new_glycan_graph.to_glycan_topology())
        }

        # TODO WIP:
        # finish in-place feature
        # probably will have and need to fix bugs related to Residue hierarchy views into new_glycan, since new_glycan keeps changing
        # the view is effectively immutable, even if AtomGroup has a mutable "+=" interface
        for path in template_topology:
            # skip empty path as first residue has been initialised already
            if len(path) == 0:
                continue

            *prev_patches, patch = path.patches
            # *prev_patches is a list, but need a tuple
            prev_res_i = new_glycan_patches_to_i[tuple(prev_patches)]
            # previous_residue = [res for res in new_glycan.iterResidues()][prev_res_i]
            old_res_i = old_glycan_patches_to_i.get(path.patches)

            # Case 1: old_residue exists and needs repair of any missing atoms
            if old_res_i is not None:
                old_residue = old_glycan.residue_graph.nodes[old_res_i]["residue"]
                new_residue = self.builder.residue_repair(old_residue)
                new_glycan += new_residue.getAtomGroup()

                new_root_atom = (
                    new_glycan.select(new_root_atom.getSelstr()).iterAtoms().__next__()
                )
                new_glycan_graph = ResidueGraph.from_AtomGroup(
                    new_glycan, new_root_atom
                )
                new_glycan_graph.identify_patches(self)
                prev_res_i = new_glycan_graph.path_to_index(prev_patches)
                new_res_i = new_glycan_graph.path_to_index(path.patches)

                previous_residue = new_glycan_graph.nodes[prev_res_i]["residue"]
                new_residue = new_glycan_graph.nodes[new_res_i]["residue"]

                new_glycan = self.builder.apply_patch(
                    patch, previous_residue, new_residue
                )

            # Case 2: residue needs to be patched onto an existing residue
            else:
                previous_residue = new_glycan_graph.nodes[prev_res_i]["residue"]
                new_glycan, new_residue = self.builder.add_residue_from_patch(
                    previous_residue, path.residue_name, patch
                )
                # update datastructures for next iterations
                new_root_atom = (
                    new_glycan.select(new_root_atom.getSelstr()).iterAtoms().__next__()
                )
                new_glycan_graph = ResidueGraph.from_AtomGroup(
                    new_glycan, new_root_atom
                )
                new_glycan_patches_to_i = {
                    path.patches: i
                    for i, path in enumerate(new_glycan_graph.to_glycan_topology())
                }

                # res_i = new_residue.getResindex()
                # # new_glycan_graph.add_node(res_i, residue=new_residue)
                # new_glycan_graph.add_edge(prev_res_i, res_i)
                # new_glycan_patches_to_i[path.patches] = res_i

        # if protein_residue is not None:
        #     # the coordinates and linkage already valid from initialisation at beginning
        #     final_glycan += new_protein_residue

        if inplace:
            ...

        final_glycan = Molecule(new_glycan, new_root_atom)
        self.assign_patches(final_glycan)
        return final_glycan

    def _initialise_new_glycan(
        self,
        glycan: Molecule | None,
        template_topology,
        protein_residue: prody.Residue | None,
        protein_link_patch,
    ):
        # handle presense/absence of amino acid root
        # four cases:
        # case both protein_residue and glycan provided
        if protein_residue is not None and glycan is not None:
            # is protein_residue already a part of glycan? e.g. from running self.find_existing_glycans(...)
            if glycan.root_residue == protein_residue:
                new_glycan = protein_residue.toAtomGroup()
                new_root_atom = new_glycan.select(
                    f"name {glycan.root_atom.getName()}"
                ).__next__()
                res = glycan.residue_graph.neighbours(
                    protein_residue.getResindex()
                ).__next__()
                new_glycan += self.builder.residue_repair(res).getAtomGroup()
            # the case of glycosylating an amino acid with a provided glycan against a template is complicated
            # because two separate things need to be done:
            #   glycan needs to be fixed against template <- no problem
            #   glycan needs to be patched onto the residue, preserving experimental ICs <- currently, can only patch using theoretical ICs
            #   - essentially, need the functionality to rotate and translate glycan in space to be close enough to protein_residue
            else:
                raise NotImplementedError(
                    "Glycosylating a residue with a template whilst preserving experimental internal coordinates is not yet supported. Either glycosylate ab initio using just the template, or glycosylate the glycan without patching onto root"
                )
                # glycosylate residue ab intio, using the 'broken' glycan as template
                # glycan = self.glycosylate(protein_residue, protein_link_patch, glycan=None, template_glycan=glycan.)
                # this glycan is now ready to be checked against true template

        # case only protein_residue provided
        elif protein_residue is not None:
            # make a new AtomGroup with only the protein_residue
            new_glycan = protein_residue.toAtomGroup()
            # have protein_residue point towards the new AtomGroup
            new_protein_residue = new_glycan.iterResidues().__next__()
            resname = template_topology.paths[0].residue_name
            # new_glycan now contains an amino acid and a single glycan
            new_glycan, _ = self.builder.add_residue_from_patch(
                new_protein_residue, resname, protein_link_patch
            )
            new_root_atom = (
                new_glycan.select(new_protein_residue.getSelstr())
                .select("name C1")
                .iterAtoms()
                .__next__()
            )

        # case only glycan provided
        elif glycan is not None:
            root_residue = glycan.root_residue
            root_atom_name = glycan.root_atom.getName()
            new_glycan = self.builder.residue_repair(root_residue).getAtomGroup()
            new_root_atom = (
                new_glycan.select(f"name {root_atom_name}").iterAtoms().__next__()
            )

        # case neither protein_residue nor glycan provided
        elif protein_residue is None and glycan is None:
            # generate ab initio first residue
            resname = template_topology.paths[0].residue_name
            # TODO: DUMMY_MAN might not be loaded into builder.Topology...
            new_glycan = self.builder.residue_ab_initio(
                resname=resname,
                resid=1,
                chain="X",
                segname="G1",
                dummy_patch="DUMMY_MAN",
            ).getAtomGroup()
            new_root_atom = new_glycan.select("name C1").iterAtoms().__next__()

        return new_glycan, new_root_atom

    def load_glycan_topologies(self, file_path: str, append: bool = True):
        try:
            new_glycan_topologies = glycan_topology.read_glycan_topology(file_path)
            if append:
                self.glycan_topologies |= new_glycan_topologies
            else:
                self.glycan_topologies = new_glycan_topologies
            return
        except UnicodeDecodeError:
            pass

        try:
            new_glycan_topologies = glycan_topology.import_connectivity_topology(
                file_path
            )
            if append:
                self.glycan_topologies |= new_glycan_topologies
            else:
                self.glycan_topologies = new_glycan_topologies
            return
        except sqlite3.DatabaseError:
            pass

        raise ValueError("File format of file not recognised.")


def _atom_id_to_index(atom_id, atom_group):
    selection_keywords = ["segment", "chain", "resid", "icode", "name"]
    sel = (
        f"{keyword} {value}"
        for keyword, value in zip(selection_keywords, atom_id.split(","))
        if value != ""
    )
    sel = " and ".join(sel)
    sel = f"({sel})"

    atom_index = atom_group.select(sel).iterAtoms().__next__().getIndex()
    return atom_index
