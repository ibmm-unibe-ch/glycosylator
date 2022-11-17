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
        glycan_paths = set(glycan.glycan_topology.paths)

        for glycan_name, template_topology in self.glycan_topologies.items():
            template_paths = set(template_topology.paths)
            if glycan_paths == template_paths:
                return glycan_name

        return "UnknownGlycan"

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

    def set_bonds(self, glycan: Molecule) -> None:
        ...

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
        # Case 2: template provided as linkage topology
        elif type(template_glycan) == GlycanTopology:
            template_topology: GlycanTopology = template_glycan
        else:
            raise TypeError(
                "Argument template_glycan should either be a str name for glycan or a GlycanTopology"
            )

        template_patches_to_i = {
            tuple(path.patches): i for i, path in enumerate(template_topology.paths)
        }

        # handle presense/absence of amino acid root
        # four cases:
        # case both protein_residue and glycan provided
        if protein_residue is not None and glycan is not None:
            # is protein_residue already a part of glycan?
            if glycan.residue_graph.graph["root_residue"] == protein_residue:
                glycan = glycan.atom_group.select(
                    f"not resname {protein_residue.getResname()}"
                ).toAtomGroup()
                root_atom = next(glycan.select("name C1").iterAtoms())
                glycan = Molecule(glycan, root_atom)
            # the case of glycosylating an amino acid with a provided glycan against a template is complicated
            # because two separate things need to be done:
            #   glycan needs to be patched onto the residue, preserving known ICs
            #   glycan needs to be fixed against template
            else:
                raise NotImplementedError(
                    "Glycosylating a residue with a template whilst preserving known ICs is not yet supported. Either glycosylate ab initio using just the template, or glycosylate the glycan without patching onto root"
                )
                # glycosylate residue ab intio, using the 'broken' glycan as template
                # glycan = self.glycosylate(protein_residue, protein_link_patch, glycan=None, template_glycan=glycan.)
                # this glycan is now ready to be checked against true template

        # case only protein_residue provided
        elif protein_residue is not None:
            # generate ab intio first residue
            resid = max(protein_residue.getAtomGroup().getResnums()) + 1
            resname = template_topology.paths[0].residue_name
            chain = protein_residue.getChid()
            segname = protein_residue.getSegname()
            glycan, dele_atoms, _ = self.builder.build_from_patch(
                protein_residue, resid, resname, chain, segname, protein_link_patch
            )
            new_glycan = self.builder.delete_atoms(glycan, dele_atoms)
            new_protein_residue = self.builder.delete_atoms(
                protein_residue.copy(), dele_atoms
            )

            root_atom = next(new_glycan.select("name C1").iterAtoms())
            glycan = Molecule(new_glycan, root_atom)

        # case only glycan provided
        elif glycan is not None:
            root_residue = glycan.root_residue
            resnum, chain, segname = (
                root_residue.getResnum(),
                root_residue.getChid(),
                root_residue.getSegname(),
            )
            # build new complete residue
            newres, missing_atoms, new_bonds = self.builder.find_missing_atoms(
                root_residue, resnum, chain, segname
            )
            # set the coordinates for missing atoms using template ICs
            ics = self.builder.Topology.topology[root_residue.getResname()]["IC"]
            self.builder.build_missing_atom_coord(newres, missing_atoms, ics)
            new_glycan = newres

        # case neither protein_residue nor glycan provided
        elif protein_residue is None and glycan is None:
            # generate ab initio first residue
            resname = template_topology.paths[0].residue_name
            # TODO: DUMMY_MAN might not be loaded into builder.Topology...
            glycan, dele_atoms, _ = self.builder.build_from_DUMMY(
                resid=1,
                resname=resname,
                chain="X",
                segname="G1",
                dummy_patch="DUMMY_MAN",
            )
            new_glycan = self.builder.delete_atoms(glycan, dele_atoms)

            root_atom = next(new_glycan.select("name C1").iterAtoms())
            glycan = Molecule(new_glycan, root_atom)

        glycan_patches_to_i = {
            tuple(path.patches): i
            for i, path in enumerate(glycan.residue_graph.to_glycan_topology())
        }
        glycan_graph = glycan.residue_graph
        atom_group = glycan.atom_group

        del_atoms = []
        # keep track of bonds -> mustn't use AtomGroup.inferBonds()
        # because accidental clashes will lead to spurious bonds
        bonds = []
        for i, path in enumerate(template_topology):
            if len(path) == 0:
                continue
            *prev_patches, patch = path.patches
            prev_template_i, template_i = (
                template_patches_to_i.get(tuple(prev_patches)),
                template_patches_to_i.get(tuple(path.patches)),
            )
            prev_res_i, res_i = (
                glycan_patches_to_i.get(tuple(prev_patches)),
                glycan_patches_to_i.get(tuple(path.patches)),
            )

            # Case 1: residue exists and may need fixing
            if res_i is not None:
                residue = glycan_graph.nodes[res_i]["residue"]
                resnum, chain, segname = (
                    residue.getResnum(),
                    residue.getChid(),
                    residue.getSegname(),
                )
                # build new complete residue
                newres, missing_atoms, new_bonds = self.builder.find_missing_atoms(
                    residue, resnum, chain, segname
                )
                bonds.extend(new_bonds)
                # set the coordinates for missing atoms using template ICs
                ics = self.builder.Topology.topology[path.residue_name]["IC"]
                self.builder.build_missing_atom_coord(newres, missing_atoms, ics)

                new_glycan += newres

                # newres.setBonds(new_bonds)

            # Case 2: second residue needs to be created ab initio
            else:
                # previous residue info
                prev_residue = glycan_graph.nodes[prev_res_i]["residue"]
                chain, segname = (
                    prev_residue.getChid(),
                    prev_residue.getSegname(),
                )
                # new res name
                resname = path.residue_name
                # new res num
                resnum = max(new_glycan.getResnums()) + 1
                res_i = max(new_glycan.getResindices()) + 1
                # make new residue ab initio
                newres, new_dele_atoms, new_bonds = self.builder.build_from_patch(
                    prev_residue, resnum, resname, chain, segname, patch
                )

                new_glycan += newres

                bonds.extend(new_bonds)
                # newres.setBonds(new_bonds)
                del_atoms.extend(new_dele_atoms)
                # atom_group += newres

                # convert newres to type Residue from AtomGroup
                newres = newres[segname, chain, resnum]

                # update datastructures for next iterations
                glycan_graph.add_edge(prev_res_i, res_i)
                glycan_graph.add_node(res_i, residue=newres)
                glycan_patches_to_i[path.patches] = res_i

        final_glycan = self.builder.delete_atoms(new_glycan, del_atoms)

        if protein_residue is not None:
            final_glycan += protein_residue

        if inplace:
            pass

        # tinker with bonds?

        # DONT DO THIS, just placeholder
        # if the ICs cause the residues to clash with themselves, or the protein
        # then spurious bonds will be inferred
        # do it properly by bookkeeping bonds!
        # final_glycan.inferBonds()
        return final_glycan

        return (
            Molecule(final_glycan, next(final_glycan.select("name C1").iterAtoms())),
            bonds,
        )

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
