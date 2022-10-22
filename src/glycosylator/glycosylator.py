import re

import networkx as nx
import prody
from prody import AtomGroup, Residue, Selection
from typing_extensions import Self

from .molecule import Molecule
from .molecule_builder import MoleculeBuilder


class Glycosylator:
    def __init__(
        self,
        glycoprotein: AtomGroup,
        charmm_topology_path: str,
        charm_parameter_path: str,
    ):
        self.glycoprotein = glycoprotein
        self.residues = [residue for residue in self.glycoprotein.iterResidues()]
        self.linking_residues = self.find_linking_residues()
        if glycoprotein.getBonds() is None:
            self._infer_glycan_bonds()

        self.charmm_topology_path = charmm_topology_path
        self.charmm_parameter_path = charm_parameter_path
        self.builder = MoleculeBuilder(charmm_topology_path, charm_parameter_path)

    @classmethod
    def fromPDB(
        cls,
        protein_file_path: str,
        charmm_topology_path: str,
        charmm_parameter_path: str,
    ):
        protein_atom_group = prody.parsePDB(protein_file_path)
        return cls(protein_atom_group, charmm_topology_path, charmm_parameter_path)

    def find_linking_residues(self) -> list[Residue]:
        sequon_pattern = re.compile("(N[^P][S|T])")
        linking_residues = []
        for chain in self.glycoprotein.iterChains():
            residues = [res for res in chain.iterResidues()]
            matches = sequon_pattern.finditer(chain.getSequence())
            linking_residues.extend([residues[match.start()] for match in matches])
        return linking_residues

    def find_existing_glycans(self):
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
            root_atom = next(selection.select("protein and name CA").iterAtoms())
            glycans.append(Molecule(selection, root_atom=root_atom))
            # glycans.append(selection)

        return glycans

    def check_protein_glycan_linkage(self, glycan: Molecule, patches: list[str]):
        residue_graph = glycan.residue_graph
        root_residue_index = residue_graph.graph["root_residue_index"]
        try:
            next_residue_index = next(residue_graph.successors(root_residue_index))
        # if no successors, return False
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
