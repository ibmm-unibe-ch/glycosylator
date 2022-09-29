from queue import Queue

import networkx as nx
import prody

from .bond_graph import build_bond_graph
from .residue_graph import build_residue_graph


class TooManyChains(Exception):
    """Raise when Molecule class has an AtomGroup with more than one chain/segment"""


class Molecule:
    def __init__(
        self,
        atom_group: prody.AtomGroup,
        root_atom: int,
    ):
        _validate_atom_group(atom_group)
        self.atom_group: prody.AtomGroup = atom_group

        self.root_atom: prody.Atom = self.atom_group[root_atom - 1]
        self.root_residue: prody.Residue = self.atom_group[
            self.root_atom.getChid(), self.root_atom.getResnum()
        ]

        self.bond_graph: nx.graph = self.build_bond_graph()
        self.residue_graph: nx.digraph = self.build_residue_graph()

        # self.guess_angles()
        # self.guess_dihedrals()
        # self.torsionals

    @classmethod
    def read_PDB(cls, pdb: str, root_atom: int = 1, **kwargs):
        atom_group = prody.parsePDB(pdb, **kwargs)
        return Molecule(atom_group, root_atom)

    def write_PDB(self, filename: str, selection: str = "all"):
        prody.writePDB(filename, self.atom_group.select(selection))

    def build_bond_graph(self):
        bonds = self.atom_group.inferBonds()
        bond_graph = nx.Graph()
        bond_graph.add_edges_from([bond.getIndices() for bond in bonds])
        return bond_graph  # bond_graph.edges() ?

    def build_residue_graph(self):
        residue_graph = nx.DiGraph()
        residue_graph.add_node(self.root_residue.getResindex())

        residues = [residue for residue in self.atom_group.iterResidues()]
        # residue index for each atom
        res_indices = self.atom_group.getResindices()

        queue = [self.root_residue.getResindex()]
        while queue:
            old_resi = queue.pop(0)
            residue = residues[old_resi]
            for atom in residue:
                # list of neighboring residues not yet in graph
                neighbor_residues = [
                    new_resi
                    for neighbor in self.bond_graph.neighbors(atom.getIndex())
                    if (new_resi := res_indices[neighbor]) not in residue_graph.nodes
                ]
                edges = [(old_resi, new_resi) for new_resi in neighbor_residues]
                residue_graph.add_edges_from(edges)
                queue.extend(neighbor_residues)

        return residue_graph

    def guess_angles(self):
        """Searches for all angles in a molecule based on the connectivity"""
        self.angles = []
        for node in self.bond_graph.nodes():
            self.angles.extend(_find_paths(self.connectivity, node, 2))

    def guess_dihedrals(self):
        """Searches for all dihedrals in a molecule based on the connectivity"""
        self.dihedrals = []
        for node in self.bond_graph.nodes():
            self.dihedrals.extend(_find_paths(self.connectivity, node, 3))

    def guess_torsionals(self, hydrogens=True):
        """Builds a list with all the torsional angles that can rotate
        Parameters:
            hydrogens: include torsional involving terminal hydrogens
        Initializes:
            torsionals: a list of serial number of atom defining a torsional angle (quadruplet)
        """
        torsionals = []
        if not hydrogens:
            elements = nx.get_node_attributes(self.bond_graph, "element")
            if not elements:
                print(
                    "Elements have not been defined (use assign_atom_type). Hydrogens cannot be excluded."
                )
                hydrogens = True

        for dihe in self.dihedrals:
            # check that quadruplet is not in cycle
            if dihe[1] in self.cycle_id and dihe[2] in self.cycle_id:
                continue

            d_dihe = []
            for a in dihe[1:-1]:
                if a in self.cycle_id:
                    a = self.cycle_id[a]
                d_dihe.append(a)

            # check direction of dihedral
            if self.directed_connectivity.has_edge(d_dihe[0], d_dihe[1]):
                pass
            elif self.directed_connectivity.has_edge(d_dihe[1], d_dihe[0]):
                dihe.reverse()
            else:
                continue

            # check if hydrogen
            if not hydrogens:
                if elements[dihe[0]] == "H" or elements[dihe[-1]] == "H":
                    continue
            # check if already in torsionals list
            exists = False
            if self.torsionals:
                for t in self.torsionals:
                    if dihe[1] == t[1] and dihe[2] == t[2]:
                        exists = True
                        break

            if exists:
                continue

            self.torsionals.append(dihe)


def _validate_atom_group(atom_group: prody.AtomGroup):
    """validates that the AtomGroup consists of a single physical molecule"""
    if atom_group.numChains() != 1 or atom_group.numSegments() != 1:
        raise TooManyChains()


def _find_paths(G, node, length, excludeSet=None):
    """Finds all paths of a given length
    Parameters:
        G: graph (netwrokx)
        node: starting node
        length: length of path
        excludedSet: set
    Returns:
        paths: list of all paths of a length starting from node
    """
    if excludeSet == None:
        excludeSet = {node}
    else:
        excludeSet.add(node)

    if length == 0:
        return [[node]]
    paths = [
        [node] + path
        for neighbor in G.neighbors(node)
        if neighbor not in excludeSet
        for path in _find_paths(G, neighbor, length - 1, excludeSet)
    ]
    excludeSet.remove(node)
    return paths
