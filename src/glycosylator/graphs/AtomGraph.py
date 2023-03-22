import networkx as nx
import Bio.PDB as bio

import glycosylator.utils as utils
from glycosylator.graphs.BaseGraph import BaseGraph


class AtomGraph(BaseGraph):
    """
    A graph representation of atoms and bonds in a contiguous molecule.
    """

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        id=None,
        apply_standard_bonds: bool = True,
        infer_bonds: bool = False,
        max_bond_length: float = None,
        restrict_residues: bool = True,
    ):
        """
        Create an AtomGraph from a PDB of a single molecule

        Parameters
        ----------
        filename : str
            Path to the PDB file
        id : str
            The ID of the molecule. By default the filename is used.
        apply_standard_bonds : bool
            Whether to apply standard bonds from known molecule connectivity.
        infer_bonds : bool
            Whether to infer bonds from the distance between atoms.
        max_bond_length : float
            The maximum distance between atoms to infer a bond.
        restrict_residues : bool
            Whether to restrict to atoms of the same residue when inferring bonds.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        id = id if id else utils.filename_to_id(filename)

        # load PDB
        structure = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        structure = structure[0].child_list[0]

        bonds = []
        if apply_standard_bonds and infer_bonds:
            raise ValueError("Cannot apply standard bonds and infer bonds at the same time.")
        elif apply_standard_bonds:
            bonds.extend(utils.apply_standard_bonds(structure))
        elif infer_bonds:
            bonds.extend(utils.infer_bonds(structure, max_bond_length, restrict_residues))

        return cls(id, bonds)


if __name__ == '__main__':
    _man = "support/examples/man9.pdb"
    man = AtomGraph.from_pdb(_man, apply_standard_bonds=True, infer_bonds=False)

    import sys

    print(sys.getsizeof(man.residues[0]))
    # print(man.nodes)
    # print(man.edges)
    # # x = utils.infer_residue_connections(man.structure)
    # import matplotlib.pyplot as plt
    # import networkx as nx

    # # nx.draw(man)
    # colors = {
    #     "C": "gray",
    #     "O": "red",
    #     "N": "blue",
    #     "S": "yellow",
    #     "P": "orange",
    #     "H": "lightgray",
    #     "F": "green",
    # }
    # _colors = [colors[i.element] for i in man.nodes]
    # g = man
    # # y = [(i.get_parent(), j.get_parent()) for i, j in x]

    # # g = nx.Graph(y)
    # # colors = {"MAN": "green", "BMA": "cyan", "MNA": "orange", "GAL": "pink", "NAG": "gray"}
    # # _colors = [colors[i.resname] for i in g.nodes]
    # nx.draw(g, node_color=_colors)
    # plt.show()
