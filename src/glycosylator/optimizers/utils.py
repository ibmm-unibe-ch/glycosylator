"""
This module contains additional utility functions for the optimizers that are Glycosylator specific.
"""


def make_scaffold_graph(
    scaffold,
    only_clashing_glycans: bool = False,
    only_clashes_between_glycans: bool = False,
) -> tuple:
    """
    Make a ResidueGraph from a glycosylated scaffold in order to optimize the conformations of the attached glycans.

    Parameters
    ----------
    scaffold: Scaffold
        The glycosylated scaffold to optimize
    only_clashing_glycans: bool
        If True, only glycans that clash with the protein or each other will will contribute rotatable edges.
    only_clashes_between_glycans: bool
        If True, only clashes between glycans will be considered. This is useful for optimizing the conformations of glycans that are already in a good position.

    Returns
    -------
    ResidueGraph
        The graph of the glycosylated scaffold
    list
        The rotatable edges of the graph
    """
    G = scaffold.make_residue_graph(True)
    rotatable_edges = []
    if only_clashing_glycans:
        glycan_residues = {
            glycan: [res for res in glycan.get_residues()]
            for glycan in scaffold.glycans.values()
        }
        _residues_flat = set(
            [res for glycan in glycan_residues.values() for res in glycan]
        )
        _added_glycans = set()
        if only_clashes_between_glycans:
            for clash in scaffold.find_clashes():
                a, b = clash
                if a.parent in _residues_flat and b.parent in _residues_flat:
                    for i in clash:
                        glycan = next(
                            glycan
                            for glycan, residues in glycan_residues.items()
                            if i.parent in residues
                        )
                        if glycan in _added_glycans:
                            continue
                        _g = glycan.make_residue_graph(True)
                        _g.add_edge(glycan.root_atom, glycan.root_residue)
                        rotatable_edges.extend(
                            _g.find_rotatable_edges(glycan.root_atom)
                        )
                        _added_glycans.add(glycan)
        else:
            for clash in scaffold.find_clashes():
                a, b = clash
                for i in clash:
                    if i.parent in _residues_flat:
                        glycan = next(
                            glycan
                            for glycan, residues in glycan_residues.items()
                            if i.parent in residues
                        )
                        if glycan in _added_glycans:
                            continue
                        _g = glycan.make_residue_graph(True)
                        _g.add_edge(glycan.root_atom, glycan.root_residue)
                        rotatable_edges.extend(
                            _g.find_rotatable_edges(glycan.root_atom)
                        )
                        _added_glycans.add(glycan)
    else:
        for glycan in scaffold.glycans.values():
            _g = glycan.make_residue_graph(True)
            _g.add_edge(glycan.root_atom, glycan.root_residue)
            rotatable_edges.extend(_g.find_rotatable_edges(glycan.root_atom))

    return G, rotatable_edges
