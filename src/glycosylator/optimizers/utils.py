"""
This module contains additional utility functions for the optimizers that are Glycosylator specific.
"""


def make_scaffold_graph(
    scaffold,
    only_clashing_glycans: bool = True,
    only_clashes_between_glycans: bool = False,
    slice: int = 0,
    include_root: bool = False,
) -> tuple["ResidueGraph", list]:
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
    slice: int
        The number of residue connections to slice from each glycan. This will sub-sample the rotatable edges of the glycans. At 0 all rotatable edges are included.
    include_root: bool
        If True, the root connection (to the scaffold) will always be included in the rotatable edges.

    Returns
    -------
    ResidueGraph
        The graph of the glycosylated scaffold
    list
        The rotatable edges of the graph
    """
    G = scaffold.make_residue_graph(True)
    _rotatable_edges = []

    if only_clashing_glycans:
        glycan_residues = {
            glycan: [res for res in glycan.get_residues()]
            for glycan in scaffold.glycans.values()
        }
        _residues_flat = set(
            [res for glycan in glycan_residues.values() for res in glycan]
        )
        _added_glycans = set()
        clashes = scaffold.find_clashes()
        if only_clashes_between_glycans:
            _check = (
                lambda a, b: a.parent in _residues_flat and b.parent in _residues_flat
            )
        else:
            _check = (
                lambda a, b: a.parent in _residues_flat or b.parent in _residues_flat
            )

            for clash in clashes:
                if _check(*clash):
                    for i in clash:
                        glycan = next(
                            (
                                glycan
                                for glycan, residues in glycan_residues.items()
                                if i.parent in residues
                            ),
                            None,
                        )

                        if glycan is None or glycan in _added_glycans:
                            continue

                        G.add_edge(glycan.root_residue, glycan.root_atom)
                        if include_root:
                            scaffold_root = next(
                                atom
                                for atom, _glycan in scaffold.glycans.items()
                                if glycan == _glycan
                            )
                            G.add_edge(scaffold_root, glycan.root_atom)
                            G.add_edge(scaffold_root.parent, scaffold_root)
                            _rotatable_edges.append((scaffold_root, glycan.root_atom))

                        if slice > 0:
                            _rotatable_edges.extend(
                                _slice_residue_connections(glycan, slice)
                            )
                        else:
                            _rotatable_edges.extend(
                                glycan.get_residue_connections(direct_by="root")
                            )

                        _added_glycans.add(glycan)

    else:
        for glycan in scaffold.glycans.values():
            if slice > 0:
                _rotatable_edges.extend(_slice_residue_connections(glycan, slice))
            else:
                _rotatable_edges.extend(
                    glycan.get_residue_connections(direct_by="root")
                )

    return G, _rotatable_edges


def _slice_residue_connections(glycan, slice: int):
    """
    Slice a number of residue connections from the glycan.

    Parameters
    ----------
    glycan: Glycan
        The glycan to slice
    slice: int
        The number of residue connections to slice from the glycan

    Returns
    -------
    list
        The sliced residue connections
    """
    bonds = glycan.get_residue_connections(direct_by="root")
    if slice > 1:
        bonds = [bonds[i] for i in range(1, len(bonds), len(bonds) // slice)]
    return bonds


if __name__ == "__main__":
    import glycosylator as gl

    scaf = gl.Scaffold.load(
        "/Users/noahhk/GIT/glycosylator/src/glycosylator/__figure_makery/scaf_glycosylated.optimized.5.pkl"
    )

    graph, rotatable_edges = make_scaffold_graph(
        scaf,
        only_clashing_glycans=True,
        only_clashes_between_glycans=True,
        slice=5,
        include_root=True,
    )

    rotatron = gl.optimizers.DistanceRotatron(
        graph, rotatable_edges, pushback=1.5, n_processes=8
    )

    # scaf_opt = gl.optimizers.optimize(scaf.copy(), rotatron)
    # scaf_opt.to_pdb("scaf_opt_new_slice_graph.pdb")

    split = gl.optimizers.split_environment(rotatron, 4)
    scaf_opt = gl.optimizers.parlallel_optimize(scaf, split)
    scaf_opt.to_pdb("scaf_opt_new_slice_parallel.pdb")
