"""
This module contains additional utility functions for the optimizers that are Glycosylator specific.
"""

import buildamol.graphs.residue_graph as residue_graph
import buildamol.structural as structural
import numpy as np


def make_scaffold_graph(
    scaffold,
    only_clashing_glycans: bool = True,
    slice: int = 0,
    include_root: bool = False,
    include_n_ancestors: int = 0,
) -> tuple["residue_graph.ResidueGraph", list]:
    """
    Make a ResidueGraph from a glycosylated scaffold in order to optimize the conformations of the attached glycans.

    Parameters
    ----------
    scaffold: Scaffold
        The glycosylated scaffold to optimize
    only_clashing_glycans: bool
        If True, only glycans that clash with the protein or each other will will contribute rotatable edges.
    slice: int
        The number of residue connections to slice from each glycan. This will sub-sample the rotatable edges of the glycans. At 0 all rotatable edges are included.
        Set to -1 to exclude all glycan-internal rotatable edges. This is useful for glycans that are already optimized. And only their placement relative to the scaffold
        needs to be optimized (i.e. you only want to optimize the root connections).
    include_root: bool
        If True, the root connection (to the scaffold) will always be included in the rotatable edges.
    include_n_ancestors: int
        If > 0, the n ancestor atoms (upstream of the scaffold anchor atom) will also be included in the rotatable edges.
        This will allow for a more versatile placement of glycans as it allows partial movement of the scaffold rather than just the glycan.
        This requires `include_root` to be True.

    Returns
    -------
    ResidueGraph
        The graph of the glycosylated scaffold
    list
        The rotatable edges of the graph
    """
    G = scaffold.make_residue_graph()

    if only_clashing_glycans:
        glycan_gen = (
            (root, glycan)
            for (root, glycan) in scaffold.get_glycans().items()
            if (glycan.count_clashes() > 0 or glycan.clashes_with_scaffold())
        )
    else:
        glycan_gen = scaffold.get_glycans().items()

    _atoms_to_fill_in = []
    _rotatable_edges = []

    _flat_residues = list(scaffold.get_residues())
    _res_coords = [res.coord for res in _flat_residues]
    dists = structural.cdist(_res_coords, _res_coords) < 8.0

    for root, glycan in glycan_gen:

        if glycan.count_residues() > 1:
            g = glycan.get_residue_graph()
            g.make_detailed(include_clashes=True)
            if slice >= 0:
                _rotatable_edges.extend(
                    glycan.get_residue_connections(direct_by="root")
                )
        else:
            g = glycan.get_residue_graph(detailed=True)

        root_residue = root.parent
        neighboring_residues = (
            _flat_residues[idx]
            for idx in np.where(dists[_flat_residues.index(root_residue)])[0]
        )
        neighboring_atoms = (
            atom for res in neighboring_residues for atom in res.child_list
        )

        _atoms_to_fill_in.extend(g.nodes)
        _atoms_to_fill_in.extend(neighboring_atoms)

        G.remove_nodes_from(g.residues)
        G.add_nodes_from(_atoms_to_fill_in)
        G.add_edges_from(g.edges)

        if slice > 0:
            _rotatable_edges = G.sample_edges(
                _rotatable_edges, len(scaffold.glycans), slice
            )

        if include_root:
            incoming = scaffold.get_glycan_root_connections(
                glycan=glycan,
                include_scaffold_ancestors=include_n_ancestors,
            )

            _rotatable_edges.extend(incoming)
            G.add_edges_from(incoming)

        _atoms_to_fill_in.clear()

    return G, _rotatable_edges


if __name__ == "__main__":
    import glycosylator as gl

    glycan = gl.glycan(
        "Gal(b1-4)GlcNAc(b1-3)[Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Gal(b1-4)Glc"
    )
    scaf = gl.Protein.from_pdb(
        "/Users/noahhk/GIT/glycosylator/__projects__/SOLF/solf.pdb"
    )
    sites = [
        scaf.get_residue(179, chain="A"),
        scaf.get_residue(141, chain="A"),
        scaf.get_residue(179, chain="C"),
        scaf.get_residue(179, chain="E"),
    ]
    scaf.apply_standard_bonds_for(*sites)
    scaf.add_hydrogens(*[i.get_atom("ND2") for i in sites])
    scaf.rename_atoms("HND22", "HD22").rename_atoms("HND21", "HD21")
    scaf.glycosylate(glycan, residues=sites)

    # scaf = gl.utils.load_pickle(
    #     "//Users/noahhk/GIT/glycosylator/__projects__/solf3/solF_plus_3G_rsr017_coot_30_man5.pdb_glycosylated_optimized.pkl"
    # )
    # scaf.apply_standard_bonds_for(*[i.parent for i in scaf.get_glycans().keys()])

    graph, rotatable_edges = make_scaffold_graph(
        scaf,
        only_clashing_glycans=True,
        slice=5,
        include_root=True,
        include_n_ancestors=2,
    )

    rotatron = gl.optimizers.DistanceRotatron(
        graph,
        rotatable_edges,
        pushback=1.5,
    )

    # scaf_opt = gl.optimizers.optimize(scaf.copy(), rotatron)
    # scaf_opt.to_pdb("scaf_opt_new_slice_graph.pdb")

    split = gl.optimizers.split_environment(rotatron, 4)
    split = [gl.optimizers.ScaffoldRotatron(i) for i in split]
    scaf_opt = gl.optimizers.parallel_optimize(scaf, split)
    scaf_opt.to_pdb("scaf_opt_new_slice_parallel.pdb")
