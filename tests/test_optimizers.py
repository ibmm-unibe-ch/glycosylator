"""
Tests for the optimizers
"""

import base
import numpy as np
import glycosylator as gl
import glycosylator.optimizers as opt


# def test_rotatron_constraints():
#     mol = gl.Molecule.from_pdb(base.MANNOSE9)
#     mol.infer_bonds(restrict_residues=False)

#     connections = mol.get_residue_connections()

#     env = opt.MultiBondRotatron(mol.make_residue_graph(detailed=True), connections)

#     initial_coords = env._coords.copy()
#     nodes = list(env.graph.nodes)

#     initial_distances = []
#     for node in nodes:
#         neighbors = env.graph.get_neighbors(node)
#         dists = np.array(
#             [np.linalg.norm(node.coord - neighbor.coord) for neighbor in neighbors]
#         )
#         initial_distances.append(dists)

#     for run in range(3):
#         v = gl.utils.visual.MoleculeViewer3D(mol.make_residue_graph(detailed=True))

#         for i in range(200):
#             action = env.action_space.sample()
#             descendant_mask = env._descendant_masks[action[0]]
#             descendant_coords = env._coords[descendant_mask].copy()
#             ancestor_coords = env._coords[~descendant_mask].copy()

#             env.step(action)

#             assert not np.allclose(descendant_coords, env._coords[descendant_mask])
#             assert np.allclose(ancestor_coords, env._coords[~descendant_mask])

#         env.apply_to_graph()
#         new_distances = []
#         for node in nodes:
#             dists = []
#             neighbors = env.graph.get_neighbors(node)
#             for neighbor in neighbors:
#                 dists.append(np.linalg.norm(node.coord - neighbor.coord))
#             new_distances.append(np.asarray(dists))

#         v.draw_edges(env.graph.edges, color="cyan", linewidth=2)
#         v.show()

#         for i in range(len(nodes)):
#             assert np.all(np.abs(initial_distances[i] - new_distances[i]) < 1e-3)

#         env.reset()
#         assert np.allclose(initial_coords, env.state)


def test_scaffold_graph_optimization():
    from pathlib import Path

    scaffold = gl.Protein.from_pdb(Path(base.HOME) / "support/examples/small.pdb")
    glycan = gl.read_iupac("test", "Neu5Ac(a2-3)Gal(b2-4)GlcNAc")

    sites = scaffold.find_n_linked_sites()
    sites = list(sites.values())[0]
    scaffold.apply_standard_bonds_for(sites)
    scaffold.add_hydrogens([atom for res in sites for atom in res.atoms])
    assert len(sites) == 1

    scaffold.attach(glycan, residues=sites, chain="new")
    assert len(scaffold.glycans) == 1
    v = scaffold.plotly(residue_graph=True)

    graph, edges = gl.make_scaffold_graph(
        scaffold,
        only_clashing_glycans=False,
        include_root=True,
        include_n_ancestors=2,
    )

    env = gl.optimizers.DistanceRotatron(graph, edges)
    env = gl.optimizers.ScaffoldRotatron(env)

    sol, val = gl.swarm_optimize(env)
    out = gl.apply_rotatron_solution(sol, env, scaffold.copy())

    v += out.move([10, 10, 10]).plotly(residue_graph=True, line_color="red")
    v.show()


def test_optimize_glycans_never_worse():
    from pathlib import Path
    from glycosylator.optimizers.glycan_optimizer import _GlycanUnit

    scaffold = gl.Protein.from_pdb(Path(base.HOME) / "support/examples/small.pdb")
    glycan = gl.read_iupac("test", "Neu5Ac(a2-3)Gal(b2-4)GlcNAc")

    sites = list(scaffold.find_n_linked_sites().values())[0]
    scaffold.apply_standard_bonds_for(sites)
    scaffold.add_hydrogens([atom for res in sites for atom in res.atoms])
    scaffold.attach(glycan, residues=sites, chain="new")

    root, glycan = next(iter(scaffold.get_glycans().items()))
    unit = _GlycanUnit(scaffold, root, glycan, include_n_ancestors=2)
    assert unit.n_edges > 3
    assert len(unit.root_edges) == 3

    # batched rotations must match the rotations of the molecule itself
    angles = np.random.default_rng(1).uniform(-np.pi, np.pi, unit.n_edges)
    X = unit.rotate(angles[None])[0]
    reference = scaffold.copy()
    for angle, (a, b) in zip(angles, unit.edges):
        reference._AtomGraph.rotate_around_edge(
            reference.get_atom(a.serial_number),
            reference.get_atom(b.serial_number),
            angle,
            descendants_only=True,
        )
    Y = np.array([reference.get_atom(a.serial_number).coord for a in unit.atom_list])
    assert np.allclose(X, Y, atol=1e-6)

    static = np.array(
        [a.coord for a in scaffold.get_atoms() if a.element != "H" and a not in unit.atoms]
    )
    from scipy.spatial import cKDTree

    before = unit.clash_cost(cKDTree(static), 2.8)

    out = gl.optimize_glycans(scaffold, seed=0, maxiter=10)
    assert out is not scaffold

    root, glycan = next(iter(out.get_glycans().items()))
    unit = _GlycanUnit(out, root, glycan, include_n_ancestors=2)
    after = unit.clash_cost(cKDTree(static), 2.8)
    assert after <= before
