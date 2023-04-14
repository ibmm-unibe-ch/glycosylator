"""
Tests to check the behaviour of the gl.AtomGraph and gl.ResidueGraph object
"""

import random
import numpy as np
import glycosylator as gl
import base


# =================================================================
# AtomGraph Tests
# =================================================================


def test_atom_graph_from_molecule():

    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol.lock_all()
    graph = gl.graphs.AtomGraph.from_molecule(mol)

    assert graph is not None, "No molecule is made"

    _received = len(list(graph.bonds))
    _expected = len(list(mol.bonds))
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    assert len(graph._locked_edges) == 0, "Molecule is not locked"

    graph = gl.graphs.AtomGraph.from_molecule(mol, locked=True)

    assert len(graph._locked_edges) == len(mol.locked_bonds), "Molecule is not locked"


def test_atom_graph_pdb_one_residue_is_non_empty():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)
    _received = len(list(mol.bonds))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_pdb_multi_residue_is_non_empty():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)
    _received = len(list(mol.bonds))
    assert _received > 0, f"Expected to find bonds, got {_received}"

    _received = len(list(mol.atoms))
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_one_residue_get_neighbors():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)

    neigs = mol.get_neighbors("C1")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_multi_residue_get_neighbors():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)

    neigs = mol.get_neighbors("C1")
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 11
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_get_descendants():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)

    o3 = next(i for i in mol.structure.get_atoms() if i.id == "O3")
    ho3 = next(i for i in mol.structure.get_atoms() if i.id == "HO3")

    _received = mol.get_descendants(ho3, o3)
    _received = len(_received)
    _expected = 22
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = mol.get_descendants(o3, ho3)
    _received = len(_received)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    # if there is only one atom, then no directionality is available and no descendants
    # can be found!
    try:
        _received = mol.get_descendants(o3, o3)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected a ValueError to be raised")

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    _received = mol.get_descendants(c5, c6)
    _received = len(_received)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"


def test_atom_graph_rotate_descendants_only():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    descendants = mol.get_descendants(c5, c6)
    others = set(i for i in mol.atoms if i not in descendants)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])

    mol.rotate_around_edge(c5, c6, np.radians(35), descendants_only=True)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])

    assert np.all(current_others == new_others), "Other atoms have also moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"


def test_atom_graph_rotate_all():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = gl.graphs.AtomGraph.from_molecule(mol)

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    descendants = mol.get_descendants(c5, c6)
    others = set(i for i in mol.atoms if i not in descendants)
    others.remove(c6)
    others.remove(c5)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])
    current_ref = np.array((c5.coord, c6.coord))

    mol.rotate_around_edge(c5, c6, np.radians(35), descendants_only=False)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])
    new_ref = np.array((c5.coord, c6.coord))

    assert not np.all(current_others == new_others), "Other atoms have not moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"
    assert np.allclose(current_ref, new_ref), "Reference atoms have moved"


# =================================================================
# ResidueGraph tests
# =================================================================


def test_residue_graph_from_molecule():

    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol.lock_all()
    graph_simple = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    v = gl.utils.visual.MoleculeViewer3D(mol)
    for edge in mol.get_residue_connections(True):
        v.draw_vector(
            f"""{edge[0].full_id[3:]} ---> {edge[1].full_id[3:]}""",
            edge[0].coord,
            1.2 * (edge[1].coord - edge[0].coord),
            color="magenta",
        )
    v.show()

    assert graph_simple is not None, "No molecule is made"

    _received = len(list(graph_simple.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    assert len(graph_simple._locked_edges) == 0, "Molecule is not locked"

    v = gl.utils.visual.MoleculeViewer3D(graph_simple)
    v.show()

    graph_detailed = gl.graphs.ResidueGraph.from_molecule(mol, detailed=True)

    assert graph_detailed is not None, "No molecule is made"

    v = gl.utils.visual.MoleculeViewer3D(graph_detailed)
    v.show()

    _received = len(list(graph_detailed.bonds))
    _expected = 40
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    assert len(graph_detailed._locked_edges) != 0, "Molecule is not locked"


def test_residue_graph_multi_residue_get_neighbors():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)

    mol = gl.graphs.ResidueGraph.from_molecule(mol, detailed=True)

    neigs = mol.get_neighbors("C1")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = mol.get_neighbors("MAN")
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = mol.get_neighbors("BMA")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = mol.get_neighbors("BMA", 2)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 6
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = mol.get_neighbors("BMA", 2, "at")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_residue_graph_residue_order():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol = mol.make_residue_graph(detailed=False)

    mol2 = gl.Molecule.from_pdb(base.MANNOSE9)
    mol2.infer_bonds(restrict_residues=False)
    mol2 = mol2.make_residue_graph(detailed=False)

    mol3 = gl.Molecule.from_pdb(base.MANNOSE9)
    mol3.infer_bonds(restrict_residues=False)
    mol3 = mol3.make_residue_graph(detailed=False)

    mol4 = gl.Molecule.from_pdb(base.MANNOSE9)
    mol4.infer_bonds(restrict_residues=False)
    mol4 = mol4.make_residue_graph(detailed=False)

    for i in range(10):
        rdx = random.randint(0, len(mol.residues) - 1)
        assert mol.residues[rdx] == mol2.residues[rdx], "Residue order is not preserved"
        assert mol.residues[rdx] == mol3.residues[rdx], "Residue order is not preserved"
        assert mol.residues[rdx] == mol4.residues[rdx], "Residue order is not preserved"


def test_residue_graph_get_descendants():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    v = gl.utils.visual.MoleculeViewer3D(graph)

    nag2 = mol.residues[0]
    nag3 = mol.residues[1]

    v.draw_point("nag2", nag2.coord, color="red")
    v.draw_point("nag3", nag3.coord, color="blue")
    v.show()

    _received = graph.get_descendants(nag2, nag3)
    _received = len(_received)
    _expected = 9
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = graph.get_descendants(nag3, nag2)
    _received = len(_received)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    bma = mol.residues[2]

    _received = graph.get_descendants(nag3, bma)
    _received = len(_received)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = graph.get_descendants(bma, nag3)
    _received = len(_received)
    _expected = 1
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"


def test_residue_graph_rotate_descendants_only():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    nag3 = mol.residues[1]
    bma = mol.residues[2]

    descendants = graph.get_descendants(nag3, bma)
    others = set(i for i in graph.residues if i not in descendants)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])

    graph.rotate_around_edge(nag3, bma, np.radians(35), descendants_only=True)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])

    assert np.all(current_others == new_others), "Other residues have also moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"


def test_residue_graph_rotate_all():
    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    nag3 = mol.residues[1]
    bma = mol.residues[2]

    descendants = graph.get_descendants(nag3, bma)
    others = set(i for i in graph.residues if i not in descendants)
    others.remove(nag3)
    others.remove(bma)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])
    current_ref = np.array((nag3.coord, bma.coord))

    graph.rotate_around_edge(nag3, bma, np.radians(35), descendants_only=False)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])
    new_ref = np.array((nag3.coord, bma.coord))

    assert not np.all(current_others == new_others), "Other residues have not moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"
    assert np.allclose(current_ref, new_ref), "Reference residues have moved"


def test_residue_graph_detailed():

    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    assert len(mol.nodes) == 11, "Wrong number of nodes"
    assert len(mol.bonds) == 10, "Wrong number of bonds"

    mol.make_detailed()

    assert len(mol.edges) == 40, "Wrong number of edges"
    assert len(mol.nodes) == 41, "Wrong number of nodes"


def test_residue_graph_detailed_get_neighbors():

    mol = gl.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol = gl.graphs.ResidueGraph.from_molecule(mol, detailed=False)
    mol.make_detailed()

    v = gl.utils.visual.MoleculeViewer3D(mol)
    v.show()

    neigs = mol.get_neighbors(3)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"
    _received = len(neigs)
    _expected = 2
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    ids = [i.id for i in neigs]
    assert "C4" in ids, "Wrong neighbor"
    assert "C1" in ids, "Wrong neighbor"

    has_residues = False
    has_atoms = False
    for node in mol.nodes:
        if node.__class__.__name__ == "Residue":
            has_residues = True
        if node.__class__.__name__ == "Atom":
            has_atoms = True

    assert has_residues, "Residues not found"
    assert has_atoms, "Atoms not found"


def test_atom_graph_lock():

    glc = gl.Molecule.from_compound("GLC")
    glc.repeat(5, "14bb")
    glc.lock_all()

    g = glc._AtomGraph

    assert len(g.nodes) == len(glc.atoms)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) > 0

    glc.unlock_all()
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 0

    glc.lock_bond(4, 5)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 1

    glc.lock_bond(4, 5, both_ways=True)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 2

    old = len(glc.bonds)
    old_edges = len(g.edges)

    glc.add_bond(5, 12)
    assert len(glc.bonds) == old + 1
    assert len(g.edges) == old_edges + 1
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 2

    glc.lock_bond(5, 12)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 3

    glc.remove_bond(5, 12)
    assert len(glc.bonds) == old
    assert len(g.edges) == old_edges
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 2
