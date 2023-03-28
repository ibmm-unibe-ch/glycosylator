"""
Tests to check the behaviour of the gl.AtomGraph and gl.ResidueGraph object
"""

import glycosylator as gl
import base


# =================================================================
# AtomGraph Tests
# =================================================================


def test_atom_graph_pdb_one_residue_is_made():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)
    assert mol is not None, "No molecule is made"


def test_atom_graph_pdb_one_residue_is_non_empty():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)
    _received = len(list(mol.bonds))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_pdb_multi_residue_is_made():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    assert mol is not None, "No molecule is made"


def test_atom_graph_pdb_multi_residue_is_non_empty():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    _received = len(list(mol.bonds))
    assert _received > 0, f"Expected to find bonds, got {_received}"

    _received = len(list(mol.atoms))
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_one_residue_get_neighbors():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)

    neigs = mol.get_neighbors("C1")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_multi_residue_get_neighbors():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)

    neigs = mol.get_neighbors("C1")
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 11
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_get_descendants():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)

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


# =================================================================
# ResidueGraph tests
# =================================================================


def test_residue_graph_pdb_is_made():
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)
    assert mol is not None, "No molecule is made"


def test_residue_graph_pdb_is_non_empty():
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)
    _received = len(list(mol.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 11
    assert _received == _expected, f"Expected {_expected} residues, got {_received}"


def test_residue_graph_atomgraph_is_made():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    mol1 = gl.graphs.ResidueGraph.from_AtomGraph(mol)
    assert mol1 is not None, "No molecule is made"
    assert mol1.to_AtomGraph() == mol, "No link to AtomGraph is made"


def test_residue_graph_atomgraph_is_non_empty():
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    mol = gl.graphs.ResidueGraph.from_AtomGraph(mol)

    _received = len(list(mol.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 11
    assert _received == _expected, f"Expected {_expected} residues, got {_received}"


def test_residue_graph_multi_residue_get_neighbors():
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)

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
    _expected = 7
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = mol.get_neighbors("BMA", 2, "at")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_residue_graph_get_descendants():
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)

    nag2 = mol.residues[0]
    nag3 = mol.residues[1]

    _received = mol.get_descendants(nag2, nag3)
    _received = len(_received)
    _expected = 9
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = mol.get_descendants(nag3, nag2)
    _received = len(_received)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    bma = mol.residues[2]

    _received = mol.get_descendants(nag3, bma)
    _received = len(_received)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = mol.get_descendants(bma, nag3)
    _received = len(_received)
    _expected = 1
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"
