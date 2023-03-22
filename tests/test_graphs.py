"""
Tests to check the behaviour of the gl.AtomGraph and gl.ResidueGraph object
"""

import glycosylator as gl
import base


def test_atom_graph_pdb_one_residue_is_made():
    """
    Test the initialisation of the AtomGraph object using a pdb file
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)
    assert mol is not None, "No molecule is made"


def test_atom_graph_pdb_one_residue_is_correct():
    """
    Test the initialisation of the AtomGraph object using a pdb file
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)
    _received = len(list(mol.bonds))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_pdb_multi_residue_is_made():
    """
    Test the initialisation of the AtomGraph object using a pdb file
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9, infer_bonds=True, apply_standard_bonds=False)
    assert mol is not None, "No molecule is made"


def test_atom_graph_pdb_multi_residue_is_non_empty():
    """
    Test the initialisation of the AtomGraph object using a pdb file
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9, infer_bonds=True, apply_standard_bonds=False)
    _received = len(list(mol.bonds))
    assert _received > 0, f"Expected to find bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_residue_graph_pdb_is_made():
    """
    Test the initialisation of the ResidueGraph object using a pdb file
    """
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)
    assert mol is not None, "No molecule is made"


def test_residue_graph_pdb_is_correct():
    """
    Test the initialisation of the ResidueGraph object using a pdb file
    """
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)
    _received = len(list(mol.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 11
    assert _received == _expected, f"Expected {_expected} residues, got {_received}"


def test_residue_graph_atomgraph_is_made():
    """
    Test the initialisation of the ResidueGraph object using an AtomGraph object
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    mol1 = gl.graphs.ResidueGraph.from_AtomGraph(mol)
    assert mol1 is not None, "No molecule is made"
    assert mol1.to_AtomGraph() == mol, "No link to AtomGraph is made"


def test_residue_graph_atomgraph_is_correct():
    """
    Test the initialisation of the ResidueGraph object using an AtomGraph object
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE9)
    mol = gl.graphs.ResidueGraph.from_AtomGraph(mol)

    _received = len(list(mol.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 11
    assert _received == _expected, f"Expected {_expected} residues, got {_received}"
