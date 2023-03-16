"""
Tests to check the behaviour of the gl.AtomGraph and gl.ResidueGraph object
"""

import glycosylator as gl
import base


def test_atom_graph_from_pdb():
    """
    Test the initialisation of the AtomGraph object using a pdb file
    """
    mol = gl.graphs.AtomGraph.from_pdb(base.MANNOSE)
    assert mol is not None, "No molecule is made"

    _received = len(list(mol.bonds))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_residue_graph_from_pdb():
    """
    Test the initialisation of the ResidueGraph object using a pdb file
    """
    mol = gl.graphs.ResidueGraph.from_pdb(base.MANNOSE9)
    assert mol is not None, "No molecule is made"

    _received = len(list(mol.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 11
    assert _received == _expected, f"Expected {_expected} residues, got {_received}"
