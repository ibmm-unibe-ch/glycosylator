"""
Tests to check the behaviour of the gl.Molecule object
"""

import numpy as np
import glycosylator as gl
import base

# # PRODY IS NOW REPLACED BY MDTRAJ SO THIS IS DEPCRECATED AND DOES NOT WORK ANYMORE
# def test_molecule_from_prody():
#     """
#     Test the initialisation of the Molecule object using a pre-made prody atom group
#     """
#     import prody

#     _raw_atom_group = prody.parsePDB(base.MANNOSE9)
#     serials = list(_raw_atom_group.getSerials())
#     _root = _raw_atom_group[serials.index(1)]

#     mol = gl.Molecule.from_prody(_raw_atom_group, _root)
#     assert mol is not None, "No molecule is made"

#     _received = len(mol.atom_group)
#     _expected = 246
#     assert _received == _expected, f"Expected {_expected} atoms, got {_received}"

#     _received = len(mol.torsionals)
#     _expected = 236
#     assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_init_from_pdb():
    """
    Test the initialisation of the Molecule object using a pdb file
    """
    _root = 0
    mol = gl.Molecule.from_pdb(base.MANNOSE9, _root)
    assert mol is not None, "No molecule is made"

    _received = mol.n_atoms
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_init_from_trajectory():
    """
    Test the initialisation of the Molecule object using a pre-made mdtraj trajectory
    """
    import mdtraj as md

    _raw_traj = md.load(base.MANNOSE9)
    _root = 0

    mol = gl.Molecule.from_trajectory(_raw_traj, _root)
    assert mol is not None, "No molecule is made"

    _received = mol.n_atoms
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"
