"""
Tests for the auxiliary structure module
"""

from copy import deepcopy
import numpy as np
import glycosylator as gl
import base
import Bio.PDB as bio

MARGIN = 1.5 * 1e-2
MANNOSE = bio.PDBParser().get_structure("MAN", base.MANNOSE)


def test_missing_proper_1():

    to_deletes = {'C1', 'C2', 'C3', 'C4', 'C5', 'O5'}

    for to_delete in to_deletes:

        _man = MANNOSE.copy()
        _man = next(_man.get_residues())
        true_coords = _man.child_dict.get(to_delete)
        assert true_coords is not None, f"No atom {to_delete} found!"

        true_coords = true_coords.coord

        _man.detach_child(to_delete)
        assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

        gl.utils.structural.fill_missing_atoms(_man)

        assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

        new_coords = _man.child_dict.get(to_delete).coord

        _diff = np.sum(np.abs(new_coords - true_coords))
        assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


def test_missing_proper_4():

    to_deletes = {'H1', 'HO1', 'HO2', 'HO3', 'HO4', 'HO6', 'O1'}

    for to_delete in to_deletes:

        _man = MANNOSE.copy()
        _man = next(_man.get_residues())

        true_coords = _man.child_dict.get(to_delete)
        assert true_coords is not None, f"No atom {to_delete} found!"

        true_coords = true_coords.coord

        _man.detach_child(to_delete)
        assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

        gl.utils.structural.fill_missing_atoms(_man)

        assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

        new_coords = _man.child_dict.get(to_delete).coord

        _diff = np.sum(np.abs(new_coords - true_coords))
        assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


def test_missing_improper_1():

    to_deletes = {'C6', 'O3', 'O4', 'O6'}

    for to_delete in to_deletes:

        _man = MANNOSE.copy()
        _man = next(_man.get_residues())

        top = deepcopy(gl.utils.defaults.get_default_topology())
        abstract = top.get_residue(_man.resname)
        for idx, i in enumerate(abstract.internal_coordinates):
            if i.is_proper:
                del abstract.internal_coordinates[idx]

        true_coords = _man.child_dict.get(to_delete)
        assert true_coords is not None, f"No atom {to_delete} found!"

        true_coords = true_coords.coord

        _man.detach_child(to_delete)
        assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

        gl.utils.structural.fill_missing_atoms(_man, top)

        assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

        new_coords = _man.child_dict.get(to_delete).coord

        _diff = np.sum(np.abs(new_coords - true_coords))
        assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


def test_missing_improper_4():

    to_deletes = {'H2', 'H3', 'H4', 'H5', 'H61', 'H62'}

    for to_delete in to_deletes:

        _man = MANNOSE.copy()
        _man = next(_man.get_residues())
        top = deepcopy(gl.utils.defaults.get_default_topology())
        abstract = top.get_residue(_man.resname)
        for idx, i in enumerate(abstract.internal_coordinates):
            if i.is_proper:
                del abstract.internal_coordinates[idx]

        true_coords = _man.child_dict.get(to_delete)
        assert true_coords is not None, f"No atom {to_delete} found!"

        true_coords = true_coords.coord

        _man.detach_child(to_delete)
        assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

        gl.utils.structural.fill_missing_atoms(_man, top)

        assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

        new_coords = _man.child_dict.get(to_delete).coord

        _diff = np.sum(np.abs(new_coords - true_coords))
        assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


def test_missing_one_random_atom():

    _man = MANNOSE.copy()
    _man = next(_man.get_residues())

    to_delete = np.random.choice(list(_man.child_dict.keys()), 1)
    to_delete = to_delete[0]
    true_coords = _man.child_dict.get(to_delete)
    assert true_coords is not None, f"No atom {to_delete} found!"

    true_coords = true_coords.coord

    _man.detach_child(to_delete)

    gl.utils.structural.fill_missing_atoms(_man)

    assert _man.child_dict.get(to_delete) is not None, f"Atom {to_delete} was not added again!"

    new_coords = _man.child_dict.get(to_delete).coord

    _diff = np.sum(np.abs(new_coords - true_coords))
    assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


def test_missing_multiple_random_atoms():

    _man = MANNOSE.copy()
    _man = next(_man.get_residues())

    to_delete = np.random.choice(list(_man.child_dict.keys()), 5, replace=False)

    true_coords = {i: _man.child_dict.get(i).coord for i in to_delete}
    for i in to_delete:
        _man.detach_child(i)

    gl.utils.structural.fill_missing_atoms(_man)

    for i in to_delete:
        assert _man.child_dict.get(i) is not None, f"Atom {i} was not added again!"

        new_coords = _man.child_dict.get(i).coord

        _diff = np.sum(np.abs(new_coords - true_coords[i]))
        assert _diff < MARGIN, f"[{i}] Difference in coordinates is {_diff=}"


def test_missing_multiple_random_atoms_galactose():

    GALACTOSE = bio.PDBParser().get_structure("GAL", base.GALACTOSE)

    _gal = GALACTOSE.copy()
    _gal = next(_gal.get_residues())

    to_delete = np.random.choice(list(_gal.child_dict.keys()), 5, replace=False)

    true_coords = {i: _gal.child_dict.get(i).coord for i in to_delete}
    for i in to_delete:
        _gal.detach_child(i)

    gl.utils.structural.fill_missing_atoms(_gal)

    for i in to_delete:
        assert _gal.child_dict.get(i) is not None, f"Atom {i} was not added again!"

        new_coords = _gal.child_dict.get(i).coord

        _diff = np.sum(np.abs(new_coords - true_coords[i]))
        assert _diff < MARGIN, f"[{i}] Difference in coordinates is {_diff=}"


def test_apply_standard_bonds():

    bonds = gl.utils.structural.apply_standard_bonds(MANNOSE)

    _recieved = len(bonds)
    _expected = 24
    _what = "bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    bonds = [set((i.id, j.id)) for i, j in bonds]

    _bond = set(("C5", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C6", "O6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C4", "C3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "O3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_apply_standard_bonds_one_atom():

    atom = {i.id: i for i in MANNOSE.get_atoms()}
    atom = atom.get("C1")

    bonds = gl.utils.structural.apply_standard_bonds(atom)
    bonds = [set((i.id, j.id)) for i, j in bonds]

    _recieved = len(bonds)
    _expected = 4
    _what = "bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "O1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_infer_bonds():

    bonds = gl.utils.structural.infer_bonds(MANNOSE)

    _recieved = len(bonds)
    _expected = 24
    _what = "bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    bonds = [set((i.id, j.id)) for i, j in bonds]

    _bond = set(("C5", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C6", "O6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C4", "C3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "O3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_infer_residue_connections():

    _man9 = bio.PDBParser().get_structure("MANNOSE9", base.MANNOSE9)
    bonds = gl.utils.structural.infer_residue_connections(_man9)

    connections = [set((i.get_parent()._id[1], j.get_parent()._id[1])) for i, j in bonds]
    bonds = [set((i.id, j.id)) for i, j in bonds]

    _received = len(bonds)
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _bond = set(("O4", "C1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("O3", "C1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((2, 3))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((5, 6))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((5, 8))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert _recieved == _expected, f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"
