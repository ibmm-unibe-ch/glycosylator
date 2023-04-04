"""
Tests to check the behaviour of the gl.Molecule object
"""

from copy import deepcopy
import numpy as np
import glycosylator as gl
import Bio.PDB as bio

import base


def test_molecule_basic():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 0

    _mol = gl.utils.defaults.__bioPDBParser__.get_structure("MAN", base.MANNOSE)
    mol = gl.Molecule(_mol)
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 0

    mol = gl.Molecule.from_smiles("OCC1OC(O)C(C(C1O)O)O", add_hydrogens=False)
    assert mol is not None

    assert len(mol.atoms) == 12
    assert len(mol.bonds) == 0

    mol = gl.Molecule.from_smiles("OCC1OC(O)C(C(C1O)O)O")
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 0

    mol = gl.Molecule(_mol)

    a = mol.get_atom(1)
    assert a is not None

    b = mol.get_atom("C1")
    assert b is not None

    a.id = "HIHIHI"
    assert b.id == "HIHIHI"
    assert a is b

    a = mol.get_residue(1)
    assert a is not None

    b = mol.get_residue(name="MAN")
    assert b is not None

    assert a is b


def test_molecule_bonds():
    mol = gl.Molecule.from_pdb(base.MANNOSE)

    mol.apply_standard_bonds()
    assert len(mol.bonds) == 24

    mol._bonds = []

    mol.infer_bonds()
    assert len(mol.bonds) == 24

    mol._bonds = []


def test_angles():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    top = gl.utils.get_default_topology()
    abstract = top.get_residue("MAN")

    for triplet, angle in mol.angles.items():

        triplet = [i.id for i in triplet]

        ics = abstract.get_internal_coordinates(*triplet, None, mode="partial")
        _angle = "bond_angle_123"
        if len(ics) == 0:
            ics = abstract.get_internal_coordinates(None, *triplet, mode="partial")
            _angle = "bond_angle_234"
            if len(ics) == 0:
                continue
        _angle = getattr(ics[0], _angle)
        assert (
            np.abs(_angle - angle) < 0.01
        ), f"Angle {angle} does not match reference {_angle}"


def test_dihedrals():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    top = gl.utils.get_default_topology()
    abstract = top.get_residue("MAN")

    for quartet, dihedral in mol.dihedrals.items():

        quartet = [i.id for i in quartet]
        ics = abstract.get_internal_coordinates(*quartet)
        if len(ics) == 0:
            ics = abstract.get_internal_coordinates(*quartet[::-1])
            if len(ics) == 0:
                continue
        _dihedral = ics[0].dihedral

        assert (
            np.abs(_dihedral - dihedral) < 0.01
        ), f"Dihedral {dihedral} does not match reference {_dihedral}"


def test_add_atoms():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    pre = len(mol.atoms)

    new = bio.Atom.Atom("C99", np.array((0.5, 1.23, -0.5)), None, 0.0, None, "C99", 1)
    mol.add_atoms(new)

    assert len(mol.atoms) == pre + 1

    _new = mol.get_atom("C99")
    assert _new is not None
    assert _new.serial_number != 1

    mol.remove_atoms(_new)
    assert len(mol.atoms) == pre

    assert "C99" not in [i.id for i in mol.atoms]


def test_add_residues():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    other = deepcopy(mol)

    residues_pre = len(mol.residues)

    new = bio.Residue.Residue((1, "H_NEW", " "), "NEW", " ")
    mol.add_residues(new)

    assert len(mol.residues) == residues_pre + 1

    _new = mol.get_residue(name="NEW")
    assert _new is not None
    assert _new.id[1] != 1

    mol.remove_residues(_new)
    assert len(mol.residues) == residues_pre

    assert "NEW" not in [i.resname for i in mol.residues]

    atoms_pre = len(mol.atoms)
    mol.add_residues(*other.residues)
    assert len(mol.residues) == residues_pre + len(other.residues)
    assert len(mol.atoms) == atoms_pre + len(other.atoms)

    _seen_serials = set()
    for atom in mol.atoms:
        assert atom.serial_number not in _seen_serials
        _seen_serials.add(atom.serial_number)


def test_get_descendants():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    descendants = mol.get_descendants("C5", "C6")
    _received = len(descendants)
    _expected = 4
    assert (
        _received == _expected
    ), f"Expected {_expected} descendants, received {_received}"

    descendants = mol.get_descendants("O3", "C3")
    _received = len(descendants)
    _expected = 21
    assert (
        _received == _expected
    ), f"Expected {_expected} descendants, received {_received}"


def test_rotate_all():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = np.radians(45)

    current_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )
    current_refs = np.array((first.coord, second.coord))
    mol.rotate_around_bond(first, second, angle)

    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )
    new_refs = np.array((first.coord, second.coord))

    assert not np.allclose(current_coords, new_coords)
    assert np.allclose(current_refs, new_refs)

    # and rotate back to revert...
    mol.rotate_around_bond(first, second, -angle)
    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )

    assert np.allclose(current_coords, new_coords)


def test_rotate_some():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = np.radians(45)

    # HO3 should not change in dicrection O3---C3
    anchor = mol.get_atom("HO3")

    current_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    current_refs = np.array((first.coord, second.coord))
    current_anchor = anchor.coord

    mol.rotate_around_bond(first, second, angle, descendants_only=True)

    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    new_refs = np.array((first.coord, second.coord))
    new_anchor = anchor.coord

    assert not np.allclose(current_coords, new_coords)
    assert np.allclose(current_refs, new_refs)
    assert np.allclose(current_anchor, new_anchor)


def test_rotate_some_inverse():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = np.radians(45)

    # HO3 should be the only change in dicrection C3---O3
    anchor = mol.get_atom("HO3")

    current_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    current_refs = np.array((first.coord, second.coord))
    current_anchor = anchor.coord

    mol.rotate_around_bond(second, first, angle, descendants_only=True)

    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    new_refs = np.array((first.coord, second.coord))
    new_anchor = anchor.coord

    assert np.allclose(current_coords, new_coords)
    assert np.allclose(current_refs, new_refs)
    assert not np.allclose(current_anchor, new_anchor)


def test_adjust_indexing():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    other = deepcopy(mol)

    mol.adjust_indexing(other)

    assert len(mol.atoms) == len(other.atoms)
    assert len(mol.residues) == len(other.residues)

    assert mol.residues[0].id[1] == 1
    assert other.residues[0].id[1] == 2

    assert mol.atoms[0].serial_number == 1
    assert other.atoms[0].serial_number == len(mol.atoms) + 1


def test_adjust_indexing_with_add_residues():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    other = deepcopy(mol)

    mol.adjust_indexing(other)

    before = len(mol.atoms)
    mol.add_residues(*other.residues)

    assert len(mol.atoms) == before * 2
    assert mol.atoms[0].serial_number == 1
    assert mol.atoms[before].serial_number == before + 1

    assert mol.residues[0].id[1] == 1
    assert len(mol.residues) == 2
