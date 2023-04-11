"""
Tests to check the behaviour of the gl.Molecule object
"""

import os
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


def test_can_write_pdb():

    import os

    glc = gl.Molecule.from_compound("GLC")
    assert glc is not None

    try:
        glc.to_pdb("test.pdb")
    except Exception as e:
        raise e

    os.remove("test.pdb")


def test_molecule_from_compound():

    glc = gl.Molecule.from_compound("GLC")
    assert glc is not None
    assert len(glc.atoms) == 24
    assert len(glc.bonds) == 24

    a = glc.atoms[0]
    assert a.full_id[0] == "GLC"
    assert a.full_id[1] == 0
    assert a.full_id[2] == "A"
    assert a.full_id[3] == ("H_GLC", 1, " ")
    assert a.full_id[4] == ("C1", " ")

    try:
        glc2 = gl.Molecule.from_compound("D-glucose")
    except ValueError:
        pass
    except Exception as e:
        raise e


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

    other = gl.Molecule.from_compound("GLC")

    residues_pre = len(mol.residues)
    atoms_pre = len(mol.atoms)

    new = bio.Residue.Residue((" ", 1, " "), "NEW", " ")
    mol.add_residues(new)

    assert len(mol.residues) == residues_pre + 1

    _new = mol.get_residue(name="NEW")
    assert _new is not None
    assert _new.id[1] != 1

    mol.remove_residues(_new)
    assert len(mol.residues) == residues_pre

    assert "NEW" not in [i.resname for i in mol.residues]

    atoms_pre = len(mol.atoms)
    mol.add_residues(*other.residues, _copy=True)
    assert len(other.residues) != 0
    assert len(mol.residues) == residues_pre + len(other.residues)
    assert len(mol.atoms) == atoms_pre + len(other.atoms)

    _seen_serials = set()
    for atom in mol.atoms:
        assert atom.serial_number not in _seen_serials
        _seen_serials.add(atom.serial_number)

    mol.remove_residues(2)
    assert len(mol.residues) == residues_pre
    assert len(mol.atoms) == atoms_pre

    new = gl.Molecule.from_compound("GLC")
    new = new.residues[0]
    mol.add_residues(new)
    assert mol.residues[-1].resname == "GLC"


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

    mol = gl.Molecule.from_compound("MAN")
    other = deepcopy(mol)

    mol.adjust_indexing(other)

    assert len(mol.atoms) == len(other.atoms)
    assert len(mol.residues) == len(other.residues)

    assert mol.residues[0].id[1] == 1
    assert other.residues[0].id[1] == 2

    assert mol.atoms[0].get_serial_number() == 1
    assert other.atoms[0].get_serial_number() == len(mol.atoms) + 1


def test_adjust_indexing_with_add_residues():

    mol = gl.Molecule.from_compound("MAN")
    other = deepcopy(mol)

    mol.adjust_indexing(other)

    before = len(mol.atoms)
    mol.add_residues(*other.residues)

    assert len(mol.atoms) == before * 2

    _seen_serials = set()
    for atom in mol.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    assert mol.residues[0].id[1] == 1
    assert len(mol.residues) == 2


def test_set_patch():

    mol = gl.Molecule.from_compound("GLC")
    mol.set_patch("14bb")

    assert mol._patch is not None
    assert not isinstance(mol._patch, str)

    mol.set_patch()
    assert mol._patch is None

    # set the patch with fancy dunder methods
    mol % "14bb"

    assert mol._patch is not None
    assert not isinstance(mol._patch, str)

    mol % None
    assert mol._patch is None

    mol %= "14bb"

    assert mol is not None
    assert mol._patch is not None
    assert not isinstance(mol._patch, str)


def test_attach():

    glc = gl.Molecule.from_compound("GLC")
    glc.set_patch("14bb")

    glc2 = deepcopy(glc)

    _current_residues = len(glc.residues)
    glc.attach(glc2)

    assert len(glc.residues) == _current_residues * 2

    # glc2 should not have been affected since it was a copy
    assert len(glc2.residues) == _current_residues

    # attach with fancy dunder methods
    glc = deepcopy(glc2)

    new = glc + glc2
    assert new is not glc
    assert len(new.residues) == _current_residues * 2
    assert len(glc.residues) == _current_residues
    assert len(glc2.residues) == _current_residues

    glc += glc2
    assert len(glc.residues) == _current_residues * 2
    assert len(glc2.residues) == _current_residues


def test_multiply():

    man = gl.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man._locked_bonds)

    n = 10
    man % "14bb"
    man = man * n

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man._locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert n * 0.75 * pre_locked < new_locked < pre_locked * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.angles.values():
        assert 100 < angle < 130

    v = gl.utils.visual.MoleculeViewer3D(man)
    v.show()


def test_repeat():

    man = gl.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man._locked_bonds)

    n = 10
    man.repeat(n, "14bb")
   
    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man._locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert n * 0.75 * pre_locked < new_locked < pre_locked * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.angles.values():
        assert 100 < angle < 130

    v = gl.utils.visual.MoleculeViewer3D(man)
    v.show()



def test_make_mannose8():

    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN

    """
    bma = gl.Molecule.from_compound("BMA")
    nag = gl.Molecule.from_compound("NAG")
    man = gl.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    nag % "14bb"
    man8 = nag * 2 + bma

    # now we attach the 13ab MAN to the BMA
    man8 % "13ab"
    man8 += man

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man % "16ab"
    man_branch = man * 2

    man_branch % "13ab"
    man_branch @ 1  # attach to residue 1 (not the default last one)
    man_branch += man

    man_branch @ None  # reset the attachment point
    man_branch % "12aa"
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    man8 @ -2
    man_branch @ 1

    man8 % "16ab"
    man8 += man_branch

    man8._bonds = []
    man8._locked_bonds = set()

    man8.infer_bonds(restrict_residues=False)

    v = gl.utils.visual.MoleculeViewer3D(man8)
    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        if not length:
            v.draw_edges([bond], color="magenta", linewidth=3)

    for angle in man8.angles.values():
        assert 100 < angle < 130

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    v.show()

    v = gl.utils.visual.MoleculeViewer3D(man8.make_residue_graph(detailed=False))
    v.show()

    v = gl.utils.visual.MoleculeViewer3D(man8.make_residue_graph(detailed=True))
    v.show()

    try:
        man8.to_pdb("man8.pdb")
    except Exception as e:
        raise e
    else:
        os.remove("man8.pdb")


def test_make_mannose8_2():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN

    """

    bma = gl.Molecule.from_compound("BMA")
    nag = gl.Molecule.from_compound("NAG")
    man = gl.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    man8 = nag % "14bb" * 2 + bma

    # now we attach the 13ab MAN to the BMA
    man8 % "13ab"
    man8 += man

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man_branch = man % "16ab" * 2

    man_branch @ 1 % "13ab"  # attach to residue 1 (not the default last one)
    man_branch += man

    man_branch % "12aa" @ None  # reset the attachment point
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    man8 @ -2 % "16ab"
    man_branch @ 1

    _man8 = deepcopy(man8)
    man8 += man_branch

    # just checkin if the one line syntax is the same as the two line syntax
    _man8 = _man8 @ -2 % "16ab" + man_branch @ 1 
    assert len(man8.atoms) == len(_man8.atoms)

    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        assert length, "Bond length is not in range 0.95 - 1.8"

    for angle in man8.angles.values():
        assert 100 < angle < 130

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    v = gl.utils.visual.MoleculeViewer3D(man8.make_residue_graph(detailed=False))
    v.show()



def test_make_mannose8_3():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN
    ```
    """

    bma = gl.Molecule.from_compound("BMA")
    nag = gl.Molecule.from_compound("NAG")
    man = gl.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    nag.set_patch("14bb")
    nag2 = nag + nag

    man8 = nag2.attach(bma)


    # now we attach the 13ab MAN to the BMA
    man8.set_patch("13ab")
    man8.attach(man)

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man.set_patch("16ab")
    man_branch = man * 2

    man_branch.set_attach_residue(1)
    man_branch.set_patch("13ab")
    man_branch.attach(man)
   
    man_branch.set_patch("12aa")
    man_branch.set_attach_residue()
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    man8.set_attach_residue(-2)
    man8.set_patch("16ab")
    man_branch.set_attach_residue(1)
    man8.attach(man_branch)

    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        assert length, "Bond length is not in range 0.95 - 1.8"

    for angle in man8.angles.values():
        assert 100 < angle < 130

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    v = gl.utils.visual.MoleculeViewer3D(man8)
    colors= ["red", "green", "blue", "magenta", "cyan", "orange", "purple", "pink", "brown", "grey", "black"]
    idx = 0
    for residue in man8.residues:

        for bond in man8.bonds:
            if bond[0].get_parent() == residue and bond[1].get_parent() == residue:
                v.draw_edges([bond], color=colors[idx], linewidth=3)
        idx += 1

    v.show()
    