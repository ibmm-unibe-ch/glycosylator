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

    b = mol.get_residue("MAN", by="name")
    assert b is not None

    assert a is b


def test_can_write_pdb():
    import os

    glc = gl.Molecule.from_compound("GLC")
    assert glc is not None

    try:
        glc.to_pdb("test.pdb")
        with open("test.pdb", "r") as f:
            lines = f.readlines()
        assert len(lines) != 0
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


def test_atomgraph_sync():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol2 = gl.Molecule.from_pdb(base.MANNOSE)

    mol.apply_standard_bonds()
    mol2.apply_standard_bonds()

    assert mol is not mol2

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 24

    assert len(mol.bonds) == len(mol2.bonds)

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph"

    for atom in mol2.get_atoms():
        assert atom in mol2._AtomGraph.nodes, f"Atom {atom} not in graph"

    # make the molecule larger
    mol % "14bb"
    mol += mol2

    assert len(mol.atoms) == 45

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph"

    mol.reindex(2, 100, 5)
    assert len(mol.atoms) == 45

    assert mol.chains[0].id == "B"
    assert mol.residues[0].id[1] == 100

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph after reindex"


def test_molecule_bonds():
    mol = gl.Molecule.from_pdb(base.MANNOSE)

    assert len(mol.bonds) == 0

    mol.apply_standard_bonds()
    assert len(mol.bonds) == 24

    mol.bonds = []
    assert len(mol.bonds) == 0

    mol.infer_bonds()
    assert len(mol.bonds) == 24


def test_molecule_get_bonds():
    glc = gl.Molecule.from_compound("GLC")
    glc.repeat(2, "14bb")

    v = glc.draw()
    v.show()
    assert len(glc.bonds) == 46

    b = glc.get_bonds("O4")
    v.draw_edges(b, color="red")
    assert len(b) == 2

    c1 = glc.get_atom("C1", residue=2)
    b = glc.get_bonds("O4", c1)
    v.draw_edges(b, color="green")
    assert len(b) == 1
    v.show()


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


def test_remove_atoms():
    glc = gl.Molecule.from_compound("GLC")
    glc.remove_atoms("C1", "O4")
    assert len(glc.atoms) == 22
    assert len(glc.bonds) == 18


def test_add_residues():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    other = gl.Molecule.from_compound("GLC")

    residues_pre = len(mol.residues)
    atoms_pre = len(mol.atoms)

    new = bio.Residue.Residue((" ", 1, " "), "NEW", " ")
    mol.add_residues(new)

    assert len(mol.residues) == residues_pre + 1

    _new = mol.get_residue("NEW")
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

    v = gl.utils.visual.MoleculeViewer3D(mol)

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = 10 # 10 degrees

    v.draw_point("first (O3)", first.coord, color="magenta")
    v.draw_point("second (C3)", second.coord, color="magenta")

    # HO3 should not change in dicrection O3---C3
    anchor = mol.get_atom("HO3")

    for i in range(5):
        current_coords = np.array(
            [
                i.coord
                for i in mol.get_atoms()
                if i != first and i != second and i != anchor
            ]
        )
        current_refs = np.array((first.coord, second.coord))
        current_anchor = anchor.coord

        mol.rotate_around_bond(first, second, angle, descendants_only=True)

        new_coords = np.array(
            [
                i.coord
                for i in mol.get_atoms()
                if i != first and i != second and i != anchor
            ]
        )
        new_refs = np.array((first.coord, second.coord))
        new_anchor = anchor.coord

        assert not np.allclose(current_coords, new_coords)
        assert np.allclose(current_refs, new_refs)
        assert np.allclose(current_anchor, new_anchor)

        v.draw_edges(mol.bonds)
        for atom in mol.atoms:
            v.draw_point(atom.id, atom.coord, opacity=0.4, showlegend=False)

    v.show()


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


def test_set_patch_or_recipe():
    mol = gl.Molecule.from_compound("GLC")
    mol.set_patch_or_recipe("14bb")

    assert mol._patch is not None
    assert not isinstance(mol._patch, str)

    mol.set_patch_or_recipe()
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


def test_attach_with_patch():
    glc = gl.Molecule.from_compound("GLC")
    glc.set_patch_or_recipe("14bb")

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


def test_attach_with_recipe():
    recipe = gl.Recipe()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    glc = gl.Molecule.from_compound("GLC")
    glc.set_patch_or_recipe(recipe)

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


def test_multiply_with_patch():
    man = gl.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)

    # set 14bb patch
    n = 10
    man % "14bb"
    man = man * n

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.angles.values():
        assert 100 < angle < 130

    v = gl.utils.visual.MoleculeViewer3D(man)
    v.show()


def test_multiply_with_recipe():
    man = gl.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man.locked_bonds)

    # make recipe for 14bb
    recipe = gl.Recipe()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    # set 14bb recipe
    n = 10
    man % recipe
    man = man * n

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man.locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert pre_locked < new_locked

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
    pre_locked = len(man.locked_bonds)

    n = 10
    man = man.repeat(n, "14bb")

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.angles.values():
        assert 100 < angle < 130

    v = gl.utils.visual.MoleculeViewer3D(man)
    v.show()


def test_repeat_with_recipe():
    man = gl.Molecule.from_compound("GLC")
    man.lock_all()

    # make recipe for 14bb
    recipe = gl.Recipe()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man.locked_bonds)

    n = 10
    man = man.repeat(n, recipe)

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man.locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert pre_locked < new_locked

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
    nag *= 2
    man8 = nag + bma

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

    v = gl.utils.visual.MoleculeViewer3D(man8)
    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        if not length:
            v.draw_edges([bond], color="magenta", linewidth=3)

    for triplet, angle in man8.angles.items():
        if not 100 < angle < 150:
            v.draw_vector(
                None,
                triplet[0].coord,
                triplet[1].coord - triplet[0].coord,
                color="red",
            )
            v.draw_vector(
                None,
                triplet[2].coord,
                triplet[1].coord - triplet[2].coord,
                color="red",
            )

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    v.draw_edges(man8.get_residue_connections(triplet=True), color="red", linewidth=3)

    v.show()

    g = man8.make_residue_graph(detailed=False)
    g2 = man8.make_residue_graph(detailed=True)

    assert g is not g2
    assert len(g.nodes) < len(g2.nodes)

    v = gl.utils.visual.MoleculeViewer3D(g)
    v.show()

    v = gl.utils.visual.MoleculeViewer3D(g2)
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
    nag.set_patch_or_recipe("14bb")
    nag2 = nag + nag

    man8 = nag2.attach(bma)

    # now we attach the 13ab MAN to the BMA
    man8.attach(man, "13ab")

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man.set_patch_or_recipe("16ab")
    man_branch = man * 2

    man_branch.set_attach_residue(1)
    man_branch.set_patch_or_recipe("13ab")
    man_branch.attach(man)

    man_branch.set_patch_or_recipe("12aa")
    man_branch.set_attach_residue()
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    # man8.set_attach_residue(-2)
    # man_branch.set_attach_residue(1)
    # man8.set_patch_or_recipe("16ab")
    man8.attach(man_branch, "16ab", at_residue=-2, other_residue=1)

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
    colors = [
        "red",
        "green",
        "blue",
        "magenta",
        "cyan",
        "orange",
        "purple",
        "pink",
        "brown",
        "grey",
        "black",
    ]
    idx = 0
    for residue in man8.residues:
        for bond in man8.bonds:
            if bond[0].get_parent() == residue and bond[1].get_parent() == residue:
                v.draw_edges([bond], color=colors[idx], linewidth=3)
        idx += 1

    v.show()


def test_make_mannose8_with_recipe():
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
    nag.repeat(2)
    man8 = nag + bma

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man % "16ab"
    man_branch = man * 2
    man_branch @ 1 % "13ab"
    man_branch += man

    man_branch % "12aa" @ None
    man_branch += man

    man8 @ -2
    man_branch @ 1

    recipe_16ab = gl.Recipe("16ab")
    recipe_16ab.add_bond(("O6", "C1"))
    recipe_16ab.add_delete("HO6", "target")
    recipe_16ab.add_delete("HO1", "source")
    recipe_16ab.add_delete("O1", "source")

    man8 % recipe_16ab
    man8 += man_branch

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
    colors = [
        "red",
        "green",
        "blue",
        "magenta",
        "cyan",
        "orange",
        "purple",
        "pink",
        "brown",
        "grey",
        "black",
    ]
    idx = 0
    for residue in man8.residues:
        for bond in man8.bonds:
            if bond[0].get_parent() == residue and bond[1].get_parent() == residue:
                v.draw_edges([bond], color=colors[idx], linewidth=3)
        idx += 1

    v.show()


def test_relabel():
    scrambled = gl.Molecule.from_compound("BMA")

    # relabel to elementwise order without specific connectivity
    counts = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}

    for atom in scrambled.atoms:
        counts[atom.element] += 1
        atom.id = atom.element + str(counts[atom.element])

    # randomly rotate the molecule
    random_vector = np.random.rand(3) * 10

    # choose some rotation axis
    axis = np.random.randint(0, len(scrambled.bonds) - 1)
    axis = scrambled.bonds[axis]

    for atom in scrambled.atoms:
        atom.coord += random_vector

    scrambled.rotate_around_bond(*axis, np.random.randint(0, 180))

    old_scrambled = set(i.id for i in scrambled.atoms)
    old_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    # now we relabel the atoms
    scrambled.autolabel()

    new_scrambled = set(i.id for i in scrambled.atoms)
    new_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    # now we check if the relabeling worked

    assert scrambled.get_bonds("C1", "C2") != []
    assert scrambled.get_bonds("C6", "H61") != []
    assert old_scrambled.difference(new_scrambled) != set()
    assert np.allclose(old_scrambled_coords, new_scrambled_coords)


def test_rotate_descendants_2():
    glc = gl.Molecule.from_compound("GLC")
    glc.lock_all()
    glc.repeat(4, "14bb")
    connections = glc.get_residue_connections()

    v = gl.utils.visual.MoleculeViewer3D(glc)

    v.draw_point(
        "root",
        glc.get_atom(1).coord,
        color="red",
        showlegend=True,
    )

    cdx = 1
    for c in connections:
        glc.unlock_bond(*c, True)

        v.draw_vector(
            f"""[{cdx}]    {c[0].full_id[3:]} -> {c[1].full_id[3:]}""",
            c[0].coord,
            1.2 * (c[1].coord - c[0].coord),
            color="cyan",
            showlegend=True,
        )
        cdx += 1

    v.draw_edges(
        glc.get_bonds(glc.residues[0]),
        color="limegreen",
        linewidth=5,
        opacity=1,
    )

    opacities = np.linspace(0.2, 0.8, len(connections))
    cdx = 0
    for c in connections:
        assert c not in glc.locked_bonds

        descendants = glc.get_descendants(*c)
        ancestors = glc.get_ancestors(*c)

        old_coords_descendants = np.array([i.coord for i in descendants])
        old_coords_ancestors = np.array([i.coord for i in ancestors])

        glc.rotate_around_bond(*c, 2, descendants_only=True)

        new_coords_descendants = np.array([i.coord for i in descendants])
        new_coords_ancestors = np.array([i.coord for i in ancestors])

        assert not np.allclose(old_coords_descendants, new_coords_descendants)
        assert np.allclose(old_coords_ancestors, new_coords_ancestors)

        v.draw_edges(glc.bonds, color="lightgreen", opacity=opacities[cdx])

        v.draw_edges(
            glc.get_bonds(glc.residues[0]),
            color="green",
            linewidth=3,
            opacity=1,  # opacities[cdx],
        )

        cdx += 1
    # for atom in glc.atoms:
    #     v.draw_point(atom.id, atom.coord, color="red", opacity=0.3, showlegend=False)

    # v.draw_edges(glc.bonds, color="red")

    v.draw_edges(glc.get_bonds(glc.residues[0]), color="teal", linewidth=5, opacity=1)
    v.show()
