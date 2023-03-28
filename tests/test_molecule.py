"""
Tests to check the behaviour of the gl.Molecule object
"""

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

    a = mol.get_atom(serial_number=1)
    assert a is not None

    b = mol.get_atom(id="C1")
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

    could_match_any = False
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
        assert np.abs(_angle - angle) < 0.01, f"Angle {angle} does not match reference {_angle}"
        could_match_any = True

    assert could_match_any, "No angles could be matched"


def test_dihedrals():
    mol = gl.Molecule.from_pdb(base.MANNOSE)
    mol.apply_standard_bonds()

    top = gl.utils.get_default_topology()
    abstract = top.get_residue("MAN")

    could_match_any = False
    for quartet, dihedral in mol.dihedrals.items():

        quartet = [i.id for i in quartet]
        ics = abstract.get_internal_coordinates(*quartet)
        if len(ics) == 0:
            ics = abstract.get_internal_coordinates(*quartet[::-1])
            if len(ics) == 0:
                continue
        _dihedral = ics[0].dihedral

        assert np.abs(_dihedral - dihedral) < 0.01, f"Dihedral {dihedral} does not match reference {_dihedral}"
        could_match_any = True

    assert could_match_any, "No dihedrals could be matched"
