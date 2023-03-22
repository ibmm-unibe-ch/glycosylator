"""
Test the behaviour of force fields and the abstract classes that store their behaviour
"""

import glycosylator as gl
import base


def test_default_settings():
    top = gl.utils.defaults.get_default_topology()
    assert top is not None, "No topology is made"

    prm = gl.utils.defaults.get_default_parameters()
    assert prm is not None, "No parameters are made"


def test_make_charmm_topology():
    top = gl.force_fields.charmm.CHARMMTopology()
    assert top is not None, "No topology is made"

    top = gl.force_fields.charmm.CHARMMTopology.from_file(base.CHARMM_TOPOLOGY_FILE)
    assert top is not None, "No topology is made"

    assert len(top.residues) == 6


def test_make_charmm_parameters():
    prm = gl.force_fields.charmm.CHARMMParameters()
    assert prm is not None, "No parameters are made"

    prm = gl.force_fields.charmm.CHARMMParameters.from_file(base.CHARMM_PARAMETERS_FILE)
    assert prm is not None, "No parameters are made"

    assert len(prm.masses) == 57
    assert len(prm.bonds) == 153
    assert len(prm.angles) == 438
    assert len(prm.dihedrals) == 839
    assert len(prm.impropers) == 14


def test_residue():
    top = gl.utils.defaults.get_default_topology()

    assert top.has_residue("MAN"), "No mannose found!"
    man = top.get_residue("MAN")

    assert len(man.atoms) == 24
    assert len(man.bonds) == 24

    assert len(man.get_internal_coordinates()) == 22
    assert len(man.get_internal_coordinates("C1", "C2", "C3", "O3")) == 1
    assert len(man.get_internal_coordinates("C1", "C2", "C3", None)) == 0
    assert len(man.get_internal_coordinates("C1", "C2", None, "O3", mode="partial")) == 1
    assert len(man.get_internal_coordinates("C1", None, None, "O3", mode="partial")) == 1
    assert len(man.get_internal_coordinates("C1", "O3", mode="anywhere_partial")) == 13
    assert len(man.get_internal_coordinates("C1", "O3", mode="anywhere")) == 1

    assert man.get_bond("C1", "C2") is not None
