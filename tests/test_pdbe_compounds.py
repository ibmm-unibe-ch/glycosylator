"""
Tests for the PDBe compounds class.
"""

import base
import Bio.PDB as bio
import glycosylator as gl
from glycosylator.resources import pdbe_compounds


def test_from_cif():

    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    assert comps is not None, "Could not load the PDBe compounds from a CIF file."
    assert len(comps.ids) == 25, "The number of compounds is not correct."


def test_getting_compounds():

    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man = comps.get("D-Mannose", by="name")
    assert man is not None

    man = comps.get("C6 H12 O6", by="formula")
    assert man is not None

    man = comps.get("OC1C(O)C(OC(O)C1O)CO", by="smiles")
    assert man is not None

    man = comps.get("MAN")
    assert man is not None

    assert len(man.atoms) == 24
    assert len(man.bonds) == 24


def test_compound_is_same():

    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man = comps.get("D-Mannose", by="name")

    ref = gl.Molecule.from_pdb(base.MANNOSE)
    ref.infer_bonds()

    man_atoms = [(i.id, i.serial_number) for i in man.atoms]
    ref_atoms = [(i.id, i.serial_number) for i in ref.atoms]
    assert man_atoms == ref_atoms

    man_bonds = set((i.id, j.id) for i, j in man.bonds)
    ref_bonds = set((i.id, j.id) for i, j in ref.bonds)

    assert man_bonds == ref_bonds


def test_compound_getting_types():

    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man_mol = comps.get("D-Mannose", by="name")
    assert isinstance(man_mol, gl.Molecule)

    man_dict = comps.get("D-Mannose", by="name", return_type="dict")
    assert isinstance(man_dict, dict)

    man_struct = comps.get("D-Mannose", by="name", return_type="structure")
    assert isinstance(man_struct, bio.Structure.Structure)
    assert len(list(man_struct.get_atoms())) == 24

    man_res = comps.get("D-Mannose", by="name", return_type="residue")
    assert isinstance(man_res, bio.Residue.Residue)
    assert len(list(man_res.get_atoms())) == 24

    glc_mol = comps.get("Glucose", by="name")
    assert isinstance(glc_mol, list)
    assert len(glc_mol) == 2
    assert isinstance(glc_mol[0], gl.Molecule)

    glc_dict = comps.get("Glucose", by="name", return_type="dict")
    assert isinstance(glc_dict, list)
    assert len(glc_dict) == 2
    assert isinstance(glc_dict[0], dict)

    glc_struct = comps.get("Glucose", by="name", return_type="structure")
    assert isinstance(glc_struct, list)
    assert len(glc_struct) == 2
    assert isinstance(glc_struct[0], bio.Structure.Structure)
    assert len(list(glc_struct[0].get_atoms())) == 24