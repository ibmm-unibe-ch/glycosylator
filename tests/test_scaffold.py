"""
Tests to check the behaviour of the Scaffold class
"""

import base
import glycosylator as gl
import os


def test_from_pdb():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    assert scaffold is not None

    assert len(scaffold.chains) == 2
    assert len(scaffold.residues) == 579
    assert len(scaffold.atoms) == 9102
    assert len(scaffold.bonds) == 0

    b = scaffold.infer_bonds()

    assert len(b) == len(scaffold.bonds) > 8500

    new_b = scaffold.infer_residue_connections(triplet=False)
    assert 570 < len(new_b) < 578
    assert len(scaffold.bonds) > len(b)


def test_attach_glycan_to_one_residue():
    glycan = gl.Molecule.from_pdb(base.MANNOSE9)
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    # make a recipe for attaching to ASN
    recipe = gl.Recipe()
    recipe.add_bond(("ND2", "C1"))
    recipe.add_delete("O1", "source")
    recipe.add_delete("HO1", "source")
    recipe.add_delete("HD22", "target")

    # set recipe as default modifier
    scaffold % recipe

    # find some ASNs
    asns = scaffold.get_residues("ASN", by="name", chain="A")
    assert len(asns) > 0

    # get the ND2 of the ASN and set as root atom for attaching
    scaffold.root_atom = scaffold.get_atom("ND2", residue=asns[0])
    glycan.root_atom = 1

    # attach to the first ASN
    current_residues = len(scaffold.residues)
    old_connections = len(scaffold._molecule_connections)
    scaffold.attach(glycan)

    assert len(scaffold.residues) == current_residues + 11
    assert len(scaffold._molecule_connections) > old_connections

    scaffold.to_pdb("test_attach_glycan_to_one_residue.pdb")
    old, new = 0, 0
    with open(base.PROTEIN) as f:
        for line in f:
            old += 1
    with open("test_attach_glycan_to_one_residue.pdb") as f:
        for line in f:
            new += 1
    assert new > old
    os.remove("test_attach_glycan_to_one_residue.pdb")


def test_attach_glycan_to_residue_list():
    glycan = gl.Molecule.from_pdb(base.MANNOSE9)
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    # make a recipe for attaching to ASN
    recipe = gl.Recipe()
    recipe.add_bond(("ND2", "C1"))
    recipe.add_delete("O1", "source")
    recipe.add_delete("HO1", "source")
    recipe.add_delete("HD22", "target")

    # set recipe as default modifier
    scaffold % recipe

    # find some ASNs
    asns = scaffold.get_residues("ASN", by="name", chain="A")
    assert len(asns) > 5

    asns = asns[:5]

    glycan.root_atom = 1

    current_residues = len(scaffold.residues)
    old_connections = len(scaffold._molecule_connections)
    scaffold.attach(glycan, residues=asns, at_atom="ND2")

    assert len(scaffold.residues) > current_residues
    assert len(scaffold._molecule_connections) > old_connections

    scaffold.to_pdb("test_attach_glycan_to_residue_list.pdb")

    with open("test_attach_glycan_to_residue_list.pdb") as f:
        lines = f.readlines()
    assert len(lines) > 10000
    os.remove("test_attach_glycan_to_residue_list.pdb")


def test_attach_glycan_with_sequon():
    glycan = gl.Molecule.from_pdb(base.MANNOSE9)
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    recipe = gl.Recipe()
    recipe.add_bond(("ND2", "C1"))
    recipe.add_delete("O1", "source")
    recipe.add_delete("HO1", "source")
    recipe.add_delete("HD22", "target")

    scaffold % recipe

    # use the N-linked sequon for ASNs
    sequon = "(N)[^P](?=S|T)"

    residues = scaffold.find(sequon)
    assert len(residues) > 0

    glycan.root_atom = 1

    old_residues = len(scaffold.residues)
    scaffold.attach(glycan, sequon=sequon, at_atom="ND2")

    new_residues = len(scaffold.residues)
    assert new_residues > old_residues

    scaffold.to_pdb("test_attach_glycan_with_sequon.pdb")

    lines = 0
    with open("test_attach_glycan_with_sequon.pdb") as f:
        for line in f:
            lines += 1
    assert lines > 14000
    os.remove("test_attach_glycan_with_sequon.pdb")
