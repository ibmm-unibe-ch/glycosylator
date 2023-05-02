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


def test_seq():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)

    seq = scaffold.seq
    assert isinstance(seq, dict)
    assert len(seq) == 2

    seq = list(seq.values())
    assert (
        seq[0]
        == "AENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMHTDIISLWDQSLKPCVKLTPLCVTLQCTNVTNNITDDMRGELKNCSFNMTTELRDKKQKVYSLFYRLDVVQINSNKEYRLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVMIRSENITNNAKNILVQFNTPVQINCTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIIRFANSSGGDLEVTTHSFNCGGEFFYCNTSGLFNSTWISNNDSITLPCRIKQIINMWQRIGQAMYAPPIQGVIRCVSNITGLILTRDGGSTNSTTETFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRCKRRV"
    )
    assert (
        seq[1]
        == "VFLGFLGAAGSTMGAASMTLTVQARNLLSGTVWGIKQLQARVLAVERYLRDQQLLGIWGCSGKLICCTNVPWNSSWSNRNLSEIWDNMTWLQWDKEISNYTQIIYGLLEESQNQQEKNEQDLLALD"
    )


def test_atomgraph_sync():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)

    for atom in scaffold.get_atoms():
        assert (
            atom in scaffold._AtomGraph.nodes
        ), f"Atom {atom} not in AtomGraph before reindex"

    scaffold.reindex()
    for atom in scaffold.get_atoms():
        assert (
            atom in scaffold._AtomGraph.nodes
        ), f"Atom {atom} not in AtomGraph after reindex"

    # add some residues
    old_residues = len(scaffold.residues)
    ser = gl.Molecule.from_compound("SER").residues[0]
    scaffold.add_residues(ser)

    assert len(scaffold.residues) == old_residues + 1

    for atom in scaffold.get_atoms():
        assert atom in scaffold._AtomGraph.nodes, f"Atom {atom} not in AtomGraph"


def test_exclude_include_chain():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    scaffold.reindex()

    scaffold.exclude_chain("B")

    sequon = "(N)(?=[A-OQ-Z][ST])"
    residues = scaffold.find(sequon)

    assert "A" in residues.keys()
    assert "B" not in residues.keys()

    residues = list(residues.values())
    assert len(residues[0]) == 20

    scaffold.include_chain("B")
    residues = scaffold.find(sequon)

    assert "A" in residues.keys()
    assert "B" in residues.keys()

    residues = list(residues.values())
    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


def test_exclude_include_residue():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    scaffold.reindex()

    asns = scaffold.find()
    asns = asns["A"][:10]

    for a in asns:
        scaffold.exclude_residue(a)

    residues = scaffold.find()
    residues = list(residues.values())

    assert len(residues[0]) == 20 - 10  # we excluded 10 of the 20 ASNs
    assert len(residues[1]) == 4

    for a in asns:
        scaffold.include_residue(a)

    residues = scaffold.find()
    residues = list(residues.values())

    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


def test_find_sequon():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    scaffold.reindex()

    sequon = "(N)(?=[A-OQ-Z][ST])"
    residues = scaffold.find(sequon)
    residues = list(residues.values())

    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


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
    sequon = "(N)(?=[^P])(?=(S|T))"

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


def test_hollow_out():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    scaffold.reindex()

    old_residues = len(scaffold.residues)

    # hollow out the protein
    scaffold.hollow_out()

    new_residues = len(scaffold.residues)

    assert new_residues < old_residues


def test_fill():
    scaffold = gl.Scaffold.from_pdb(base.PROTEIN)
    scaffold.reindex()

    old_residues = len(scaffold.residues)

    old_residues_dict = {}
    for chain in scaffold.chains:
        old_residues_dict[chain] = set(chain.child_list)

    # hollow out the protein
    scaffold.hollow_out()

    new_residues = len(scaffold.residues)

    assert new_residues < old_residues

    new_residues_dict = {}
    for chain in scaffold.chains:
        assert old_residues_dict[chain].issuperset(set(chain.child_list))

    # fill the protein back in
    scaffold.fill()

    new_residues = len(scaffold.residues)
    assert new_residues == old_residues

    new_residues_dict = {}
    for chain in scaffold.chains:
        new_residues_dict[chain] = set(chain.child_list)
    assert new_residues_dict == old_residues_dict
