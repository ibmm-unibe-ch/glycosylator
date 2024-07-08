"""
Tests to check the behaviour of the Scaffold class
"""

import base
import glycosylator as gl
import os


def test_from_pdb():
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
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
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    seq = scaffold.get_sequence()
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
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

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
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
    scaffold.reindex()

    scaffold.exclude_chain("B")

    sequon = "(N)(?=[A-OQ-Z][ST])"
    residues = scaffold.find_glycosylation_sites(sequon)

    assert "A" in residues.keys()
    assert "B" not in residues.keys()

    residues = list(residues.values())
    assert len(residues[0]) == 20

    scaffold.include_chain("B")
    residues = scaffold.find_glycosylation_sites(sequon)

    assert "A" in residues.keys()
    assert "B" in residues.keys()

    residues = list(residues.values())
    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


def test_exclude_include_residue():
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
    scaffold.reindex()

    asns = scaffold.find_glycosylation_sites()
    asns = asns["A"][:10]

    for a in asns:
        scaffold.exclude_residue(a)

    residues = scaffold.find_glycosylation_sites()
    residues = list(residues.values())

    assert len(residues[0]) == 20 - 10  # we excluded 10 of the 20 ASNs
    assert len(residues[1]) == 4

    for a in asns:
        scaffold.include_residue(a)

    residues = scaffold.find_glycosylation_sites()
    residues = list(residues.values())

    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


def test_find_sequon():
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
    scaffold.reindex()

    sequon = "(N)(?=[A-OQ-Z][ST])"
    residues = scaffold.find_glycosylation_sites(sequon)
    residues = list(residues.values())

    assert len(residues[0]) == 20
    assert len(residues[1]) == 4


def test_attach_glycan_to_one_residue_simple():
    glycan = gl.glycan(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    # glycan.infer_bonds(restrict_residues=False)
    scaffold.apply_standard_bonds()

    assert len(glycan.residues) != 0
    assert len(scaffold.glycans) == 0

    scaffold.attach(
        glycan, residues=scaffold.get_residues("ASN", by="name", chain="A")[0]
    )

    assert len(scaffold.glycans) == 1


def test_attach_glycan_to_residue_list():
    glycan = gl.glycan(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    # glycan.infer_bonds(restrict_residues=False)
    scaffold.apply_standard_bonds()
    glycan.infer_glycan_tree()

    assert len(glycan.residues) != 0
    assert len(scaffold.glycans) == 0

    scaffold.attach(
        glycan, residues=scaffold.get_residues("ASN", by="name", chain="A")[:5]
    )
    assert len(scaffold.glycans) == 5
    scaffold.to_pdb("scaf.glycv2.pdb")


def test_attach_glycan_with_sequon():
    glycan = gl.Glycan.from_pdb(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    scaffold.attach(glycan, sequon="O-linked")

    assert len(scaffold.glycans) > 1


def test_attach_with_operator_inplace():
    glycan = gl.Glycan.from_pdb(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    link = gl.linkage("ND2", "C1", None, ["O1", "HO1"])
    scaffold % link @ scaffold.get_residues("ASN", by="name", chain="A")[0]
    scaffold += glycan

    assert len(scaffold.glycans) == 1


def test_attach_with_operator_inplace_twice():
    glycan = gl.Glycan.from_pdb(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    link = gl.linkage("ND2", "C1", None, ["O1", "HO1"])
    scaffold % link @ scaffold.get_residues("ASN", by="name", chain="A")[0]
    scaffold += glycan

    scaffold @ scaffold.get_residues("ASN", by="name", chain="A")[1]
    scaffold += glycan

    assert len(scaffold.glycans) == 2


def test_attach_with_operator_copy():
    glycan = gl.Glycan.from_pdb(base.MANNOSE9)
    scaffold = gl.Protein.from_pdb(base.PROTEIN)

    scaffold.reindex()
    glycan.reindex()
    glycan.infer_bonds(restrict_residues=False)
    scaffold.infer_bonds()

    assert len(glycan.residues) == 11

    scaffold @ scaffold.get_residues("ASN", by="name", chain="A")[0]
    out = scaffold + glycan

    assert out is not scaffold
    assert len(out.glycans) == 1
    assert len(scaffold.glycans) == 0


def test_hollow_out():
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
    scaffold.reindex()

    old_residues = len(scaffold.residues)

    # hollow out the protein
    scaffold.hollow_out()

    new_residues = len(scaffold.residues)

    assert new_residues < old_residues


def test_fill():
    scaffold = gl.Protein.from_pdb(base.PROTEIN)
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


def test_find_glycans():
    scaffold = gl.Protein.from_pdb(base.PDB_GLYCOSYLATED_PROTEIN)
    assert len(scaffold.glycans) == 0, "glycans found in scaffold without searching"

    glycans = scaffold.find_glycans()
    assert len(glycans) == 3, "glycans not found in scaffold"

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    for i, glycan in enumerate(glycans.values()):
        glycan.draw2d(ax=axs[i])
        axs[i].set_title(f"glycan {i+1}")
    plt.show()
