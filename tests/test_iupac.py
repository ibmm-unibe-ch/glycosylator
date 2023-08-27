import base as base
import glycosylator as gl


def test_from_iupac():
    s = "Gal(b1-4)GlcNAc(b1-3)[Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Gal(b1-4)Glc"
    glycan = gl.glycan(s)

    res_graph = glycan.make_residue_graph()

    assert len(res_graph.nodes) == 7
    r = glycan.get_residue(2)
    assert len(res_graph.get_neighbors(r)) == 3
    r = glycan.get_residue("FCA")
    assert len(res_graph.get_neighbors(r)) == 1
    assert res_graph.get_neighbors(r).pop().id == "NAG"


def test_from_iupac2():
    s = "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    glycan = gl.glycan(s)

    res_graph = glycan.make_residue_graph()

    assert len(res_graph.nodes) == 7
    r = glycan.get_residue(2)
    assert len(res_graph.get_neighbors(r)) == 2
    r = glycan.get_residue(3)
    assert len(res_graph.get_neighbors(r)) == 3

    r = glycan.get_residues("BMA")
    assert len(r) == 1
    r = glycan.get_residues("MAN")
    assert len(r) == 4


def test_to_iupac():
    s = "Fuc(a1-2)Gal(b1-3)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)Glc"
    glycan = gl.glycan(s)

    out = glycan.to_iupac(add_terminal_conformation=False)
    assert out == s

    out = glycan.to_iupac(add_terminal_conformation=True)
    assert out == s + "(a1-"


def test_to_iupac2():
    s = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)Glc(b1-"
    glycan = gl.glycan(s)

    out = glycan.to_iupac(add_terminal_conformation=False)
    assert out == s[:-7] + "b-Glc"

    out = glycan.to_iupac(add_terminal_conformation=True)
    assert out == s
