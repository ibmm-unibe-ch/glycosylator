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
    assert res_graph.get_neighbors(r).pop().resname == "NAG"


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


def test_from_iupac3():
    s = "Glc(b1-4)Glc(b1-3)[Gal(b1-4)]GlcNAc"
    glycan = gl.Glycan.from_iupac(None, s)

    res_graph = glycan.make_residue_graph()

    assert len(res_graph.nodes) == 4
    r = glycan.get_residue(2)
    assert len(res_graph.get_neighbors(r)) == 1
    r = glycan.get_residue(3)
    assert len(res_graph.get_neighbors(r)) == 2


def test_from_iupac_with_pdbids():
    s = "GLC(b1-4)BGC(b1-3)[Gal(b1-4)]GlcNAc"
    glycan = gl.Glycan.from_iupac(None, s)

    s2 = "Glc(b1-4)BGC(b1-3)[Gal(b1-4)]GlcNAc"
    glycan2 = gl.Glycan.from_iupac(None, s2)

    res_graph = glycan.make_residue_graph()

    assert len(res_graph.nodes) == 4
    r = glycan.get_residue(2)
    assert len(res_graph.get_neighbors(r)) == 1
    r = glycan.get_residue(3)
    assert len(res_graph.get_neighbors(r)) == 2

    assert glycan.get_residue(-1).resname == "GLC"
    assert glycan.get_residue(-2).resname == "BGC"
    assert glycan2.get_residue(-1).resname == "BGC"


def test_from_iupac_with_DL():
    s = "Gal(b1-4)GlcNAc(b1-3)[Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Gal(b1-4)Glc"
    glycan = gl.Glycan.from_iupac(None, s)

    s2 = "l-Galp(b1-4)GlcNAc(b1-3)[l-Galp(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Gal(b1-4)Glc"
    glycan2 = gl.Glycan.from_iupac(None, s2)

    assert glycan.get_residue(-1).resname != glycan2.get_residue(-1).resname
    assert glycan.get_residue(-3).resname != glycan2.get_residue(-3).resname


def test_from_iupac_with_other_linkages():
    my_link = gl.linkage("C1", "C1", ["O1", "HO1"], ["O1", "HO1"], id="MYLINK")
    gl.add_linkage(my_link)

    s = "Glc(%MYLINK)Glc"
    _s = gl.utils.IUPACParser()(s)

    glycan = gl.Glycan.from_iupac(None, s)

    assert glycan.get_residue_connections(triplet=False)[0].atom1.id == "C1"
    assert glycan.get_residue_connections(triplet=False)[0].atom2.id == "C1"


def test_to_iupac():
    s = "Fuc(a1-2)Gal(b1-3)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)Glc"
    glycan = gl.glycan(s)

    out = glycan.to_iupac(add_terminal_conformation=False)
    assert out == s

    out = glycan.to_iupac(add_terminal_conformation=True)
    assert out == s + "(a1-"


def test_to_iupac2():
    s = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)Glc(b1-"
    _s = gl.utils.IUPACParser()(s)
    glycan = gl.glycan(s)

    out = glycan.to_iupac(add_terminal_conformation=False)
    assert out == s[:-7] + "b-Glc"

    out = glycan.to_iupac(add_terminal_conformation=True)
    assert out == s
