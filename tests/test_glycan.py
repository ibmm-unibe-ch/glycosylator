import base as base
import glycosylator as gl


def test_extend_simple_linear():
    iupac_partial = "Gal(b1-4)GlcNAc"
    glycan = gl.glycan(iupac_partial)

    iupac_full = "Gal(b1-4)GlcNAc(b1-3)" + iupac_partial
    glycan.extend_to(iupac_full)

    glycan.show3d()
    assert glycan.to_iupac(False) == iupac_full


def test_extend_branched():
    iupac_partial = "Gal(b1-4)GlcNAc"
    glycan = gl.glycan(iupac_partial)

    iupac_full = "Glc(b1-4)Glc(b1-4)[Gal(b1-4)]GlcNAc"
    try:
        full = gl.glycan(iupac_full)
    except Exception as e:
        raise ValueError("The full glycan is not valid")
    glycan.extend_to(iupac_full)

    glycan.show3d(residue_graph=True)


def test_extend_branched2():
    full_iupac = "Gal(b1-3)[Neu5Ac(a2-6)[Neu5Ac(a2-8)Neu5Ac(a2-3)[GalNAc(b1-4)]Glc(b1-4)]Gal(b1-4)GlcNAc(b1-6)]GalNAc"
    partial_iupac = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc"
    glycan = gl.glycan(partial_iupac)
    glycan.extend_to(full_iupac)
    glycan.optimize()

    ref_glycan = gl.glycan(full_iupac)
    ref_glycan.optimize()

    res_graph = glycan.make_residue_graph()
    ref_res_graph = ref_glycan.make_residue_graph()
    assert len(res_graph.nodes) == len(ref_res_graph.nodes)
    assert len(res_graph.edges) == len(ref_res_graph.edges)

    assert glycan.to_biopython() == ref_glycan.to_biopython()

    graph1 = glycan.make_residue_graph()
    graph2 = ref_glycan.make_residue_graph()

    verdicts = []
    for residue in glycan.get_residues():
        matches = [res for res in ref_glycan.get_residues() if res.matches(residue)]
        found_match = False
        for match in matches:
            neighs1 = graph1.get_neighbors(residue)
            neighs2 = graph2.get_neighbors(match)
            if len(neighs1) == len(neighs2):
                verdicts.append(True)
                found_match = True
                break
        if not found_match:
            verdicts.append(False)

    assert all(verdicts)

    reconstruct_full = gl.Glycan.from_iupac(None, glycan.to_iupac(False))
    graph3 = reconstruct_full.make_residue_graph()

    graph_hists = {}
    for graph in [graph1, graph2, graph3]:
        hist = {}
        for node in graph.nodes:
            n = len(graph.get_neighbors(node))
            hist[n] = hist.get(n, 0) + 1
        graph_hists[graph] = hist

    assert graph_hists[graph1] == graph_hists[graph2] == graph_hists[graph3]

    glycan.to_pdb("glycan.pdb")
    ref_glycan.to_pdb("ref_glycan.pdb")

    # v = gl.utils.visual.MoleculeViewer3D()
    # v.draw_edges(*res_graph.edges, color="black")
    # v.draw_edges(*ref_res_graph.edges, color="red")
    # v.draw_edges(*reconstruct_full.make_residue_graph().edges, color="blue")
    # v.show()


def test_extend_big():
    iupac = "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
    glycan = gl.glycan("b-GlcNAc")
    glycan.extend_to(iupac)
    glycan.to_pdb("glycan_big.pdb")
    
def test_draw2d():
    iupac = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc"
    glycan = gl.glycan(iupac)

    glycan.show2d()
