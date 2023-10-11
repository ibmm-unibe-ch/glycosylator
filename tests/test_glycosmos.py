import base as base
import glycosylator as gl

REF_ID = "G78791QP"
REF_IUPAC = "Gal(b1-3)[Fuc(a1-2)]Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(a1-"


def test_get_id_from_iupac():
    id = gl.resources.get_glytoucan_id_from_iupac(REF_IUPAC)
    assert id == REF_ID


def test_get_iupac_from_id():
    iupac = gl.resources.get_iupac_from_glycosmos(REF_ID)
    assert iupac == REF_IUPAC


def test_make_glycan_from_id():
    g = gl.Glycan.from_glycosmos(REF_ID)
    assert g.to_iupac() == REF_IUPAC
    assert g.get_glytoucan_id() == REF_ID
    g.show3d()
    g.show2d()


def test_find_matching_glycans():
    ids = gl.resources.find_glytoucan_ids_from_iupac(REF_IUPAC)
    assert len(ids) > 0
    for i, seq in ids:
        assert REF_IUPAC in seq
