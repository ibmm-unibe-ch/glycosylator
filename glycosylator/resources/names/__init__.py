"""
These are name mappings for glycan residues between their PDB identifiers and IUPAC names.
"""

__pdb_id_to_iupac_name_mapping__ = {
    "TNO": "TalN",
    "KDO": "Kdo",
    "KDM": "Kdn",
    "KDN": "b-Kdn",
    "NDG": "GlcNAc",
    "NGZ": "GlcNAc",
    "NAG": "b-GlcNAc",
    "MUB": "MurNAc",
    "AMU": "b-MurNAc",
    "PSE": "Pse",
    "Z0F": "Ido",
    "ZCD": "Ido",
    "4N2": "b-Ido",
    "AL1": "AltN",
    "6AN": "6dAltNAc",
    "ADA": "GalA",
    "GTK": "b-GalA",
    "GTR": "b-GalA",
    "TLO": "TalNAc",
    "66O": "6dGul",
    "4LG": "4eLeg",
    "A0N": "AllN",
    "QUN": "QuiNAc",
    "Z9W": "QuiNAc",
    "4LC": "MurNGc",
    "AFD": "All",
    "Z2D": "All",
    "ALL": "b-All",
    "WOO": "b-All",
    "BM3": "ManNAc",
    "BM7": "b-ManNAc",
    "GLA": "Gal",
    "GXL": "Gal",
    "GAL": "b-Gal",
    "GIV": "b-Gal",
    "GZL": "b-Gal",
    "DHA": "Dha",
    "Z3U": "Dig",
    "HSY": "Xyl",
    "XYS": "Xyl",
    "LXC": "b-Xyl",
    "XYP": "b-Xyl",
    "XYZ": "b-Xyl",
    "T6T": "Tag",
    "LXB": "b-GulNAc",
    "TYV": "Tyv",
    "95Z": "ManN",
    "ABE": "Abe",
    "GCU": "GlcA",
    "BDP": "b-GlcA",
    "MUR": "b-Mur",
    "1S4": "b-Mur",
    "A2G": "GalNAc",
    "YYQ": "GalNAc",
    "NGA": "b-GalNAc",
    "IIO": "IdoN",
    "PZU": "Par",
    "4GL": "Gul",
    "GUP": "Gul",
    "GL0": "b-Gul",
    "Z8H": "b-Gul",
    "GUN": "GulN",
    "PA1": "GlcN",
    "GCS": "b-GlcN",
    "GLC": "Glc",
    "BGC": "b-Glc",
    "Z8T": "b-Glc",
    "IDR": "IdoA",
    "LEG": "Leg",
    "289": "DDmanHep",
    "A5C": "Tal",
    "SDY": "b-Tal",
    "ZEE": "b-Tal",
    "MAN": "Man",
    "BMA": "b-Man",
    "DDA": "Oli",
    "RAE": "Oli",
    "Z5J": "b-Oli",
    "ACI": "Aci",
    "LGU": "GulA",
    "SHD": "Alt",
    "Z6H": "Alt",
    "3MK": "b-Alt",
    "SIA": "Neu5Ac",
    "SLB": "b-Neu5Ac",
    "GMH": "LDmanHep",
    "NEX": "Neu",
    "Z9N": "Fru",
    "BDF": "b-Fru",
    "FRU": "b-Fru",
    "LFR": "b-Fru",
    "COL": "Col",
    "MAV": "ManA",
    "BEM": "b-ManA",
    "6TT": "6dTal",
    "BAC": "Bac",
    "B6D": "Bac",
    "X6X": "GalN",
    "1GN": "b-GalN",
    "RIB": "Rib",
    "YYM": "Rib",
    "Z6J": "Rib",
    "0MK": "b-Rib",
    "32O": "b-Rib",
    "BDR": "b-Rib",
    "RIP": "b-Rib",
    "API": "Api",
    "XXM": "b-Api",
    "64K": "Ara",
    "AHR": "Ara",
    "ARA": "Ara",
    "BXY": "Ara",
    "ARB": "b-Ara",
    "BXX": "b-Ara",
    "FUB": "b-Ara",
    "SEJ": "b-Ara",
    "FCA": "Fuc",
    "FUC": "Fuc",
    "FCB": "b-Fuc",
    "FUL": "b-Fuc",
    "GYE": "b-Fuc",
    "ATN": "AltNAc",
    "AUA": "AllA",
    "NAA": "b-AllNAc",
    "6TN": "6dTalNAc",
    "SOE": "Sor",
    "UEA": "b-Sor",
    "G6D": "Qui",
    "YYK": "b-Qui",
    "X1X": "TalA",
    "X0X": "b-TalA",
    "PSV": "Psi",
    "SF6": "Psi",
    "SF9": "b-Psi",
    "TTV": "b-Psi",
    "RAM": "Rha",
    "XXR": "Rha",
    "RM4": "b-Rha",
    "6AT": "6dAlt",
    "NGC": "Neu5Gc",
    "NGE": "b-Neu5Gc",
    "LDY": "Lyx",
    "Z4W": "b-Lyx",
    "HSQ": "IdoNAc",
    "LXZ": "IdoNAc",
    "UNF": "FucNAc",
    "49T": "b-FucNAc",
    "ATA": "AltA",
}

# just in case someone external has already incorporated this
__pdb_id_to_name_mapping__ = __pdb_id_to_iupac_name_mapping__
"""
WARNING: This is deprecated, use `__pdb_id_to_iupac_name_mapping__` instead.
"""

__iupac_name_to_pdb_id_mapping__ = {
    value: key for key, value in __pdb_id_to_iupac_name_mapping__.items()
}

__beta_compounds_pdb_ids__ = set(
    [
        key
        for key, value in __pdb_id_to_iupac_name_mapping__.items()
        if value.startswith("b-")
    ]
)
"""
A set of PDB ids for monosaccharides that are beta-configured.
"""


def id_to_name(id: str) -> str:
    """
    Translate a monosaccharide PDB id to a name

    Parameters
    ----------
    id: str
        The PDB id of the monosaccharide

    Returns
    -------
    str
        The name of the monosaccharide.
        If the id is not found, it is returned as is.
    """
    return __pdb_id_to_iupac_name_mapping__.get(id, id)


def is_beta(id: str) -> bool:
    """
    Check if a monosaccharide PDB id is beta-configured.

    Parameters
    ----------
    id: str
        The PDB id of the monosaccharide

    Returns
    -------
    bool
        True if the monosaccharide is beta-configured, False otherwise.
    """
    return id in __beta_compounds_pdb_ids__
