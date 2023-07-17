"""
Icons for drawing glycan structures in 2D
"""

import os
import matplotlib.pyplot as plt

DIR = os.path.dirname(__file__)

__loaded_icons__ = {}

__translation_table__ = {
    "GLC": "Glc",
    "NDG": "GlcNAc",
    "NGZ": "GlcNAc",
    "PA1": "GlcN",
    "GCU": "GlcA",
    "G6D": "Qui",
    "Z9W": "QuiNAc",
    "DDA": "Oli",
    "RAE": "Oli",
    "BAC": "Bac",
    "B6D": "Bac",
    "API": "Api",
    "XXM": "Api",
    "MAN": "Man",
    "BM3": "ManNAc",
    "95Z": "ManN",
    "MAV": "ManA",
    "RAM": "Rha",
    "XXR": "Rha",
    "RHN": "RhaNAc",
    "TYV": "Tyv",
    "64K": "Ara",
    "AHR": "Ara",
    "ARA": "Ara",
    "BXY": "Ara",
    "KDM": "Kdn",
    "PSE": "Pse",
    "GMH": "LDmanHep",
    "Z9N": "Fru",
    "GLA": "Gal",
    "GXL": "Gal",
    "A2G": "GalNAc",
    "YYQ": "GalNAc",
    "X6X": "GalN",
    "ADA": "GalA",
    "LDY": "Lyx",
    "LEG": "Leg",
    "KDO": "Kdo",
    "T6T": "Tag",
    "4GL": "Gul",
    "GUP": "Gul",
    "LXB": "GulNAc",
    "GUN": "GulN",
    "LGU": "GulA",
    "66O": "6dGul",
    "ABE": "Abe",
    "HSY": "Xyl",
    "XYS": "Xyl",
    "DHA": "Dha",
    "SOE": "Sor",
    "SHD": "Alt",
    "Z6H": "Alt",
    "ATN": "AltNAc",
    "AL1": "AltN",
    "AL2": "AltA",
    "6AT": "6dAlt",
    "6AN": "6dAltNAc",
    "PZU": "Par",
    "RIB": "Rib",
    "YYM": "Rib",
    "Z6J": "Rib",
    "ACI": "Aci",
    "289": "DDmanHep",
    "PSV": "Psi",
    "SF6": "Psi",
    "AFD": "All",
    "Z2D": "All",
    "NAA": "AllNAc",
    "AUA": "AllA",
    "A0N": "AllN",
    "Z3U": "Dig",
    "SIA": "Neu5Ac",
    "MUB": "MurNAc",
    "A5C": "Tal",
    "TLO": "TalNAc",
    "TNO": "TalN",
    "X1X": "TalA",
    "6TT": "6dTal",
    "6TN": "6dTalNAc",
    "COL": "Col",
    "NGC": "Neu5Gc",
    "4LG": "4eLeg",
    "4LC": "MurNGc",
    "Z0F": "Ido",
    "ZCD": "Ido",
    "HSQ": "IdoNAc",
    "LXZ": "IdoNAc",
    "IIO": "IdoN",
    "IDR": "IdoA",
    "NEX": "Neu",
    "MUR": "Mur",
    "FCA": "Fuc",
    "FUC": "Fuc",
    "UNF": "FucNAc",
}

__base_color_mapping__ = {
    "GLC": "blue",
    "NDG": "blue",
    "NGZ": "blue",
    "PA1": "blue",
    "GCU": "blue",
    "G6D": "blue",
    "Z9W": "blue",
    "DDA": "blue",
    "RAE": "blue",
    "BAC": "blue",
    "B6D": "blue",
    "API": "blue",
    "XXM": "blue",
    "MAN": "green",
    "BM3": "green",
    "95Z": "green",
    "MAV": "green",
    "RAM": "green",
    "XXR": "green",
    "RHN": "green",
    "TYV": "green",
    "64K": "green",
    "AHR": "green",
    "ARA": "green",
    "BXY": "green",
    "KDM": "green",
    "PSE": "green",
    "GMH": "green",
    "Z9N": "green",
    "GLA": "yellow",
    "GXL": "yellow",
    "A2G": "yellow",
    "YYQ": "yellow",
    "X6X": "yellow",
    "ADA": "yellow",
    "LDY": "yellow",
    "LEG": "yellow",
    "KDO": "yellow",
    "T6T": "yellow",
    "4GL": "orange",
    "GUP": "orange",
    "LXB": "orange",
    "GUN": "orange",
    "LGU": "orange",
    "66O": "orange",
    "ABE": "orange",
    "HSY": "orange",
    "XYS": "orange",
    "DHA": "orange",
    "SOE": "orange",
    "SHD": "pink",
    "Z6H": "pink",
    "ATN": "pink",
    "AL1": "pink",
    "AL2": "pink",
    "6AT": "pink",
    "6AN": "pink",
    "PZU": "pink",
    "RIB": "pink",
    "YYM": "pink",
    "Z6J": "pink",
    "ACI": "pink",
    "289": "pink",
    "PSV": "pink",
    "SF6": "pink",
    "AFD": "purple",
    "Z2D": "purple",
    "NAA": "purple",
    "AUA": "purple",
    "A0N": "purple",
    "Z3U": "purple",
    "SIA": "purple",
    "MUB": "purple",
    "A5C": "lightblue",
    "TLO": "lightblue",
    "TNO": "lightblue",
    "X1X": "lightblue",
    "6TT": "lightblue",
    "6TN": "lightblue",
    "COL": "lightblue",
    "NGC": "lightblue",
    "4LG": "lightblue",
    "4LC": "lightblue",
    "Z0F": "brown",
    "ZCD": "brown",
    "HSQ": "brown",
    "LXZ": "brown",
    "IIO": "brown",
    "IDR": "brown",
    "NEX": "brown",
    "MUR": "brown",
    "FCA": "red",
    "FUC": "red",
    "UNF": "red",
}


def get_icon_path(name: str) -> str:
    """
    Get the path to an icon PNG file

    Parameters
    ----------
    name : str

    Returns
    -------
    str
    """
    return os.path.join(DIR, name + ".png")


def get_icon(name: str) -> "ndarray":
    """
    Get an icon by name

    Parameters
    ----------
    name : str
        The name of the icon

    Returns
    -------
    "ndarray"
        The icon image
    """
    if name.startswith("b-"):
        name = name[2:]
    if not name in __translation_table__.values():
        name = __translation_table__.get(name, "unknown")
    if name not in __loaded_icons__:
        __loaded_icons__[name] = plt.imread(get_icon_path(name))
    return __loaded_icons__[name]


def get_base_color(name: str) -> str:
    """
    Get the only base color of a monosaccharide

    Parameters
    ----------
    name : str
        The name of the monosaccharide

    Returns
    -------
    str
        The base color of the monosaccharide
    """
    if name.startswith("b-"):
        name = name[2:]
    if not name in __translation_table__.values():
        name = __translation_table__.get(name, "unknown")
    return __base_color_mapping__.get(name, "gray")
