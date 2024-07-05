"""
Pre-defined sequons for protein glycosylation sites
"""

from typing import Union
from buildamol.core import Linkage
import glycosylator.resources as resources

import os
import json

DIR = os.path.dirname(__file__)

GLYCOSYLATOR_SEQUONS_FILE = os.path.join(DIR, "sequons.json")
"""
File for default sequons for protein glycosylation sites
"""


class Sequon:
    """
    A class representing a sequon pattern for glycosylation
    """

    def __init__(self, id: str, pattern: str):
        self.id = id
        self.pattern = pattern
        self.linkages = {}

    def add_linkage(self, resname: str, linkage: str):
        """
        Add a linkage id that is associated with a given residue

        Parameters
        ----------
        resname : str
            Residue name
        linkage : str
            Linkage id
        """
        self.linkages[resname] = linkage

    def get_linkage(self, resname: str, _topology=None) -> "Linkage":
        """
        Get the linkage id associated with a given residue

        Parameters
        ----------
        resname : str
            Residue name
        _topology : Topology
            A particular topology to use to get the linkage from.

        Returns
        -------
        Linkage
            The Linkage object
        """
        if _topology is None:
            _topology = resources.get_default_topology()
        l = self.linkages[resname]
        return _topology.get_patch(l)

    def to_dict(self):
        """
        Convert the sequon to a dictionary

        Returns
        -------
        dict
            Sequon dictionary
        """
        return {
            "id": self.id,
            "pattern": self.pattern,
            "linkages": self.linkages,
        }

    @classmethod
    def _from_dict(cls, _dict):
        """
        Create a sequon from a dictionary

        Parameters
        ----------
        _dict : dict
            Sequon dictionary

        Returns
        -------
        Sequon
            Sequon
        """
        sequon = cls(_dict["id"], _dict["pattern"])
        sequon.linkages = _dict["linkages"]
        return sequon


def load_sequons(filename: str) -> dict:
    """
    Load a json file with sequons

    Parameters
    ----------
    filename : str
        Path to the json file

    Returns
    -------
    dict
        Sequons
    """
    with open(filename, "r") as f:
        sequons = json.load(f)

    _sequons = {}
    for key, val in sequons.items():
        _sequons[key] = Sequon._from_dict(val)
    return _sequons


SEQUONS = load_sequons(GLYCOSYLATOR_SEQUONS_FILE)
"""
Default sequons for protein glycosylation sites
"""


def save_sequons(filename: str):
    """
    Save the current sequon regex patterns to a json file

    Parameters
    ----------
    filename : str
        Path to the json file
    """
    with open(filename, "w") as f:
        sequons = {k: v.to_dict() for k, v in SEQUONS.items()}
        json.dump(sequons, f, indent=4)


def get_sequon(key: str) -> "Sequon":
    """
    Get a sequon regex pattern

    Parameters
    ----------
    key : str
        Sequon key

    Returns
    -------
    Sequon or None
    """
    return SEQUONS.get(key, None)


def add_sequon(
    sequon_or_pattern: Union[str, Sequon], key: str = None, overwrite: bool = False
):
    """
    Add a sequon to the default sequons

    Parameters
    ----------
    sequon_or_pattern : str
        Sequon object or regex pattern
    key : str
        Sequon key, if a regex pattern is provided.
    overwrite : bool
        If True, the new sequon is permanently included in the default sequons file.
    """
    if isinstance(sequon_or_pattern, str):
        sequon_or_pattern = Sequon(key, sequon_or_pattern)
    elif key is None:
        key = sequon_or_pattern.id

    if overwrite:
        save_sequons(GLYCOSYLATOR_SEQUONS_FILE + ".bak")

    SEQUONS[key] = sequon_or_pattern
    if overwrite:
        save_sequons(GLYCOSYLATOR_SEQUONS_FILE)


def restore_sequons(overwrite: bool = True):
    """
    Restore the default sequon regex patterns

    Parameters
    ----------
    overwrite : bool
        If True, the backup is permanently set as the default again.
    """
    if not os.path.isfile(GLYCOSYLATOR_SEQUONS_FILE + ".bak"):
        raise FileNotFoundError("No backup file found")
    SEQUONS.update(load_sequons(GLYCOSYLATOR_SEQUONS_FILE + ".bak"))
    if overwrite:
        os.rename(GLYCOSYLATOR_SEQUONS_FILE + ".bak", GLYCOSYLATOR_SEQUONS_FILE)


def set_default_sequons(sequons: dict):
    """
    Set the default sequons

    Parameters
    ----------
    sequons : dict
        Sequons
    """
    global SEQUONS
    SEQUONS = sequons


def get_default_sequons() -> dict:
    """
    Get the default sequons

    Returns
    -------
    dict
        Sequons
    """
    return SEQUONS


def available_sequons() -> list:
    """
    Get the available sequons

    Returns
    -------
    list
        Sequons
    """
    return list(SEQUONS.keys())


__all__ = [
    "SEQUONS",
    "get_sequon",
    "add_sequon",
    "load_sequons",
    "save_sequons",
    "restore_sequons",
    "Sequon",
    "set_default_sequons",
    "get_default_sequons",
    "available_sequons",
]
