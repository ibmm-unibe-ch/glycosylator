"""
Pre-defined sequons for protein glycosylation sites
"""

import os
import json

DIR = os.path.dirname(__file__)

GLYCOSYLATOR_SEQUONS_FILE = os.path.join(DIR, "sequons.json")
"""
File for default sequons for protein glycosylation sites
"""


def get_sequon(key: str):
    """
    Get a sequon regex pattern

    Parameters
    ----------
    key : str
        Sequon key

    Returns
    -------
    str
        Sequon regex pattern, or None if not found
    """
    return SEQUONS.get(key, None)


def set_sequon(key: str, pattern: str, overwrite: bool = False):
    """
    Set a sequon regex pattern

    Parameters
    ----------
    key : str
        Sequon key
    pattern : str
        Sequon regex pattern
    overwrite : bool
        If True, the new sequon is permanently included in the default sequons file.
    """
    if overwrite:
        save_sequons(GLYCOSYLATOR_SEQUONS_FILE + ".bak")
    SEQUONS[key] = pattern
    if overwrite:
        save_sequons(GLYCOSYLATOR_SEQUONS_FILE)


def load_sequons(filename: str) -> dict:
    """
    Load a json file with sequon regex patterns

    Parameters
    ----------
    filename : str
        Path to the json file

    Returns
    -------
    dict
        Sequon regex patterns
    """
    with open(filename, "r") as f:
        sequons = json.load(f)
    return sequons


def save_sequons(filename: str):
    """
    Save the current sequon regex patterns to a json file

    Parameters
    ----------
    filename : str
        Path to the json file
    """
    with open(filename, "w") as f:
        json.dump(SEQUONS, f, indent=4)


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
    SEQUONS = load_sequons(GLYCOSYLATOR_SEQUONS_FILE + ".bak")
    if overwrite:
        os.rename(GLYCOSYLATOR_SEQUONS_FILE + ".bak", GLYCOSYLATOR_SEQUONS_FILE)


SEQUONS = load_sequons(GLYCOSYLATOR_SEQUONS_FILE)
"""
Default sequons for protein glycosylation sites
"""

__all__ = [
    "SEQUONS",
    "get_sequon",
    "set_sequon",
    "load_sequons",
    "save_sequons",
    "restore_sequons",
]
