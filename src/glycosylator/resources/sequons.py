"""
Pre-defined sequons for protein glycosylation sites
"""

import json

SEQUONS = {
    "N-linked": "(N)(?=[A-OQ-Z][ST])",
    "O-linked": "(S|T)(?=[A-OQ-Z][ST])",
}
"""
Default sequons for protein glycosylation sites
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


def set_sequon(key: str, pattern: str):
    """
    Set a sequon regex pattern

    Parameters
    ----------
    key : str
        Sequon key
    pattern : str
        Sequon regex pattern
    """
    SEQUONS[key] = pattern


def load_sequons(filename: str):
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
        SEQUONS.update(json.load(f))
    return SEQUONS


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
