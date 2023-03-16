import os
import pickle
import Bio.PDB as bio

import glycosylator.utils.constants as constants
import glycosylator.force_fields.charmm as charmm

# =================================================================
# Default settings of values
# =================================================================

DEFAULT_BOND_LENGTH = 1.6
"""
The default length of a bond in Angstrom
"""

DEFAULT_CHARMM_TOPOLOGY_FILE = os.path.join(constants.RESOURCES, "CHARMM.top.pkl")
"""
The path to the default CHARMM topology file
"""

DEFAULT_CHARMM_PARAMETERS_FILE = os.path.join(constants.RESOURCES, "CHARMM.prm.pkl")
"""
The path to the default CHARMM parameters file
"""


# =================================================================
# Default instances of auxiliary classes
# =================================================================


__default_instances__ = dict(
    bioPDBParser=bio.PDBParser(),
    Topology=charmm.CHARMMTopology.load(DEFAULT_CHARMM_TOPOLOGY_FILE),
    Parameters=charmm.CHARMMParameters.load(DEFAULT_CHARMM_PARAMETERS_FILE),
)
"""
Default instance dictionary
"""

__bioPDBParser__ = __default_instances__["bioPDBParser"]
"""
The default instance of Bio.PDB.PDBParser
"""

# =================================================================
# Auxiliary functions
# =================================================================


def get_default_instance(key):
    """
    Get the default instance of a class

    Parameters
    ----------
    key : str
        The key of the default instance

    Returns
    -------
    obj
        The default instance of the class
    """
    return __default_instances__[key]


def set_default_instance(key, obj):
    """
    Set the default instance of a class

    Parameters
    ----------
    key : str
        The key of the default instance
    obj
        The new default instance
    """
    __default_instances__[key] = obj


def set_default_topology(obj, overwrite: bool = False):
    """
    Set a CHARMMTopology object as the new default.

    Parameters
    ----------
    obj : CHARMMTopology
        The CHARMMTopology object to set as the new default
    overwrite : bool
        If set to `True`, the new object will be permanently saved as
        the default. Otherwise, the new object will only be used for
        the current session.
    """
    if not isinstance(obj, charmm.CHARMMTopology):
        raise TypeError("The object must be a CHARMMTopology instance.")
    __default_instances__["Topology"] = obj
    if overwrite:
        obj.save(DEFAULT_CHARMM_TOPOLOGY_FILE)


def set_default_parameters(obj, overwrite: bool = False):
    """
    Set a CHARMMParameters object as the new default.

    Parameters
    ----------
    obj : CHARMMParameters
        The CHARMMParameters object to set as the new default
    overwrite : bool
        If set to `True`, the new object will be permanently saved as
        the default. Otherwise, the new object will only be used for
        the current session.
    """
    if not isinstance(obj, charmm.CHARMMParameters):
        raise TypeError("The object must be a CHARMMParameters instance.")
    __default_instances__["Parameters"] = obj
    if overwrite:
        obj.save(DEFAULT_CHARMM_PARAMETERS_FILE)


def get_default_parameters():
    """
    Get the default CHARMMParameters object

    Returns
    -------
    CHARMMParameters
        The default CHARMMParameters object
    """
    return __default_instances__["Parameters"]


def get_default_topology():
    """
    Get the default CHARMMTopology object

    Returns
    -------
    CHARMMTopology
        The default CHARMMTopology object
    """
    return __default_instances__["Topology"]
