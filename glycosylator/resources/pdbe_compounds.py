"""
The glycosylator PDBE compounds library with additional curated sugars
"""

import os
from buildamol.resources.pdbe_compounds import *
import buildamol.resources.pdbe_compounds as core

_DIR = os.path.dirname(__file__)

GLYCOSYLATOR_COMPOUNDS_FILE = os.path.join(_DIR, "compounds.pkl")
"""
The path to the glycosylator-curated compounds file (pickle)
"""

core.defaults.DEFAULT_PDBE_COMPONENT_FILES["glycosylator"] = GLYCOSYLATOR_COMPOUNDS_FILE


def load_glycosylator_compounds():
    """
    Load the glycosylator compounds library
    """
    compounds = PDBECompounds.load(GLYCOSYLATOR_COMPOUNDS_FILE)
    core.set_default_compounds(compounds)


def set_default_compounds(compounds: PDBECompounds, overwrite: bool = False):
    """
    Set the default compounds library

    Parameters
    ----------
    compounds : buildamol.PDBECompounds
        The compounds library
    overwrite : bool
        Whether to overwrite the existing default library
    """
    if overwrite:
        current = core.get_default_compounds()
    core.set_default_compounds(compounds, overwrite=False)
    if overwrite:
        current.save(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak")
        compounds.save(GLYCOSYLATOR_COMPOUNDS_FILE)


def restore_default_compounds(overwrite: bool = True):
    """
    Restore the default compounds library

    Parameters
    ----------
    overwrite : bool
        If True, the backup is permanently set as the default again.
    """
    if not os.path.isfile(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak"):
        raise FileNotFoundError("No backup file found")
    set_default_compounds(PDBECompounds.load(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak"))
    if overwrite:
        os.rename(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak", GLYCOSYLATOR_COMPOUNDS_FILE)


def reference_glycan_residue_ids(acceptable_types: list = None):
    """
    A set of compound ids (residue names), which are considered to be glycans.

    Note
    ----
    This assumes that the compounds and their residue share the same name as is the default case for the PDBE compounds library.

    Parameters
    ----------
    acceptable_types : list
        A list of acceptable compound types. If None, "SACCHARIDE" is used.
    """
    if acceptable_types is None:
        acceptable_types = ["SACCHARIDE"]
    compounds = core.get_default_compounds()
    ids = [
        c for c in compounds.ids if compounds._compounds[c]["type"] in acceptable_types
    ]
    return set(ids)


__all__ = [
    "load_glycosylator_compounds",
    "reference_glycan_residue_ids",
] + core.__all__
