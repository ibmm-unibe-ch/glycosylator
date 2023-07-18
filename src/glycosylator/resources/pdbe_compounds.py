"""
The glycosylator PDBE compounds library with additional curated sugars
"""

import os
import biobuild.resources.pdbe_compounds as core

DIR = os.path.dirname(__file__)

GLYCOSYLATOR_COMPOUNDS_FILE = os.path.join(DIR, "compounds.json")
"""
The path to the glycosylator-curated compounds file (JSON)
"""

core.defaults.DEFAULT_PDBE_COMPONENT_FILES["glycosylator"] = GLYCOSYLATOR_COMPOUNDS_FILE


def load_glycosylator_compounds():
    """
    Load the glycosylator compounds library
    """
    compounds = core.PDBECompounds.from_json(GLYCOSYLATOR_COMPOUNDS_FILE)
    core.set_default_compounds(compounds)


def set_default_compounds(compounds: core.PDBECompounds, overwrite: bool = False):
    """
    Set the default compounds library

    Parameters
    ----------
    compounds : biobuild.PDBECompounds
        The compounds library
    overwrite : bool
        Whether to overwrite the existing default library
    """
    if overwrite:
        current = core.get_default_compounds()
    core.set_default_compounds(compounds, overwrite=False)
    if overwrite:
        current.to_json(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak")
        compounds.to_json(GLYCOSYLATOR_COMPOUNDS_FILE)


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
    set_default_compounds(
        core.PDBECompounds.from_json(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak")
    )
    if overwrite:
        os.rename(GLYCOSYLATOR_COMPOUNDS_FILE + ".bak", GLYCOSYLATOR_COMPOUNDS_FILE)


__all__ = [
    "set_default_compounds",
    "restore_default_compounds",
    "load_glycosylator_compounds",
]
