"""
This module contains functions for interacting with PubChem.
"""

import pubchempy as pcp

_grammar_dict = {
    1: "st",
    2: "nd",
    3: "rd",
}


def query(_query: str, by: str = "name", idx: int = 0):
    """
    Query PubChem for a compound.

    Parameters
    ----------
    _query : str
        The query string
    by : str, optional
        The type of query to perform. One of "name", "cid", "smiles", "inchi", "inchikey".
        Defaults to "name".
    idx : int, optional
        The index of the compound to return if multiple compounds are found. Defaults to 0.
        Only one compound will be returned at a time.

    Returns
    -------
    tuple
        The 2D and 3D representations of the compound
    """
    _2d_comp = pcp.get_compounds(_query, by)
    _3d_comp = pcp.get_compounds(_query, by, record_type="3d")
    if len(_2d_comp) == 0 or len(_3d_comp) == 0:
        raise ValueError(f"Could not find molecule {_query} in PubChem")
    elif len(_2d_comp) > 1 or len(_3d_comp) > 1:
        Warning(
            f"Found multiple molecules for {_query} in PubChem, using the {idx+1}{_grammar_dict.get(idx+1, 'th')} one!"
        )
    return _2d_comp[idx], _3d_comp[idx]
