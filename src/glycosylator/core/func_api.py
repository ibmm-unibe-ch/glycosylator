"""
Some basic functions to interact with the core of glycosylator.
"""
import os
import Bio.PDB as bio

import glycosylator.core.molecule as molecule
import glycosylator.core.scaffold as scaffold
import glycosylator.core.linkage as linkage
import glycosylator.resources as resources


def read_pdb(
    filename: str,
    kind: str = "molecule",
    model: int = 0,
    chain: str = None,
    root_atom: int = None,
):
    """
    Read a PDB file into a Bio.PDB.Structure

    Parameters
    ----------
    filename : str
        The filename of the PDB file to read
    kind : str
        The kind of object to make. Can be either "molecule" or "scaffold"
    model : int
        The model to read from the PDB file (if multiple models are present, only applies to "molecules")
    chain : str
        The chain to read from the PDB file (if multiple chains are present and only one is desired, only applies to "molecules")
    root_atom : int
        Set a specific atom as the root atom (in case a "molecule" is made). If None, the first atom is used.

    Returns
    -------
    object
        The object of type `kind` that was read from the PDB file
    """
    kind = kind.lower()
    if kind == "molecule":
        return molecule.Molecule.from_pdb(
            filename, root_atom=root_atom, model=model, chain=chain
        )
    elif kind == "scaffold":
        return scaffold.Scaffold.from_pdb(filename)
    else:
        raise ValueError(f"Unknown kind: {kind}")


def make_molecule(mol: str):
    """
    Generate a molecule from an input string. This string can be a PDB id, filename, SMILES or InChI string, IUPAC name or abbreviation.
    This function will try its best to automatically generate the molecule with minimal user effort. However, using a dedicated class method is
    recommended for more efficient and predictable results.

    Parameters
    ----------
    mol : str
        The input string
    """
    if isinstance(mol, bio.Structure.Structure):
        return molecule.Molecule(mol)

    if not isinstance(mol, str):
        raise ValueError("input must be a string")

    # ------------------
    # mol may be a PDB id
    # ------------------
    if len(mol) == 3:
        if resources.has_compound(mol):
            return resources.get_compound(mol)
    try:
        return molecule.Molecule.from_pubchem(mol)
    except:
        pass

    if os.path.isfile(mol):
        return molecule.Molecule.from_pdb(mol)

    try:
        return molecule.Molecule.from_smiles(mol)
    except:
        pass

    raise ValueError(f"Could not generate molecule from input: {mol}")


def polymerize(
    molecule: molecule.Molecule, n: int, patch_or_recipe=None, inplace: bool = False
) -> molecule.Molecule:
    """
    Polymerize a molecule

    Parameters
    ----------
    molecule : Molecule
        The molecule to polymerize
    n : int
        The number of monomers to add
    patch_or_recipe : str or Recipe
        The patch or recipe to use for polymerization. If None, the default patch is used.
    inplace : bool
        Whether to polymerize the molecule in place or return a new molecule

    Returns
    -------
    Molecule
        The polymerized molecule
    """
    if patch_or_recipe is None and molecule._patch is None:
        raise ValueError(
            "No patch or recipe provided and no default is set on the molecule"
        )
    return molecule.repeat(n, patch_or_recipe, inplace=inplace)
