"""
Some basic functions to interact with the core of glycosylator.
"""
import os
from typing import Union
import Bio.PDB as bio

import glycosylator.structural as structural
import glycosylator.utils as utils
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


def glycan(id: str, g: Union[str, list], _topology=None):
    """
    Make a molecule from an IUPAC glycan string in condensed format or a glycan structure graph.

    Parameters
    ----------
    id : str
        The id of the molecule to create
    g : str or list
        The glycan string or graph.
        If a string, it must be in IUPAC condensed format - currently, neither extended nor short formats are supported.
        If a list, it must be a list of tuples of the form (parent, child, linkage), where parent and child are strings
        representing the residue names and an "@{id}" suffix to distinguish individual residues. The linkage must be
        a valid id of any defined linkage in the provided or default topology, e.g. "14bb" or "16ab" for the default
        CHARMM topology.
    _topology : dict
        The topology to use. If None, the default topology is used.

    Returns
    -------
    molecule : Molecule
        The created molecule

    Examples
    --------
    To generate a small glycan of the structure:
    ```
    ~ --- NAG                  MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN
    ```
    the IUPAC string would be:
    >>> iupac = "Man(a1-6)[Man(a1-3)]Man(b1-4)Nag(b1-4)Nag(b1-" # notice the final "b1-" to indicate where the glycan attaches to a scaffold

    Which can be parsed into a molecule with:
    >>> mol = glycan("my_glycan", iupac)

    Alternatively, the same molecule can be generated from a graph:
    >>> graph = [
    ("NAG@1", "NAG@2", "14bb"),
    ("NAG@2", "BMA@1", "14bb"),
    ("BMA@1", "MAN@1", "13ab"),
    ("BMA@1", "MAN@2", "16ab"),
    ] # notice the @{id} after each residue name 
    >>> mol = glycan("my_glycan", graph)

    The `@{id}` suffix is used to distinguish between residues with the same name, for example the two mannoses (`MAN@1` and `MAN@2`).
    The in the example above the ids reflect the number of residues with the same name, hence `NAG@2` connects to `BMA@1` "the second NAG connecting to the first BMA".
    However, this is not a strict requirement. Any numeric or string value that will mark each residue as a unique node is acceptable.
    Hence, also the following graph is valid where the index simply reflects the order of the residues in the molecule: 
    
    >>> graph = [
    ("NAG@1", "NAG@2", "14bb"),
    ("NAG@2", "BMA@3", "14bb"),
    ("BMA@3", "MAN@4", "13ab"),
    ("BMA@3", "MAN@5", "16ab"),
    ]
    >>> mol = glycan("my_glycan", graph)
    """
    if isinstance(g, str):
        g = structural.IUPACParser().parse(g)
    elif not isinstance(g, list):
        raise ValueError("g must be either a string or a list")
    if not isinstance(g[0], (tuple, list)) or len(g[0]) != 3:
        raise ValueError("g must be a list of tuples of length 3")
    mol = _parse_iupac_graph(id, g, _topology)
    return mol


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


def make_scaffold(scaf: str):
    """
    Generate a scaffold from an input string. This string can be a PDB id, filename, a SMILES or InChI string, or a IUPAC name or abbreviation.
    """
    if isinstance(scaf, bio.Structure.Structure):
        return scaffold.Scaffold(scaf)

    if not isinstance(scaf, str):
        raise ValueError("input must be a string")

    try:
        mol = make_molecule(scaf)
        return scaffold.Scaffold(mol.structure)
    except:
        raise ValueError(f"Could not generate scaffold from input: {scaf}")


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


def _parse_iupac_graph(id, glycan_segments, _topology=None):
    """
    Make a molecule from a list of glycan segments that were generated by the IUPACParser class

    Parameters
    ----------
    id : str
        The id of the molecule
    glycan_segments : list
        A list of glycan segments

    Returns
    -------
    Molecule
        The molecule
    """
    if not _topology:
        _topology = utils.defaults.get_default_topology()

    # Check that all segments have a known patch
    for segment in glycan_segments:
        link = segment[-1]
        if not _topology.has_patch(link):
            raise ValueError(
                f"No patch/recipe available for linkage: {link}. Try adding a patch or recipe to the topology."
            )

    mol = None
    first_mol = None
    second_mol = None
    at_residue = None
    other_residue = None
    residue_id_mapping = {}
    for i, segment in enumerate(glycan_segments):
        first, second, link = segment

        first_name = first.split("@")[0]
        second_name = second.split("@")[0]

        if first in residue_id_mapping:
            at_residue = residue_id_mapping[first]
            first_mol = mol
        else:
            first_mol = make_molecule(first_name)
            residue_id_mapping[first] = len(residue_id_mapping) + 1
            at_residue = None

        if second in residue_id_mapping:
            other_residue = residue_id_mapping[second]
            second_mol = mol
        else:
            second_mol = make_molecule(second_name)
            residue_id_mapping[second] = len(residue_id_mapping) + 1
            other_residue = None

        if not mol:
            mol = first_mol
        first_mol.attach(
            second_mol,
            link,
            at_residue=at_residue,
            other_residue=other_residue,
            _topology=_topology,
        )

    mol.id = id
    return mol


__all__ = [
    "read_pdb",
    "make_molecule",
    "polymerize",
    "make_scaffold",
    "glycan",
]

if __name__ == "__main__":
    # glycans = [
    #     ("NAG@1", "NAG@2", "14bb"),
    #     ("NAG@2", "BMA@1", "14bb"),
    #     ("BMA@1", "MAN@1", "13ab"),
    #     ("BMA@1", "MAN@2", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    # glycans = [
    #     ("NAG@1", "NAG@2", "14bb"),
    #     ("NAG@2", "BMA@3", "14bb"),
    #     ("BMA@3", "MAN@4", "13ab"),
    #     ("BMA@3", "MAN@5", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    iupac = "Man(a1-6)[Man(a1-3)]Man(b1-4)Nag(b1-4)Nag(b1-"
    mol = glycan("test", iupac)
    mol.show()

    # glycans = [
    #     ("NAG@A", "NAG@B", "14bb"),
    #     ("NAG@B", "BMA@C", "14bb"),
    #     ("BMA@C", "MAN@D", "13ab"),
    #     ("BMA@C", "MAN@E", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    pass
