"""
Auxiliary functions related to structure handling
"""

import os
import pickle
import warnings

import Bio.PDB as bio

import glycosylator.utils.defaults as defaults
import glycosylator.utils.constants as constants

DEFAULT_BOND_LENGTH = defaults.DEFAULT_BOND_LENGTH


# =================================================================
# Structure related classes
# =================================================================


class Bond(tuple):
    """
    A bond between two atoms
    """

    def __new__(cls, atom1, atom2):
        return super(Bond, cls).__new__(cls, (atom1, atom2))

    def __init__(self, atom1, atom2):
        super().__init__()
        self._frozen = False

    @property
    def is_frozen(self):
        return self._frozen

    @is_frozen.setter
    def is_frozen(self, value):
        self._frozen = value

    def freeze(self):
        self._frozen = True

    def unfreeze(self):
        self._frozen = False

    def __hash__(self):
        return hash(tuple(sorted(self)))

    def __eq__(self, other):
        return set(self) == set(other)

    def __repr__(self):
        return f"Bond({self[0]}, {self[1]}, frozen={self.is_frozen})"


# =================================================================
# Structure related functions
# =================================================================


def infer_residue_connections(structure, bond_length: float = None):
    """
    Infer the connectivity graph of residues from the distances between atoms of residue pairs.
    This will establish only bonds between close-by atoms from non-identical residues.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host at least two Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
    bond_length : float
        The maximum distance between two atoms to be considered a bond.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms from different Residues that are bonded.
    """
    if bond_length is None:
        bond_length = DEFAULT_BOND_LENGTH

    bonds = []
    _seen_residues = set()
    for residue1 in structure.get_residues():
        for residue2 in structure.get_residues():
            if residue1 == residue2:
                continue
            elif residue2 in _seen_residues:
                continue

            atoms1 = list(residue1.get_atoms())
            atoms2 = list(residue2.get_atoms())

            _neighbors = bio.NeighborSearch(atoms1 + atoms2)
            _neighbors = _neighbors.search_all(radius=bond_length)

            bonds.extend(
                [n for n in _neighbors if n[0].get_parent() != n[1].get_parent()],
            )

        _seen_residues.add(residue1)

    return bonds


def infer_bonds(structure, bond_length: float = None, restrict_residues: bool = True):
    """
    Generate a connectivity graph by inferring bonds from the distance between atoms.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
    bond_length : float
        The maximum distance between two atoms to be considered a bond.
    restrict_residues : bool
        If set to `True`, only bonds between atoms of the same residue will be considered.

    Returns
    -------
    bonds : list
        The connectivity graph of the molecule, storing tuples of `Bio.PDB.Atom` objects.
    """
    if bond_length is None:
        bond_length = DEFAULT_BOND_LENGTH

    if restrict_residues:
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            _neighbors = bio.NeighborSearch(atoms)
            bonds = _neighbors.search_all(radius=bond_length)
    else:
        atoms = list(structure.get_atoms())
        _neighbors = bio.NeighborSearch(atoms)
        bonds = _neighbors.search_all(radius=bond_length)

    return bonds
