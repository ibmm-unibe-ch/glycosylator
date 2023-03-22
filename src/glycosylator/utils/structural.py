"""
Auxiliary functions related to structure handling
"""

import os
import pickle
import warnings

import numpy as np
import Bio.PDB as bio
import numba

import glycosylator.utils.defaults as defaults
import glycosylator.utils.constants as constants


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
        bond_length = defaults.DEFAULT_BOND_LENGTH

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

    # bonds = [Bond(i) for i in bonds]

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
        bond_length = defaults.DEFAULT_BOND_LENGTH

    if restrict_residues:
        bonds = []
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            _neighbors = bio.NeighborSearch(atoms)
            bonds.extend(
                _neighbors.search_all(radius=bond_length),
            )

    else:
        atoms = list(structure.get_atoms())
        _neighbors = bio.NeighborSearch(atoms)
        bonds = _neighbors.search_all(radius=bond_length)

    # # convert to Bond-tuples
    # bonds = [Bond(i) for i in bonds]

    return bonds


def apply_standard_bonds(structure, _topology=None):
    """
    Apply bonds according to loaded standard topologies. This will compute a list of tuples with bonded
    atoms from the same residue. It will not infer residue-residue connections! It is possible to provide a single
    atom as input structure, in which case all bonds that involve said atom are returned.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which hosts Atoms:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
        - `Bio.PDB.Atom`

    _topology : Topology
        The topology to use for applying the bonds.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms that are bonded.
    """
    if structure.level == "R":
        residue = structure

        if _topology is None:
            _topology = defaults.get_default_topology()

        if not _topology.has_residue(residue.resname):
            warnings.warn(f"[ignoring] No abstract residue found in Topology for {residue.resname}!")
            return

        atoms = residue.child_dict
        abstract_residue = _topology.get_residue(residue.resname)
        bonds = [(atoms.get(bond.atom1.id), atoms.get(bond.atom2.id)) for bond in abstract_residue.bonds]
        bonds = [bond for bond in bonds if bond[0] and bond[1]] # make sure to have no None entries...
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_standard_bonds(residue, _topology)
        bonds = [bond for bond in bonds if bond[0].id == structure.id or bond[1].id == structure.id]
        return bonds

    else:
        bonds = []
        for residue in structure.get_residues():
            bonds.extend(apply_standard_bonds(residue, _topology))
        return bonds


def fill_missing_atoms(structure, _topology=None):
    """
    Fill missing atoms in residues.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`

    _topology
        A specific topology to use for references.
    """
    if structure.level in set(("S", "M", "C")):
        for residue in structure.get_residues():
            fill_missing_atoms(residue, _topology=_topology)
        return
    elif structure.level == "A":
        raise ValueError("Cannot fill missing atoms in an Atom object!")
    elif structure.level != "R":
        raise ValueError(f"Invalid level for structure: {structure.level}")

    if not _topology:
        _topology = defaults.get_default_topology()

    if not _topology.has_residue(structure.resname):
        warnings.warn(f"[ignoring] No abstract residue found in Topology for {structure.resname}!")
        return
    _abstract = _topology.get_residue(structure.resname)

    # first get all atoms that need to be added
    _atoms_to_add = _abstract.get_missing_atoms(structure)

    # now iterate through the atoms and impute them one by one
    # try to keep atoms that could not be imputed in a queue
    # in hopes of imputing them later, but abort if it's no good...
    _fail_timer = 0
    while len(_atoms_to_add) != 0:
        _atom = _atoms_to_add.pop(0)
        _success = compute_atom_coordinate(_atom, structure, _abstract)
        if not _success:
            _atoms_to_add.append(_atom)
            _fail_timer += 1
            if _fail_timer == len(_atoms_to_add) + 1:
                warnings.warn(f"Could not impute atoms: {[i.id for i in _atoms_to_add]}! Aborting...", RuntimeWarning)
                return
            continue
        _fail_timer = 0


def compute_atom_coordinate(atom, residue, abstract):
    """
    Compute the coordinates of an atom from internal coordinates from an AbstractResidue into a structure.

    Parameters
    ----------
    atom : AbstractAtom
        The atom to compute the coordinates for.
    residue : Bio.PDB.Residue
        The target residue wherein to add the atom.
    abstract : AbstractResidue
        The abstract residue to reference data from.

    Returns
    -------
    success: bool
        Whether the atom could successfully be integrated.
        It may be that a specific atom could not be computed due
        to lacking reference atoms. In this case it may be necessary
        to first impute other atoms.
    """

    # get the internal coordinates to impute the atom positions from
    # since we can either compute atoms at position 1 or 4, we will specifically
    # search for ICs with these properties...
    _coord_func = compute_atom1_from_others
    _ics = abstract.get_internal_coordinates(atom, None, None, None, mode="partial")
    if len(_ics) == 0:
        _coord_func = compute_atom4_from_others
        _ics = abstract.get_internal_coordinates(None, None, None, atom, mode="partial")
        if len(_ics) == 0:
            # warnings.warn(f"Cannot find suitable IC for {atom}!")
            return False

    # try finding an IC that has available reference atoms
    for ic in _ics:
        _ref = ic.get_reference_atoms(residue)
        if len(_ref) == 3:
            break
    else:
        # warnings.warn(f"Cannot find suitable IC for {atom}!")
        return False

    # get the coordinates of the reference atoms
    _ref = [i.coord for i in _ref]

    # now use the internal coordinates to compute the new atom coordinate
    _coord = _coord_func(*_ref, ic)

    # create the new atom
    _new_atom = atom.to_biopython()
    _new_atom.set_coord(_coord)

    # add the new atom to the residue
    residue.add(_new_atom)
    return True


def compute_atom1_from_others(coords2, coords3, coords4, ic):
    """
    Compute the coordinates of the first atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : AbstractInternalCoordinates
        The internal coordinates

    Returns
    -------
    coords1 : array-like
        The coordinates of the first atom
    """
    if ic.is_proper:
        return _IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=coords2,
            r=-ic.bond_length_12,
            theta=-np.radians(ic.bond_angle_123),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        _vec = _IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=np.full(3, 0),
            r=1,
            theta=np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
        if ic.bond_length_13:
            _vec *= ic.bond_length_13
        else:
            BC = np.linalg.norm(coords2 - coords3)
            AC = ic.bond_length_12
            AB = np.sqrt(BC**2 + AC**2 - 2 * BC * AC * np.cos(np.radians(ic.bond_angle_123)))
            _vec *= AB
        final = coords3 + _vec
        return final


def compute_atom4_from_others(coords1, coords2, coords3, ic):
    """
    Compute the coordinates of the fourth atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : AbstractInternalCoordinates
        The internal coordinates

    Returns
    -------
    coords4 : array-like
        The coordinates of the fourth atom
    """
    if ic.is_proper:
        return _IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        return _IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )


@numba.njit
def _rotation_matrix(axis, angle):
    """
    Compute the rotation matrix about an arbitrary axis in 3D

    Source
    ------
    Stackoverflow thread: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    Parameters
    ----------
    axis : array-like
        The axis to rotate around
    angle : float
        The angle to rotate by (in radians)

    Returns
    -------
    rotation_matrix : array-like
        The rotation matrix
    """
    # axis = np.asarray(axis)
    # angle = np.asarray(angle)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def _IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    compute the coordinates of a fourth atom from a proper internal coordinate
    system and the other three atom coordinates.

    Parameters
    ----------
    a, b, c : np.ndarray
        coordinates of the other three atoms
    anchor : np.ndarray
        coordinates of the anchor atom relative to which the new coordinate should be calculated
    r : float
        bond length of the new atom relative to the anchor
    theta : float
        bond angle between the new atom and its plane partners
    dihedral : float
        dihedral angle of the internal coordinates
    """
    ab = b - a
    bc = c - b

    # compute normalized bond vectors for available atoms
    ab /= np.linalg.norm(ab)
    bc /= np.linalg.norm(bc)

    # compute plane vector for atoms 1-2-3
    plane_abc = np.cross(ab, bc)
    plane_abc /= np.linalg.norm(plane_abc)

    # rotate the plane vector around the middle bond (2-3) to get the plane 2-3-4
    _rot = _rotation_matrix(bc, dihedral)
    plane_bcd = np.dot(_rot, plane_abc)
    plane_bcd /= np.linalg.norm(plane_bcd)

    # rotate the middle bond around the new plane
    _rot = _rotation_matrix(plane_bcd, theta)
    cd = np.dot(_rot, bc)
    cd /= np.linalg.norm(cd)

    # compute the coordinates of the fourth atom
    d = anchor + r * cd
    return d


if __name__ == '__main__':

    from copy import deepcopy

    to_deletes = {'H2', 'H3', 'H4', 'H5', 'H61', 'H62'}

    MANNOSE = bio.PDBParser().get_structure("MAN", "./support/examples/MAN.pdb")

    for to_delete in to_deletes:

        _man = MANNOSE.copy()
        _man = next(_man.get_residues())
        top = deepcopy(defaults.get_default_topology())
        abstract = top.get_residue(_man.resname)
        for idx, i in enumerate(abstract.internal_coordinates):
            if i.is_proper:
                del abstract.internal_coordinates[idx]

        true_coords = _man.child_dict.get(to_delete)
        assert true_coords is not None, f"No atom {to_delete} found!"

        true_coords = true_coords.coord

        _man.detach_child(to_delete)
        assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

        fill_missing_atoms(_man, top)

        assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

        # io = bio.PDBIO()
        # io.set_structure(_man)
        # io.save("test.pdb")

        new_coords = _man.child_dict.get(to_delete).coord

        _diff = np.sum(np.abs(new_coords - true_coords))
        assert _diff < 1, f"[{to_delete}] Difference in coordinates is {new_coords - true_coords} ({_diff=})"
