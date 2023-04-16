"""
Functions to infer structural data such as missing atom coordinates or bond connectivity
"""

import warnings
import numpy as np

from Bio.PDB import NeighborSearch

import glycosylator.utils.abstract as _abstract
import glycosylator.utils.defaults as defaults
import glycosylator.structural.base as base
import glycosylator.structural.neighbors as neighbors


def compute_residue_radius(residue):
    """
    Compute the radius of a residue by computing the distance between its center of mass
    and the furthest atom.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        The residue to compute the radius for.

    Returns
    -------
    radius : float
        The radius of the residue.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    radius = np.max(distances)
    return radius


def compute_outlier_atoms(residue, f: float = 1.5):
    """
    Compute which atoms of a residue are especially far away from the residues center of mass.
    This function compute the distances between the center of mass of the residue and its atoms
    and returns all atoms are further than `f * p75` away from the center of mass, where `p75`
    is the 75th percentile of the distances.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        The residue to compute the outlier atoms for.
    f : float
        The factor to multiply the 75th percentile with.

    Returns
    -------
    outlier_atoms : list
        A list of atoms that are considered outliers.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    p75 = np.percentile(distances, 75)
    f = f * p75

    outlier_atoms = [atom for atom, distance in zip(atoms, distances) if distance > f]
    return outlier_atoms


def infer_surface_residues(
    structure,
    cutoff: int = 75,
    fraction: float = None,
):
    """
    Infer residues that are likely to be on the surface of the structure
    using the Solvent accessible surface area (SASA) of the structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the surface residues from.
    n_points : int
        The number of points to sample on the surface of the structure.
    cutoff : int
        The cutoff to use for classifying residues as surface residues.
    fraction : float
        The fraction of residues to classify as surface residues. In this case,
        the cutoff is adjusted to match the fraction of residues.

    Returns
    -------
    surface_residues : list
        A list of residues that are likely to be on the surface of the structure.
    """

    sasa = defaults.get_default_instance("bioSASA")
    sasa.compute(structure, level="R")

    sasa_values = np.array([residue.sasa for residue in structure.get_residues()])
    sasa_values = sasa_values / sasa_values.max() * 100

    if fraction is not None:
        cutoff = np.percentile(sasa_values, 100 - fraction * 100)

    surface_residues = [
        residue
        for residue, sasa in zip(structure.get_residues(), sasa_values)
        if sasa > cutoff
    ]
    return surface_residues


def infer_residue_connections(
    structure,
    bond_length: float = None,
    triplet: bool = False,
):
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
    triplet : bool
        If True, bonds between atoms of the same residue are also added, if one
        of the atoms is considered bonded to another residue. Like this residue connections
        are described not by a single bond with a pair of atoms, but two bonds with a triplet of atoms.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms from different Residues that are bonded.
    """
    if not triplet:
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

                _neighbors = NeighborSearch(atoms1 + atoms2)
                _neighbors = _neighbors.search_all(radius=bond_length)

                _neighbors = (
                    i for i in _neighbors if i[0].element != "H" and i[1].element != "H"
                )

                bonds.extend(
                    [n for n in _neighbors if n[0].get_parent() != n[1].get_parent()],
                )

            _seen_residues.add(residue1)
    else:
        connections = infer_residue_connections(
            structure, bond_length=bond_length, triplet=False
        )
        base_bonds = infer_bonds(
            structure, bond_length=bond_length, restrict_residues=True
        )
        _additional_bonds = []
        for atom1, atom2 in connections:
            _new = [bond for bond in base_bonds if atom1 in bond]
            _additional_bonds.extend(_new)
        bonds = connections + _additional_bonds

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
            _neighbors = NeighborSearch(atoms)
            bonds.extend(
                _neighbors.search_all(radius=bond_length),
            )

    else:
        atoms = list(structure.get_atoms())
        _neighbors = NeighborSearch(atoms)
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
            warnings.warn(
                f"[ignoring] No abstract residue found in Topology for {residue.resname}!"
            )
            return

        atoms = residue.child_dict
        abstract_residue = _topology.get_residue(residue.resname)
        bonds = [
            (atoms.get(bond.atom1.id), atoms.get(bond.atom2.id))
            for bond in abstract_residue.bonds
        ]
        bonds = [
            bond for bond in bonds if bond[0] and bond[1]
        ]  # make sure to have no None entries...
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_standard_bonds(residue, _topology)
        bonds = [
            bond
            for bond in bonds
            if bond[0].id == structure.id or bond[1].id == structure.id
        ]
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
        warnings.warn(
            f"[ignoring] No abstract residue found in Topology for {structure.resname}!"
        )
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
        _success = impute_atom_from_IC(_atom, structure, _abstract)
        if not _success:
            _atoms_to_add.append(_atom)
            _fail_timer += 1
            if _fail_timer == len(_atoms_to_add) + 1:
                warnings.warn(
                    f"Could not impute atoms: {[i.id for i in _atoms_to_add]}! Aborting...",
                    RuntimeWarning,
                )
                return
            continue
        _fail_timer = 0


def compute_internal_coordinates(bonds: list):
    """
    Compute internal coordinates for a structure.

    Parameters
    ----------
    bonds : list
        A list of tuples of atoms that are bonded.
        The atoms must be `Bio.PDB.Atom` objects with coordinates.

    Returns
    -------
    list
        A list of AbstractInternalCoordinates
    """
    quartets = neighbors.compute_quartets(bonds)
    ics = []
    for quartet in quartets:
        angle_123 = base.compute_angle(quartet[0], quartet[1], quartet[2])
        angle_234 = base.compute_angle(quartet[1], quartet[2], quartet[3])
        dihedral = base.compute_dihedral(quartet[0], quartet[1], quartet[2], quartet[3])
        l_12 = base.compute_distance(quartet[0], quartet[1])
        # l_23 = base.compute_distance(quartet[1], quartet[2])
        l_34 = base.compute_distance(quartet[2], quartet[3])

        ic = _abstract.AbstractInternalCoordinates(
            quartet[0].id,
            quartet[1].id,
            quartet[2].id,
            quartet[3].id,
            l_12,
            l_34,
            angle_123,
            angle_234,
            dihedral,
            improper=quartet.improper,
        )
        ics.append(ic)
    return ics


def impute_atom_from_IC(atom, residue, abstract):
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
    success, _coord = compute_atom_coordinate_from_IC(atom, residue, abstract)
    if not success:
        return False

    # create the new atom
    _new_atom = atom.to_biopython()
    _new_atom.set_coord(_coord)

    # add the new atom to the residue
    residue.add(_new_atom)
    return True


def compute_atom_coordinate_from_IC(atom, residue, abstract):
    """
    Compute the coordinates of an atom from internal coordinates from an AbstractResidue into a structure.

    Parameters
    ----------
    atom : AbstractAtom
        The atom to compute the coordinates for.
    residue : Bio.PDB.Residue
        The target residue where other reference atoms are located.
    abstract : AbstractResidue
        The abstract residue to reference internal coordinates from.

    Returns
    -------
    success: bool
        Whether the atom could successfully be integrated.
    coords : np.ndarray
        The coordinates of the atom
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
            return False, None

    # try finding an IC that has available reference atoms
    for ic in _ics:
        _ref = ic.get_reference_atoms(residue)
        if len(_ref) == 3:
            break
    else:
        # warnings.warn(f"Cannot find suitable IC for {atom}!")
        return False, None

    # get the coordinates of the reference atoms
    _ref = [i.coord for i in _ref]

    # now use the internal coordinates to compute the new atom coordinate
    _coord = _coord_func(*_ref, ic)

    return True, _coord


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
        return base._IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=coords2,
            r=-ic.bond_length_12,
            theta=-np.radians(ic.bond_angle_123),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        _vec = base._IC_to_xyz(
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
            AB = np.sqrt(
                BC**2 + AC**2 - 2 * BC * AC * np.cos(np.radians(ic.bond_angle_123))
            )
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
        return base._IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        return base._IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
