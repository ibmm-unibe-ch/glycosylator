"""
Abstract classes for storing force field data from CHARMM topology and parameter files
"""

from typing import Union, NamedTuple
import attr
import Bio.PDB as bio
import numpy as np


class AbstractEntity:
    """
    A representation of a single collection entity (Residue, Patch, etc.)
    """

    def __init__(self, id=None):
        self.id = id
        self._atoms = {}
        self.bonds = []
        self.internal_coordinates = []

    @property
    def atoms(self):
        return list(self._atoms.values())

    def get_atom(self, id) -> "AbstractAtom":
        """
        Get an atom by its ID

        Parameters
        ----------
        id: str
            The ID or type of the atom

        Returns
        -------
        atom: AbstractAtom
            The atom with the given ID.
            If no atom is found, None is returned.
        """
        if isinstance(id, (list, tuple)):
            return [self.get_atom(i) for i in id]
        atom = self._atoms.get(id, None)
        return atom

    def get_atoms_by_type(self, _type):
        """
        Get a list of atoms by their type

        Parameters
        ----------
        _type: str
            The type of the atom

        Returns
        -------
        atoms: list
            A list of atoms with the given type.
        """
        atoms = [i for i in self._atoms.values() if i.type == _type]
        return atoms

    def has_atom(self, _atom):
        """
        Check if the entity has an atom

        Parameters
        ----------
        _atom: AbstractAtom
            The atom to check for

        Returns
        -------
        has_atom: bool
            True if the atom is in the entity, False otherwise
        """
        if isinstance(_atom, AbstractAtom):
            _atom = _atom.id
        return _atom in self._atoms

    def get_missing_atoms(self, obj):
        """
        Get a list of atoms that are missing from the entity

        Parameters
        ----------
        obj
            The object to compare to which hosts atoms and has a `get_atoms()` method.
            E.g. a `bio.PDB.Residue`.

        Returns
        -------
        missing_atoms: list
            A list of atoms that are missing from the entity
        """
        _atoms = set(i.id for i in obj.get_atoms())
        missing_atoms = [atom for atom in self.atoms if atom.id not in _atoms]
        return missing_atoms

    def get_bond(self, id1, id2) -> "AbstractBond":
        """
        Get a bond by its atom IDs
        """
        if isinstance(id1, AbstractAtom):
            id1 = id1.id
        if isinstance(id2, AbstractAtom):
            id2 = id2.id
        for bond in self.bonds:
            if bond.atom1.id == id1 and bond.atom2.id == id2:
                return bond
            elif bond.atom1.id == id2 and bond.atom2.id == id1:
                return bond
        return None

    def get_internal_coordinates(self, *ids, mode: str = "exact"):
        """
        Get internal coordinates by their atom IDs

        Parameters
        ----------
        ids: str
            The IDs of the atoms in the internal coordinate (this can also be any data object that has an `id` attribute)
        mode: str
            The mode to use for matching the internal coordinates.
            Supported modes are:
            - `exact`: The internal coordinate must match the given atom IDs exactly (requires that four ids are given)
            - `partial`: The internal coordinate must match the given atom IDs where they are provided, but wildcards can be set to None (requires that four ids or None are given in order)
            - `multi_partial`: The internal coordinate must match any given atom IDs where they are provided, but wildcards can be set to None (requires that four ids or None are given in order)
            - `anywhere`: The internal coordinate must contain all of the given atom IDs in any order (requires at least one id)
            - `anywhere_partial`: The internal coordinate must contain any of the given atom IDs in any order (requires at least one id)

        Returns
        -------
        ics: list
            A list of internal coordinates
        """

        if len(ids) == 0:
            return self.internal_coordinates

        ids = tuple(i.id if hasattr(i, "id") else i for i in ids)

        if mode == "exact":

            if len(ids) != 4:
                raise ValueError(
                    "Exact mode requires that four ids are given to match the internal coordinates"
                )

            for ic in self.internal_coordinates:
                if ids == ic.ids or ids[::-1] == ic.ids:
                    return [ic]
            return []

        elif mode == "partial":

            if len(ids) != 4:
                raise ValueError(
                    "Partial mode requires that four ids or None are given to match the internal coordinates"
                )

            ids = np.array(ids)
            mask = ids != None

            ics = [
                ic
                for ic in self.internal_coordinates
                if np.all(ids[mask] == np.array(ic.ids)[mask])
            ]
            return ics

        elif mode == "anywhere":

            ids = set(ids)
            ics = [ic for ic in self.internal_coordinates if ids.issubset(set(ic.ids))]
            return ics

        elif mode == "anywhere_partial":

            ids = set(ids)
            ics = [
                ic
                for ic in self.internal_coordinates
                if len(ids.intersection(set(ic.ids))) != 0
            ]
            return ics

        elif mode == "multi_partial":

            if len(ids) != 4:
                raise ValueError(
                    "Multi partial mode requires that four ids or None are given to match the internal coordinates"
                )

            ids = np.array(ids)
            mask = ids != None

            ics = [
                ic
                for ic in self.internal_coordinates
                if np.any(ids[mask] == np.array(ic.ids)[mask])
            ]
            return ics

        else:
            raise ValueError(f"Unknown mode {mode}")

    def add_atom(self, atom):
        """
        Add an atom to the residue
        """
        self._atoms[atom.id] = atom

    def add_bond(self, bond):
        """
        Add a bond to the residue
        """
        self.bonds.append(bond)

    def add_internal_coordinates(self, ic):
        """
        Add an internal coordinate to the residue
        """
        self.internal_coordinates.append(ic)

    def add_ic(self, ic):
        """
        Add an internal coordinate to the residue
        """
        self.add_internal_coordinates(ic)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.id})"


class AbstractResidue(AbstractEntity):
    """
    A representation of a single Residue (compound)
    """

    def __init__(self, id=None):
        super().__init__(id)
        self.PDB_id = None
        self.CHARMM_id = None

    def to_biopython(self, seqid: int = 1):
        """
        Convert to a biopython residue.

        Parameters
        ----------
        seqid: int
            The sequence ID of the residue

        Returns
        -------
        residue: bio.PDB.Residue
            The biopython residue
        """
        residue = bio.PDB.Residue.Residue(("H_" + self.id, seqid, " "), self.id, " ")

        for atom in self.atoms:
            residue.add(atom.to_biopython())

        return residue


@attr.s(hash=True)
class AbstractAtom:
    """
    A representation of a single Atom
    """

    id = attr.ib(type=str, hash=True)
    type = attr.ib(default=None, type=str, hash=True)
    charge = attr.ib(default=None, type=float, repr=False)
    mass = attr.ib(default=None, type=float, repr=False)
    _element = attr.ib(default=None, type=str, repr=False)
    _parent = attr.ib(default=None, repr=False)
    is_wildcard = attr.ib(default=False, type=bool, repr=False)
    coord = attr.ib(default=None, type=np.ndarray, repr=False)

    @property
    def element(self):
        if self._element is None:
            return self.id[0]
        return self._element

    def to_biopython(self, serial_number: int = 1):
        """
        Returns a Bio.PDB.Atom object

        Parameters
        ----------
        serial_number: int
            The serial number of the atom
        """
        if self.coord:
            coord = self.coord
        else:
            coord = np.full(3, np.nan)
        new = bio.Atom.Atom(
            self.id,
            coord=coord,
            serial_number=serial_number,
            bfactor=0.0,
            occupancy=0.0,
            altloc="",
            fullname=self.id,
            element=self.element,
            pqr_charge=self.charge,
        )
        return new

    def get_parent(self):
        """
        Get the parent of the atom
        """
        return self._parent

    def set_parent(self, obj):
        """
        Set the parent of the atom
        """
        self._parent = obj

    def __repr__(self):
        return f"AbstractAtom({self.id})"


@attr.s
class AbstractBond:
    """
    A representation of a single Bond between two atoms (or atom types)
    optionally, a bond spring constant (K) and length can be supplied.
    """

    atom1 = attr.ib(type=Union[AbstractAtom, str])
    atom2 = attr.ib(type=Union[AbstractAtom, str])

    K = attr.ib(type=float, default=None)
    length = attr.ib(type=float, default=None)

    @property
    def atoms(self):
        return self.atom1, self.atom2

    def get_parent(self):
        """
        Get the parent of the bond (i.e. it's residue)
        """
        if isinstance(self.atom1, str):
            return None
        return self.atom1.get_parent()

    def find_atoms(self, obj):
        """
        Find the atoms of the bond in an object.
        This will return a tuple of identified atoms with the
        same id if they exist in the object, None for any atom
        that was not found.
        """
        if isinstance(self.atom1, str):
            atom1 = obj.get_atom(self.atom1)
        else:
            atom1 = obj.get_atom(self.atom1.id)

        if isinstance(self.atom2, str):
            atom2 = obj.get_atom(self.atom2)
        else:
            atom2 = obj.get_atom(self.atom2.id)

        return atom1, atom2

    def migrate_atoms(self, obj):
        """
        Migrate the atoms of the bond to a new object.
        This will update the atom1 and atom2 attributes of the bond
        to point to the atoms in the new object if they can be found.
        """
        atom1, atom2 = self.find_atoms(obj)
        self.atom1 = atom1 if atom1 else self.atom1
        self.atom2 = atom2 if atom2 else self.atom2

    def __get_item__(self, key):
        return self.atoms[key]

    def __repr__(self):
        return f"AbstractBond({self.atom1.id}, {self.atom2.id})"


@attr.s
class AbstractInternalCoordinates:
    """
    A representation of a single Internal Coordinate

    Parameters
    ----------
    atom1
        The first atom in the internal coordinate
    atom2
        The second atom in the internal coordinate
    atom3
        The third atom in the internal coordinate
    atom4
        The fourth atom in the internal coordinate
    bond_length_12
        The bond length between atoms 1 and 2
    bond_length_34
        The bond length between atoms 3 and 4
    bond_angle_123
        The bond angle between atoms 1, 2 and 3
    bond_angle_234
        The bond angle between atoms 2, 3 and 4
    dihedral
        The dihedral angle between atoms 1, 2, 3 and 4
    bond_length_13
        The bond length between atoms 1 and 3 (optional, for impropers)
    improper
        Whether the internal coordinate is an improper (optional)
    """

    # With this set of data, the position of atom 1 may be determined based on the
    # positions of atoms 2-4, and the position of atom 4 may be determined from the
    # positions of atoms 1-3, allowing the recursive generation of coordinates for
    # all atoms in the structure based on a three-atom seed.

    atom1 = attr.ib(type=Union[AbstractAtom, str])
    atom2 = attr.ib(type=Union[AbstractAtom, str])
    atom3 = attr.ib(type=Union[AbstractAtom, str])
    atom4 = attr.ib(type=Union[AbstractAtom, str])
    bond_length_12 = attr.ib(type=float, repr=False)
    bond_length_34 = attr.ib(type=float, repr=False)
    bond_angle_123 = attr.ib(type=float, repr=False)
    bond_angle_234 = attr.ib(type=float, repr=False)
    dihedral = attr.ib(type=float, repr=False)
    bond_length_13 = attr.ib(default=None, type=float, repr=False)
    improper = attr.ib(default=False, type=bool, repr=False)

    @property
    def angles(self):
        """
        Returns the bond angles that constitute the internal coordinate
        """
        return (self.bond_angle_123, self.bond_angle_234)

    @property
    def lengths(self):
        """
        Returns the bond lengths that constitute the internal coordinate
        """
        if self.improper:
            return (self.bond_length_13, self.bond_length_34)
        else:
            return (self.bond_length_12, self.bond_length_34)

    @property
    def atoms(self):
        """
        Returns the atoms that constitute the internal coordinate
        """
        return (self.atom1, self.atom2, self.atom3, self.atom4)

    @property
    def ids(self):
        """
        Returns the ids of the atoms that constitute the internal coordinate
        """
        if hasattr(self.atom1, "id"):
            return (self.atom1.id, self.atom2.id, self.atom3.id, self.atom4.id)
        else:
            return (self.atom1, self.atom2, self.atom3, self.atom4)

    @property
    def is_improper(self):
        """
        Returns True if the internal coordinate is improper
        """
        return self.improper

    @property
    def is_proper(self):
        """
        Returns True if the internal coordinate is proper
        """
        return not self.improper

    def get_reference_atoms(self, _src):
        """
        Get reference atoms with available coordinates from a source.
        Note, that this method does not check if reference coordinates
        are actually available, so make sure to only provide source objects
        that have cartesian coordinates contained.

        Parameters
        ----------
        _src : object
            The source object to get the reference atom from. This can be:
            - a list or tuple of Atoms
            - an object with a `get_atoms` method (e.g. a biopython Residue)

        Returns
        -------
        list
            A list of reference atoms with available coordinates
        """

        if hasattr(_src, "get_atoms"):
            _src = _src.get_atoms()
        elif not isinstance(_src, (list, tuple, set)):
            raise ValueError("Invalid source object")

        _src = {a.id: a for a in _src}
        _set_src = set(_src.keys())

        _intersection = set(self.ids).intersection(_set_src)
        if len(_intersection) == 0:
            return []

        _ref_atoms = []
        for _id in self.ids:
            if _id in _intersection:
                _ref_atoms.append(_src[_id])

        return _ref_atoms

    def __repr__(self):
        return f"AbstractInternalCoordinates({self.atom1.id}, {self.atom2.id}, {self.atom3.id}, {self.atom4.id})"


class AbstractPatch(AbstractEntity):
    """
    A representation of a single Patch
    """

    def __init__(self, id=None) -> None:
        super().__init__(id)
        self._delete_ids = []

    @property
    def deletes(self):
        """
        Returns the atom IDs to delete
        in a tuple of lists where the first list
        contains the atom IDs to delete from the
        first structure and the second one from the second structure
        """
        deletes = (
            [i[1:] for i in self._delete_ids if i[0] == "1"],
            [i[1:] for i in self._delete_ids if i[0] == "2"],
        )
        return deletes

    @property
    def IC_atom_ids(self):
        """
        Returns a set of all atom IDs of all atoms for which the patch also stores
        internal coordinates.
        """
        ids = set()
        for ic in self.internal_coordinates:
            if isinstance(ic.atom1, str):
                ids.update(ic.atoms)
                continue
            ids.update(ic.ids)
        return ids

    def add_delete(self, id):
        """
        Add an atom ID to delete
        """
        self._delete_ids.append(id)

    add_id_to_delete = add_delete


@attr.s
class AbstractAngle:

    atom1 = attr.ib(type=str)
    atom2 = attr.ib(type=str)
    atom3 = attr.ib(type=str)
    angle = attr.ib(type=float, repr=False)
    K = attr.ib(type=float, default=None, repr=False)
    urey_bradley_k = attr.ib(type=float, default=None, repr=False)
    urey_bradley_length = attr.ib(type=float, default=None, repr=False)

    @property
    def atoms(self):
        return self.atom1, self.atom2, self.atom3

    @property
    def Kub(self):
        """A synonym of 'urey_bradley_k'"""
        return self.urey_bradley_k


@attr.s
class AbstractDihedral:

    atom1 = attr.ib(type=str)
    atom2 = attr.ib(type=str)
    atom3 = attr.ib(type=str)
    atom4 = attr.ib(type=str)
    angle = attr.ib(type=float)
    K = attr.ib(type=float, default=None, repr=False)
    multiplicity = attr.ib(type=int, default=1, repr=False)

    @property
    def atoms(self):
        return self.atom1, self.atom2, self.atom3, self.atom4

    @property
    def n(self):
        """A synonym of 'multiplicity'"""
        return self.multiplicity


@attr.s
class AbstractImproper:

    atom1 = attr.ib(type=str)
    atom2 = attr.ib(type=str)
    atom3 = attr.ib(type=str)
    atom4 = attr.ib(type=str)
    angle = attr.ib(type=float, repr=False)
    K = attr.ib(type=float, default=None, repr=False)

    @property
    def atoms(self):
        return self.atom1, self.atom2, self.atom3, self.atom4


class AbstractNonBonded(NamedTuple):

    atom: str
    epsilon: float
    min_radius: float
    epsilon_14: float = None
    min_radius_14: float = None

    @property
    def sigma(self):
        """
        Estimate the atom size.
        """
        return self.min_radius * 2 ** (1 / 6)
