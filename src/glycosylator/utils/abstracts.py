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

    def get_atom(self, id):
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

    def get_bond(self, id1, id2):
        """
        Get a bond by its atom IDs
        """
        if isinstance(id1, AbstractAtom):
            id1 = id1.id
        if isinstance(id2, AbstractAtom):
            id2 = id2.id
        for bond in self.bonds:
            if bond[0].id == id1 and bond[1].id == id2:
                return bond
            elif bond[0].id == id2 and bond[1].id == id1:
                return bond
        return None

    def get_internal_coordinates(self, id1, id2, id3, id4):
        """
        Get an internal coordinate by its atom IDs
        """
        if isinstance(id1, AbstractAtom):
            id1 = id1.id
        if isinstance(id2, AbstractAtom):
            id2 = id2.id
        if isinstance(id3, AbstractAtom):
            id3 = id3.id
        if isinstance(id4, AbstractAtom):
            id4 = id4.id
        for ic in self.internal_coordinates:
            if ic[0].id == id1 and ic[1].id == id2 and ic[2].id == id3 and ic[3].id == id4:
                return ic
            elif ic[0].id == id4 and ic[1].id == id3 and ic[2].id == id2 and ic[3].id == id1:
                return ic
        return None

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


@attr.s
class AbstractAtom:
    """
    A representation of a single Atom
    """

    id = attr.ib(type=str)
    type = attr.ib(default=None, type=str)
    charge = attr.ib(default=None, type=float)
    mass = attr.ib(default=None, type=float)
    element = attr.ib(default=None, type=str)
    _parent = attr.ib(default=None)
    is_wildcard = attr.ib(default=False, type=bool)

    def to_biopython(self, parent=None):
        """
        Returns a Bio.PDB.Atom object

        Parameters
        ----------
        parent: Bio.PDB.Residue
            The parent residue of the atom
        """
        new = bio.Atom.Atom(
            self.id,
            coord=np.full(3, np.nan),
            serial_number=np.nan,
            bfactor=0.0,
            occupancy=0.0,
            altloc="",
            fullname=self.id,
            element=self.element,
            pqr_charge=self.charge,
        )
        if parent:
            new.set_parent(parent)
        elif self._parent:
            new.set_parent(self._parent)
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

    def __repr__(self):
        return f"AbstractBond({self.atom1.id}, {self.atom2.id})"


@attr.s
class AbstractInternalCoordinates:
    """
    A representation of a single Internal Coordinate
    """

    atom1 = attr.ib(type=Union[AbstractAtom, str])
    atom2 = attr.ib(type=Union[AbstractAtom, str])
    atom3 = attr.ib(type=Union[AbstractAtom, str])
    atom4 = attr.ib(type=Union[AbstractAtom, str])
    bond_length_12 = attr.ib(type=float)
    bond_length_34 = attr.ib(type=float)
    bond_angle_123 = attr.ib(type=float)
    bond_angle_234 = attr.ib(type=float)
    dihedral = attr.ib(type=float)
    bond_length_13 = attr.ib(default=None, type=float)
    improper = attr.ib(default=False, type=bool)

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

    def get_parent(self):
        """
        Get the parent of the Internal Coordinate (i.e. it's residue)
        """
        return self.atom1.get_parent()

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
        """
        return self._delete_ids

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
    angle = attr.ib(type=float)
    K = attr.ib(type=float, default=None)
    urey_bradley_k = attr.ib(type=float, default=None)
    urey_bradley_length = attr.ib(type=float, default=None)

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
    K = attr.ib(type=float, default=None)
    multiplicity = attr.ib(type=int, default=1)

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
    angle = attr.ib(type=float)
    K = attr.ib(type=float, default=None)

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


if __name__ == '__main__':

    x = AtomTypeMass("x", 1.23)
