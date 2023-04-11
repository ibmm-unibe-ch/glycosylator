"""
The base class for classes storing biopython structures
"""

from typing import Union
import warnings

import Bio.PDB as bio

import glycosylator.structural as structural
import glycosylator.utils as utils
from glycosylator.graphs import AtomGraph, ResidueGraph


class BaseEntity:
    def __init__(self, structure):
        self._base_struct = structure
        self._id = structure.id

        self._AtomGraph = None
        self._ResidueGraph = None

        self._bonds = []
        self._locked_bonds = set()

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        root_atom: Union[str, int] = None,
        id: str = None,
    ):
        """
        Read a Molecule from a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        return cls(struct, root_atom)

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        self._id = new_id
        self._base_struct.id = new_id

    @property
    def structure(self):
        return self._base_struct

    @property
    def bonds(self):
        return self._bonds

    @property
    def atom_triplets(self):
        """
        Compute triplets of three consequtively bonded atoms
        """
        if len(self._bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_triplets(self._bonds)

    @property
    def atom_quartets(self):
        """
        Compute quartets of four consequtively bonded atoms
        """
        if len(self._bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_quartets(self._bonds)

    @property
    def angles(self):
        """
        Compute all angles of consecutively bonded atom triplets within the molecule.

        Returns
        -------
        angles : dict
            A dictionary of the form {atom_triplet: angle}
        """
        return {
            triplet: structural.compute_angle(*triplet)
            for triplet in self.atom_triplets
        }

    @property
    def dihedrals(self):
        """
        Compute all dihedrals of consecutively bonded atom quartets within the molecule.

        Returns
        -------
        dihedrals : dict
            A dictionary of the form {atom_quartet: dihedral}
        """
        return {
            quartet: structural.compute_dihedral(*quartet)
            for quartet in self.atom_quartets
        }

    @property
    def AtomGraph(self):
        if not self._AtomGraph:
            self.make_atom_graph()
        return self._AtomGraph

    @property
    def ResidueGraph(self):
        if not self._ResidueGraph:
            self.make_residue_graph()
        return self._ResidueGraph

    def get_atoms(self, *atoms: Union[int, str, tuple], by: str = None) -> list:
        """
        Get one or more atoms from the structure either based on their
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues all of them are returned in a list.
        It is a safer option to use the full_id or serial number to retrieve a specific atom.
        If no search parameters are provided, the underlying iterator of the structure is returned.

        Parameters
        ----------
        atoms
            The atom id, serial number or full_id tuple,
            this supports multiple atoms to search for. However,
            only one type of parameter is supported per call.
        by : str
            The type of parameter to search for. Can be either 'id', 'serial' or 'full_id'
            If None is given, the parameter is inferred from the datatype of the atoms argument
            'serial' in case of `int`, 'id' in case of `str`, `full_id` in case of a tuple.

        Returns
        -------
        atom : list
            The atom(s)
        """
        if len(atoms) == 0:
            return self._base_struct.get_atoms()

        if by is None:
            if isinstance(atoms[0], int):
                by = "serial"
            elif isinstance(atoms[0], str):
                by = "id"
            elif isinstance(atoms[0], tuple):
                by = "full_id"
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
                )

        if by == "id":
            atoms = [i for i in self._base_struct.get_atoms() if i.id in atoms]
        elif by == "serial":
            atoms = [
                i for i in self._base_struct.get_atoms() if i.serial_number in atoms
            ]
        elif by == "full_id":
            atoms = [i for i in self._base_struct.get_atoms() if i.full_id in atoms]
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
            )

        return atoms

    def get_atom(self, atom: Union[int, str, tuple], by: str = None):
        """
        Get an atom from the structure either based on its
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues, only the first one is returned.
        To get all atoms matching the criteria, use the `get_atoms` method.

        Parameters
        ----------
        atom
            The atom id, serial number or full_id tuple
        by : str
            The type of parameter to search for. Can be either 'id', 'serial' or 'full_id'
            Because this looks for one specific atom, this parameter can be inferred from the datatype
            of the atom parameter. If it is an integer, it is assumed to be the serial number,
            if it is a string, it is assumed to be the atom id and if it is a tuple, it is assumed
            to be the full_id.

        Returns
        -------
        atom : bio.Atom.Atom
            The atom
        """
        if isinstance(atom, bio.Atom.Atom):
            return atom

        if by is None:
            if isinstance(atom, int):
                by = "serial"
            elif isinstance(atom, str):
                by = "id"
            elif isinstance(atom, tuple):
                by = "full_id"
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
                )

        if by == "id":
            atom = next(i for i in self._base_struct.get_atoms() if i.id == atom)
        elif by == "serial":
            atom = next(
                i for i in self._base_struct.get_atoms() if i.serial_number == atom
            )
        elif by == "full_id":
            atom = next(i for i in self._base_struct.get_atoms() if i.full_id == atom)
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
            )

        return atom

    def get_bonds(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom] = None,
        either_way: bool = True,
    ):
        """
        Get one or multiple bonds from the molecule. If only one atom is provided, all bonds
        that are connected to that atom are returned.

        Parameters
        ----------
        atom1
            The atom id, serial number or full_id tuple of the first atom
        atom2
            The atom id, serial number or full_id tuple of the second atom
        either_way : bool
            If True, the order of the atoms does not matter, if False, the order of the atoms
            does matter. By setting this to false, it is possible to also search for bonds that have
            a specific atom in position 1 or 2 depending on which argument was set, while leaving the other input as none.

        Returns
        -------
        bond : list
            The bond(s)
        """
        if atom1:
            atom1 = self.get_atom(atom1)
        if atom2:
            atom2 = self.get_atom(atom2)

        if atom1 and atom2:
            if either_way:
                return [i for i in self._bonds if atom1 in i and atom2 in i]
            else:
                return [
                    i
                    for i in self._bonds
                    if atom1 in i and atom2 in i and i[0] == atom1 and i[1] == atom2
                ]
        elif atom1:
            if either_way:
                return [i for i in self._bonds if atom1 in i]
            else:
                return [i for i in self._bonds if atom1 in i and i[0] == atom1]
        elif atom2:
            if either_way:
                return [i for i in self._bonds if atom2 in i]
            else:
                return [i for i in self._bonds if atom2 in i and i[1] == atom2]
        else:
            raise ValueError("No atom provided")

    def add_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
    ):
        """
        Add a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        self._bonds.append(bond)

    def remove_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        either_way: bool = True,
    ):
        """
        Remove a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        either_way : bool
            If True, the bond will be removed in both directions, i.e. if the bond is (1, 2)
            it will be removed if either (1, 2) or (2, 1) is provided.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        if either_way:
            if bond in self._bonds:
                self._bonds.remove(bond)
            if bond[::-1] in self._bonds:
                self._bonds.remove(bond[::-1])
        else:
            self._bonds.remove(bond)
        if bond in self._locked_bonds:
            self._locked_bonds.remove(bond)
        if either_way and bond[::-1] in self._locked_bonds:
            self._locked_bonds.remove(bond[::-1])

    def purge_bonds(self, atom: Union[int, str, bio.Atom.Atom]):
        """
        Remove all bonds connected to an atom

        Parameters
        ----------
        atom
            The atom to remove the bonds from, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        atom = self.get_atom(atom)
        bonds = [i for i in self._bonds if atom in i]
        for bond in bonds:
            if bond in self._bonds:
                self._bonds.remove(bond)
            elif bond[::-1] in self._bonds:
                self._bonds.remove(bond[::-1])

        bonds = [i for i in self._locked_bonds if atom in i]
        for bond in bonds:
            if bond in self._locked_bonds:
                self._locked_bonds.remove(bond)
            elif bond[::-1] in self._locked_bonds:
                self._locked_bonds.remove(bond[::-1])

    def lock_all(self, both_ways: bool = True):
        """
        Lock all bonds in the structure

        Parameters
        ----------
        both_ways : bool
            If True, the bond is locked in both directions
            i.e. atom1 --- atom2 direction will be unavailable for
            rotation as well as atom2 --- atom1 direction as well.
        """
        self._locked_bonds = set(self._bonds)
        if both_ways:
            self._locked_bonds.update(b[::-1] for b in self._bonds)

    def unlock_all(self):
        """
        Unlock all bonds in the structure
        """
        self._locked_bonds = set()

    def lock_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        both_ways: bool = False,
    ):
        """
        Lock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        both_ways : bool
            If True, the bond is locked in both directions. By default the bond is only locked in the specified direction.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        self._locked_bonds.add(bond)
        if both_ways:
            self._locked_bonds.add(bond[::-1])

    def unlock_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        both_ways: bool = False,
    ):
        """
        Unlock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        both_ways : bool
            If True, the bond is unlocked in both directions. By default the bond is only unlocked in the specified direction.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        if both_ways:
            if bond in self._locked_bonds:
                self._locked_bonds.remove(bond)
            if bond[::-1] in self._locked_bonds:
                self._locked_bonds.remove(bond[::-1])
            return
        self._locked_bonds.remove(bond)

    def is_locked(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
    ):
        """
        Check if a bond is locked

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        Returns
        -------
        bool
            True if the bond is locked, False otherwise
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        return bond in self._locked_bonds

    def infer_bonds(
        self, max_bond_length: float = None, restrict_residues: bool = True
    ) -> list:
        """
        Infer bonds between atoms in the structure

        Parameters
        ----------
        max_bond_length : float
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.
        restrict_residues : bool
            Whether to restrict bonds to only those in the same residue.
            If False, bonds between atoms in different residues are also inferred.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.infer_bonds(
            self._base_struct, max_bond_length, restrict_residues
        )
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

    def infer_residue_connections(
        self, max_bond_length: float = None, triplet: bool = True
    ) -> list:
        """
        Infer bonds between atoms that connect different residues in the structure

        Parameters
        ----------
        max_bond_length : float
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.
        triplet : bool
            Whether to include bonds between atoms that are in the same residue
            but neighboring a bond that connects different residues. This is useful
            for residues that have a side chain that is connected to the main chain.
            This is mostly useful if you intend to use the returned list for some purpose,
            because the additionally returned bonds are already present in the structure 
            from inference or standard-bond applying and therefore do not actually add any 
            particular information to the Molecule object itself.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded and considered residue connections.

        Examples
        --------
        For a molecule with the following structure:
        ```     
             connection -->  OA     OB --- H
                            /  \\   /
              (1)CA --- (2)CA   (1)CB 
               /         \\        \\
           (6)CA         (3)CA    (2)CB --- (3)CB
               \\         /
             (5)CA --- (4)CA
        ``` 
        The circular residue A and linear residue B are connected by a bond between
        (1)CA and the oxygen OA and (1)CB. By default, because OA originally is associated
        with residue A, only the bond OA --- (1)CB is returned. However, if `triplet=True`,
        the bond OA --- (1)CA is also returned, because the entire connecting "bridge" between residues
        A and B spans either bond around OA.
        >>> mol.infer_residue_connections(triplet=False)
        [("OA", "(1)CB")]
        >>> mol.infer_residue_connections(triplet=True)
        [("OA", "(1)CB"), ("OA", "(2)CA")]
        """
        bonds = structural.infer_residue_connections(
            self._base_struct, max_bond_length, triplet
        )
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

    def compute_angle(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        atom3: Union[str, int, bio.Atom.Atom],
    ):
        """
        Compute the angle between three atoms where atom2 is the middle atom.

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        atom3
            The third atom

        Returns
        -------
        float
            The angle in degrees
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        atom3 = self.get_atom(atom3)
        return structural.compute_angle(atom1, atom2, atom3)

    def compute_dihedral(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        atom3: Union[str, int, bio.Atom.Atom],
        atom4: Union[str, int, bio.Atom.Atom],
    ):
        """
        Compute the dihedral angle between four atoms

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        atom3
            The third atom
        atom4
            The fourth atom

        Returns
        -------
        float
            The dihedral angle in degrees
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        atom3 = self.get_atom(atom3)
        atom4 = self.get_atom(atom4)
        return structural.compute_dihedral(atom1, atom2, atom3, atom4)

    def make_atom_graph(self, locked: bool = True):
        """
        Generate an AtomGraph for the Molecule

        Parameters
        ----------
        locked : bool
            If True, the graph will also migrate the information on any locked bonds into the graph.
        """
        if len(self._bonds) == 0:
            warnings.warn(
                "No bonds found (yet), be sure to first apply or infer bonds to generate a graph"
            )
            return
        self._AtomGraph = AtomGraph(self._base_struct.id, self._bonds)
        if locked:
            for bond in self._locked_bonds:
                self._AtomGraph.lock_edge(*bond)
        return self._AtomGraph

    def make_residue_graph(self, detailed: bool = False, locked: bool = True):
        """
        Generate a ResidueGraph for the Molecule

        Parameters
        ----------
        detailed : bool
            If True the graph will include the residues and all atoms that form bonds
            connecting different residues. If False, the graph will only include the residues
            and their connections without factual bonds between any existing atoms.
        locked : bool
            If True, the graph will also migrate the information on any locked bonds into the graph.
            This is only relevant if detailed is True.
        """
        if len(self._bonds) == 0:
            warnings.warn(
                "No bonds found (yet), be sure to first apply or infer bonds to generate a graph"
            )
            return
        self._ResidueGraph = ResidueGraph.from_molecule(self, detailed, locked)
        return self._ResidueGraph

    def to_pdb(self, filename: str):
        """
        Write the molecule to a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        io = bio.PDBIO()
        io.set_structure(self._base_struct)
        io.save(filename)
