from collections import defaultdict
import warnings
import networkx as nx
import numpy as np

from copy import deepcopy

from typing import Union

import Bio.PDB as bio
from glycosylator.graphs import AtomGraph, ResidueGraph
import glycosylator.utils as utils
import glycosylator.structural as structural


class Molecule:
    """
    A molecule to add onto a protein.

    Parameters
    ----------
    structure : bio.PDB.Structure
        A biopython structure that stores the atomic structure
    root_atom : str or int or bio.PDB.Atom
        The id or the serial number of the root atom
        at which the molecule would be attached to a another
        structure such as protein scaffold or another Molecule.
    """

    def __init__(
        self,
        structure,
        root_atom: Union[str, int, bio.Atom.Atom] = None,
    ):
        self._id = structure.id
        self._base_struct = structure

        if len(structure.child_list) == 1:
            self._chain = structure.child_list[0]
        else:
            raise ValueError("Molecule class only supports structures with one model")
        if len(self._chain.child_list) == 1:
            self._chain = self._chain.child_list[0]
        else:
            raise ValueError(
                "Molecule class only supports structures with one chain/segment"
            )

        if root_atom:
            if not root_atom in self._chain.get_atoms():
                raise ValueError("The root atom is not in the structure")
        self._root_atom = root_atom

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

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        id: str = None,
        root_atom: Union[str, int] = None,
        add_hydrogens: bool = True,
    ):
        """
        Read a Molecule from a SMILES string

        Parameters
        ----------
        smiles : str
            The SMILES string
        id : str
            The id of the Molecule. By default the provided smiles string is used.
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        add_hydrogens : bool
            Whether to add hydrogens to the molecule
        """
        struct = structural.read_smiles(smiles, add_hydrogens)
        struct.id = id if id else smiles
        return cls(struct, root_atom)

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value
        self._base_struct.id = value

    @property
    def root_atom(self):
        return self._root_atom

    @root_atom.setter
    def root_atom(self, value):
        if value is None:
            self._root_atom = None
        else:
            self._root_atom = self.get_atom(value)

    @property
    def root_residue(self):
        if self._root_atom:
            return self.root_atom.get_parent()

    @property
    def structure(self):
        return self._base_struct

    @property
    def chain(self):
        return self._chain

    @property
    def residues(self):
        return list(self._chain.get_residues())

    @property
    def atoms(self):
        return list(self._chain.get_atoms())

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

    def get_residue(self, seqid: int = None, name: str = None):
        """
        Get a residue from the structure either based on its
        seqid or its residue name. Note, if there are multiple
        residues with the same name, all of them are returned
        in a list!

        Parameters
        ----------
        seqid : int
            The sequence id of the residue
        name : str
            The name of the residue

        Returns
        -------
        residue : bio.Residue.Residue or list
            The residue(s)
        """
        if seqid:
            return self._chain.child_list[seqid - 1]
        elif name:
            res = [i for i in self._base_struct.get_residues() if i.resname == name]
            if len(res) == 1:
                res = res[0]
            return res
        else:
            raise ValueError("Either seqid or name must be provided")

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
        for instance there are multiple 'C1' from different residues, the first one is returned.

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

    def get_root(self):
        return self.root_atom

    def set_root(self, atom):
        self.root_atom = atom

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
        return self._AtomGraph

    def make_residue_graph(self, detailed: bool = True, locked: bool = True):
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

    def reindex(self):
        """
        Reindex the atoms and residues in the molecule.
        You can use this method if you made substantial changes
        to the molecule and want to be sure that there are no gaps in the
        atom and residue numbering.
        """
        j = 0
        idx = 0

        residues = list(self.residues)
        parents = set([i.get_parent() for i in residues])
        for residue, parent in zip(residues, parents):
            parent.detach_child(residue.id)

        for residue, parent in zip(residues, parents):
            idx += 1
            residue.id = (residue.id[0], idx, *residue.id[2:])
            for atom in residue.get_atoms():
                j += 1
                atom.serial_number = j
                atom.set_parent(residue)
            parent.add(residue)

        self._AtomGraph = None
        self._ResidueGraph = None

    def adjust_indexing(self, other: "Molecule"):
        """
        Adjust another Molecule's residue and atom indexing to continue this Molecule's indexing.
        E.g. any residue with seqid 1 in the other molecule will be adjusted to have seqid 3 if this molecule
        already has two residues. And correspondingly, the atom numbering will be adjusted.

        Parameters
        ----------
        other : Molecule
            The other molecule to adjust
        """
        residues = other.residues
        parents = [i.get_parent() for i in residues]
        for residue, parent in zip(residues, parents):
            parent.detach_child(residue.id)

        rdx = len(self.residues)
        adx = len(self.atoms)
        for residue, parent in zip(residues, parents):
            rdx += 1
            residue.id = (residue.id[0], rdx, *residue.id[2:])
            for atom in residue.get_atoms():
                adx += 1
                atom.serial_number = adx
                atom.set_parent(residue)
            parent.add(residue)

    def add_residues(
        self,
        *residues: bio.Residue.Residue,
        adjust_seqid: bool = True,
        _copy: bool = False
    ):
        """
        Add residues to the structure

        Parameters
        ----------
        residues : bio.Residue.Residue
            The residues to add
        adjust_seqid : bool
            If True, the seqid of the residues is adjusted to
            match the current number of residues in the structure
            (i.e. a new residue can be given seqid 1, and it will be adjusted
            to the correct value of 3 if there are already two other residues in the molecule).
        _copy : bool
            If True, the residues are copied before adding them to the molecule.
            This is useful if you want to add the same residue to multiple molecules, while leaving
            them and their original parent structures intakt.
        """
        rdx = len(self.residues)
        adx = len(self.atoms)
        for residue in residues:
            if _copy:
                residue = residue.copy()
                residue.detach_parent()
            else:
                p = residue.get_parent()
                if p:
                    p.detach_child(residue.id)
            if adjust_seqid:
                rdx += 1
                residue.id = (residue.id[0], rdx, *residue.id[2:])
            for atom in residue.get_atoms():
                adx += 1
                atom.serial_number = adx
                atom.set_parent(residue)

            self._chain.add(residue)

    def remove_residues(self, *residues: Union[int, bio.Residue.Residue]):
        """
        Remove residues from the structure

        Parameters
        ----------
        residues : int or bio.Residue.Residue
            The residues to remove, either the object itself or its seqid
        """
        for residue in residues:
            if isinstance(residue, int):
                residue = self._chain.child_list[residue - 1]
            self.remove_atoms(*residue.child_list)
            self._chain.detach_child(residue.id)

    def add_atoms(self, *atoms: bio.Atom.Atom, residue=None, _copy: bool = False):
        """
        Add atoms to the structure. This will automatically adjust the atom's serial number to
        fit into the structure.

        Parameters
        ----------
        atoms : bio.Atom.Atom
            The atoms to add
        residue : int or str
            The residue to which the atoms should be added,
            this may be either the seqid or the residue name,
            if None the atoms are added to the last residue.
            Note, that if multiple identically named residues
            are present, the first one is chosen, so using the
            seqid is a safer option!
        _copy : bool
            If True, the atoms are copied and then added to the structure.
            This will leave the original atoms (and their parent structures) untouched.
        """
        if residue is not None:
            if isinstance(residue, int):
                residue = self._chain.child_list[residue - 1]
            elif isinstance(residue, str):
                residue = next(i for i in self._chain.child_list if i.resname == i)
            target = residue
        else:
            target = self._chain.child_list[-1]

        _max_serial = len(self.atoms)
        for atom in atoms:
            if _copy:
                atom = deepcopy(atom)
            _max_serial += 1
            atom.serial_number = _max_serial
            target.add(atom)

    def remove_atoms(self, *atoms: Union[int, str, tuple, bio.Atom.Atom]):
        """
        Remove one or more atoms from the structure

        Parameters
        ----------
        atoms
            The atoms to remove, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        for atom in atoms:
            atom = self.get_atom(atom)

            bonds = [i for i in self._bonds if atom in i]
            for bond in bonds:
                self._bonds.remove(bond)

            bonds = [i for i in self._locked_bonds if atom in i]
            for bond in bonds:
                self._locked_bonds.remove(bond)

            atom.get_parent().detach_child(atom.id)

        # reindex the atoms
        adx = 0
        for atom in self.atoms:
            adx += 1
            atom.serial_number = adx

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
    ):
        """
        Remove a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        self._bonds.remove(bond)
        if bond in self._locked_bonds:
            self._locked_bonds.remove(bond)

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
        self._locked_bonds.remove(bond)
        if both_ways and bond[::-1] in self._locked_bonds:
            self._locked_bonds.remove(bond[::-1])

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

    def infer_missing_atoms(self, _topology=None):
        """
        Infer missing atoms in the structure

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        """
        structural.fill_missing_atoms(self._base_struct, _topology)

    def apply_standard_bonds(self, _topology=None) -> list:
        """
        Get the standard bonds for the structure

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.apply_standard_bonds(self._base_struct, _topology)
        self._bonds.extend([b for b in bonds if b not in self._bonds])
        return bonds

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

    def attach(
        self,
        other: "Molecule",
        patch: Union[utils.abstract.AbstractPatch, str] = None,
        at: Union[int, str, tuple, bio.Atom.Atom] = None,
        other_at: Union[int, str, tuple, bio.Atom.Atom] = None,
        _topology=None,
    ):
        """
        Attach another structure to this one.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        patch : str or AbstractResidue
            A patch to apply when attaching. If none is given, an appropriate patch
            is searched based on the provided bonds at which to attach. However,
            if the patch is not found, an error is raised. Hence, it is safer to provide this parameter,
            in which case no `at`-bonds need to be specified as the patch holds this information. Either a
            patch object can be supplied or its id by which it can be retrieved from the used topology.
        at : tuple
            This molecule's bond that should be replaced by the bond attaching this and the other molecule.
        other_at : tuple
            The other molecule's bond that should be replaced by the bond attaching this and the other molecule.
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        """
        if at is None and patch is None:
            at = other.root_atom
        if other_at is None and patch is None:
            other_at = self.root_atom

        if not at or not other_at:
            raise ValueError("Cannot attach: no atoms or patch specified")

        if patch is not None:
            if isinstance(patch, str):
                if not _topology:
                    _topology = utils.get_default_topology()
                patch = _topology.get_patch(patch)

        structural.apply_patch(
            self._base_struct, other._base_struct, at, other_at, patch, _topology
        )

    def rotate_around_bond(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        angle: float,
        descendants_only: bool = False,
    ):
        """
        Rotate the structure around a bond

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        angle
            The angle to rotate by in degrees
        descendants_only
            Whether to only rotate the descendants of the bond
            (sensible only for linear molecules, or bonds that are not part of a circular structure).
        
        Examples
        --------
        For a molecule starting as:
        ```
                     OH
                    /
        (1)CH3 --- CH 
                    \\
                    CH2 --- (2)CH3
        ```
        we can rotate around the bond `(1)CH3 --- CH` by 180° using

        >>> import numpy as np
        >>> angle = np.radians(180) # rotation angle needs to be in radians
        >>> mol.rotate_around_bond("(1)CH3", "CH", angle)

        and thus achieve the following:
        ```
                     CH2 --- (2)CH3
                    /
        (1)CH3 --- CH 
                    \\
                     OH
        ```

        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if (atom1, atom2) in self._locked_bonds:
            raise RuntimeError("Cannot rotate around a locked bond")

        self.AtomGraph.rotate_around_edge(atom1, atom2, angle, descendants_only)

    def get_descendants(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
    ):
        """
        Get the atoms downstream of a bond. This will return the set
        of all atoms that are connected after the bond atom1-atom2 in the direction of atom2,
        the selection can be reversed by reversing the order of atoms (atom2-atom1).

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom

        Returns
        -------
        set
            A set of atoms

        Examples
        --------
        For a molecule
        ```
                     OH
                    /
        (1)CH3 --- CH 
                    \\
                    CH2 --- (2)CH3
        ```
        >>> mol.get_descendants("(1)CH3", "CH")
        {"OH", "CH2", "(2)CH3"}
        >>> mol.get_descendants("CH", "CH2")
        {"(2)CH3"}
        >>> mol.get_descendants("CH2", "CH")
        {"OH", "(1)CH3"}
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        return self.AtomGraph.get_descendants(atom1, atom2)

    def get_neighbors(
        self,
        atom: Union[int, str, tuple, bio.Atom.Atom],
        n: int = 1,
        mode: str = "upto",
    ):
        """
        Get the neighbors of an atom.

        Parameters
        ----------
        atom
            The atom
        n
            The number of bonds that may separate the atom from its neighbors.
        mode
            The mode to use. Can be "upto" or "at". If `upto`, all neighbors that are at most `n` bonds away
            are returned. If `at`, only neighbors that are exactly `n` bonds away are returned.

        Returns
        -------
        set
            A set of atoms


        Examples
        --------
        For a molecule
        ```
                     O --- (2)CH2
                    /         \\
        (1)CH3 --- CH          OH
                    \\
                    (1)CH2 --- (2)CH3
        ```
        >>> mol.get_neighbors("(2)CH2", n=1)
        {"O", "OH"}
        >>> mol.get_neighbors("(2)CH2", n=2, mode="upto")
        {"O", "OH", "CH"}
        >>> mol.get_neighbors("(2)CH2", n=2, mode="at")
        {"CH"}
        """
        atom = self.get_atom(atom)
        return self.AtomGraph.get_neighbors(atom, n, mode)

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

    # def _get_atom(self, _atom):
    #     """
    #     Get an atom based on its id, serial number or check if it is available...
    #     """
    #     if isinstance(_atom, str):
    #         _atom = self.get_atoms(_atom, by="id")
    #         if len(_atom) == 0:
    #             raise ValueError(f"Atom {_atom} not found in structure")

    #     elif isinstance(_atom, int):

    #     elif isinstance(_atom, tuple):
    #         if not _atom in self._full_id_atoms:
    #             raise ValueError(f"Atom {_atom} not found in structure")
    #         return self._full_id_atoms[_atom]

    #     elif isinstance(_atom, bio.Atom.Atom):
    #         if _atom in self._base_struct.get_atoms():
    #             return _atom
    #         else:
    #             warnings.warn(f"Atom {_atom} not found in structure")
    #             return _atom

    #     else:
    #         raise TypeError(f"Invalid type for atom: {type(_atom)}")


#         _validate_atom_group(atom_group)
#         self.atom_group: prody.AtomGroup = atom_group

#         self.root_atom: prody.Atom = root_atom
#         self.root_residue: prody.Residue = self.atom_group[self.root_atom.getChid(), self.root_atom.getResnum()]

#         self.atom_graph: AtomGraph = AtomGraph(self.atom_group, self.root_atom)
#         self.residue_graph: ResidueGraph = ResidueGraph.from_AtomGraph(self.atom_graph)

#         self.guess_angles()
#         self.guess_dihedrals()
#         self.guess_torsionals()

#     @classmethod
#     def from_PDB(cls, pdb: str, root_atom_serial: int, **kwargs):
#         atom_group = prody.parsePDB(pdb, **kwargs)
#         serials = list(atom_group.getSerials())
#         root_atom_index = serials.index(root_atom_serial)
#         return Molecule(atom_group, atom_group[root_atom_index])

#     @property
#     def glycan_topology(self):
#         return self.residue_graph.to_glycan_topology()

#     def write_PDB(self, filename: str, selection: str = "all"):
#         prody.writePDB(filename, self.atom_group.select(selection).toAtomGroup())

#     def guess_angles(self):
#         """Searches for all angles in a molecule based on the connectivity"""
#         self.angles = []
#         for node in self.atom_graph.nodes():
#             self.angles.extend(_find_paths(self.atom_graph, node, 2))

#     def guess_dihedrals(self):
#         """Searches for all dihedrals in a molecule based on the connectivity"""
#         self.dihedrals = []
#         for node in self.atom_graph.nodes():
#             self.dihedrals.extend(_find_paths(self.atom_graph, node, 3))

#     def guess_torsionals(self, hydrogens=True):
#         """Builds a list with all the torsional angles that can rotate
#         Parameters:
#             hydrogens: include torsional involving terminal hydrogens
#         Initializes:
#             torsionals: a list of serial number of atom defining a torsional angle (quadruplet)
#         """
#         cycles = nx.cycle_basis(self.atom_graph.to_undirected(as_view=True))
#         # TODO: can cycle_id just be a flat set of all atoms in any cycles?
#         cycle_id = {atom: i for i, cycle in enumerate(cycles) for atom in cycle}

#         elements = self.atom_group.getElements()
#         torsionals = []
#         for dihedral in self.dihedrals:
#             if dihedral[1] in cycle_id and dihedral[2] in cycle_id:
#                 # skip dihedral if both middle atoms in a cycle
#                 continue
#             if not hydrogens:
#                 if elements[dihedral[0]] == "H" or elements[dihedral[-1]] == "H":
#                     # skip dihedral if either of outer atoms are hydrogen
#                     continue
#             # use the direction of bond_graph to orient dihedral
#             # to then check for uniqueness
#             # i.e. ABCD and DCBA are same and should not be repeated
#             if self.atom_graph.has_edge(dihedral[2], dihedral[1]):
#                 dihedral.reverse()
#             if dihedral not in torsionals:
#                 torsionals.append(dihedral)

#         self.torsionals = torsionals

#     def measure_dihedral_angle(self, dihedral_atoms):
#         """Calculates dihedral angle for 4 atoms.
#         Parameters:
#             torsional: list of atom serial numbers
#         Returns:
#             angle: dihedral angle in degrees
#         """

#         idx = np.argsort(dihedral_atoms)
#         vec_sel = self.atom_group.select(f"index {' '.join(map(str, dihedral_atoms))}")
#         c0, c1, c2, c3 = vec_sel.getCoords()[idx, :]

#         q1 = c1 - c0
#         q2 = c2 - c1
#         q3 = c3 - c2

#         q1xq2 = np.cross(q1, q2)
#         q2xq3 = np.cross(q2, q3)

#         n1 = q1xq2 / np.sqrt(np.dot(q1xq2, q1xq2))
#         n2 = q2xq3 / np.sqrt(np.dot(q2xq3, q2xq3))

#         u1 = n2
#         u3 = q2 / (np.sqrt(np.dot(q2, q2)))
#         u2 = np.cross(u3, u1)

#         cos_theta = np.dot(n1, u1)
#         sin_theta = np.dot(n1, u2)
#         angle = np.degrees(-np.arctan2(sin_theta, cos_theta))

#         return angle

#     def get_all_torsional_angles(self):
#         """Computes all the torsional angles of the molecule
#         Return:
#             angles: angles in degrees
#         """
#         return [self.measure_dihedral_angle(torsional) for torsional in self.torsionals]

#     def rotate_bond(self, torsional, theta, absolute=False):
#         """Rotate the molecule around a torsional angle. Atom affected are in direct graph.
#         Parameters:
#             torsional: index of torsional angle (in torsionals) or list of indices of atoms defining the torsional angle.
#             theta: amount (degrees)
#             absolute: defines if theta is the increment or the absolute value of the angle
#         Returns:
#             c_angle: angle before the rotation in degrees
#         """
#         if type(torsional) == int:
#             torsional = self.torsionals[torsional]
#         elif torsional not in self.torsionals:
#             raise ValueError("Invalid Torsional")

#         atoms = []
#         a1 = torsional[-2]

#         atoms = nx.descendants(self.atom_graph, a1)
#         sel = self.atom_group.select(f"index {' '.join(map(str, atoms))}")
#         t = torsional[1:-1]
#         v1, v2 = self.atom_group.select(f"index {' '.join(map(str, t))}").getCoords()[np.argsort(t), :]
#         axis = v2 - v1
#         c_angle = 0.0
#         if absolute:
#             c_angle = self.measure_dihedral_angle(torsional)
#             theta = theta - c_angle

#         coords = sel.getCoords()
#         M = _rotation_matrix(axis, np.radians(theta))
#         coords = M.dot(coords.transpose())
#         sel.setCoords(coords.transpose() + v2 - np.dot(M, v2))


# def _validate_atom_group(atom_group: prody.AtomGroup):
#     """validates that the AtomGroup consists of a single physical molecule"""
#     pass
#     # if atom_group.numChains() != 1 or atom_group.numSegments() != 1:
#     #     raise TooManyChains()


# def _find_paths(G, node, length, excludeSet=None):
#     """Finds all paths of a given length
#     Parameters:
#         G: graph (netwrokx)
#         node: starting node
#         length: length of path
#         excludedSet: set
#     Returns:
#         paths: list of all paths of a length starting from node
#     """
#     if excludeSet == None:
#         excludeSet = {node}
#     else:
#         excludeSet.add(node)

#     if length == 0:
#         return [[node]]
#     paths = [[node] + path for neighbor in G.neighbors(node) if neighbor not in excludeSet for path in _find_paths(G, neighbor, length - 1, excludeSet)]
#     excludeSet.remove(node)
#     return paths


if __name__ == "__main__":

    mol = Molecule.from_pdb("support/examples/MAN.pdb")
    mol.infer_bonds()

    # d = mol._get_descendents(mol.get_atoms(id="HO3"), mol.get_atoms(id="O3"))
    # assert len(d) == 22

    # d = mol._get_descendents(mol.get_atoms(id="C5"), mol.get_atoms(id="C6"))
    # assert len(d) == 4

    # d = mol._get_descendents(mol.get_atoms(id="C4"), mol.get_atoms(id="C5"))

    # # mol.rotate_around_bond("C5", "C6", np.radians(25), True)

    mol.remove_atoms("C5", "HO4", "C1", "H61")
