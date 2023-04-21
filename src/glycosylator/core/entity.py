"""
The base class for classes storing biopython structures
"""

from copy import deepcopy
from typing import Union
import warnings

import Bio.PDB as bio

import glycosylator.structural as structural
import glycosylator.utils as utils
import glycosylator.graphs as graphs


class BaseEntity:
    def __init__(self, structure):
        self._base_struct = structure
        self._id = structure.id

        self._bonds = []
        # self._locked_bonds = set()

        self._AtomGraph = graphs.AtomGraph(self.id, [])

        # let the molecule store a patch to use for attaching other
        # molecules to it
        self._patch = None

        self._working_chain = None
        self._root_atom = None
        self._attach_residue = None

    @classmethod
    def from_pdb(
        cls,
        filename: str,
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
        return cls(struct)

    @classmethod
    def load(cls, filename: str):
        """
        Load a Molecule from a pickle file

        Parameters
        ----------
        filename : str
            Path to the file
        """
        obj = utils.load_pickle(filename)
        if obj.__class__.__name__ != cls.__name__:
            raise TypeError(
                f"Object loaded from {filename} is not a {cls.__name__} but a {type(obj)}"
            )
        return obj

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
    def patch(self):
        return self._patch

    @patch.setter
    def patch(self, value):
        if value is None:
            self._patch = None
        else:
            self._patch = self.get_residue(value)

    @property
    def bonds(self):
        return self._bonds
        # return list(self._AtomGraph.edges)

    @bonds.setter
    def bonds(self, value):
        if value is None or len(value) == 0:
            self._bonds = []
            self._AtomGraph.clear_edges()
        self._bonds = value
        self._AtomGraph.clear_edges()
        self._AtomGraph.add_edges_from(value)

    @property
    def locked_bonds(self):
        return self._AtomGraph._locked_edges

    @locked_bonds.setter
    def locked_bonds(self, value):
        if value is None or len(value) == 0:
            self._AtomGraph._locked_edges = set()
        elif isinstance(value, set):
            self._AtomGraph._locked_edges = value
        else:
            raise TypeError("locked_bonds must be a set")

    @property
    def chains(self):
        return sorted(list(self._base_struct.get_chains()))

    @property
    def residues(self):
        return sorted(list(self._base_struct.get_residues()))

    @property
    def atoms(self):
        return sorted(self._AtomGraph.nodes)
        # return list(self._base_struct.get_atoms())

    @property
    def atom_triplets(self):
        """
        Compute triplets of three consequtively bonded atoms
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_triplets(self.bonds)

    @property
    def atom_quartets(self):
        """
        Compute quartets of four consequtively bonded atoms
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_quartets(self.bonds)

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
    def _chain(self):
        if self._working_chain is None:
            return self.chains[-1]
        else:
            return self._working_chain

    @_chain.setter
    def _chain(self, value):
        self._working_chain = value

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
    def attach_residue(self):
        return self._attach_residue

    @attach_residue.setter
    def attach_residue(self, value):
        if value is None:
            self._attach_residue = None
        else:
            self._attach_residue = self.get_residue(value)

    def get_root(self):
        return self.root_atom

    def set_root(self, atom):
        self.root_atom = atom

    def save(self, filename: str):
        """
        Save the molecule object to a pickle file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        utils.save_pickle(self, filename)

    def view(self):
        """
        Open a browser window to view the molecule in 3D
        """
        utils.visual.MoleculeViewer3D(self).show()

    def get_attach_residue(self):
        return self._attach_residue

    def set_attach_residue(self, residue: Union[int, bio.Residue.Residue] = None):
        """
        Set the residue that is used for attaching other molecules to this one.

        Parameters
        ----------
        residue
            The residue to be used for attaching other molecules to this one
        """
        if isinstance(residue, int):
            residue = self.get_residue(residue)
        self._attach_residue = residue

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
        if (atom1, atom2) in self.locked_bonds:
            raise RuntimeError("Cannot rotate around a locked bond")

        self._AtomGraph.rotate_around_edge(atom1, atom2, angle, descendants_only)

    def get_ancestors(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
    ):
        """
        Get the atoms upstream of a bond. This will return the set
        of all atoms that are connected before the bond atom1-atom2 in the direction of atom1,
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
        >>> mol.get_ancestors("(1)CH3", "CH")
        set()
        >>> mol.get_ancestors("CH", "CH2")
        {"(1)CH3", "OH"}
        >>> mol.get_ancestors("CH2", "CH")
        {"(2)CH3"}
        """
        return self.get_descendants(atom2, atom1)

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

        return self._AtomGraph.get_descendants(atom1, atom2)

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
        return self._AtomGraph.get_neighbors(atom, n, mode)

    def reindex(self, start_resid: int = 1, start_atomid: int = 1):
        """
        Reindex the atoms and residues in the molecule.
        You can use this method if you made substantial changes
        to the molecule and want to be sure that there are no gaps in the
        atom and residue numbering.

        Parameters
        ----------
        start_resid : int
            The starting residue id
        start_atomid : int
            The starting atom id
        """
        j = start_atomid
        idx = start_resid

        residues = list(self.residues)
        parents = [i.get_parent() for i in residues]
        for residue, parent in zip(residues, parents):
            parent.detach_child(residue.id)

        for residue, parent in zip(residues, parents):
            # for residue in self.residues:
            residue.id = (residue.id[0], idx, *residue.id[2:])
            idx += 1
            for atom in residue.get_atoms():
                atom.serial_number = j
                j += 1
                atom.set_parent(residue)
            parent.add(residue)

    def get_residues(self):
        return self._base_struct.get_residues()

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
            elif isinstance(atoms[0], bio.Atom.Atom):
                return atoms
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
            if atom in self.atoms:
                return atom
            else:
                return self.get_atom(atom.serial_number, by="serial")

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
        atom1: Union[int, str, tuple, bio.Atom.Atom, bio.Residue.Residue],
        atom2: Union[int, str, tuple, bio.Atom.Atom] = None,
        either_way: bool = True,
        residue_internal: bool = True,
    ):
        """
        Get one or multiple bonds from the molecule. If only one atom is provided, all bonds
        that are connected to that atom are returned.

        Parameters
        ----------
        atom1
            The atom id, serial number or full_id tuple of the first atom.
            This may also be a residue, in which case all bonds between atoms in that residue are returned.
        atom2
            The atom id, serial number or full_id tuple of the second atom
        either_way : bool
            If True, the order of the atoms does not matter, if False, the order of the atoms
            does matter. By setting this to false, it is possible to also search for bonds that have
            a specific atom in position 1 or 2 depending on which argument was set, while leaving the other input as none.
        residue_internal : bool
            If True, only bonds where both atoms are in the given residue (if `atom1` is a residue) are returned.
            If False, all bonds where either atom is in the given residue are returned.

        Returns
        -------
        bond : list
            The bond(s)
        """
        if isinstance(atom1, bio.Residue.Residue):
            if residue_internal:
                return [
                    i
                    for i in self.bonds
                    if i[0].get_parent() == atom1 and i[1].get_parent() == atom1
                ]
            else:
                return [
                    i
                    for i in self.bonds
                    if i[0].get_parent() == atom1 or i[1].get_parent() == atom1
                ]

        if atom1:
            atom1 = self.get_atom(atom1)
        if atom2:
            atom2 = self.get_atom(atom2)

        if atom1 and atom2:
            if either_way:
                return [i for i in self.bonds if atom1 in i and atom2 in i]
            else:
                return [
                    i
                    for i in self.bonds
                    if atom1 in i and atom2 in i and i[0] == atom1 and i[1] == atom2
                ]
        elif atom1:
            if either_way:
                return [i for i in self.bonds if atom1 in i]
            else:
                return [i for i in self.bonds if atom1 in i and i[0] == atom1]
        elif atom2:
            if either_way:
                return [i for i in self.bonds if atom2 in i]
            else:
                return [i for i in self.bonds if atom2 in i and i[1] == atom2]
        else:
            raise ValueError("No atom provided")

    def get_residue(
        self,
        residue: Union[int, str, tuple, bio.Residue.Residue],
        by: str = None,
        chain=None,
    ):
        """
        Get a residue from the structure either based on its
        name, serial number or full_id. Note, if multiple residues match the requested criteria,
        for instance there are multiple 'MAN' from different chains, only the first one is returned.

        Parameters
        ----------
        residue
            The residue id, seqid or full_id tuple
        by : str
            The type of parameter to search for. Can be either 'name', 'seqid' or 'full_id'
            By default, this is inferred from the datatype of the residue parameter.
            If it is an integer, it is assumed to be the sequence identifying number,
            if it is a string, it is assumed to be the residue name and if it is a tuple, it is assumed
            to be the full_id.
        chain : str
            Further restrict to a residue from a specific chain.
        """
        if isinstance(residue, bio.Residue.Residue):
            if residue in self.residues:
                return residue
            else:
                return self.get_residue(residue.id[1], by="seqid", chain=chain)

        if by is None:
            if isinstance(residue, int):
                by = "seqid"
            elif isinstance(residue, str):
                by = "name"
            elif isinstance(residue, tuple):
                by = "full_id"
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'name', 'seqid' or 'full_id'"
                )

        if by == "name":
            residue = (
                i for i in self._base_struct.get_residues() if i.resname == residue
            )
        elif by == "seqid":
            residue = (
                i for i in self._base_struct.get_residues() if i.id[1] == residue
            )
        elif by == "full_id":
            residue = (
                i for i in self._base_struct.get_residues() if i.full_id == residue
            )
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'name', 'seqid' or 'full_id'"
            )
        if chain is not None:
            chain = self.get_chain(chain)
            residue = (i for i in residue if i.get_parent() == chain)
        return next(residue)

    def get_chain(self, chain: str):
        """
        Get a chain from the structure either based on its
        name.

        Parameters
        ----------
        chain
            The chain id
        Returns
        -------
        chain : bio.Chain.Chain
            The chain
        """
        if isinstance(chain, bio.Chain.Chain):
            if chain in self.chains:
                return chain
            else:
                return self.get_chain(chain.id)
        return next(i for i in self.chains if i.id == chain)

    def add_residues(
        self,
        *residues: bio.Residue.Residue,
        adjust_seqid: bool = True,
        _copy: bool = False,
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
            p = residue.get_parent()
            if p:
                p.detach_child(residue.id)
            if _copy and p is not None:
                r = deepcopy(residue)
                p.add(residue)
                residue = r

            rdx += 1
            if adjust_seqid and residue.id[1] != rdx:
                residue.id = (residue.id[0], rdx, *residue.id[2:])

            self._chain.add(residue)
            residue._generate_full_id()

            for atom in residue.child_list:
                if adjust_seqid:
                    adx += 1
                    atom.set_serial_number(adx)
                atom.set_parent(residue)
                self._AtomGraph.add_node(atom)

    def remove_residues(self, *residues: Union[int, bio.Residue.Residue]) -> list:
        """
        Remove residues from the structure

        Parameters
        ----------
        residues : int or bio.Residue.Residue
            The residues to remove, either the object itself or its seqid

        Returns
        -------
        list
            The removed residues
        """
        _residues = []
        for residue in residues:
            if isinstance(residue, int):
                residue = self._chain.child_list[residue - 1]

            for atom in residue.child_list:
                self.purge_bonds(atom)
                self._AtomGraph.remove_node(atom)

            self._chain.detach_child(residue.id)
            _residues.append(residue)
        return _residues

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
            atom.set_serial_number(_max_serial)
            target.add(atom)
            self._AtomGraph.add_node(atom)

    def remove_atoms(self, *atoms: Union[int, str, tuple, bio.Atom.Atom]) -> list:
        """
        Remove one or more atoms from the structure

        Parameters
        ----------
        atoms
            The atoms to remove, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        Returns
        -------
        list
            The removed atoms
        """
        _atoms = []
        for atom in atoms:
            atom = self.get_atom(atom)
            self._AtomGraph.remove_node(atom)
            self.purge_bonds(atom)
            atom.get_parent().detach_child(atom.id)
            _atoms.append(atom)

        # reindex the atoms
        adx = 0
        for atom in self.atoms:
            adx += 1
            atom.serial_number = adx

        return _atoms

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

        if (atom1, atom2) not in self._bonds:
            self._bonds.append((atom1, atom2))
        if not self._AtomGraph.has_edge(atom1, atom2):
            self._AtomGraph.add_edge(atom1, atom2)

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
        if either_way:
            self.remove_bond(atom1, atom2, either_way=False)
            self.remove_bond(atom2, atom1, either_way=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        if (atom1, atom2) in self._bonds:
            self._bonds.remove((atom1, atom2))
        if self._AtomGraph.has_edge(atom1, atom2):
            self._AtomGraph.remove_edge(atom1, atom2)
        if (atom1, atom2) in self.locked_bonds:
            self.locked_bonds.remove((atom1, atom2))

    def purge_bonds(self, atom: Union[int, str, bio.Atom.Atom] = None):
        """
        Remove all bonds connected to an atom

        Parameters
        ----------
        atom
            The atom to remove the bonds from, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms. If None, all bonds
            are removed.
        """
        if atom is None:
            self.bonds = []
            return

        atom = self.get_atom(atom)
        bonds = (i for i in self.bonds if atom in i)
        for bond in bonds:
            self.remove_bond(*bond)

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
        self.locked_bonds = set(self.bonds)
        if both_ways:
            self.locked_bonds.update(b[::-1] for b in self.bonds)

    def unlock_all(self):
        """
        Unlock all bonds in the structure
        """
        self._AtomGraph.unlock_all()

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
        if both_ways:
            self.lock_bond(atom2, atom1, both_ways=False)
            self.lock_bond(atom2, atom1, both_ways=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._AtomGraph.lock_edge(atom1, atom2)

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
        if both_ways:
            self.unlock_bond(atom1, atom2, both_ways=False)
            self.unlock_bond(atom2, atom1, both_ways=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if (atom1, atom2) in self._AtomGraph.edges:
            self._AtomGraph.unlock_edge(atom1, atom2)
        if both_ways and (atom2, atom1) in self._AtomGraph.edges:
            self._AtomGraph.unlock_edge(atom2, atom1)

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
        return bond in self.locked_bonds

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
        to_add = [b for b in bonds if b not in self.bonds]
        self._bonds.extend(to_add)
        self._AtomGraph.add_edges_from(to_add)
        return bonds

    def get_residue_connections(self, triplet: bool = True):
        """
        Get bonds between atoms that connect different residues in the structure
        This method is different from `infer_residue_connections` in that it works
        with the already present bonds in the molecule instead of computing new ones.

        Parameters
        ----------
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
        set
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = set(i for i in self.bonds if i[0].get_parent() != i[1].get_parent())

        if triplet:
            _new = set()
            for atom1, atom2 in bonds:
                neighs = self.get_neighbors(atom1)
                neighs.remove(atom2)
                neighs -= set(i for i in neighs if i.element == "H")
                if len(neighs) == 1:
                    neigh = neighs.pop()
                    if neigh.get_parent() == atom1.get_parent():
                        _new.add((atom1, neigh))
                    continue

                neighs = self.get_neighbors(atom2)
                neighs.remove(atom1)
                neighs -= set(i for i in neighs if i.element == "H")
                if len(neighs) == 1:
                    neigh = neighs.pop()
                    if neigh.get_parent() == atom2.get_parent():
                        _new.add((atom2, neigh))
                    continue

            bonds.update(_new)
        return bonds

    def infer_residue_connections(
        self, bond_length: Union[float, tuple] = None, triplet: bool = True
    ) -> list:
        """
        Infer bonds between atoms that connect different residues in the structure

        Parameters
        ----------
        bond_length : float or tuple
            If a float is given, the maximum distance between atoms to consider them bonded.
            If a tuple, the minimal and maximal distance between atoms. If None, the default value is min 0.8 Angstrom, max 1.6 Angstroms.
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
        `(1)CA` and the oxygen `OA` and `(1)CB`. By default, because OA originally is associated
        with residue A, only the bond `OA --- (1)CB` is returned. However, if `triplet=True`,
        the bond `OA --- (1)CA` is also returned, because the entire connecting "bridge" between residues
        A and B spans either bond around `OA`.
        >>> mol.infer_residue_connections(triplet=False)
        [("OA", "(1)CB")]
        >>> mol.infer_residue_connections(triplet=True)
        [("OA", "(1)CB"), ("OA", "(2)CA")]
        """
        bonds = structural.infer_residue_connections(
            self._base_struct, bond_length, triplet
        )
        self._bonds.extend(bonds)
        self._AtomGraph.add_edges_from(bonds)
        return bonds

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
        to_add = [b for b in bonds if b not in self.bonds]
        self._bonds.extend(to_add)
        self._AtomGraph.add_edges_from(to_add)
        return bonds

    def autolabel(self, _compounds=None):
        """
        Relabel the atom ids in the current structure to match the standard labelling
        according to available compounds. This is useful if you want to use some pre-generated
        PDB file that may have used a different labelling scheme for atoms. Using this method,
        the labels are adjusted which makes it possible to use standard bond inference and patching.

        Parameters
        ----------
        _compounds : PDBECompounds
            The compounds object to use for relabelling. If None, the default compounds object is used.
        """
        if _compounds is None:
            _compounds = utils.defaults.get_default_compounds()
        _compounds.relabel_atoms(self._base_struct)

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

        Returns
        -------
        AtomGraph
            The generated graph
        """
        graph = deepcopy(self._AtomGraph)
        if not locked:
            graph.unlock_all()
        return graph

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
        graph = graphs.ResidueGraph.from_molecule(self, detailed, locked)
        return graph

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

    def infer_missing_atoms(self, _topology=None, _compounds=None):
        """
        Infer missing atoms in the structure based on a reference topology or compounds database.
        By default, if a residue is not available in the topology, the compounds database is used.

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        _compounds
            A specific compounds object to use for referencing.
            If None, the default compounds object is used.
        """
        structural.fill_missing_atoms(self._base_struct, _topology, _compounds)

    def _direct_bonds(self, bonds: list, by: str = "serial", save: bool = True):
        """
        Sort a given list of bonds such that the first atom is the "earlier" atom
        in the bond.

        Parameters
        ----------
        bonds : list
            A list of tuples of atom pairs that are bonded and connect different residues
        by : str
            The attribute to sort by. Can be either "serial", "resid" or "root".
            In the case of "serial", the bonds are sorted by the serial number of the first atom.
            In the case of "resid", the bonds are sorted by the residue number of the first atom.
            In this case, bonds connecting atoms from the same residue are sorted by the serial number of the first atom.
            In the case of "root" the bonds are sorted based on the graph distances to the root atom,
            provided that the root atom is set (otherwise the atom with serial 1 is used).
        save : bool
            Whether to save the sorted bonds as the new bonds in the Molecule.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded and connect different residues
        """
        if by == "serial":
            directed = [
                bond if bond[0].serial_number < bond[1].serial_number else bond[::-1]
                for bond in bonds
            ]
        elif by == "root":
            root = self._root_atom
            if root is None:
                root = self.get_atom(1)
            directed = self._AtomGraph.direct_edges(root, bonds)
        elif by == "resid":
            directed = [
                bond if not should_revert(bond) else bond[::-1] for bond in bonds
            ]
        if save:
            for old, new in zip(bonds, directed):
                self.remove_bond(*old)
                self.add_bond(*new)
        return directed

    def __mod__(self, patch):
        """
        Add a patch to the molecule using the % operator (i.e. mol % patch)
        """
        self.set_patch(patch)
        return self

    def __xor__(self, atom):
        """
        Set the root atom for the molecule using the ^ operator (i.e. mol ^ atom)
        """
        self.set_root(atom)
        return self

    def __matmul__(self, residue) -> "Molecule":
        """
        Set the residue at which the molecule should be attached to another molecule
        using the @ operator (i.e. mol @ 1, for residue 1)
        """
        self.set_attach_residue(residue)
        return self


should_revert = lambda bond: (
    bond[0].parent.id[1] > bond[1].parent.id[1]
    or (
        bond[0].parent.id[1] == bond[1].parent.id[1]
        and bond[0].serial_number > bond[1].serial_number
    )
)
"""
A helper function to sort bonds and determine whether they should be reversed
"""

if __name__ == "__main__":
    # f = "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    # e = BaseEntity.from_pdb(f)
    # e.infer_bonds(restrict_residues=False)
    # e.get_residue_connections()
    pass