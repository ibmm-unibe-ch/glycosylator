"""
The `Molecule` class is a wrapper around a biopython structure and a core part 
of glycosylator functionality. It provides a convenient interface to available
structures in PDB format, as well as an easy way to create new molecules from
SMILES strings or reference compounds from the PDBECompounds database.

To make a new molecule, the easiest way is to use the `Molecule.from_compound` classmethod.
This method is linked to the PDBECompounds database, which contains a large number of reference
compounds (currently restricted to sugars, due to the nature of "glycosylator", however). The database
is queried using the `by` parameter, which can be one of the following:
- "id" for the PDB id
- "name" for the name of the compound (must match any known synonym of the iupac name)
- "formula" for the chemical formula (usually ambiguous and will therefore often raise an error)
- "smiles" for the SMILES string (also accepts InChI)

.. code-block:: python

    from glycosylator import Molecule

    # create a new galactose molecule
    gal = Molecule.from_compound("GAL") # use the PDB id

Alternatively, a molecule can be created from a PDB file using the `Molecule.from_pdb` classmethod. 
In this case the molecule may be as large as the PDB format supports (currently 99999 atoms).

.. code-block:: python

    my_glycan = Molecule.from_pdb("my_glycan.pdb")

However, unless loaded from a PDB file, each `Molecule` starts of with a single chain and residue,
specifying a single compound, such as Galactose or Glucose. Starting from there, multiple molecules
can be patched together to form complex structures. 

.. code-block:: python

    # create a glucose molecule
    glc = Molecule.from_compound("GLC")

    # create a galactose molecule
    gal = Molecule.from_compound("GAL")

    # now attach the galactose to the glucose
    # using a 1-4 beta linkage
    glc.attach(gal, "14bb")

    # Now the glc object is a disaccharide holding two residues, one glucose and one galactose

In the above example, the `attach` method is used to attach the galactose molecule to the glucose, but for those
among us who prefer a more shorty syntax, the `+` operator will do the same thing. In this case we need
to specify the linkage type (known as `patch`) as an attribute to the receiving molecule using `set_patch` first.

.. code-block:: python

    # specify that incoming molecules shall be
    # attached using a 1-4 beta linkage
    glc.set_patch("14bb")

    # now attach the galactose to the glucose
    glc += gal

For those who still want a shorter syntax, `set_patch` is represeted by the `%` opreator, so the above example
can be written as

.. code-block:: python

    glc = glc % "14bb" + gal

    # or (equivalently)
    # glc % "14bb" 
    # glc += gal 

.. note::

    The `attach` method and the `+` operator are not commutative, i.e. `glc + gal` is not the same as `gal + glc`!
    Also, the `attach` method (and `+=` operator) will modify the accepting molecule object inplace, while the `+` operator will
    return a new object containing both molecules.
    Hence, if you are working with very large structures, it is recommended to use the `attach` method instead of the `+` operator,
    to avoid copying the structures again and again. 
    
Naturally, if we can _add_  different molecules, we may also add the _same_ molecule multiple times. Which is exactly what
what the `repeat` method does, or its shorthand-form, the `*` operator.

.. code-block:: python

    # create cellulose from glucose
    cel = glc.repeat(10, "14bb")

    # or (equivalently)
    cel = glc % "14bb" * 10

Molecules are mutable objects that allow not only easy incorporation of new `biopython` atoms and residues,
or removal of old ones, but also spacial transfomations such as selective rotations around specific bonds.

The Molecule class is also the main interface to the `AtomGraph` and `ResidueGraph` classes, which
provide a convenient way to access the underlying structure as a graph. Especially, the ResidueGraph
class provides a convenient way to access the overall structure of larger molecules such as glycans.
"""

from copy import deepcopy

from typing import Union

import Bio.PDB as bio

import glycosylator.core.entity as entity
import glycosylator.utils as utils
import glycosylator.structural as structural


class Molecule(entity.BaseEntity):
    """
    A molecule to add onto a scaffold.

    Parameters
    ----------
    structure : bio.PDB.Structure
        A biopython structure that stores the atomic structure
    root_atom : str or int or bio.PDB.Atom
        The id or the serial number of the root atom
        at which the molecule would be attached to a another
        structure such as protein scaffold or another Molecule.
    model : int
        The model to use from the structure. Defaults to 0. This may be any
        valid identifier for a model in the structure, such as an integer or string.
    chain : str
        The chain to use from the structure. Defaults to the first chain in the structure.
    """

    def __init__(
        self,
        structure,
        root_atom: Union[str, int, bio.Atom.Atom] = None,
        model: int = 0,
        chain: str = None,
    ):
        super().__init__(structure)

        model = structure.child_dict[model]

        if len(model.child_list) == 0:
            raise ValueError("The model is empty")
        elif len(model.child_list) == 1:
            self._chain = model.child_list[0]
        else:
            if chain:
                self._chain = model.child_dict[chain]
            else:
                self._chain = model.child_list[0]

        if root_atom:
            if not root_atom in self._chain.get_atoms():
                raise ValueError("The root atom is not in the structure")
        self._root_atom = root_atom

        # let the molecule also store the residue at which it should be attached to
        # another molecule
        self._attach_residue = None

    @classmethod
    def empty(cls, id: str = None):
        """
        Create an empty Molecule object

        Parameters
        ----------
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        structure = structural.make_empty_structure(id)
        return cls(structure)

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        root_atom: Union[str, int] = None,
        id: str = None,
        model: int = 0,
        chain: str = None,
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
        model : int
            The model to use from the structure. Defaults to 0. This may be any
            valid identifier for a model in the structure, such as an integer or string.
        chain : str
            The chain to use from the structure. Defaults to the first chain in the structure.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        return cls(struct, root_atom, model=model, chain=chain)

    @classmethod
    def from_compound(
        cls,
        compound: str,
        by: str = "id",
        root_atom: Union[str, int] = None,
    ) -> "Molecule":
        """
        Create a Molecule from a reference compound from the PDBECompounds database

        Parameters
        ----------
        compound : str
            The compound to search for
        by : str
            The field to search by. This can be
            - "id" for the PDB id
            - "name" for the name of the compound (must match any known synonym of the iupac name)
            - "formula" for the chemical formula
            - "smiles" for the SMILES string (also accepts InChI)
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        """
        mol = utils.defaults.get_default_compounds().get(compound, by=by)
        if isinstance(mol, list):
            raise ValueError(
                f"Multiple compounds found using '{by}={compound}', choose any of these ids specifically {[i.id for i in mol]}"
            )
        mol.set_root(root_atom)
        return mol

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

    @property
    def chain(self):
        return self._chain

    @property
    def residues(self):
        return self._chain.child_list

    @property
    def atoms(self):
        _atoms = [atom for residue in self.residues for atom in residue]
        return _atoms

    def get_root(self):
        return self.root_atom

    def set_root(self, atom):
        self.root_atom = atom

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

    def get_patch(self):
        return self._patch

    def set_patch(
        self, patch: Union[str, utils.abstract.AbstractPatch] = None, _topology=None
    ):
        """
        Set a patch to be used for attaching other molecules to this one

        Parameters
        ----------
        patch : str or utils.abstract.AbstractPatch
            The patch to be used. Can be either a string with the name of the patch or an instance of the patch
            If None is given, the current patch is removed.
        _topology
            The topology to use for referencing the patch.
        """
        if not _topology:
            _topology = utils.get_default_topology()
        if patch is None:
            self._patch = None
        elif isinstance(patch, str):
            self._patch = _topology.get_patch(patch)
        elif isinstance(patch, utils.abstract.AbstractPatch):
            self._patch = patch
        else:
            raise ValueError(f"Unknown patch type {type(patch)}")

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
        rdx = len(self.residues)
        adx = len(self.atoms)
        other.reindex(rdx + 1, adx + 1)

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

    def attach(
        self,
        other: "Molecule",
        patch: Union[utils.abstract.AbstractPatch, str] = None,
        at_residue: Union[int, bio.Residue.Residue] = None,
        other_residue: Union[int, bio.Residue.Residue] = None,
        _topology=None,
    ):
        """
        Attach another structure to this one.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        patch : str or AbstractResidue
            A patch to apply when attaching. If none is given, the default patch that was set earlier
            on the molecule is used. If no patch was set, an AttributeError is raised. If a string is
            given, it is interpreted as the name of a patch in the topology.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the last residue of the molecule.
        other_residue : int or Residue
            The residue of the other molecule to attach. If None, the first residue of the other molecule.
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        """
        if not _topology:
            _topology = utils.get_default_topology()

        if not patch:
            patch = self._patch
        if not patch:
            raise AttributeError(
                "No patch was set for this molecule. Either set a default patch or provide a patch when attaching."
            )
        if isinstance(patch, str):
            patch = _topology.get_patch(patch)
        if not patch:
            raise ValueError(
                "No patch was found with the given name. Either set a default patch or provide a patch when attaching."
            )

        p = structural.__default_keep_copy_patcher__
        p.apply(patch, self, other, at_residue, other_residue)
        p.merge()
        return self

    def repeat(self, n: int, patch=None, inplace: bool = True):
        """
        Repeat the molecule n times into a homo-polymer.

        Parameters
        ----------
        n : int
            The number or units of the final polymer.
        patch : str or AbstractPatch
            The patch to use when patching individual units together.
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.

        Returns
        -------
        molecule
            The modified molecule (either the original object or a copy)
        """
        if not isinstance(n, int):
            raise TypeError("Can only multiply a molecule by an integer")
        if n <= 0:
            raise ValueError("Can only multiply a molecule by a positive integer")
        if not patch and not self._patch:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        _other = deepcopy(self)
        if not inplace:
            obj = deepcopy(self)
        else:
            obj = self

        _patch = obj._patch
        if patch:
            obj % patch

        for i in range(n - 1):
            obj += _other

        if _patch:
            obj % _patch

        return obj

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

    def __add__(self, other) -> "Molecule":
        """
        Add two molecules together. This will return a new molecule.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only add two molecules together")

        patch = self._patch
        if not patch:
            patch = other._patch_to_use
        if not patch:
            raise RuntimeError(
                "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 + mol2)"
            )

        p = structural.__default_copy_copy_patcher__
        p.apply(patch, self, other)
        new = p.merge()
        return new

    def __iadd__(self, other) -> "Molecule":
        """
        Attach another molecule to this one
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only add two molecules together")

        patch = self._patch
        if not patch:
            patch = other._patch_to_use
        if not patch:
            raise RuntimeError(
                "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 += mol2)"
            )

        self.attach(other, patch)
        return self

    def __mul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the * operator (i.e. mol * 3)
        This requires that the molecule has a patch defined
        """
        if not isinstance(n, int):
            raise TypeError("Can only multiply a molecule by an integer")
        if n <= 0:
            raise ValueError("Can only multiply a molecule by a positive integer")
        if not self._patch:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        new = self + self
        for i in range(n - 2):
            new += self

        return new

    def __imul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the *= operator (i.e. mol *= 3)
        This requires that the molecule has a patch defined
        """
        if not self._patch:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        self.repeat(n)

        return self

    def __matmul__(self, residue) -> "Molecule":
        """
        Set the residue at which the molecule should be attached to another molecule
        using the @ operator (i.e. mol @ 1, for residue 1)
        """
        self.set_attach_residue(residue)
        return self

    def __repr__(self):
        return f"Molecule({self.id})"


if __name__ == "__main__":

    # from timeit import timeit

    # # man = Molecule.from_compound("MAN")
    # glc = Molecule.from_compound("GLC")

    # # man.adjust_indexing(glc)
    # # assert glc.atoms[0].serial_number == len(man.atoms) + 1

    # t1 = timeit()

    # glc.repeat(5, "14bb")
    # # glc % "14bb"
    # # glc *= 10

    # t2 = timeit()

    # print(t2 - t1)

    from glycosylator.utils import visual

    # v = visual.MoleculeViewer3D(glc)
    # v.show()

    man = Molecule.from_pdb("support/examples/membrane.pdb", model=4)
    man.infer_bonds()

    v = visual.MoleculeViewer3D(man.make_atom_graph())
    v.show()
