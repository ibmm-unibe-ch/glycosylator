"""
The `Molecule` class is a wrapper around a biopython structure and a core part 
of glycosylator functionality. It provides a convenient interface to molecular structures
and their properties, such as atoms, bonds, residues, chains, etc. The purpose of a `Molecule` is to
be attached onto a `Scaffold`. 


Making Molecules
================

To make a new molecule, the easiest way is to use the `Molecule.from_compound` classmethod.
This method is linked to the PDBECompounds database, which contains a large number of reference
compounds (currently restricted to sugars and amino acids due to the nature of "glycosylator"). The database
is queried using the `by` parameter, which can be one of the following:
- "id" for the PDB id (default)
- "name" for the name of the compound (must match any known synonym of the iupac name)
- "formula" for the chemical formula (usually ambiguous and will therefore often raise an error)
- "smiles" for the SMILES string (also accepts InChI)

.. code-block:: python

    from glycosylator import Molecule

    # create a new galactose molecule
    gal = Molecule.from_compound("GAL") # use the PDB id

    # create a new glucose molecule
    glc = Molecule.from_compound("alpha-d-glucose", by="name") # use the name

Alternatively, an existing structure can be loaded to a molecule from a PDB file using the `Molecule.from_pdb` classmethod. 
In this case the molecule may be as large as the PDB format supports (currently 99999 atoms).

.. code-block:: python

    my_glycan = Molecule.from_pdb("my_glycan.pdb")

Unless loaded from a PDB file, each `Molecule` starts with a single chain and residue,
specifying a single compound, such as Galactose or Glucose. Starting from there, multiple molecules
can be assembled together to form complex structures. 

Of course, an existing `biopython.Structure` object can also be used to create a `Molecule`:

.. code-block:: python

    from Bio.PDB import PDBParser

    parser = PDBParser()
    structure = parser.get_structure("my_structure", "my_structure.pdb")

    # do any modifications to the structure here
    # ...
    
    my_molecule = Molecule(structure)


Connecting Molecules
====================

Since most modifications are not simply single residues but rather complex structures, the second main purpose of a `Molecule` is to be easily 
connected to other Molecules to form a larger structure. To this end, the `Molecule` class provides a number of methods to easily assemble complex structures from small single residue molecules.


Forming Polymers
----------------

The simplest way to generate a large structure is probably the `repeat` method, which will repeat the given molecule `n` times
to form a homo-polymer.

.. code-block:: python

    # create a glucose molecule
    glc = Molecule.from_compound("GLC")

    # create cellulose from glucose
    # using a 1-4 beta-beta glycosidic linkeage
    glc.repeat(10, "14bb")

    # Now we have a cellulose of 10 glucoses

In the above example we used the `repeat` method explicitly, but we could also achieve the same with the short-hand `*=`. For this to work, we need to specify 
the linkeage type beforehand.

.. code-block:: python

    # specify the "default" linkeage type for connecting 
    # other molecules to this glucose
    glc.patch = "14bb"

    # now make a cellulose by multiplying glucoses
    glc *= 20

    # Now we have a cellulose of 20 glucoses

    
If we wish to keep `glc` as a single residue Glucose and still get our desired cellulose, 
we can set `inplace=False` when calling `repeat` or simply use the `*` operator, both of which 
will have the same effect of creating a new copy that houses the appropriate residues.

.. code-block:: python

    cel = glc.repeat(10, "14bb", inplace=False)

    # or (equivalently)
    glc.patch = "14bb"
    cel = glc * 10


    
Connecting different Molecules
-------------------------------

What if we want to connect different molecules? For example, we may want to connect a Galactose to a Glucose to form Lactose.
This can be achieved using the `attach` method, which will attach a given molecule to to another molecule.

.. code-block:: python

    glc = Molecule.from_compound("GLC")
    gal = Molecule.from_compound("GAL")

    # attach the galactose to the glucose
    # (we want a copy, so we set inplace=False just like with 'repeat')
    lac = glc.attach(gal, "14bb", inplace=False)

    # Now we have a lactose molecule


    
In the above example, the `attach` method is used to attach the galactose molecule to the glucose, but for those
among us who prefer a more shorty syntax, the `+` operator will do the same thing. In this case we need
to specify the linkage type (known as `patch`) as an attribute to the receiving molecule using `set_patch_or_recipe` first.

.. code-block:: python

    # specify that incoming molecules shall be
    # attached using a 1-4 beta linkage
    glc.set_patch_or_recipe("14bb")

    # now attach the galactose to the glucose
    glc += gal

For those who still want a shorter syntax, `set_patc_or_recipeh` is represeted by the `%` opreator, so the above example
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
    A molecule consists of a single chain.

    Parameters
    ----------
    structure : bio.PDB.Structure
        A biopython structure object
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
        super().__init__(structure, model)

        if not chain or len(self._model.child_list) == 1:
            self._working_chain = self._model.child_list[0]
        else:
            self._working_chain = self._model.child_dict.get(chain)
            if not self._working_chain:
                raise ValueError("The chain {} is not in the structure".format(chain))

        if root_atom:
            self.set_root(root_atom)
            if not root_atom in self._working_chain.get_atoms():
                raise ValueError("The root atom is not in the structure")
        self._root_atom = root_atom

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

    def get_residue_connections(self, triplet: bool = True, direct_by: str = "resid"):
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
        direct_by : str
            The attribute to sort by. Can be either "serial", "resid" or "root".
            In the case of "serial", the bonds are sorted by the serial number of the first atom.
            In the case of "resid", the bonds are sorted by the residue id of the first atom.
            In the case of "root", the bonds are sorted by the root atom of the first atom.
            Set to None to not sort the bonds.

        Returns
        -------
        set
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = super().get_residue_connections(triplet)
        if direct_by is not None:
            direct_connections = None
            if direct_by == "resid":
                direct_connections = self.get_residue_connections(
                    triplet=False, direct_by=None
                )
                direct_connections1, direct_connections2 = set(
                    a for a, b in direct_connections
                ), set(b for a, b in direct_connections)
                direct_connections = direct_connections1.union(direct_connections2)
            bonds = self._direct_bonds(bonds, direct_by, direct_connections)
        return set(bonds)

    def repeat(self, n: int, patch_or_recipe=None, inplace: bool = True):
        """
        Repeat the molecule n times into a homo-polymer.

        Parameters
        ----------
        n : int
            The number or units of the final polymer.
        patch_or_recipe : str or Patch or Recipe
            The patch or recipe to use when patching individual units together.
            If noe is given, the default patch or recipe is used (if defined).
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
        if not patch_or_recipe and not self._patch:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        _other = deepcopy(self)
        if not inplace:
            obj = deepcopy(self)
        else:
            obj = self

        _patch = obj._patch
        if patch_or_recipe:
            obj % patch_or_recipe

        for i in range(n - 1):
            obj += _other

        if _patch:
            obj % _patch

        return obj

    def attach(
        self,
        other: "Molecule",
        patch_or_recipe: Union[
            utils.abstract.AbstractPatch, str, utils.abstract.AbstractRecipe
        ] = None,
        at_residue: Union[int, bio.Residue.Residue] = None,
        other_residue: Union[int, bio.Residue.Residue] = None,
        inplace: bool = True,
        _topology=None,
    ):
        """
        Attach another structure to this one using a Patch or a Recipe.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        patch_or_recipe : str or Patch or Recipe
            Either a Patch to apply when attaching or a Recipe to use when stitching.
            If None is defined, the default patch or recipe that was set earlier on the molecule is used.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the defined `attach_residue` is used.
        other_residue : int or Residue
            The residue in the other molecule to attach this molecule to. If None, the defined `attach_residue` of the other molecule is used.
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.
        _topology : Topology
            The topology to use when attaching. If None, the topology of the molecule is used. Only used if the patch is a string.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only attach a Molecule to another Molecule")

        if not inplace:
            obj = deepcopy(self)
        else:
            obj = self

        if not patch_or_recipe:
            patch_or_recipe = obj._patch
            if not patch_or_recipe:
                raise ValueError("Cannot attach a molecule without a patch defined")

        if isinstance(patch_or_recipe, str):
            if not _topology:
                _topology = utils.get_default_topology()
            patch_or_recipe = _topology.get_patch(patch_or_recipe)

        if isinstance(patch_or_recipe, utils.abstract.AbstractPatch):
            obj.patch(
                other,
                patch_or_recipe,
                at_residue=at_residue,
                other_residue=other_residue,
                _topology=_topology,
            )
        elif isinstance(patch_or_recipe, utils.abstract.AbstractRecipe):
            obj.stitch(
                other,
                patch_or_recipe,
                at_residue=at_residue,
                other_residue=other_residue,
            )
        return obj

    def patch(
        self,
        other: "Molecule",
        patch: Union[utils.abstract.AbstractPatch, str] = None,
        at_residue: Union[int, bio.Residue.Residue] = None,
        other_residue: Union[int, bio.Residue.Residue] = None,
        _topology=None,
    ):
        """
        Attach another structure to this one using a Patch.

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

        if not patch:
            patch = self._patch
        if not patch:
            raise AttributeError(
                "No patch was set for this molecule. Either set a default patch or provide a patch when attaching."
            )
        if isinstance(patch, str):
            if not _topology:
                _topology = utils.get_default_topology()
            patch = _topology.get_patch(patch)
        if not patch:
            raise ValueError(
                "No patch was found with the given name. Either set a default patch or provide a patch when attaching."
            )

        p = structural.__default_keep_copy_patcher__
        p.apply(patch, self, other, at_residue, other_residue)
        p.merge()
        return self

    def stitch(
        self,
        other: "Molecule",
        recipe: utils.abstract.AbstractRecipe = None,
        remove_atoms=None,
        other_remove_atoms=None,
        at_atom=None,
        other_at_atom=None,
        at_residue=None,
        other_residue=None,
    ):
        """
        Stitch two molecules together by removing atoms and connecting them with a bond. This works without a pre-defined patch.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        recipe : Recipe
            The recipe to use when stitching. If None, the default recipe that was set earlier on the molecule is used (if defined).
        remove_atoms : list of int
            The atoms to remove from this molecule. Only used if no recipe is provided.
        other_remove_atoms : list of int
            The atoms to remove from the other molecule. Only used if no recipe is provided.
        at_atom : int or str or Bio.PDB.Atom
            The atom forming the bond in this molecule. If a string is provided, an `at_residue` needs to be defined from which to get the atom. If None is provided, the root atom is used (if defined). Only used if no recipe is provided.
        other_at_atom : int or str or Bio.PDB.Atom
            The atom to attach to in the other molecule. If a string is provided, an `other_residue` needs to be defined from which to get the atom. If None is provided, the root atom is used (if defined). Only used if no recipe is provided.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the `attach_residue` is used. Only used if a recipe is provided and the atoms
        """
        if not recipe and not remove_atoms:
            if not self._patch:
                raise AttributeError(
                    "No recipe was set for this molecule and no manual instructions were found. Either set a default recipe, provide a recipe when stitching, or provide the information about removed and bonded atoms directly."
                )
            recipe = self._patch

        if recipe:
            return self.stitch(
                other,
                remove_atoms=recipe.deletes[0],
                other_remove_atoms=recipe.deletes[1],
                at_atom=recipe.bonds[0][0],
                other_at_atom=recipe.bonds[0][1],
                at_residue=at_residue,
                other_residue=other_residue,
            )

        if not isinstance(other, Molecule):
            raise TypeError("Can only stitch two molecules together")

        if not at_atom:
            at_atom = self._root_atom
            if not at_atom:
                raise ValueError(
                    "No atom to attach to was provided and no root atom was defined in this molecule"
                )

        if not other_at_atom:
            other_at_atom = other._root_atom
            if not other_at_atom:
                raise ValueError(
                    "No atom to attach to was provided and no root atom was defined in the other molecule"
                )

        # at_atom = self.get_atom(at_atom, residue=at_residue)
        # other_at_atom = other.get_atom(other_at_atom, residue=other_residue)

        p = structural.__default_keep_copy_stitcher__
        p.apply(self, other, remove_atoms, other_remove_atoms, at_atom, other_at_atom)
        self = p.merge()
        return self

    def __add__(self, other) -> "Molecule":
        """
        Add two molecules together. This will return a new molecule.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only add two molecules together")

        patch = self._patch
        if not patch:
            patch = other._patch
        if not patch:
            raise RuntimeError(
                "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 + mol2)"
            )

        if isinstance(patch, utils.abstract.AbstractPatch):
            p = structural.__default_copy_copy_patcher__
            p.apply(patch, self, other)
        elif isinstance(patch, utils.abstract.AbstractRecipe):
            p = structural.__default_copy_copy_stitcher__
            p.apply(
                self,
                other,
                patch.deletes[0],
                patch.deletes[1],
                patch.bonds[0][0],
                patch.bonds[0][1],
            )
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
            patch = other._patch
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

    def __repr__(self):
        return f"Molecule({self.id})"


if __name__ == "__main__":
    from glycosylator import Recipe

    recipe = Recipe()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    glc = Molecule.from_compound("GLC")
    glc.set_patch_or_recipe(recipe)

    glc2 = deepcopy(glc)

    _current_residues = len(glc.residues)
    glc.attach(glc2)
    pass
    # import pickle

    # tmp = pickle.load(open("/Users/noahhk/GIT/glycosylator/tmp.pickle", "rb"))
    # tmp.get_residue_connections()
    # exit()
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

    # man = Molecule.from_pdb("support/examples/membrane.pdb", model=4)
    # man.infer_bonds()
    # man.

    man = Molecule.from_pdb("/Users/noahhk/GIT/glycosylator/support/examples/MAN9.pdb")
    man.infer_bonds(restrict_residues=False)

    g = man.make_residue_graph(True)
    v = visual.MoleculeViewer3D(g)

    for c in man.get_residue_connections():
        v.draw_vector(
            None, c[0].coord, 1.3 * (c[1].coord - c[0].coord), color="magenta"
        )
    v.show()
