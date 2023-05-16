"""
This module contains classes for default linkage instructions to attach molecules together and to scaffolds.


Types of Linkages 
==================

There are two kinds of linkages: `Patch` and `Recipe`. They specify in essence the same information, namely:
- which atoms are to be connected (one bond may be created)
- which atoms are to be deleted in the process from the two connecting structures

The only difference between the two is that the `Patch` additionally also specifies the geometry of the resulting structure after 
the new bond has been applied. This is to say, it also specifies the internal coordinates of the atoms in the immediate vicinity of
the newly created bond. This is especially useful, because it allows for molecules to be connected purely based on geometry and
is thus faster than the `Recipe` approach. The latter requires a step of optimization where the connected residues are rotated around the
newly formed bond in order to find a geometry that is energetically favorable.

Therefore, the `Patch` approach is preferred, but it is not always possible to use it - especially if it is not possible to establish 
the internal coordinates of bond atom neighborhood _a priori_ with confidence. 

Creating Linkages
==================

Linkages are created by instantiating the `Patch` or `Recipe` classes and then adding the desired instructions to them.
They share most of their methods for adding instructions, namely:

- `add_bond`: adds a bond between two atoms
- `add_delete`: deletes an atom from the source or target molecule
- `add_internal_coordinates`: (Patch only) adds internal coordinates for the atoms in the immediate vicinity of the newly formed bond


Example
-------

For instance, if we wished to make a new `Recipe` for attaching a glycan to a protein at an Asparagine residue,
we could do the following:

.. code-block:: python

    import glycosylator as gl

    # make a recipe for attaching to ASN
    my_recipe = gl.Recipe(id="glycan_to_asn")
    
    # add the bond to be created
    my_recipe.add_bond(("ND2", "C1"))

    # now specify which atoms are to be deleted from the source (glycan) and target (protein)
    my_recipe.add_delete("O1", "source")
    my_recipe.add_delete("HO1", "source")
    my_recipe.add_delete("HD22", "target")


To now use the recipe in order to attach a glycan to a protein, we could do the following:

.. code-block:: python

    # load some molecules from pickle files    
    my_prot = gl.Scaffold.load("my_protein.pkl")
    my_glycan = gl.Molecule.load("my_glycan.pkl")


    # find some ASN residue the protein (just the first ASN we can find)
    asn = my_prot.get_residue("ASN", by="name", chain="A")

    # now attach the glycan to the ASN residue
    my_prot.attach(my_glycan, my_recipe, residues=[asn])


.. note::

    On a technical note, the `Recipe` and `Patch` classes are just for convenience forwarded version of the `abstract.AbstractPatch` and `abstract.AbstractRecipe` classes.
    Hence, if you find functions or methods within _glycosylator_  that specifically ask for `AbstractPatch` or `AbstractRecipe` instances, you can use the `Patch` and `Recipe` classes.
    The other way around, you can also use the `AbstractPatch` and `AbstractRecipe` classes where `Patch` and `Recipe` are expected.

"""

import glycosylator.utils.abstract as abstract


class Patch(abstract.AbstractPatch):
    """
    Using the `Patch` class, a template reaction instruction is stored for attaching molecules to one another.
    """

    pass


class Recipe(abstract.AbstractRecipe):
    """
    Using the `Recipe` class, a template reaction instruction is stored for stitching molecules together.
    """

    pass


def write_json(patch_or_recipe, filename: str):
    """
    Write a `Patch` or `Recipe` to a JSON file.

    Parameters
    ----------
    patch_or_recipe : Patch or Recipe
        The patch or recipe to write to file
    filename : str
        The filename to write to
    """
    raise NotImplementedError
