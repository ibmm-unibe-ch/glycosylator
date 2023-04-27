"""
This module contains classes for default linkeage instructions to attach molecules together and to scaffolds.
These classes are essentially forwarded versions of the `abstract` classes from the utilities and are offered for convenience.
"""

import glycosylator.utils.abstract as abstract


class Recipe(abstract.AbstractRecipe):
    """
    Using the `Recipe` class, a template reaction instruction is stored for stitching molecules together.
    """

    pass


class Patch(abstract.AbstractPatch):
    """
    Using the `Patch` class, a template reaction instruction is stored for attaching molecules to one another.
    """

    pass
