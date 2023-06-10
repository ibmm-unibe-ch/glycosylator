"""
The `structural` module contains classes and functions to work with molecular structures.
It is at the heart of glycosylator functionality and provides most of the useful features.

Almost all functions and classes in this module are integrated into the main glycosylator API
through methods of the `Molecule` and `Scaffold` classes - so it is usually not necessary to 
use this module directly. However there are some useful features that are not directly integrated
into the API, and in some cases users may want to use this module directly to access them.
"""
from glycosylator.structural.infer import *
from glycosylator.structural.smiles import *
from glycosylator.structural.base import *
from glycosylator.structural.patch import (
    Patcher,
    patch,
    __default_copy_copy_patcher__,
    __default_keep_copy_patcher__,
)

from glycosylator.structural.stitch import (
    Stitcher,
    stitch,
    __default_copy_copy_stitcher__,
    __default_keep_copy_stitcher__,
)

from glycosylator.structural.neighbors import (
    AtomNeighborhood,
    ResidueNeighborhood,
    compute_quartets,
    compute_triplets,
    generate_triplets,
)

from glycosylator.structural.iupac import (
    IUPACParser,
)
