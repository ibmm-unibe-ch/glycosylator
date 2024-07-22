import glycosylator.utils as __utils
import glycosylator.utils.visual as visual
from glycosylator.core import *

# from buildamol import load_sugars
from buildamol.core import *
import buildamol.graphs as graphs
import buildamol.structural as structural

# this module overwrites some of the imports from the buildamol resources module
from glycosylator.resources import *

import glycosylator.resources as resources
import glycosylator.optimizers as optimizers
from glycosylator.optimizers import *

# little hack to make sure the global utils refers to the glycosylator.utils module
# and not the glycosylator.optimizers.utils module
utils = __utils

# from the resources module
load_glycosylator_compounds()

__version__ = "5.6.3"
__author__ = "Noah Kleinschmidt"
