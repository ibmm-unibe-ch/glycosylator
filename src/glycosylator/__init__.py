from glycosylator.core import *
import glycosylator.utils.visual as visual

# from biobuild import load_sugars
from biobuild.core import *

# this module overwrites some of the imports from the biobuild resources module
from glycosylator.resources import *

# from the resources module
load_glycosylator_compounds()

__version__ = "4.1.15"
__author__ = "Noah Kleinschmidt"
