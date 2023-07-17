# from biobuild import load_sugars
from biobuild.core import *
from biobuild.resources import *


from glycosylator.core import *

# this module overwrites some of the imports from the biobuild resources module
from glycosylator.resources import *

# from the resources module
load_glycosylator_compounds()

__version__ = "4.0.12"
__author__ = "Noah Kleinschmidt"
