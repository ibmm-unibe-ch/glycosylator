import os
import pickle
import Bio.PDB as bio

import glycosylator.utils.constants as constants


# =================================================================
# Default settings of values
# =================================================================

DEFAULT_BOND_LENGTH = 1.6
"""
The default length of a bond in Angstrom
"""

DEFAULT_CHARMM_TOPOLOGY_FILE = os.path.join(constants.RESOURCES, "CHARMM.top.pkl")
"""
The path to the default CHARMM topology file
"""

DEFAULT_CHARMM_PARAMETERS_FILE = os.path.join(constants.RESOURCES, "CHARMM.prm.pkl")
"""
The path to the default CHARMM parameters file
"""


# =================================================================
# Default instances of auxiliary classes
# =================================================================

__bioPDBParser__ = bio.PDBParser()
"""
Default PDBParser using BioPython
"""
