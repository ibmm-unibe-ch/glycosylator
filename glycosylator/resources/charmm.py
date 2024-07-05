from buildamol.resources.charmm import *
import buildamol.utils.defaults

import os

DIR = os.path.dirname(os.path.abspath(__file__))
buildamol.utils.defaults.DEFAULT_CHARM_TOPOLOGY_FILE = os.path.join(DIR, "topology.xml")

buildamol.utils.defaults.set_default_instance(
    "Topology", read_topology(buildamol.utils.defaults.DEFAULT_CHARM_TOPOLOGY_FILE)
)

glyco_link_patch_name_mapping = {
    "ASN-glyco": "NGLB",
    "SER-glyco": "OGLB",
    "THR-glyco": "SGPB",
    "CER-glyco": "CERB",
}
"""
Mapping between the standard (IC-free) patches and the CHARMM patches.
"""

available_charmm_patches = {}
"""
Mapping of all patches that are available from CHARMM
(only those that have ICs are included). 
"""
for patch in buildamol.get_default_topology().linkages:
    prefix = patch.id[:2] if patch.id[:1].isdigit() else patch.id[:3]
    available_charmm_patches.setdefault(prefix, [])
    if patch.has_IC:
        available_charmm_patches[prefix].append(patch.id)
