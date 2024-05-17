from buildamol.resources.charmm import *
import buildamol.utils.defaults

import os

DIR = os.path.dirname(os.path.abspath(__file__))
buildamol.utils.defaults.DEFAULT_CHARM_TOPOLOGY_FILE = os.path.join(
    DIR, "topology.json"
)

buildamol.utils.defaults.set_default_instance(
    "Topology", read_topology(buildamol.utils.defaults.DEFAULT_CHARM_TOPOLOGY_FILE)
)
