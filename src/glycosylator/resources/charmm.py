from biobuild.resources.charmm import *
import biobuild.utils.defaults as defaults

import os

DIR = os.path.dirname(os.path.abspath(__file__))
defaults.DEFAULT_CHARM_TOPOLOGY_FILE = os.path.join(DIR, "topology.json")

defaults.set_default_instance(
    "Topology", read_topology(defaults.DEFAULT_CHARM_TOPOLOGY_FILE)
)
