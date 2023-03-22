"""
This module defines classes to represent PDB structures as graph objects.
At the base is the `AtomGraph` object which directly translates into a PDB structure and can be saved as one. 
It stores all atoms and their connectivity. As an abstraction, the `AtomGraph` can produce a `ResidueGraph` object which represents each residue as a sphere of a certain volume
instead of all individually contained atoms.
"""

from glycosylator.graphs.AtomGraph import AtomGraph
from glycosylator.graphs.ResidueGraph import ResidueGraph

__all__ = ["AtomGraph", "ResidueGraph"]
