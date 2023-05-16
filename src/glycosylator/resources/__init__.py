"""
External data resource integrations to _glycosylator_. 

PDBe Compounds
==============

The Protein Data Bank in Europe (PDBe) provides a database of small molecules that are found in PDB structures.
They are available via the `PDBe Component Library <https://www.ebi.ac.uk/pdbe/pdb-component-library/>`. _Glycosylator_
integrates a subset of the database directly for Molecule creation.

CHARMM
======

CHARMM is a molecular simulation software that is widely used in the field of computational chemistry.
It is developed and maintained by the `CHARMM Development Project <http://www.charmm.org/>`_.
_Glycosylator_ integrates the CHARMM force field for use in conformational optimization.

Sequons
=======

Sequons are short amino acid sequences that are recognized by the cell's glycosylation machinery.
They are used to identify potential glycosylation sites in proteins. _Glycosylator_ integrates a small database of sequons
that can be used to identify potential glycosylation sites in proteins.
"""
import glycosylator.resources.pdbe_compounds as pdbe_compounds
import glycosylator.resources.sequons as sequons
import glycosylator.resources.charmm as charmm


PDBECompounds = pdbe_compounds.PDBECompounds
CHARMMTopology = charmm.CHARMMTopology
CHARMMParameters = charmm.CHARMMParameters
SEQUONS = sequons.SEQUONS
