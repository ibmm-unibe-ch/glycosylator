"""
External data resource integrations to _glycosylator_. 

PDBe Compounds
==============

The Protein Data Bank in Europe (PDBe) provides a database of small molecules that are found in PDB structures.
They are available via the `PDBe Component Library <https://www.ebi.ac.uk/pdbe/pdb-component-library/>`. _Glycosylator_
integrates a subset of the database directly for Molecule creation. Since this database is integrated directly in glycosylator it can be used offline.

PubChem
=======

PubChem maintains an extensive database of small molecules. _Glycosylator_ integrates the PubChem API for Molecule creation.
Since this database is queried online and not integrated directly in glycosylator it requires an internet connection.
PubChem follows a different atom labelling convention than the CHARMM force field. Hence, many compounds may not be compatible
with the pre-set `Patches` in _glycosylator_ and may thus not work at all or produce erroneous results. To assist in importing 
PubChem compounds anyway, _glycosylator_ provides a function `pubchem_to_cif` that can be used to convert PubChem compound
data into a CIF file which can be more easily edited and subsequently loaded into _glycosylator_ using the `PDBECompounds` class
or using `Molecule.from_cif()`.

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
import glycosylator.resources.pubchem as pubchem
import glycosylator.resources.sequons as sequons
import glycosylator.resources.charmm as charmm

PDBECompounds = pdbe_compounds.PDBECompounds
CHARMMTopology = charmm.CHARMMTopology
CHARMMParameters = charmm.CHARMMParameters
SEQUONS = sequons.SEQUONS

pubchem_query = pubchem.query

read_topology = charmm.read_topology
read_parameters = charmm.read_parameters
load_topology = charmm.load_topology
load_parameters = charmm.load_parameters
save_topology = charmm.save_topology
save_parameters = charmm.save_parameters

read_compounds = pdbe_compounds.read_compounds
load_compounds = pdbe_compounds.load_compounds
save_compounds = pdbe_compounds.save_compounds
has_compound = pdbe_compounds.has_compound
get_compound = pdbe_compounds.get_compound


has_patch = charmm.has_patch
get_patch = charmm.get_patch
add_patch = charmm.add_patch
