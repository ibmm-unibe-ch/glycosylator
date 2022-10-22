import prody

from .molecule_builder import MoleculeBuilder


class Glycosylator:
    def __init__(self, protein_atom_group, charmm_topology_path, charm_parameter_path):
        self.charmm_topology_path = charmm_topology_path
        self.charmm_parameter_path = charm_parameter_path
        self.builder = MoleculeBuilder(charmm_topology_path, charm_parameter_path)

    @classmethod
    def fromPDB(cls, protein_file_path, charmm_topology_path, charmm_parameter_path):
        protein_atom_group = prody.parsePDB(protein_file_path)
        return cls(protein_atom_group, charmm_topology_path, charmm_parameter_path)
