"""
The base class for the Patcher and Stitcher classes.
"""


class Connector:
    def __init__(self, copy_target: bool = False, copy_source: bool = False):
        self.copy_target = copy_target
        self.copy_source = copy_source

        self.target = None
        self.source = None

        self._anchors = (None, None)
        self._target_residue = None
        self._source_residue = None

    def set_target(self, target: "Molecule"):
        """
        Set the target molecule

        Parameters
        ----------
        target : Molecule
            The target molecule
        """
        self.target = target

    def set_source(self, source: "Molecule"):
        """
        Set the source molecule

        Parameters
        ----------
        source : Molecule
            The source molecule
        """
        self.source = source

    def get_anchors(self, _ref_atoms, target_residue=None, source_residue=None):
        """
        Find appropriate anchor atoms in the source and target molecules

        Parameters
        ----------
        _ref_atoms : tuple
            The two atoms that form the bond between the two molecules.
            These may be specified by name or by index.
        target_residue : int or str or Residue
            The residue in the target molecule to which the source molecule will be patched.
            By default, the last residue in the target molecule will be used if it contains
            appropriate anchor atoms.

        source_residue : int or str or Residue
            The residue in the source molecule that will be patched into the target molecule.
            By default, the first residue in the source molecule will be used if it contains
            appropriate anchor atoms.

        Returns
        -------
        tuple
            The two atoms that form the bond between the two molecules
        """
        if not self.target:
            raise AttributeError("No target set")
        if not self.source:
            raise AttributeError("No source set")

        ref_atom_1 = self.target.get_atoms(_ref_atoms[0])
        ref_atom_2 = self.source.get_atoms(_ref_atoms[1])

        if target_residue:
            if isinstance(target_residue, int):
                target_residue = self.target.get_residue(id=target_residue)
            elif isinstance(target_residue, str):
                target_residue = self.target.get_residue(name=target_residue)
            ref_atom_1 = [i for i in ref_atom_1 if i.get_parent() == target_residue]

        if source_residue:
            if isinstance(source_residue, int):
                source_residue = self.source.get_residue(id=source_residue)
            elif isinstance(source_residue, str):
                source_residue = self.source.get_residue(name=source_residue)
            ref_atom_2 = [i for i in ref_atom_2 if i.get_parent() == source_residue]

        if len(ref_atom_1) == 0:
            raise ValueError("No anchor atom found in target molecule")
        if len(ref_atom_2) == 0:
            raise ValueError("No anchor atom found in source molecule")

        ref_atom_1 = ref_atom_1[-1]
        ref_atom_2 = ref_atom_2[0]

        self._anchors = (ref_atom_1, ref_atom_2)
        self._target_residue = ref_atom_1.get_parent()
        self._source_residue = ref_atom_2.get_parent()
        return ref_atom_1, ref_atom_2