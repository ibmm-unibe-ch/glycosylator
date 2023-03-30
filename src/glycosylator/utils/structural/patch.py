"""
Functions to patch molecules together
"""

import Bio.PDB as bio
import numpy as np
from copy import deepcopy


class PatchError(Exception):
    pass


class Patcher:
    """
    This class is responsible for patching molecules together
    based on an `AbstractPatch` object and two `Molecule` objects,
    of which one is the target and one is the source. Residues from the source
    are integrated into the target molecule.

    Parameters
    ----------
    copy_target : bool
        Whether to copy the target molecule before patching
    copy_source : bool
        Whether to copy the source molecule before patching
    """
    
    def __init__(self, copy_target: bool = False, copy_source: bool = False):
        self.patch = None
        self.target = None
        self.source = None
        self.copy_target = copy_target
        self.copy_source = copy_source

        self._anchors = (None, None)
    
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
    
    def set_patch(self, patch: "AbstractPatch"):
        """
        Set the patch

        Parameters
        ----------
        patch : AbstractPatch
            The patch
        """
        self.patch = patch

    def patch_molecules(self, patch: "AbstractPatch", target: "Molecule", source: "Molecule"):
        """
        Parameters
        ----------
        patch : AbstractPatch
            The patch to apply
        target : Molecule
            The target molecule
        source : Molecule
            The source molecule
        """
        self.patch = patch
        self.target = target
        self.source = source

        if self.copy_target:
            self.target = deepcopy(self.target)
        if self.copy_source:
            self.source = deepcopy(self.source)
    
    def get_anchor_atoms(self):
        """
        Returns the two atoms that will be used for anchoring structure alignment.
        These atoms form a bond between the two molecules.

        Returns
        -------
        tuple
            The two atoms that form the bond between the two molecules
            the first atom is from molecule1 and the second is from molecule2.
        """
        if not self.patch:
            raise PatchError("No patch set")
        if not self.target:
            raise PatchError("No target set")
        if not self.source:
            raise PatchError("No source set")
        
        _ref_atoms = {int(i[0])-1: i[1:] for i in self.patch.bonds[0].atoms}

        ref_atom_1 = self.target.get_atom(_ref_atoms[0])
        ref_atom_2 = self.source.get_atom(_ref_atoms[1])

        self._anchors = (ref_atom_1, ref_atom_2)
        return ref_atom_1, ref_atom_2

    def _align_anchors(self):
        """
        Align the source to the target molecule such that the anchors are in the correct distance 
        """

if __name__ == "__main__":
    
    import glycosylator as gl

    man = "support/examples/MAN.pdb"
    man1 = gl.Molecule.from_pdb(man)
    man1.infer_bonds()

    man2 = deepcopy(man1)
    # now make sure that man2 has some different coordinates
    man2.rotate_around_bond(1, 2, 35)
    r = np.random.rand(3) * 0.1
    for i in man2.atoms:
        i.coord += r    

    top = gl.get_default_topology()
    abstract = top.get_residue(man1.id)
    
    patch = top.get_patch("12aa")

    pather = Patcher()
    pather.patch_molecules(patch, man1, man2)

