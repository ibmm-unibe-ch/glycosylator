"""
Functions to patch molecules together
"""

import Bio.PDB as bio
import numpy as np
from scipy.spatial.transform import Rotation
from copy import deepcopy


import glycosylator.utils.abstract as abstract
import glycosylator.utils.structural.base as base
import glycosylator.utils.structural.infer as infer


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

    def patch_molecules(
        self,
        patch: "AbstractPatch" = None,
        target: "Molecule" = None,
        source: "Molecule" = None,
    ):
        """
        Patch two molecules together.
        This will merge the source molecule's residues into the target molecule.
        If no patch, target, and source have been set already, they can be provided
        as arguments to this method.

        Parameters
        ----------
        patch : AbstractPatch
            The patch to apply
        target : Molecule
            The target molecule
        source : Molecule
            The source molecule

        Returns
        -------
        Molecule
            The patched target molecule
        """
        if patch:
            self.patch = patch
        if target:
            self.target = target
        if source:
            self.source = source

        if self.copy_target:
            self.target = deepcopy(self.target)
        if self.copy_source:
            self.source = deepcopy(self.source)

        # a dictionary to store the anchor atoms in the source
        # molecule and their original coordinates
        self._source_computed_anchors = {}

        self.get_anchor_atoms()

        # compute internal coordinates to rebuild the source molecule
        _ics = infer.compute_internal_coordinates(self.source.bonds)
        self._source_ICs = abstract.AbstractEntity()
        self._source_ICs.internal_coordinates = _ics

        # self._v.draw_edges(self.source.bonds, color="red", opacity=0.5)

        # align the source molecule to the target molecule
        self._align_anchors()

        # now start imputing the source molecule into the target molecule
        # starting with the atoms that are described in the patch's internal coordinates
        self._infer_patch_IC_neighbors()

        # and now trans-rotate the source molecule such that
        # the inferred neighboring atoms overlay their imputed counterparts
        self._transpose_source()

        self._delete_atoms()
        self._merge()

        return self.target

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

        _ref_atoms = {int(i[0]) - 1: i[1:] for i in self.patch.bonds[0].atoms}

        ref_atom_1 = self.target.get_atom(_ref_atoms[0])
        ref_atom_2 = self.source.get_atom(_ref_atoms[1])

        self._anchors = (ref_atom_1, ref_atom_2)
        return ref_atom_1, ref_atom_2

    def _align_anchors(self):
        """
        Align the source to the target molecule such that the anchors are in the correct distance
        """
        ic_34 = self._match_IC_34()
        if len(ic_34) == 0:
            raise PatchError("No matching IC found")
        ic = ic_34[0]

        atom1 = self.target.get_atom(ic.atom1[1:])
        atom2 = self.target.get_atom(ic.atom2[1:])
        new_anchor_source = infer.compute_atom4_from_others(
            atom1.coord,
            atom2.coord,
            self._anchors[0].coord,
            ic,
        )

        # store the old coordinates and set the new coordinates
        self._source_computed_anchors[self._anchors[1]] = self._anchors[1].coord
        self._anchors[1].set_coord(new_anchor_source)

        # self._v.draw_point("anchor_source (new)", self._anchors[1].coord, color="green")

    def _infer_patch_IC_neighbors(self):
        """
        Infer the neighbors of the internal coordinates in the patch.
        This is required to rebuild the source molecule.
        """
        _obj = {
            "1": self.target,
            "2": self.source,
        }
        for ic in self.patch.internal_coordinates:
            if ic.atom4[0] == "1":
                continue

            atom4 = _obj[ic.atom4[0]].get_atom(ic.atom4[1:], by="id")
            if atom4 in self._source_computed_anchors.keys():
                continue

            atom1 = _obj[ic.atom1[0]].get_atom(ic.atom1[1:], by="id")
            atom2 = _obj[ic.atom2[0]].get_atom(ic.atom2[1:], by="id")
            atom3 = _obj[ic.atom3[0]].get_atom(ic.atom3[1:], by="id")

            _new_coords = infer.compute_atom4_from_others(
                atom1.coord, atom2.coord, atom3.coord, ic
            )

            # self._v.draw_point(
            #     atom4.id + " (old)",
            #     atom4.coord,
            #     color="brown",
            #     opacity=0.6,
            # )

            self._source_computed_anchors[atom4] = atom4.coord
            atom4.set_coord(_new_coords)

            # self._v.draw_point(
            #     atom4.id
            #     + " (new should-be)"
            #     + f"tb-deleted={atom4.id in self.patch.deletes[1]}",
            #     atom4.coord,
            #     color="blue",
            #     opacity=0.6,
            # )

    def _transpose_source(self):
        """
        Transpose the source molecule to overlay the
        """

        if len(self._source_computed_anchors) < 3:
            raise PatchError("Not enough anchors to transpose")

        _old_coords = np.stack(list(self._source_computed_anchors.values()))
        _new_coords = np.array([i.coord for i in self._source_computed_anchors.keys()])

        # for n in _old_coords:
        #     self._v.draw_point("old", n, color="red")

        # for n in _new_coords:
        #     self._v.draw_point("new", n, color="limegreen")

        # compute translation vector
        old_centroid = _old_coords.mean(axis=0)
        new_centroid = _new_coords.mean(axis=0)
        translation_vector = new_centroid - old_centroid

        # self._v.draw_vector(
        #     "translation vector",
        #     old_centroid,
        #     translation_vector,
        #     color="teal",
        # )

        _relative_old_coords = _old_coords - old_centroid
        _relative_new_coords = _new_coords - new_centroid

        H = (_relative_old_coords).T.dot(_relative_new_coords)
        U, S, VT = np.linalg.svd(H)
        R = VT.T @ U.T

        for atom in self.source.get_atoms():

            # self._v.draw_point(
            #     atom.id + " (old)",
            #     atom.coord,
            #     color="brown",
            #     opacity=0.3,
            #     showlegend=False,
            # )

            if atom in self._source_computed_anchors.keys():
                # vec = self._source_computed_anchors[atom] - old_centroid
                # new_coord = (R @ vec) + old_centroid + translation_vector

                # self._v.draw_point(
                #     atom.id + " (new computed from old)",
                #     new_coord,
                #     color="purple",
                #     opacity=0.6,
                # )
                continue

            vec = atom.coord - old_centroid
            new_coord = (R @ vec) + old_centroid + translation_vector
            atom.set_coord(new_coord)

        #     self._v.draw_point(
        #         atom.id + " (new)",
        #         atom.coord,
        #         color="lightblue",
        #         opacity=0.6,
        #     )

        # self._v.draw_edges(self.source.bonds, color="blue", opacity=0.6)

    def _match_IC(self, n_target: int, n_source: int):
        """
        Get appropriate internal coordinates from the patch
        that include the anchor atoms of the target and source molecule
        at specific positions (starting at 1).

        Parameters
        ----------
        n_target : int
            The position of the target anchor atom in the internal coordinate
        n_source : int
            The position of the source anchor atom in the internal coordinate

        Returns
        -------
        list
            All available internal coordinates that include the anchor atoms at the specified positions
        """

        ids = [None, None, None, None]
        ids[n_target - 1] = "1" + self._anchors[0].id
        ids[n_source - 1] = "2" + self._anchors[1].id

        ics = self.patch.get_internal_coordinates(
            *ids,
            mode="partial",
        )

        if len(ics) == 0:
            raise PatchError(
                "No internal coordinates found for anchor atoms at positions {} and {}".format(
                    n_target, n_source
                )
            )

        return ics

    def _match_IC_34(self, ic=None):
        """
        Get appropriate internal coordinates from the patch
        that include the anchor atoms of the target and source molecule
        at positions 3 and 4 and have three atoms from the target molecule.
        Or check if a given internal coordinate matches the criteria.

        Returns
        -------
        list or bool
            All available internal coordinates that include the anchor atoms at the specified positions
            or True if the given internal coordinate matches the criteria.
        """
        if ic:
            return ic.atom1[0] == "1" and ic.atom2[0] == "1" and ic.atom3[0] == "1"

        ids = ["1" + i.id for i in self.target.get_atoms()]
        ics = self._match_IC(3, 4)
        ics = [i for i in ics if i.atom1 in ids and i.atom2 in ids and i.atom3 in ids]
        if len(ics) == 0:
            raise PatchError(
                "No internal coordinates found for anchor atoms at positions 3 and 4"
            )
        return ics

    def _delete_atoms(self):
        """
        Delete all atoms that need to be deleted according to the patch
        """
        delete_from_target, delete_from_source = self.patch.deletes
        self.target.remove_atoms(*delete_from_target)
        self.source.remove_atoms(*delete_from_source)

    def _merge(self):
        """
        Merge the source molecule into the target molecule
        """
        self.target.adjust_indexing(self.source)
        self.target.add_residues(*self.source.residues)
        self.target._bonds.extend(self.source._bonds)
        self.target._bonds.append(self._anchors)
        self.target._locked_bonds.update(self.source._locked_bonds)


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

    man1.lock_all()
    man2.lock_all()

    top = gl.get_default_topology()

    patcher = Patcher(True, True)
    for i in ("12ab", "14bb", "12ab"):
        patch = top.get_patch(i)
        new = patcher.patch_molecules(patch, man1, man2)

        seen_atoms = set()
        for atom in new.atoms:
            assert atom.serial_number not in seen_atoms
            seen_atoms.add(atom.serial_number)

        res_con = gl.utils.structural.infer_residue_connections(new.chain, triplet=True)

        v2 = gl.utils.visual.MoleculeViewer3D(new)
        v2.draw_edges(edges=list(new.bonds), color="blue", opacity=1)
        v2.draw_edges(edges=list(new._locked_bonds), color="red", linewidth=3)
        v2.draw_edges(edges=res_con, color="limegreen", linewidth=4)
        v2.show()
