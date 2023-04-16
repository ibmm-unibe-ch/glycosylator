"""
Functions to stitch molecules together even in the absence of a patch
"""

from typing import Union

import numpy as np
from scipy.spatial.transform import Rotation
from copy import deepcopy

import Bio.PDB as bio

import glycosylator.core.molecule as molecule
import glycosylator.structural as structural
import glycosylator.structural.connector as base
import glycosylator.optimizers as optimizers


class Stitcher(base.Connector):
    """
    This class is responsible for stitching molecules together
    when there is not a patch available which specifies the immediate
    geometry of their connection. The only required input is a pair of atoms
    that should form a new bond, and two tuples of atoms (at least one from each molecule)
    that are being removed in the process.

    Parameters
    ----------
    copy_target : bool
        Whether to copy the target molecule before stitching
    copy_source : bool
        Whether to copy the source molecule before stitching
    """

    def __init__(self, copy_target: bool = False, copy_source: bool = False):
        super().__init__(copy_target, copy_source)
        self._removals = (None, None)
        self._policy = None, None, None

    def stitch(
        self,
        target: "Molecule",
        source: "Molecule",
        target_removals: tuple,
        source_removals: tuple,
        target_atom: Union[int, bio.Atom.Atom] = None,
        source_atom: Union[int, bio.Atom.Atom] = None,
    ) -> "Molecule":
        """
        Stitch the source and target molecules together

        Parameters
        ----------
        target : Molecule
            The target molecule
        source : Molecule
            The source molecule. This will be attached to the target molecule.
        target_removals : tuple
            A tuple of atoms to be removed from the target molecule. These may either be the atom objects directly or their serial numbers.
        source_removals : tuple
            A tuple of atoms to be removed from the source molecule. These may either be the atom objects directly or their serial numbers.
        target_atom : int or bio.Atom.Atom
            The atom on the target molecule to which the source molecule will be attached. This may either be the atom object directly or its serial number.
            If none is provided, the molecule's "root atom" is used (if defined).
        source_atom : int or bio.Atom.Atom
            The atom on the source molecule to which the target molecule will be attached. This may either be the atom object directly or its serial number.


        Returns
        -------
        Molecule
            The stitched molecule
        """

        if self.copy_target:
            target = deepcopy(target)
        if self.copy_source:
            source = deepcopy(source)

        self.target = target
        self.source = source

        self._anchors = self._find_anchors(target_atom, source_atom)
        self._removals = self._find_removals(target_removals, source_removals)

        self._align_anchors()

        self._remove_atoms()
        self._optimize()

        return self.target

    def _align_anchors(self):
        """
        Align the anchor atoms such that they are in the right positions to form a bond
        This will move the source molecule such that its anchor atom is in the position
        of one of the target molecule's removed atoms
        """

        ref_target_atom = self._anchors[0]
        neighs = self.target.get_neighbors(ref_target_atom)
        ref_target_atom = next(i for i in neighs if i in self._removals[0])

        ref_source_atom = self._anchors[1]
        neighs = self.source.get_neighbors(ref_source_atom)
        ref_source_atom = next(i for i in neighs if i in self._removals[1])

        # self._v.draw_point("anchor target", self._anchors[0].coord, color="blue")
        # self._v.draw_point("anchor source", self._anchors[1].coord, color="blue")
        # self._v.draw_point("ref_target_atom", ref_target_atom.coord, color="red")
        # self._v.draw_point("ref_source_atom", ref_source_atom.coord, color="red")

        _old_coords = np.stack(
            [
                self._anchors[1].coord,
                ref_source_atom.coord,
            ]
        )
        _new_coords = np.array(
            [
                ref_target_atom.coord,
                self._anchors[0].coord,
            ]
        )

        # compute translation vector
        old_centroid = _old_coords.mean(axis=0)
        new_centroid = _new_coords.mean(axis=0)
        translation_vector = new_centroid - old_centroid

        _relative_old_coords = _old_coords - old_centroid
        _relative_new_coords = _new_coords - new_centroid

        H = (_relative_old_coords).T.dot(_relative_new_coords)
        U, S, VT = np.linalg.svd(H)
        R = VT.T @ U.T

        # self._v.draw_edges(self.source.bonds, color="black", opacity=0.5)

        for atom in self.source.get_atoms():
            vec = atom.coord - old_centroid
            new_coord = (R @ vec) + old_centroid + translation_vector
            atom.set_coord(new_coord)

        # self._v.draw_edges(self.source.bonds, color="orange", opacity=0.5)
        # self._v.draw_point("anchor source (new)", self._anchors[1].coord, color="teal")

    def _find_anchors(
        self,
        target_atom=None,
        source_atom=None,
        target_residue=None,
        source_residue=None,
    ):
        """
        Find the anchor atoms for the target and source molecules

        Parameters
        ----------
        target_atom : int or bio.Atom.Atom
            The atom on the target molecule to which the source molecule will be attached. This may either be the atom object directly or its serial number.
            If none is provided, the molecule's "root atom" is used (if defined).
        source_atom : int or bio.Atom.Atom
            The atom on the source molecule to which the target molecule will be attached. This may either be the atom object directly or its serial number.
        target_residue : int or bio.Residue.Residue
            The residue on the target molecule to which the source molecule will be attached. This may either be the residue object directly or its serial number.
            If none is provided, the molecule's "attach residue" is used (if defined).
        source_residue : int or bio.Residue.Residue
            The residue on the source molecule to which the target molecule will be attached. This may either be the residue object directly or its serial number.
            If none is provided, the molecule's "attach residue" is used (if defined).
        """
        if target_atom is None:
            target_atom = self.target.root_atom
        if source_atom is None:
            source_atom = self.source.root_atom
        if target_residue is None:
            target_residue = self.target.attach_residue
        if source_residue is None:
            source_residue = self.source.attach_residue

        return self._get_anchors(
            (target_atom, source_atom), target_residue, source_residue
        )

    def _find_removals(self, target_removals, source_removals):
        """
        Find the atoms to remove from the target and source molecules
        while stitching them together
        """
        target_removals = [self.target.get_atom(atom) for atom in target_removals]
        source_removals = [self.source.get_atom(atom) for atom in source_removals]

        return (set(target_removals), set(source_removals))

    def _remove_atoms(self):
        """
        Remove the atoms specified in the removals list
        """
        self.target.remove_atoms(*self._removals[0])
        self.source.remove_atoms(*self._removals[1])

    def _optimize(self):
        """
        Optimize the geometry of the source molecule
        """

        tmp = molecule.Molecule.empty(self.target.id)
        self.target.adjust_indexing(self.source)

        anchor_1_serial = self._anchors[0].serial_number
        anchor_2_serial = self._anchors[1].serial_number

        _index_mappings = {
            i.full_id[3:]: i.serial_number for i in self.target.get_atoms()
        }
        _index_mappings.update(
            {i.full_id[3:]: i.serial_number for i in self.source.get_atoms()}
        )

        tmp.add_residues(
            self._target_residue, self._source_residue, adjust_seqid=False, _copy=True
        )

        bonds = list(
            i
            for i in self.target.bonds
            if i[0] in self._target_residue.child_list
            and i[1] in self._target_residue.child_list
        )
        for b in bonds:
            b = (b[0].full_id, b[1].full_id)
            tmp.add_bond(*b)
        bonds = list(
            i
            for i in self.source.bonds
            if i[0] in self._source_residue.child_list
            and i[1] in self._source_residue.child_list
        )
        for b in bonds:
            b = (b[0].full_id, b[1].full_id)
            tmp.add_bond(*b)

        tmp.add_bond(anchor_1_serial, anchor_2_serial)
        tmp.lock_all()

        graph = tmp.make_residue_graph()
        graph.make_detailed(include_heteroatoms=True)

        edges = list(sorted(tmp.get_residue_connections()))
        env = gl.optimizers.MultiBondRotatron(graph, edges)

        agent = optimizers.EvoRotator(env, verbose=True, mutation_stdev=0.3, popsize=50)
        best = agent(100)
        self._policy = edges, *best

    def merge(self):
        """
        Merge the source molecule into the target molecule
        """
        self.target.add_residues(*self.source.residues)

        for bond in self.source.bonds:
            b = (bond[0].full_id, bond[1].full_id)
            self.target.add_bond(*b)
            if self.source.is_locked(*bond):
                self.target.lock_bond(*b)

        self.target.add_bond(self._anchors[0], self._anchors[1])
        self.target.reindex()

        # v = gl.utils.visual.MoleculeViewer3D(self.target)

        if self._policy[0] is not None:

            edges, angles, _ = self._policy

            for bond, angle in zip(edges, angles):

                bond = (bond[0].full_id, bond[1].full_id)
                if len(self.target.get_bonds(*bond, either_way=False)) == 0:
                    self.target.add_bond(*bond)

                bond = self.target.get_bonds(*bond, either_way=False)[0]

                # v.draw_vector(
                #     "rotation",
                #     bond[0].coord,
                #     bond[1].coord - bond[0].coord,
                #     color="orange",
                # )

                self.target.rotate_around_bond(
                    *bond, np.degrees(angle), descendants_only=True
                )

            self._policy = None, None, None

        # v.draw_edges(self.target.bonds, color="red")
        # v.show()
        return self.target


__default_keep_copy_stitcher__ = Stitcher(False, True)
__default_copy_copy_stitcher__ = Stitcher(True, True)


if __name__ == "__main__":

    import glycosylator as gl

    glc = gl.Molecule.from_compound("GLC")
    glc2 = gl.Molecule.from_compound("GLC")

    glc2.rotate_around_bond(6, 5, 68)
    glc2.rotate_around_bond(3, 4, 41)
    s = Stitcher(True, True)

    v = None
    for i in range(20):

        s.stitch(glc, glc2, ("O1", "HO1"), ("HO4",), "C1", "O4")
        new_glc = s.merge()

        if v is None:
            v = gl.utils.visual.MoleculeViewer3D(new_glc)
        else:
            v.draw_edges(new_glc.bonds, color="red", opacity=0.5)
    v.show()
    pass
