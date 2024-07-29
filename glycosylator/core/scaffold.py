"""
The `Scaffold` class is used to represent a scaffold structure such as a protein or membrane on to which a modification in form of one or more `Glycan`s can be added.
The Scaffold also inherets from BuildAMol and can be used in many ways similar to Molecules and Glycans. 
However, it lacks some features such as `repeat` method.


Attaching Glycans
-----------------

The primary purpose of the `Scaffold` is to have `Glycan`s attached to it.
This is done by using the `attach` method, which offers a versatile interface to attach one or more `Glycan`s to the `Scaffold` structure.
Alternatively, one may use the `glycosylate` method (which is just an alias for `attach`) or the stand-alone `glycosylate` function (the difference is
that the function by default makes a copy of the scaffold rather than working inplace).

We can attach Glycans to scaffolds via the `attach` method to:
- one residue (the attach_residue)
- a list of specified residues
- all residues that match a specific sequence pattern (sequon, only for protein scaffolds)


Attaching to a specific residue
-------------------------------
Below are two examples of attaching a glycan to a protein scaffold at a specific residue.
In the first we use default "pythonic" methods to properly attach the glycan to the protein:

.. code-block:: python

    import glycosylator as gl

    # load some Glycans from pickle files 
    # (if we have already some processed ones...)
    my_prot = gl.Scaffold.load("my_protein.pkl")
    my_glycan = gl.Glycan.load("my_glycan.pkl")

    # find some ASN residue the protein to set
    # as the attach-residue where the glycan should be attached
    # (just the first ASN we can find)
    asn.attach_residue = my_prot.get_residue("ASN")

    # get a linkage
    # (glycosylator already has a linkage for default N-glycosylation,
    # but we could define our own linkage if we have some special case)
    link = gl.get_linkage("ASN-glyco")

    # now attach the glycan to the protein at the ASN residue
    my_prot.attach(my_glycan, link)

In this second version we change we achieve the same but with an operator "short-hand" syntax because we feel "extra cocky":

.. code-block:: python

    # set the linkage as default modifier and specify the ASN as attach-residue
    my_prot % link @ my_prot.get_residue("ASN")

    # now attach the glycan to the ASN residue
    my_prot += my_glycan

    # or in one line (which is NOT inplace, however, but makes a copy)
    glycosylated = my_prot % link @ my_prot.get_residue("ASN") + my_glycan


Attaching to a list of residues
-------------------------------

Often we want to modify a number of residues in the same way, for example when we want to attach a glycan to all
residues in a protein that match a specific sequence pattern (e.g. all N-glycosylation sequons). We can do this by providing
a list of `residues` to the `attach` method. 

.. code-block:: python

    # find all ASN residues in the protein
    residues_to_modify = my_prot.get_residues("ASN", by="name")
    
    # now attach the glycan to all ASN residues
    my_prot.attach(my_glycan, link, residues=residues_to_modify)
  
This code will automatically add the glycan to all ASN residues in the list `residues_to_modify`.
"""

from collections import defaultdict
from typing import Union
import warnings
import numpy as np
from scipy.spatial.distance import cdist

from glycosylator.core.Glycan import Glycan
import glycosylator.resources as resources

from buildamol.core import Molecule, Linkage, molecule
import buildamol.core.entity as entity
import buildamol.utils.defaults as defaults
import buildamol.base_classes as base_classes
import buildamol.utils.auxiliary as aux
import buildamol.structural as structural


__all__ = ["Scaffold", "scaffold", "glycosylate"]

_acceptable_glycan_anchor_atoms = ("O", "N")
"""
The elements to restrict to when searching for anchor atoms for glycans that are directly attached to the scaffold.
"""

_acceptable_glycan_root_atoms = ("C",)
"""
The elements to restrict to when searching for the root atom of a glycan.
"""


def scaffold(scaf=None) -> "Scaffold":
    """
    Make a Scaffold object from a structure or a file.

    Parameters
    ----------
    scaf : str or object
        This can be a path to a file or a structure-like object like a BuildAMol or RDKit molecule.

    Returns
    -------
    Scaffold
        A Scaffold object.
    """
    try:
        mol = molecule(scaf)
        scaf = Scaffold(mol.structure)
        scaf.add_bonds(*mol.get_bonds())
        return scaf
    except Exception as e:
        raise ValueError("Could not make a Scaffold from the input.") from e


def glycosylate(
    scaffold: "Scaffold",
    glycan: "Glycan",
    link: "Linkage" = None,
    residues: list = None,
    sequon: str = None,
    inplace: bool = False,
) -> "Scaffold":
    """
    Glycosylate a scaffold with a glycan.

    Parameters
    ----------
    scaffold : Scaffold
        The Scaffold to glycosylate.
    glycan : Glycan
        The Glycan to attach.
    link : Linkage
        The Linkage to use. If a sequon is used, the sequon can be used to determine the Linkage automatically. If a default glycosylation linkage is available for the specified residues no linkage needs to be defined.
    residues : list, optional
        A list of Residue objects to glycosylate.
    sequon : str, optional
        A regex pattern to search the sequence for matching residues.
    inplace : bool, optional
        Whether to modify the Scaffold inplace.

    Returns
    -------
    Scaffold
        The glycosylated Scaffold.
    """
    if not sequon and (not residues and not link):
        raise ValueError(
            "Either a list of residues or a sequon pattern must be provided."
        )

    return scaffold.attach(
        glycan,
        link,
        residues=residues,
        sequon=sequon,
        chain="new",
        inplace=inplace,
    )


class Scaffold(entity.BaseEntity):
    """
    The :class:`Scaffold` class is used to represent a scaffold structure such as a protein or membrane onto which a
    modification in form of one or more :class:`glycosylator.core.Glycan.Glycan`s can be added.
    """

    def __init__(self, structure, model: int = 0) -> None:
        super().__init__(structure, model)

        self._excluded_chains = set()
        self._excluded_residues = set()

        self._attached_glycans = defaultdict(dict)
        self._glycan_chain_mapping = {}
        self._internal_residues = set()
        self._internal_bonds = set()
        self.type = "scaffold"

    def index(self, residue: Union[int, base_classes.Residue]) -> int:
        """
        Returns the index of an amino acid residue within the sequence
        of a protein scaffold.

        Parameters
        ----------
        residue : int or Residue
            The residue to get the index for

        Returns
        -------
        str
            The chain id
        int
            The index of the residue within the sequence of the chain
        """
        res = self.get_residue(residue)
        if res is None:
            raise ValueError(f"Residue '{residue}' not found")
        chain = res.parent
        index = chain.child_list.index(res)
        return chain.id, index

    @property
    def glycans(self) -> list:
        """
        Returns a list of all glycans that have been attached to the scaffold.
        """
        return list(self._attached_glycans.values())

    def get_glycans(self) -> dict:
        """
        Returns a dictionary of all glycan Glycans that have been attached to the scaffold
        The dictionary contains the atoms which are connected to the glycan Glycans as keys and the glycan Glycans as values.
        """
        return self._attached_glycans

    @property
    def glycan_residues(self) -> set:
        """
        Returns a set of all residues that belong to an attached glycan.
        """
        return set(
            [
                res
                for glycan in self._attached_glycans.values()
                for res in glycan.get_residues()
            ]
        )

    @property
    def scaffold_residues(self) -> set:
        """
        Returns a set of all residues that belong to the scaffold structure.
        """
        return set(self.get_residues()) - self.glycan_residues

    @property
    def internal_residues(self) -> set:
        """
        Returns a set of all residues that have been labelled as "internal" based on
        Solvent Accessible Surface Area (SASA) calculations. These residues are excluded from the main
        scaffold structure by default.
        """
        return self._internal_residues

    @classmethod
    def from_molecule(cls, mol) -> "Scaffold":
        """
        Create a Scaffold from a Molecule

        Parameters
        ----------
        mol : Molecule
            The molecule to create the scaffold from

        Returns
        -------
        Scaffold
            A scaffold created from the molecule
        """
        new = cls(mol._base_struct, model=mol._model.id)
        new._add_bonds(*mol.bonds)
        return new

    def copy(self):
        internal_residues = set(self._internal_residues)
        internal_bonds = set(self._internal_bonds)

        self.fill()
        new = entity.BaseEntity.copy(self)
        for bond in internal_bonds:
            bond = new.get_atom(bond[0].serial_number), new.get_atom(
                bond[1].serial_number
            )
            bond = new.get_bond(*bond)
            new._internal_bonds.add(bond)
            new._bonds.remove(bond)
            new._AtomGraph.remove_edge(*bond)

        for res in internal_residues:
            res = new.get_residue(res.serial_number)
            p = res.parent
            p.detach_child(res.get_id())
            res.set_parent(p)
            new._internal_residues.add(res)
            for atom in res.get_atoms():
                new._AtomGraph.remove_node(atom)

        for root, glycan in new._attached_glycans.items():
            _chain = new._glycan_chain_mapping[glycan]

            for chain in glycan.get_chains():
                glycan._model.child_dict.pop(chain.get_id(), None)
                chain._new_id()
                glycan._model.child_dict[chain.get_id()] = chain

            for res in glycan.get_residues():
                chain = res.parent
                chain.child_dict.pop(res.get_id(), None)
                _chain.child_dict.pop(res.get_id(), None)
                res._new_id()
                _chain.child_dict[res.get_id()] = res
                chain.child_dict[res.get_id()] = res

            for atom in glycan.get_atoms():
                res = atom.parent
                res.child_dict.pop(atom.get_id(), None)
                atom._new_id()
                res.child_dict[atom.get_id()] = atom

            glycan._AtomGraph.clear()
            bonds = glycan.get_bonds()
            for bond in bonds:
                atom1, atom2 = bond
                order = bond.order
                glycan._AtomGraph.add_edge(atom1, atom2, order=order, bond_obj=bond)
                new._AtomGraph.add_edge(atom1, atom2, order=order, bond_obj=bond)
            new._AtomGraph.add_edge(
                root,
                glycan.root_atom,
                order=1,
                bond_obj=new.get_bond(root, glycan.root_atom),
            )

        return new

    def reindex(
        self, start_chainid: int = 1, start_resid: int = 1, start_atomid: int = 1
    ):
        super().reindex(start_chainid, start_resid, start_atomid)
        count = self.count_atoms() + 1
        for glycan in self._attached_glycans.values():
            glycan.reindex(start_atomid=count)
            count += glycan.count_atoms()
        return self

    def exclude_chain(self, chain: Union[str, base_classes.Chain]):
        """
        Exclude a chain from the scaffold from letting any
        Glycans attach to any of its residues. This will also
        prevent any sequon matches to be found by the `find` method.

        Parameters
        ----------
        chain : str or Chain
            The chain to exclude
        """
        if isinstance(chain, str):
            _chain = next((i for i in self._model.get_chains() if i.id == chain), None)
        else:
            _chain = chain
        if _chain is not None:
            self._excluded_chains.add(_chain)
        else:
            raise ValueError(f"Chain '{chain}' not found")

    def exclude_chains(self, *chains: Union[str, base_classes.Chain]):
        """
        Exclude multiple chains from the scaffold from letting any
        Glycans attach to any of its residues. This will also
        prevent any sequon matches to be found by the `find` method.

        Parameters
        ----------
        chains : str or Chain
            The chains to exclude
        """
        for chain in chains:
            self.exclude_chain(chain)

    def include_chain(self, chain: Union[str, base_classes.Chain]):
        """
        Include a previously excluded chain in the scaffold again to let Glycans
        attach to its residues. This will also
        allow sequon matches to be found by the `find` method.

        Parameters
        ----------
        chain : str or Chain
            The chain to include
        """
        if isinstance(chain, str):
            _chain = next((i for i in self._model.get_chains() if i.id == chain), None)
        else:
            _chain = chain
        if _chain is None:
            raise ValueError(f"Chain '{chain}' not found")
        elif _chain in self._excluded_chains:
            self._excluded_chains.remove(_chain)

    def include_chains(self, *chains: Union[str, base_classes.Chain]):
        """
        Include multiple previously excluded chains in the scaffold again to let Glycans
        attach to their residues. This will also
        allow sequon matches to be found by the `find` method.

        Parameters
        ----------
        chains : str or Chain
            The chains to include
        """
        for chain in chains:
            self.include_chain(chain)

    def exclude_residue(self, residue: Union[int, base_classes.Residue]):
        """
        Exclude a residue of the scaffold from modification
        by a Glycan

        Parameters
        ----------
        residue : int or Residue
            The residue to exclude
        """
        res = self.get_residue(residue)
        if res is not None:
            self._excluded_residues.add(res)
        else:
            raise ValueError(f"Residue '{residue}' not found")

    def exclude_residues(self, *residues: Union[int, base_classes.Residue]):
        """
        Exclude multiple residues of the scaffold from modification
        by a Glycan

        Parameters
        ----------
        residues : int or Residue
            The residues to exclude
        """
        for residue in residues:
            self.exclude_residue(residue)

    def include_residue(self, residue: Union[int, base_classes.Residue]):
        """
        Include a previously excluded residue of the scaffold again for modification
        by a Glycan

        Parameters
        ----------
        residue : int or Residue
            The residue to include
        """
        res = self.get_residue(residue)
        if res is None:
            raise ValueError(f"Residue '{residue}' not found")
        elif res in self._excluded_residues:
            self._excluded_residues.remove(res)

    def include_residues(self, *residues: Union[int, base_classes.Residue]):
        """
        Include multiple previously excluded residues of the scaffold again for modification
        by a Glycan

        Parameters
        ----------
        residues : int or Residue
            The residues to include
        """
        for residue in residues:
            self.include_residue(residue)

    def is_excluded_chain(self, chain: Union[str, base_classes.Chain]) -> bool:
        """
        Check if a chain is excluded from the scaffold

        Parameters
        ----------
        chain : str or Chain
            The chain to check

        Returns
        -------
        bool
            True if the chain is excluded
        """
        if isinstance(chain, str):
            _chain = next((i for i in self._model.get_chains() if i.id == chain), None)
        else:
            _chain = chain
        if _chain is None:
            raise ValueError(f"Chain '{chain}' not found")
        return _chain in self._excluded_chains

    def is_excluded_residue(self, residue: Union[int, base_classes.Residue]) -> bool:
        """
        Check if a residue is excluded from the scaffold

        Parameters
        ----------
        residue : int or Residue
            The residue to check

        Returns
        -------
        bool
            True if the residue is excluded
        """
        res = self.get_residue(residue)
        if res is None:
            raise ValueError(f"Residue '{residue}' not found")
        return res in self._excluded_residues

    def get_glycan(self, anchor_or_index) -> Glycan:
        """
        Get a glycan that is attached to the scaffold

        Parameters
        ----------
        anchor_or_index : int or str or Atom or Residue
            This can be either the anchor atom to which the glycan attached
            or the parent residue of the anchor atom. Alternatively it can be an index
            from the glycan id (this requires that the glycan id string ends in this number; a string is also acceptible).

        Returns
        -------
        Glycan
            The glycan that is attached to the scaffold
        """
        if isinstance(anchor_or_index, (int, str)):
            for root, glycan in self._attached_glycans.items():
                if glycan.id.endswith(str(anchor_or_index)):
                    return glycan
            return None
        elif isinstance(anchor_or_index, base_classes.Atom):
            return self._attached_glycans.get(anchor_or_index, None)
        elif isinstance(anchor_or_index, base_classes.Residue):
            for root, glycan in self._attached_glycans.items():
                if root.parent == anchor_or_index:
                    return glycan
            return None

    def extend_glycan(self, glycan: Glycan, new_iupac_string: str):
        """
        Extend a glycan by adding new residues onto it based on an IUPAC string.
        This method will forward to `Glycan.extend_to`.

        Parameters
        ----------
        glycan : Glycan
            The glycan to extend
        new_iupac_string : str
            The IUPAC string of the new residue to add
        """
        # first identify the glycan
        for root, g in self._attached_glycans.items():
            if g == glycan:
                break
        else:
            raise ValueError(f"Glycan '{glycan}' not found")

        glycan = g
        glycan.extend_to(new_iupac_string)

        adx = self.count_atoms()
        for atom in glycan.get_atoms():
            atom.serial_number = adx
            adx += 1

        # now make sure all the residues and bonds are assosciated with the scaffold
        chain = self._glycan_chain_mapping[glycan]
        for res in glycan.get_residues():
            if chain.child_dict.get(res.get_id(), None):
                continue
            chain.child_list.append(res)
            chain.child_dict[res.get_id()] = res
        # self.add_residues(
        #     *(i for i in glycan.get_residues() if i not in chain.child_list),
        #     adjust_seqid=False,
        #     chain=chain,
        # )
        self._AtomGraph.migrate_bonds(glycan._AtomGraph)
        self._set_bonds(*glycan.get_bonds())

    def find_glycans(
        self, chain: str = "same", infer_bonds: bool = False, autolabel: bool = False
    ):
        """
        Find existing glycans in the scaffold structure.

        Parameters
        ----------
        chain : str
            The chain to which the glycans should be added. If 'same', the glycans are added to the
            same chain as the scaffold residue onto which they are attached. If 'new' all glycans are added
            to a single new chain. If 'each', each glycan is added to a new chain of it's own.
        infer_bonds : bool
            If True, all atomic bonds are inferred for the glycans. This is useful if you do not want to infer the bonds for the entire structure and only need the data for the glycans.
        autolabel : bool
            Whether to autolabel the atoms in the glycans according to the IUPAC naming scheme. This requires that all bonds are available for the glycans. If you are unsure, also set `infer_bonds` to True. Note however, that due to aberrant bonding patterns in larger glycans the automatic labelling might not always work as expected!
        """
        ref_residues = resources.reference_glycan_residue_ids()
        _glycan_residues = [r for r in self.get_residues() if r.resname in ref_residues]
        glycan_residues = np.empty(len(_glycan_residues), dtype=object)
        glycan_residues[:] = _glycan_residues
        glycan_residue_centers = np.array([r.coord for r in _glycan_residues])

        # due to the automatic unpacking of numpy this should be fine
        # (if not for some reason, change to the version below...)
        glycan_residue_atoms = np.array(_glycan_residues, dtype=object)
        # glycan_residue_atoms = np.array(
        #     [list(r.child_list) for r in glycan_residues], dtype=object
        # )

        glycan_residue_atom_coords = [
            np.array([a.coord for a in atoms]) for atoms in glycan_residue_atoms
        ]

        yet_to_add = set(glycan_residues)

        scaffold_residues = np.array(
            list((set(self.get_residues()).difference(yet_to_add))), dtype=object
        )
        scaffold_residue_centers = np.array([r.coord for r in scaffold_residues])

        glycans = {}
        glycan_chain_map = {}
        gdx = 1
        min_length, max_length = (
            defaults.DEFAULT_BOND_LENGTH / 2,
            defaults.DEFAULT_BOND_LENGTH * 1.3,
        )

        residue_mapping = {}

        # first parse through all glycan_residues to find the ones
        # that are attached directly to the scaffold
        for idx, glycan_residue in enumerate(glycan_residues):
            if not glycan_residue in yet_to_add:
                continue

            # get close-by residues
            mask = (
                np.linalg.norm(scaffold_residue_centers - glycan_residue.coord, axis=1)
                < 10
            )
            if not np.any(mask):
                continue

            close_residues = scaffold_residues[mask]
            close_by_atoms = list((i for r in close_residues for i in r.child_list))
            close_by_atom_coords = np.array([a.coord for a in close_by_atoms])

            residue_atom_coords = glycan_residue_atom_coords[idx]

            # find atoms that would be in the right distance range to be bonding
            dists = cdist(residue_atom_coords, close_by_atom_coords)
            _match_mask = (dists < max_length) & (dists > min_length)
            if not np.any(_match_mask):
                continue

            # now get the bonding atoms
            _matching_x_indices, _matching_y_indices = np.where(_match_mask)

            for i in _matching_x_indices:
                incoming_root = glycan_residue.child_list[i]
                if incoming_root.element in _acceptable_glycan_root_atoms:
                    break
            else:
                continue

            for j in _matching_y_indices:
                scaffold_root = close_by_atoms[j]
                if scaffold_root.element in _acceptable_glycan_anchor_atoms:
                    break
            else:
                continue

            # create the glycan
            glycan_residue.parent.detach_child(glycan_residue.get_id())
            glycan = Glycan.empty("glycan" + str(gdx))
            glycan.add_residues(glycan_residue, adjust_seqid=False)
            glycan.set_root(incoming_root)
            glycans[scaffold_root] = glycan
            residue_mapping[glycan_residue] = glycan
            yet_to_add.remove(glycan_residue)
            gdx += 1

        # then parse through all glycans to find the ones that are attached
        # to other glycan residues
        if len(yet_to_add) > 0:
            runs = 0
            ref = len(yet_to_add)
            preemptive_mapping = {}
            while len(yet_to_add) > 0:
                glycan_residue = yet_to_add.pop()

                # get close-by residues
                d = np.linalg.norm(
                    glycan_residue_centers - glycan_residue.coord, axis=1
                )
                mask = (0 < d) & (d < 10)
                if not np.any(mask):
                    yet_to_add.add(glycan_residue)
                    runs += 1
                    if runs > ref * 3:
                        raise RuntimeError(
                            f"Cannot infer glycan association for at least one of {yet_to_add}"
                        )
                    continue
                close_residues = glycan_residues[mask]
                close_by_atoms = list((i for r in close_residues for i in r.child_list))
                close_by_atom_coords = np.array([a.coord for a in close_by_atoms])

                residue_atoms = glycan_residue.child_list
                residue_atom_coords = np.array([a.coord for a in residue_atoms])

                # find atoms that would be in the right distance range to be bonding
                dists = cdist(residue_atom_coords, close_by_atom_coords)
                _match_mask = (dists < max_length) & (dists > min_length)

                if not np.any(_match_mask):
                    yet_to_add.add(glycan_residue)
                    runs += 1
                    if runs > ref * 3:
                        raise RuntimeError(
                            f"Cannot infer glycan association for at least one of {yet_to_add}"
                        )
                    continue

                # now get the bonding atoms
                _matching_x_indices, _matching_y_indices = np.where(_match_mask)

                for i in _matching_x_indices:
                    incoming_root = residue_atoms[i]
                    if incoming_root.element in _acceptable_glycan_root_atoms:
                        break
                else:
                    yet_to_add.add(glycan_residue)
                    runs += 1
                    if runs > ref * 3:
                        raise RuntimeError(
                            f"Cannot infer glycan association for at least one of {yet_to_add}"
                        )
                    continue

                for j in _matching_y_indices:
                    close_by_root = close_by_atoms[j]
                    if close_by_root.element in _acceptable_glycan_anchor_atoms:
                        break

                else:
                    yet_to_add.add(glycan_residue)
                    runs += 1
                    if runs > ref * 3:
                        raise RuntimeError(
                            f"Cannot infer glycan association for at least one of {yet_to_add}"
                        )
                    continue

                # get the glycan that the close_by_residue is part of
                close_by_residue = close_by_root.parent
                glycan = residue_mapping.get(close_by_residue, None)
                if glycan is None:
                    preemptive_mapping[glycan_residue] = (
                        close_by_residue,
                        incoming_root,
                        close_by_root,
                    )
                    continue

                glycan.add_residues(glycan_residue, adjust_seqid=False)
                glycan._set_bond(close_by_root, incoming_root)
                self._set_bond(close_by_root, incoming_root)
                residue_mapping[glycan_residue] = glycan
                yet_to_add.discard(glycan_residue)

            # now add any residues which could not be added directly
            # to the glycans
            for glycan_residue, (
                close_by_residue,
                glycan_root,
                close_by_root,
            ) in preemptive_mapping.items():
                glycan = residue_mapping.get(close_by_residue, None)
                if glycan is None:
                    warnings.warn(
                        RuntimeWarning(
                            "[weird residue] No close-by glycan found for: "
                            + str(close_by_residue)
                            + " @ "
                            + str(close_by_residue.full_id)
                        )
                    )
                    continue

                glycan.add_residues(glycan_residue, adjust_seqid=False)
                glycan._set_bond(close_by_root, glycan_root)
                self._set_bond(close_by_root, glycan_root)
                residue_mapping[glycan_residue] = glycan
                yet_to_add.discard(glycan_residue)

        # postprocessing
        if chain == "new":
            _chain = base_classes.Chain("")
            self.add_chains(_chain)

        for root, glycan in glycans.items():

            for res in glycan.get_residues():
                glycan._set_bonds(*self.get_bonds(res))

            if infer_bonds:
                glycan.infer_bonds(restrict_residues=True)
            if autolabel:
                glycan.autolabel()

            glycan.infer_glycan_tree()
            self._set_bond(root, glycan.root_atom)

            if chain == "each":
                _chain = base_classes.Chain("")
                self.add_chains(_chain)
            elif chain == "same":
                _chain = root.parent.parent

            glycan.rename_chain(glycan._working_chain, _chain.id)

            for res in glycan.get_residues():
                if _chain.child_dict.get(res.get_id(), None):
                    continue
                _chain.child_list.append(res)
                _chain.child_dict[res.get_id()] = res

            glycan_chain_map[glycan] = _chain
            self._AtomGraph.migrate_bonds(glycan._AtomGraph)
            self._set_bonds(*glycan.get_bonds())
            glycan._scaffold = self

        self._attached_glycans.update(glycans)
        self._glycan_chain_mapping.update(glycan_chain_map)
        return glycans

    def add_glycan(
        self, glycan: "Glycan", root: Union[int, base_classes.Atom], chain: str = "same"
    ):
        """
        Manually add a glycan to the scaffold structure. This is useful if you have a scaffold
        where you know which glycans are attached to which root atoms, so you do not have to
        search for them, which may take some time, depending on the size of the scaffold.

        Note
        ----
        This method simply adds an entry to the glycan registry of the scaffold, it does **not**
        attach any external glycans to the scaffold structure! Use the `attach` or `glycosylate` methods
        for doing so!

        Parameters
        ----------
        glycan : Glycan
            The glycan to add. All residues of the glycan should be part of the scaffold structure.

        root : int or Atom
            The root atom of the glycan. This is the atom to which the glycan is attached.

        chain : str
            The chain to which the glycan should be added. If 'same', the glycan is added to the
            same chain as the scaffold residue onto which it is attached. If 'new', the glycan is added
            to a new chain. Alternatively, a specific chain can be provided.
        """
        root = self.get_atom(root)
        if root is None:
            raise ValueError(f"Root Atom '{root}' not found")
        if not glycan.root_atom:
            raise ValueError("Glycan has no root atom defined")

        if chain == "same":
            chain = root.parent.parent
        elif chain == "new":
            chain = base_classes.Chain("")
            self.add_chains(chain)
        else:
            _chain = self.get_chain(chain)
            if _chain is None:
                raise ValueError(f"Chain '{chain}' not found")
            chain = _chain

        self._attached_glycans[root] = glycan
        self._glycan_chain_mapping[glycan] = chain
        glycan.rename_chain(glycan._working_chain, chain.id)
        glycan._scaffold = self
        return self

    def remove_glycan(self, glycan: "Glycan"):
        """
        Remove a glycan from the scaffold structure

        Parameters
        ----------
        glycan : Glycan
            The glycan to remove
        """
        if glycan not in self._attached_glycans.values():
            raise ValueError(f"Glycan '{glycan}' not found")
        chain = self._glycan_chain_mapping[glycan]
        for res in glycan.get_residues():
            chain.detach_child(res.get_id())
        root = next(k for k, v in self._attached_glycans.items() if v == glycan)
        self._attached_glycans.pop(root)
        self._glycan_chain_mapping.pop(glycan)
        self._remove_bond(root, glycan.root_atom)
        glycan._scaffold = None

        # re add a hydrogen to the root atom
        hydrogenator = structural.Hydrogenator()
        hydrogenator.add_hydrogens(root, self)
        return self

    def count_glycans(self) -> int:
        """
        Returns the number of glycans that are currently attached to the scaffold

        Returns
        -------
        int
            The number of glycans
        """
        return len(self.glycans)

    def get_glycan_root_connections(
        self, glycan: Union["Glycan", int] = None, include_scaffold_ancestors: int = 0
    ) -> list:
        """
        Get the bonds that connect the glycans to the scaffold.

        Parameters
        ----------
        glycan : Glycan or int
            A specific glycan to get the connections for. If None, all glycans are considered.
        include_scaffold_ancestors : int
            By default (0) only the direct bonds between scaffold-anchor atom and glycan-root atom are included.
            This may cause some glycans to be virtually un-optimizable due to bad orientation of the scaffold residue.
            To mitigate this, one can include a number of "upstream" bonds from the scaffold-anchor atoms up to n neighbors
            to ensure a better orientation of the glycan by also allowing adjustments in the scaffold structure.
            Set this value to the number of ancestor bonds to include.

        Returns
        -------
        list
            A list of tuples specifying the connections between the scaffold and the glycans. As well as any scaffold ancestors that are included.
        """
        _filter = lambda x: x.element != "H" and x not in glycan.get_atoms()
        connections = []
        if glycan is not None:
            if not isinstance(glycan, Glycan):
                g = self.get_glycan(glycan)
                if g is None:
                    raise ValueError(
                        f"Glycan '{glycan}' not found. Either provide the glycan directly or use the index of the glycan."
                    )
                glycan = g
            root = next(k for k, v in self._attached_glycans.items() if v == glycan)
            glycan_gen = ((root, glycan),)
        else:
            glycan_gen = self.get_glycans().items()

        for root, glycan in glycan_gen:
            connections.append((root, glycan.root_atom))
            if include_scaffold_ancestors > 0:
                neighs = self.get_neighbors(
                    root, include_scaffold_ancestors, filter=_filter
                )
                for n in neighs:
                    neighs2 = self.get_neighbors(n, 1, filter=_filter)
                    if len(neighs2) < 2:
                        continue
                    for n2 in neighs2:
                        if len(self.get_neighbors(n2, 1, filter=_filter)) < 2:
                            continue
                        if glycan.root_atom in self.get_descendants(n2, n):
                            connections.append((n2, n))
                            break
        return connections

    def add_residues(
        self,
        *residues: base_classes.Residue,
        chain: str = None,
        adjust_seqid: bool = True,
        _copy: bool = False,
    ):
        """
        Add residues to the structure

        Parameters
        ----------
        residues : Residue
            The residues to add
        chain : str
            The chain to add the residues to. If None, the residues are added to the last chain
        adjust_seqid : bool
            If True, the seqid of the residues is adjusted to
            match the current number of residues in the structure
            (i.e. a new residue can be given seqid 1, and it will be adjusted
            to the correct value of 3 if there are already two other residues in the Glycan).
        _copy : bool
            If True, the residues are copied before adding them to the structure.
            This is useful if you want to add the same residue to multiple Glycans, while leaving
            them and their original parent structures intact.
        """
        if chain is None:
            chain = self._chain
        elif chain not in self.chains:
            chain = self._model.child_dict.get(chain)
        if chain is None:
            raise ValueError(f"Chain '{chain}' not found")
        else:
            _current_chain = self._chain
            self._chain = chain
            super().add_residues(*residues, adjust_seqid=adjust_seqid, _copy=_copy)
            self._chain = _current_chain

    def hollow_out(self, cutoff: Union[int, float] = 0.75):
        """
        Remove all residues and atoms from the structure that are not accessible on the surface.
        This is useful for creating a simplified representation of the structure, which will ease
        the computational load for structure optimization of attached glycan molecules.

        Parameters
        ----------
        cutoff : int or float
            The cutoff SASA value for residues to be classfied as "surface residues".
            If an integer is provided, the value is interpreted as the absolute threshold
            for the SASA value. If a float is provided, the value is interpreted as a
            percentage of the maximal SASA value in the structure.
        """
        # get surface residues
        _attached_residues = self.glycan_residues
        structure = aux.DummyStructure(
            [res for res in self.residues if res not in _attached_residues]
        )
        residues = structural.infer_surface_residues(structure, cutoff=cutoff)

        # works fine with a list-comprehension...
        to_remove = (
            res
            for res in self.residues
            if res not in residues and res not in _attached_residues
        )

        for residue in to_remove:
            chain = residue.parent
            chain.detach_child(residue.get_id())
            residue.set_parent(chain)

            for atom in residue.get_atoms():
                atom_bonds = self._get_bonds(atom, None)
                self._AtomGraph.remove_node(atom)
                self._internal_bonds.update(atom_bonds)
                for b in atom_bonds:
                    self._bonds.remove(b)

            self._internal_residues.add(residue)

    def fill(self):
        """
        Fill the structure with the residues that were removed by the `hollow_out` method.
        """
        for residue in self._internal_residues:
            chain = residue.parent
            chain.add(residue)
            for atom in residue.get_atoms():
                self._AtomGraph.add_node(atom)

        for bond in self._internal_bonds:
            self._AtomGraph.add_edge(*bond)
            self._bonds.append(bond)

        self._internal_bonds.clear()
        self._internal_residues.clear()

    def to_pdb(self, filename: str, symmetric: bool = False):
        self.fill()  # fill the structure before writing it to a file
        super().to_pdb(filename, symmetric=symmetric)

    def to_molfile(self, filename: str):
        self.fill()
        super().to_molfile(filename)

    def to_json(self, filename: str):
        self.fill()
        super().to_json(filename)

    def to_biopython(self):
        self.fill()
        return super().to_biopython()

    def to_rdkit(self):
        self.fill()
        return super().to_rdkit()

    def to_pybel(self):
        self.fill()
        return super().to_pybel()

    def to_stk(self):
        self.fill()
        return super().to_stk()

    def to_smiles(self):
        self.fill()
        return super().to_smiles()

    def make_residue_graph(self, detailed: bool = False, locked: bool = True):
        graph = super().make_residue_graph(detailed, locked)
        # also add all remaining non connected (scaffold residues)
        for residue in self.residues:
            if residue not in self._internal_residues and residue not in graph:
                graph.add_node(residue)
        return graph

    get_residue_graph = make_residue_graph

    def attach(
        self,
        mol: "Glycan",
        link: "Linkage" = None,
        residues: list = None,
        chain="each",
        inplace: bool = True,
        _topology=None,
    ):
        """
        Attach a Glycan to the scaffold.

        Parameters
        ----------
        mol : Glycan
            The Glycan to attach. This glycan will be attached via the attach_residue so be sure to specify it if not the first residue should be used!
        link : Linkage
            The linkage to use when attaching the glycan. If None, the default linkage that was set earlier on the scaffold is used (if defined) or a linkage is identified based on the residue names to which the glycan should be attached.
        residues : list
            A list of residues at which to attach copies of the glycan. This is an alternative to using a sequon.
        chain : str
            The chain to which the glycan should be attached. If "same", the glycan is attached to the same chain as the respective scaffold residue.
            This can be set to "new" to create a new chain for all incoming glycans, or "each" to create a new chain for each glycan that is attached.
        inplace : bool
            Whether to modify the scaffold inplace or copy it before attaching the glycan(s) at the specified residue(s).
        _topology
            A particular topology to use for referning linkages from (if a sequon is provided).

        Returns
        -------
        scaffold
            The scaffold with the Glycan(s) attached, which is either a copy of the original scaffold or the original.
        """
        # ==========================================================
        #                    check input parameters
        # ==========================================================
        if not isinstance(mol, Glycan):
            raise ValueError(f"Expected a Glycan object, got {type(mol)} instead.")

        if mol.attach_residue is None:
            mol.attach_residue = 1

        # first check if all residues are present in the scaffold and all have a linkage defined to attach the glycans
        # either from provided input or default linkages
        if residues is not None:
            if isinstance(residues, (str, int, base_classes.Residue)):
                residues = [residues]
            must_have_glyco_linkage = self._linkage is None and link is None
            _residues = []
            for res in residues:
                res = self.get_residue(res)
                if not res:
                    raise ValueError(f"Residue '{res}' not found")
                _residues.append(res)
                if must_have_glyco_linkage:
                    _link1 = resources.get_linkage(res.resname + "-glyco")
                    if _link1 is None:
                        raise ValueError(
                            f"No glycosylation linkage found for residue '{res.resname}'. Either provide a linkage or specify a default using 'add_linkage(linkage(..., id='{res.resname}-glyco'))'."
                        )
            residues = _residues

        elif self.attach_residue is None:
            raise ValueError(
                "Either a list of residues or an attach_residue must be available to glycosylate anything."
            )

        elif self.attach_residue is not None:
            residues = [self.attach_residue]

        _topology = _topology or resources.get_default_topology()

        if not inplace:
            scaffold = self.copy()
        else:
            scaffold = self

        adx = self.count_atoms()

        # ==========================================================
        #          attach Glycan to multiple residues
        # ==========================================================
        if residues:

            one_new_chain = chain == "new"
            each_new_chain = chain == "each"

            if one_new_chain:
                _chain = base_classes.Chain(chr(65 + len(scaffold.chains)))
                scaffold._model.add(_chain)

            if link is None and self._linkage is not None:
                link = self._linkage

            s = structural.__default_keep_keep_stitcher__

            for residue in residues:

                # get the linkage (this should work since we already checked for linkages above)
                if link is None:
                    _link = resources.get_linkage(residue.resname + "-glyco")
                else:
                    _link = link

                if each_new_chain:
                    _chain = base_classes.Chain(chr(65 + len(scaffold.chains)))
                    scaffold._model.add(_chain)
                elif chain == "same" or chain is None:
                    _chain = residue.parent

                # copy the Glycan to avoid changing the original
                _mol = mol.copy()
                _mol.root_atom = _mol.get_atom(
                    _link.bond[1], residue=_mol.attach_residue
                )

                scaffold, _mol = s.apply(
                    scaffold,
                    _mol,
                    *_link.deletes,
                    *_link.bond,
                    residue,
                    _mol.attach_residue,
                )

                # adjust the glycan-id if it is not given
                if not _mol.id:
                    _mol.id = f"glycan-{len(scaffold.glycans) + 1}"

                # new policy for adding residues
                for res in _mol.get_residues():
                    _chain.link(res)
                    for atom in res.child_list:
                        atom.serial_number = adx
                        adx += 1

                # scaffold.add_residues(*_mol.get_residues(), chain=_chain)
                scaffold._bonds.extend(_mol.get_bonds())
                scaffold._AtomGraph.migrate_bonds(_mol._AtomGraph)
                scaffold._add_bond(s._anchors[0], s._anchors[1])

                root = scaffold.get_atom(s._anchors[0], residue=residue)
                scaffold._attached_glycans[root] = _mol
                scaffold._glycan_chain_mapping[_mol] = _chain
                if hasattr(_mol, "_scaffold"):
                    _mol._scaffold = self
            return scaffold

    # alias
    glycosylate = attach

    def __add__(self, mol):
        """
        Attach a Glycan to the scaffold
        """
        return self.attach(mol, inplace=False)

    def __iadd__(self, mol):
        """
        Attach a Glycan to the scaffold
        """
        return self.attach(mol, inplace=True)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.id})"


if __name__ == "__main__":
    import glycosylator as gl

    s = Scaffold.from_pdb(
        "/Users/noahhk/GIT/glycosylator/__projects__/solf2/solf2_man5_glycosylated_raw.pdb"
    )
    g = s.find_glycans(infer_bonds=True)
    assert len(g) > 0
    o = s.copy()
    o.reindex()
    o.to_pdb("o.pdb")
    # g = next(iter(g.values()))
    # s.extend_glycan(g, "ManNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-")
    # s.snfg().show()

    # s = Scaffold.from_pdb(
    #     "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    # )
    # s.reindex()
    # s.infer_bonds()
    # g = gl.glycan(
    #     "GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)[GalNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
    # )

    # s.attach(g, residues=[*s.get_residues("ASN", by="name")[:3]])
    # out = s.copy()
    # s.show3d(residue_graph=True)

    # original_glycosylated = "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.pdb"
    # s = Scaffold.from_pdb(original_glycosylated)
    # s.reindex()
    # glycans = s.find_glycans(True, False)
    # glycan1 = list(glycans.values())[0]
    # glycan1.autolabel()
    # pass

    # DO NOT DELETE ANYTHING BELOW THIS LINE
    # f1 = "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    # s = Scaffold.from_pdb(f1)
    # s.reindex()
    # s.infer_bonds()
    # s1 = s.copy()
    # assert s1 is not s
    # assert s1.count_bonds() == s.count_bonds()

    # exit()
    # # s.hollow_out()
    # # sequon = "(N)(?=[A-OQ-Z][ST])"
    # # residues = s.find(sequon)

    # # f = "/Users/noahhk/GIT/glycosylator/scaffold.pkl"
    # # s.save(f)

    # # s = Scaffold.load("scaffold.pkl")
    # # _s = Scaffold(s.structure)
    # # _s.bonds = s.bonds
    # # s = _s

    # mol = gl.Glycan.from_json(
    #     "/Users/noahhk/GIT/glycosylator/support/examples/man8.json"
    # )
    # mol.root_atom = 1

    # # link = gl.linkage("ND2", "C1", ["HD22"], ["O1", "HO1"])

    # # out = glycosylate(s, mol, link, sequon="N-linked")
    # # out.make_residue_graph().show()
    # # exit()

    # # mol = gl.Glycan.from_pdb("/Users/noahhk/GIT/glycosylator/support/examples/man9.pdb")
    # # mol.infer_bonds(restrict_residues=False)
    # # mol.reindex()
    # # mol = gl.Glycan.load("glycan.pkl")

    # # mol.root_atom = 1

    # import time

    # print("Start attaching...")
    # t1 = time.time()
    # residues = s.find("N-linked")
    # residues = list(residues.values())[0]
    # link = gl.linkage("ND2", "C1", ["HD22"], ["O1", "HO1"])
    # s.attach(
    #     mol,
    #     link,
    #     # remove_atoms=("HD22",),
    #     # mol_remove_atoms=("HO1", "O1"),
    #     # at_atom="ND2",
    #     # sequon="N-linked",
    #     residues=residues,
    # )
    # t2 = time.time()
    # print("Oh yeah - it works!")
    # from datetime import datetime

    # s.to_pdb(f"final_scaffold_superduper.{datetime.now()}.pdb")
    # print(f"Time: {t2 - t1}")

    # g = s.find_glycans()
    # _g = list(g.values())[0]
    # pass
