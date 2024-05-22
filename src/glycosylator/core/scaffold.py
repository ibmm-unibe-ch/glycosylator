"""
The `Scaffold` class is used to represent a scaffold structure such as a protein or membrane on to which a modification in form of one or more `Glycan`s can be added.
The Scaffold also inherets from Biobuild and can be used in many ways similar to Molecules and Glycans. 
However, it lacks some features such as `repeat` method.


Attaching Glycans
-----------------

The primary purpose of the `Scaffold` is to have `Glycan`s attached to it.
This is done by using the `attach` method, which offers a versatile interface to attach one or more `Glycan`s to the `Scaffold` structure.

.. note:: 

    Connecting a `Glycan` and `Scaffold` is based on `Recipe`s (not on `Patch`es!).
    Furthermore, "attaching" a `Glycan` to a `Scaffold` is defined as "connecting the Glycan's and Scaffold's root atoms".
    In practice, (re-)setting the objects' root atoms is handled behind the scenes, but is is important to keep this in
    mind when using the `attach` method and when writing `Recipe`s because they need to create bonds between the root atoms
    and cannot specify bonds between other atoms!


We can attach Glycans to scaffolds via the `attach` method to:
- one specific residue
- a list of residues
- all residues that match a specific sequence pattern (sequon, only for protein scaffolds)


Attaching to a specific residue
-------------------------------
Below are two examples of attaching a glycan to a protein scaffold at a specific residue.
In the first we use default "pythonistic" methods to properly attach the glycan to the protein:

.. code-block:: python

    import glycosylator as gl

    # load some Glycans from pickle files 
    # (if we have already some processed ones...)
    my_prot = gl.Scaffold.load("my_protein.pkl")
    my_glycan = gl.Glycan.load("my_glycan.pkl")

    # set the root atom of the glycan (this is required for attaching to a scaffold)
    my_glycan.root_atom = my_glycan.get_atom("C1", residue=1) # get the C1 of the first atom

    # find some ASN residue the protein 
    # (just the first ASN we can find)
    asn = my_prot.get_residue("ASN")

    # and set the root atom to be the ASN's ND2 atom
    my_prot.root_atom = my_prot.get_atom("ND2", residue=asn)

    # get a linkage
    # (glycosylator already has a linkage for N-glycosylation)
    link = gl.get_linkage("ASN-glyco")

    # now attach the glycan to the protein
    my_prot.attach(my_glycan, link)

In this second version we change we achieve the same but with an operator "short-hand" syntax because we feel "extra cocky":

.. code-block:: python

    # set the root atom on the glycan
    my_glycan ^ my_glycan.get_atom("C1", residue=1)
    
    # find an ASN residue in the protein and set its ND2 atom as root
    asn = my_prot.get_residue("ASN")
    my_prot ^ my_prot.get_atom("ND2", residue=asn)

    # set the linkage as default modifier
    my_prot % link

    # now attach the glycan to the ASN residue
    my_prot += my_glycan


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
  
This code will automatically add the glycan to all ASN residues in the list `residues_to_modify` and will use the `at_atom` argument to
connect the glycan's root atom to the ND2 atom of each ASN residue.


Attaching to residues matching a sequence pattern
-------------------------------------------------

The `Scaffold` class offers the `find` method which accepts a `regex` pattern to search 
the protein sequence for matching residues. This is useful for attaching to all residues that match a specific sequence pattern,
such as all N-glycosylation sequons. The `find` method returns a dictionary keyed by chain-id containing lists of matching residues.
For convenience, the `attach` method accepts a `sequon` argument where a regex pattern can directly be provided to do the searching automatically:

.. code-block:: python

    # attach to all N-glycosylation sequons           * sequon pattern *
    my_prot.attach(my_glycan, recipe=my_recipe, sequon="(N)(?=[^P][ST])", at_atom="ND2")

"""

from collections import defaultdict
from typing import Union


from glycosylator.core.Glycan import Glycan
import glycosylator.resources as resources
import glycosylator.utils.iupac as iupac
import glycosylator.utils.visual as visual

from buildamol.core import Molecule, Linkage, molecule
import buildamol.core.entity as entity
import buildamol.core.base_classes as base_classes
import buildamol.utils.auxiliary as aux
import buildamol.structural as structural


import re

__all__ = ["Scaffold", "scaffold", "glycosylate"]


def scaffold(scaf=None) -> "Scaffold":
    """
    Make a Scaffold object from a structure or a file.

    Parameters
    ----------
    scaf : str or object
        This can be a path to a file or a structure-like object like a Biobuild or RDKit molecule.

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
    except:
        raise ValueError("Could not make a Scaffold from the input.")


def glycosylate(
    scaffold: "Scaffold",
    glycan: "Glycan",
    link: "Linkage" = None,
    residues: list = None,
    sequon: str = None,
    inplace: bool = False,
) -> "Scaffold":
    """
    Glycosylate a Scaffold with a Glycan.

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
        self._internal_residues = set()
        self._internal_bonds = set()

    @property
    def seq(self) -> str:
        """
        Returns the sequences of the scaffold (if it is a protein)
        for each chain in the structure

        Returns
        -------
        dict
            A dictionary of the sequences of the scaffold for each chain in the structure
        """
        compounds = resources.get_default_compounds()
        seqs = {}
        for chain in self._base_struct.get_chains():
            ids = [residue.resname for residue in chain.get_residues()]
            ids = compounds.translate_ids_3_to_1(ids)
            seqs[chain] = "".join(ids)

        return seqs

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
    def glycans(self) -> dict:
        """
        Returns a dictionary of all glycan Glycans that have been attached to the scaffold
        The dictionary contains the atoms which are connected to the glycan Glycans as keys and the glycan Glycans as values.


        Returns
        -------
        dict
        """
        # double call to dict to make a copy
        return dict(self._attached_glycans)

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
        new = cls(mol._base_struct, model=mol._model)
        for bond in mol.bonds:
            new.add_bond(*bond)
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

        new._attached_glycans = defaultdict(dict)
        for root, glycan in self._attached_glycans.items():
            new_glycan = Glycan.empty(glycan.id)
            for res in glycan.get_residues():
                _res = new.get_residue(res.serial_number)
                new_glycan._chain.add(_res)
                new_glycan._AtomGraph.add_nodes_from(_res.get_atoms())
                _res.set_parent(new.get_chain(res.parent.id))

            new_glycan.set_root(new.get_atom(glycan.root_atom.serial_number))

            new_glycan._add_bonds(
                *(
                    (new.get_atom(b[0].serial_number), new.get_atom(b[1].serial_number))
                    for b in glycan.bonds
                )
            )
            new._attached_glycans[new.get_atom(root.serial_number)] = new_glycan
        return new

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

    def find_glycans(self, infer_bonds: bool = False, autolabel: bool = False) -> dict:
        """
        Find all glycan residues that are present in the scaffold structure.

        Parameters
        ----------
        infer_bonds : bool
            If True, all atomic bonds are inferred for the glycans. This is useful if not all atomic
            bonds should be inferred for an entire scaffold structure, but only for the glycans.
        autolabel : bool
            If True, the atoms within the glycans are automatically labelled according to the
            IUPAC naming scheme. This requires that bonds are available for all atoms in the glycans.
            If you are unsure, also set `infer_bonds` to True. Note however, that that due to aberrant
            bonding patterns in larger glycans the automatic labelling might not always work as expected!

        Returns
        -------
        dict
            A dictionary of all glycans that are present in the scaffold structure
        """
        residue_connections = self.get_residue_connections(triplet=False)
        if len(residue_connections) == 0:
            residue_connections = self.infer_residue_connections(triplet=False)
        residue_connections = [
            (i, i.parent, j, j.parent) for i, j in residue_connections
        ]

        ref_residues = resources.reference_glycan_residue_ids()
        glycan_residues = [r for r in self.get_residues() if r.resname in ref_residues]

        glycans = {}
        gdx = 1
        yet_to_add = set(glycan_residues)
        for a, res_a, b, res_b in residue_connections:
            # first find find all glycan "starts" by searching for residues
            # that are attached to a non-glycan residue
            if res_a not in glycan_residues and res_b in glycan_residues:
                glycan = Glycan.empty("glycan" + str(gdx))
                glycan.add_residues(res_b, adjust_seqid=False)

                glycan._add_bonds(*self.get_bonds(res_b))
                glycan.set_root(b)
                glycans[a] = glycan
                yet_to_add.remove(res_b)
                gdx += 1
                continue

            # then add all residues that are attached to a glycan residue
            # and are glycan residues themselves
            elif res_a in glycan_residues and res_b not in glycan_residues:
                glycan = Glycan.empty("glycan" + str(gdx))
                glycan.add_residues(res_a, adjust_seqid=False)
                glycan._add_bonds(*self.get_bonds(res_a))
                glycan.set_root(a)
                glycans[b] = glycan
                yet_to_add.remove(res_a)
                gdx += 1
                continue

        # then add all remaining residues to the glycans
        _counter = 0
        while len(yet_to_add) > 0 or _counter > 200:
            for a, res_a, b, res_b in residue_connections:
                if res_a in yet_to_add and res_b in yet_to_add:
                    continue
                elif res_a not in glycan_residues and res_b in glycan_residues:
                    continue
                for glycan in glycans.values():
                    if (
                        res_a in glycan.get_residues()
                        and res_b not in glycan.get_residues()
                    ):
                        glycan.add_residues(res_b, adjust_seqid=False)
                        glycan._add_bonds(*self.get_bonds(res_b))
                        glycan._add_bond(a, b)
                        yet_to_add.discard(res_b)
                        break
                    elif (
                        res_b in glycan.get_residues()
                        and res_a not in glycan.get_residues()
                    ):
                        glycan.add_residues(res_a, adjust_seqid=False)
                        glycan._add_bonds(*self.get_bonds(res_a))
                        glycan._add_bond(a, b)
                        yet_to_add.discard(res_a)
                        break
            _counter += 1

        glycans.update(self._attached_glycans)

        for root, glycan in glycans.items():
            if infer_bonds:
                glycan.infer_bonds(restrict_residues=True)
            if autolabel:
                glycan.autolabel()
            if infer_bonds:
                glycan.infer_residue_connections(triplet=False)

            # make sure to fill the glycan trees
            for bond in glycan.get_residue_connections(triplet=False, direct_by="root"):
                a, b = bond
                res_a, res_b = a.parent, b.parent
                if res_a not in glycan_residues or res_b not in glycan_residues:
                    continue
                pass
                id = iupac.make_link_id(bond)
                _link = Linkage(
                    id=id, description="<auto-generated from scaffold: labels only>"
                )
                _link.atom1 = a
                _link.atom2 = b
                glycan._glycan_tree.add(res_a, res_b, _link)

            # make sure that all glycan residues still officially belong to the scaffold
            for res in glycan.get_residues():
                res.set_parent(root.parent.parent)

        self._attached_glycans = glycans
        return glycans

    def count_glycans(self) -> int:
        """
        Returns the number of glycans that are currently attached to the scaffold

        Returns
        -------
        int
            The number of glycans
        """
        return len(self.glycans)

    def find_glycosylation_sites(self, sequon: str = "N-linked") -> list:
        """
        Find the residues at which Glycans can be attached to the scaffold

        Parameters
        ----------
        sequon : str
            A regex pattern sequon to search for, which must contain a single capturing group
            specifying the residue at which the Glycan should be attached. Defaults to the "N-linked"
            sequon which matches Asn-X-Ser/Thr sites. Alternatively, "O-linked" can be used, which matches
            Ser/Thr-X-Ser/Thr sites.

        Returns
        -------
        dict
            A dictionary of lists of matching residues in each chain.
        """

        if sequon in resources.SEQUONS:
            sequon = resources.SEQUONS[sequon]
        elif isinstance(sequon, str):
            sequon = resources.Sequon("new", sequon)

        seqs = self.seq
        sequons = {
            chain: re.finditer(sequon.pattern, seqs[chain], re.IGNORECASE)
            for chain in self.chains
            if chain not in self._excluded_chains
        }
        sequons = {
            key: [
                # adjust indexing
                self._model.child_dict[key.get_id()].child_list[0].id[1]
                + m.start(0)
                - 1
                for m in value
            ]
            for key, value in sequons.items()
            if value is not None
        }

        _residues = {}
        for chain, indices in sequons.items():
            chain = self._model.child_dict[chain.get_id()]
            if not chain.child_list:
                continue
            cdx = chain.child_list[0].id[1]
            _residues[chain] = [
                chain.child_list[idx - (cdx - 1)]
                for idx in indices
                if chain.child_list[idx - (cdx - 1)] not in self._excluded_residues
            ]
        return _residues

    def find_n_linked_sites(self) -> dict:
        """
        Find all N-linked glycosylation sites in the protein scaffold.

        Returns
        -------
        dict
            A dictionary of lists of matching residues in each chain.
        """
        return self.find_glycosylation_sites(sequon="N-linked")

    def find_o_linked_sites(self) -> dict:
        """
        Find all O-linked glycosylation sites in the protein scaffold.

        Returns
        -------
        dict
            A dictionary of lists of matching residues in each chain.
        """
        return self.find_glycosylation_sites(sequon="O-linked")

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
        sequon: Union[str, resources.Sequon] = None,
        chain=None,
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
        sequon : str or Sequon
            Instead of providing one or multiple residues at which to attach the glycan, you can also provide a
            sequon which will be used to find the residues at which to attach the glycan.
        chain : str
            The chain to which the glycan should be attached. If None, the glycan is attached to the same chain as the respective scaffold residue.
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

        if isinstance(sequon, str):
            sequon = resources.get_sequon(sequon)

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
                    _link2 = None
                    if sequon is not None:
                        _link2 = sequon.get_linkage(res.resname, _topology)
                    if _link1 is None and _link2 is None:
                        raise ValueError(
                            f"No glycosylation linkage found for residue '{res.resname}'. Either provide a linkage, provide a sequon with defined linkages, or specify a default using 'add_linkage('{res.resname}-glyco')'."
                        )
            residues = _residues

        elif not sequon and self.attach_residue is None:
            raise ValueError(
                "Either a list of residues or a sequon pattern must be provided if no attach_residue is defined."
            )
        elif self.attach_residue is not None:
            residues = [self.attach_residue]

        _topology = _topology or resources.get_default_topology()

        if not inplace:
            scaffold = self.copy()
        else:
            scaffold = self

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
                    _link = None
                    if sequon is not None:
                        _link = sequon.get_linkage(residue.resname, _topology)
                    if not _link:
                        _link = resources.get_linkage(residue.resname + "-glyco")
                else:
                    _link = link

                if each_new_chain:
                    _chain = base_classes.Chain(chr(65 + len(scaffold.chains)))
                    scaffold._model.add(_chain)
                elif not chain:
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

                scaffold.add_residues(*_mol.get_residues(), chain=chain)
                scaffold._bonds.extend(_mol.get_bonds())
                scaffold._AtomGraph.migrate_bonds(_mol._AtomGraph)
                scaffold._add_bond(s._anchors[0], s._anchors[1])

                root = scaffold.get_atom(s._anchors[0], residue=residue)
                scaffold._attached_glycans[root] = _mol

            return scaffold

        # ==========================================================
        #          attach Glycan using a sequon
        # ==========================================================
        elif sequon is not None:
            res_dict = scaffold.find_glycosylation_sites(sequon)

            for chain_id in res_dict:
                residues = res_dict[chain_id]

                scaffold = scaffold.attach(
                    mol,
                    residues=residues,
                    chain=chain,
                    sequon=sequon,
                    inplace=True,
                )
            return scaffold

    # alias
    glycosylate = attach

    def draw_snfg(
        self, glycans: bool = True, glyco_sites: bool = True, legend: bool = True
    ) -> visual.ScaffoldViewer2D:
        """
        Draw the scaffold with attached glycans and glycosylation sites in a 2D representation.

        Parameters
        ----------
        glycans : bool
            If True, the glycans are drawn
        glyco_sites : bool
            If True, the glycosylation sites are drawn
        legend : bool
            If True, a legend is drawn

        Returns
        -------
        visual.ScaffoldViewer2D
            A 2D viewer of the scaffold
        """
        viewer = visual.ScaffoldViewer2D(self)
        if glycans:
            viewer.draw_glycans()
        if glyco_sites:
            viewer.draw_glycosylation_sites(draw_legend=legend)
        return viewer

    draw2d = draw_snfg

    def show_snfg(
        self, glycans: bool = True, glyco_sites: bool = True, legend: bool = True
    ):
        """
        Show the scaffold with attached glycans and glycosylation sites in a 2D representation.

        Parameters
        ----------
        glycans : bool
            If True, the glycans are drawn
        glyco_sites : bool
            If True, the glycosylation sites are drawn
        legend : bool
            If True, a legend is drawn
        """
        viewer = self.draw_snfg(glycans, glyco_sites, legend)
        viewer.show()

    show2d = show_snfg

    draw3d = Molecule.draw
    show3d = Molecule.show

    # def py3dmol(self):
    #     import py3Dmol
    #     import os

    #     self.to_pdb(".tmp.pdb")
    #     view = py3Dmol.view(width=800, height=600)
    #     view.addModel(open(".tmp.pdb").read(), "pdb")
    #     view.setStyle({"cartoon": {"color": "spectrum"}})
    #     os.remove(".tmp.pdb")
    #     return view

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
        return f"Scaffold({self.id})"


if __name__ == "__main__":
    import glycosylator as gl

    s = Scaffold.from_pdb(
        "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    )
    s.reindex()
    s.infer_bonds()
    g = gl.glycan(
        "GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)[GalNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
    )

    s.attach(g, sequon="N-linked")
    out = s.copy()
    s.show_snfg()

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
