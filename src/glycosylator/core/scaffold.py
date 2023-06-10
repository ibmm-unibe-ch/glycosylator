"""
The `Scaffold` class is used to represent a scaffold structure such as a protein or membrane on to which a modification in form of one or more `Molecule`s can be added.
As a class the `Scaffold` shares many features with the `Molecule` class, but it is not a subclass of it but a sister class. 

As such it lacks certain features such as the `repeat` or `from_compound` methods. 
A Scaffold can be edited just like a Molecule, so refer to the `Molecule` documentation for more information on how to do that.


Attaching Molecules
===================

The primary purpose of the `Scaffold` is to have `Molecule`s attached to it.
This is done by using the `attach` method, which offers a versatile interface to attach one or more `Molecule`s to the `Scaffold` structure.

.. note:: 

    Connecting a `Molecule` and `Scaffold` is based on `Recipe`s (not on `Patch`es!).
    Furthermore, "attaching" a `Molecule` to a `Scaffold` is defined as "connecting the Molecule's and Scaffold's root atoms".
    In practice, (re-)setting the objects' root atoms is handled behind the scenes, but is is important to keep this in
    mind when using the `attach` method and when writing `Recipe`s because they need to create bonds between the root atoms
    and cannot specify bonds between other atoms!


We can attach molecules to scaffolds via the `attach` method to:
- one specific residue
- a list of residues
- all residues that match a specific sequence pattern (sequon, only for protein scaffolds)


Attaching to a specific residue
-------------------------------
Below are two examples of attaching a glycan to a protein scaffold at a specific residue.
In the first we use default "pythonistic" methods to properly attach the glycan to the protein:

.. code-block:: python

    import glycosylator as gl

    # load some molecules from pickle files
    my_prot = gl.Scaffold.load("my_protein.pkl")
    my_glycan = gl.Molecule.load("my_glycan.pkl")

    # set the root atom of the glycan (this is required for attaching to a scaffold)
    my_glycan.root_atom = my_glycan.get_atom("C1", residue=1) # get the C1 of the first atom

    # find some ASN residue the protein 
    # (just the first ASN we can find)
    asn = my_prot.get_residue("ASN")

    # make a recipe for attaching to ASN
    my_recipe = gl.Recipe(id="glycan_to_asn")
    # ... add intstructions here ... 
    
    # now attach the glycan to the ASN residue
    # (we use a list of the one residue here)
    my_prot.attach(my_glycan, recipe=my_recipe, residues=[asn], at_atom="ND2")


In the second example we make use of the fact that "attaching" always connects the 
root atoms, which allows us to use the `attach` method without specifying the `at_atom` argument:

.. code-block:: python

    my_glycan.root_atom = my_glycan.get_atom("C1", residue=1)

    # set the root atom on the protein as ND2 of the ASN residue
    asn = my_prot.get_residue("ASN")
    my_prot.root_atom = my_prot.get_atom("ND2", residue=asn)

    # now we can directly attach
    my_prot.attach(my_glycan, my_recipe)


And in this final version we change the second example to use the operator "short-hand" syntax because we feel "extra cocky":

.. code-block:: python

    # set the root atom on the glycan
    my_glycan ^ my_glycan.get_atom("C1", residue=1)
    
    # find an ASN residue in the protein and set its ND2 atom as root
    asn = my_prot.get_residue("ASN")
    my_prot ^ my_prot.get_atom("ND2", residue=asn)

    # set the recipe as default modifier
    my_prot % my_recipe

    # now attach the glycan to the ASN residue
    my_prot += my_glycan


Attaching to a list of residues
-------------------------------

Often we want to modify a number of residues in the same way, for example when we want to attach a glycan to all
residues in a protein that match a specific sequence pattern (e.g. all N-glycosylation sequons). We can do this by providing
a list of `residues` to the `attach` method:

.. code-block:: python

    # find all ASN residues in the protein
    residues_to_modify = [res for res in my_prot.get_residues() if res.resname == "ASN"]
    
    # now attach the glycan to all ASN residues
    my_prot.attach(my_glycan, recipe=my_recipe, residues=residues_to_modify, at_atom="ND2")
  
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
from copy import deepcopy
from typing import Union

import Bio.PDB as bio

import glycosylator.resources as resources
import glycosylator.utils.defaults as defaults
import glycosylator.structural as structural
import glycosylator.utils.auxiliary as aux
import glycosylator.core.entity as entity

import re


class Scaffold(entity.BaseEntity):
    """
    The :class:`Scaffold` class is used to represent a scaffold structure such as a protein or membrane on to which a
    modification in form of one or more :class:`glycosylator.core.molecule.Molecule`s can be added.
    """

    def __init__(self, structure, model: int = 0) -> None:
        super().__init__(structure, model)

        self._excluded_chains = set()
        self._excluded_residues = set()

        self._attached_molecules = defaultdict(dict)
        self._internal_residues = set()

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
        compounds = defaults.get_default_compounds()
        seqs = {}
        for chain in self._base_struct.get_chains():
            ids = [residue.resname for residue in chain.get_residues()]
            ids = compounds.translate_ids_3_to_1(ids)
            seqs[chain.get_id()] = "".join(ids)

        return seqs

    @property
    def attached_molecules(self) -> dict:
        """
        Returns a dictionary of all molecules that have been attached to the scaffold
        The dictionary contains the ids of molecules as keys and dictionaries of the associated
        scaffold residues and the molecule's residues and bonds as values.

        Returns
        -------
        dict
            A dictionary of all molecules that have been attached to the scaffold
        """
        return self._attached_molecules

    @property
    def internal_residues(self) -> set:
        """
        Returns a set of all residues that have been labelled as "internal" based on
        Solvent Accessible Surface Area (SASA) calculations. These residues are excluded from the main
        scaffold structure by default.

        Returns
        -------
        set
            A set of all residues that have been added to the scaffold
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

    def exclude_chain(self, chain: Union[str, bio.Chain.Chain]):
        """
        Exclude a chain from the scaffold from letting any
        molecules attach to any of its residues. This will also
        prevent any sequon matches to be found by the `find` method.

        Parameters
        ----------
        chain : str or Chain
            The chain to exclude
        """
        if isinstance(chain, str):
            _chain = self._model.child_dict.get(chain)
        if _chain is not None:
            self._excluded_chains.add(_chain)
        else:
            raise ValueError(f"Chain '{chain}' not found")

    def include_chain(self, chain: Union[str, bio.Chain.Chain]):
        """
        Include a previously excluded chain in the scaffold again to let molecules
        attach to its residues. This will also
        allow sequon matches to be found by the `find` method.

        Parameters
        ----------
        chain : str or Chain
            The chain to include
        """
        if isinstance(chain, str):
            _chain = self._model.child_dict.get(chain)
        if _chain is None:
            raise ValueError(f"Chain '{chain}' not found")
        elif _chain in self._excluded_chains:
            self._excluded_chains.remove(_chain)

    def exclude_residue(self, residue: Union[int, bio.Residue.Residue]):
        """
        Exclude a residue of the scaffold from modification
        by a molecule

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

    def include_residue(self, residue: Union[int, bio.Residue.Residue]):
        """
        Include a previously excluded residue of the scaffold again for modification
        by a molecule

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

    def find(self, sequon: str = "N-linked") -> list:
        """
        Find the residues at which molecules can be attached to the scaffold

        Parameters
        ----------
        sequon : str
            A regex pattern sequon to search for, which must contain a single capturing group
            specifying the residue at which the molecule should be attached. Defaults to the "N-linked"
            sequon which matches Asn-X-Ser/Thr sites. Alternatively, "O-linked" can be used, which matches
            Ser/Thr-X-Ser/Thr sites.

        Returns
        -------
        list
            A list of all residues that match the sequon
        """

        if sequon in resources.SEQUONS:
            sequon = resources.SEQUONS[sequon]

        seqs = self.seq
        sequons = {
            chain.get_id(): re.finditer(sequon, seqs[chain.get_id()], re.IGNORECASE)
            for chain in self.chains
            if chain not in self._excluded_chains
        }
        sequons = {
            key: [
                # adjust indexing
                self._model.child_dict[key].child_list[0].id[1] + m.start(0) - 1
                for m in value
            ]
            for key, value in sequons.items()
            if value is not None
        }

        _residues = {}
        for id, indices in sequons.items():
            chain = self._model.child_dict[id]
            cdx = chain.child_list[0].id[1]
            _residues[id] = [
                chain.child_list[idx - (cdx - 1)]
                for idx in indices
                if chain.child_list[idx - (cdx - 1)] not in self._excluded_residues
            ]
        return _residues

    def add_residues(
        self,
        *residues: bio.Residue.Residue,
        chain: str = None,
        adjust_seqid: bool = True,
        _copy: bool = False,
    ):
        """
        Add residues to the structure

        Parameters
        ----------
        residues : bio.Residue.Residue
            The residues to add
        chain : str
            The chain to add the residues to. If None, the residues are added to the last chain
        adjust_seqid : bool
            If True, the seqid of the residues is adjusted to
            match the current number of residues in the structure
            (i.e. a new residue can be given seqid 1, and it will be adjusted
            to the correct value of 3 if there are already two other residues in the molecule).
        _copy : bool
            If True, the residues are copied before adding them to the molecule.
            This is useful if you want to add the same residue to multiple molecules, while leaving
            them and their original parent structures intakt.
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
        the computational load for structure optimization of attached molecules.

        Parameters
        ----------
        cutoff : int or float
            The cutoff SASA value for residues to be classfied as "surface residues".
            If an integer is provided, the value is interpreted as the absolute threshold
            for the SASA value. If a float is provided, the value is interpreted as a
            percentage of the maximal SASA value in the structure.
        """
        # get surface residues
        structure = aux.DummyStructure(
            [res for res in self.residues if res not in self._molecule_residues]
        )
        residues = structural.infer_surface_residues(structure, cutoff=cutoff)

        # works fine with a list-comprehension...
        to_remove = (
            res
            for res in self.residues
            if res not in residues and res not in self._molecule_residues
        )

        for residue in to_remove:
            chain = residue.get_parent()
            chain.detach_child(residue.id)
            residue.set_parent(chain)

            for atom in residue.get_atoms():
                self._AtomGraph.remove_node(atom)
                self._purge_bonds(atom)

            self._internal_residues.add(residue)

    def fill(self):
        """
        Fill the structure with the residues that were removed by the `hollow_out` method.
        """
        for residue in self._internal_residues:
            chain = residue.get_parent()
            chain.add(residue)
            for atom in residue.get_atoms():
                self._AtomGraph.add_node(atom)
                bonds = self._get_bonds(atom, None)
                self._AtomGraph.add_edges_from(bonds)

        self._internal_residues.clear()

    def to_pdb(self, filename: str):
        self.fill()  # fill the structure before writing it to a file
        return super().to_pdb(filename)

    def make_residue_graph(self, detailed: bool = False, locked: bool = True):
        graph = super().make_residue_graph(detailed, locked)
        # also add all remaining non connected (scaffold residues)
        for residue in self.residues:
            residue.coord = residue.center_of_mass()
            if residue not in graph:
                graph.add_node(residue)
        return graph

    def attach(
        self,
        mol: "Molecule",
        recipe: "AbstractRecipe" = None,
        remove_atoms: tuple = None,
        mol_remove_atoms: tuple = None,
        residues: list = None,
        sequon: str = None,
        at_atom: str = None,
        chain=None,
        _copy: bool = False,
    ):
        """
        Attach a molecule to the scaffold via their root atoms.

        Parameters
        ----------
        mol : Molecule
            The molecule to attach
        recipe : Recipe
            The recipe to use when stitching. If None, the default recipe that was set earlier on the molecule is used (if defined).
        remove_atoms : tuple
            The atoms to remove from the scaffold while stitching the molecule to it. These must be the atom ids (e.g. "HO4")
            and they must be part of the scaffold's root residue.
        mol_remove_atoms : tuple
            The atoms to remove from the molecule while stitching it to the scaffold. These must be the atom ids (e.g. "HO4")
            and they must be part of the molecule's root residue.
        residues : list
            A list of residues at which to attach copies of the molecule. This is an alternative to using a sequon.
        sequon : str
            Instead of providing one or multiple residues at which to attach the molecule, you can also provide a
            sequon, which will be used to find the residues at which to attach the molecule.
        at_atom: str
            In case a sequon or a list of residues are provided, this specifies the atom of the matching residues that shall be used for stitching.
            This must be the atom's id (e.g. "CA"). If no sequon and no residue list is provided, this parameter is ignored and the root atom of the scaffold is used for attachment.
        chain : str or int or bio.PDB.Chain.Chain
            The chain to which the molecule should be attached. If None, the molecule is attached to the same chain as the root residue.
            This can be set to "new" to create a new chain for the molecule, or "each" to create a new chain for each molecule that is attached
            if a sequon is provided. If a sequon is provided and this is set to "new" all molecules will be attached to the same new chain.

        _copy : bool
            Whether to copy the scaffold before attaching the molecule(s) at the specified residue(s).

        Returns
        -------
        scaffold
            The scaffold with the molecule(s) attached, which is either a copy of the original scaffold or the original.
        """
        # ==========================================================
        #                    check input parameters
        # ==========================================================
        if not recipe and not remove_atoms:
            if not self._patch:
                raise AttributeError(
                    "No recipe was set for this molecule and no manual instructions were found. Either set a default recipe, provide a recipe when stitching, or provide the information about removed and bonded atoms directly."
                )
            recipe = self._patch

        if recipe:
            return self.attach(
                mol,
                remove_atoms=recipe.deletes[0],
                mol_remove_atoms=recipe.deletes[1],
                residues=residues,
                sequon=sequon,
                at_atom=at_atom,
                chain=chain,
                _copy=_copy,
            )

        if _copy:
            scaffold = deepcopy(self)
        else:
            scaffold = self

        # ==========================================================
        #          attach molecule to multiple residues
        # ==========================================================

        if sequon is not None or residues is not None:
            _root = scaffold.root_atom
            _attach_residue = scaffold.attach_residue

            if not isinstance(at_atom, str):
                raise ValueError(
                    f"When providing a sequon, you must also provide an 'at_atom' of type str (got {type(at_atom)})"
                )

            _chain = chain
            if chain == "new":
                _chain = bio.Chain.Chain(chr(65 + len(scaffold.chains)))
                self._model.add(_chain)
            if chain == "each":
                _chain = "new"

            # ----------------------------------------------------------
            #                attach molecule using a sequon
            # ----------------------------------------------------------
            if sequon is not None:
                res_dict = scaffold.find(sequon)
                for chain_id in res_dict:
                    residues = res_dict[chain_id]
                    scaffold = scaffold.attach(
                        mol,
                        remove_atoms=remove_atoms,
                        mol_remove_atoms=mol_remove_atoms,
                        residues=residues,
                        at_atom=at_atom,
                        chain=_chain,
                        _copy=False,
                    )
                scaffold.root_atom = _root
                scaffold.attach_residue = _attach_residue
                return scaffold

            # ----------------------------------------------------------
            #           attach molecule using a list of residues
            # ----------------------------------------------------------
            elif residues is not None:
                for res in residues:
                    if res not in scaffold.get_residues():
                        res = scaffold.get_residue(res)
                    root = res.child_dict.get(at_atom)
                    if root is None:
                        raise ValueError(
                            f"No atom with id '{at_atom}' found in residue {res}"
                        )

                    scaffold.root_atom = root
                    scaffold.attach_residue = scaffold.root_residue
                    scaffold = scaffold.attach(
                        mol,
                        remove_atoms=remove_atoms,
                        mol_remove_atoms=mol_remove_atoms,
                        chain=_chain,
                        _copy=False,
                    )

            scaffold.root_atom = _root
            scaffold.attach_residue = _attach_residue
            return scaffold

        # ==========================================================
        #          attach molecule to a single residue
        # ==========================================================

        if self.root_atom is None:
            raise ValueError("No root atom defined on the scaffold")
        elif mol.root_atom is None:
            raise ValueError("No root atom defined on the molecule")

        if not isinstance(remove_atoms, (list, tuple)):
            raise ValueError(
                f"'remove_atoms' must be a list or tuple (got {type(remove_atoms)})"
            )
        if not isinstance(mol_remove_atoms, (list, tuple)):
            raise ValueError(
                f"'mol_remove_atoms' must be a list or tuple (got {type(mol_remove_atoms)})"
            )

        s = structural.__default_keep_copy_stitcher__

        if chain is None:
            chain = scaffold.root_residue.get_parent()
        if chain == "new":
            chain = bio.Chain.Chain(chr(65 + len(scaffold.chains)))
            self._model.add(chain)

        scaffold, _mol = s.apply(
            scaffold,
            mol,
            remove_atoms,
            mol_remove_atoms,
            scaffold.root_atom,
            mol.root_atom,
        )

        connections = _mol.infer_residue_connections(triplet=True)
        connections.append(s._anchors)
        scaffold._attached_molecules[_mol.id][scaffold.root_atom.get_parent()] = (
            _mol.residues,
            _mol.bonds,
            connections,
        )

        scaffold.add_residues(*_mol.residues, chain=chain)
        scaffold._bonds.extend(_mol.bonds)
        scaffold._AtomGraph.add_edges_from(_mol._AtomGraph.edges)
        scaffold.add_bond(s._anchors[0].serial_number, s._anchors[1].serial_number)

        return scaffold

    def __add__(self, mol):
        """
        Attach a molecule to the scaffold via their root atoms.
        """
        return self.attach(mol, _copy=True)

    def __iadd__(self, mol):
        """
        Attach a molecule to the scaffold via their root atoms.
        """
        return self.attach(mol, _copy=False)

    def __repr__(self) -> str:
        return f"Scaffold({self.id})"


if __name__ == "__main__":
    import glycosylator as gl

    f1 = "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    s = Scaffold.from_pdb(f1)
    s.reindex()
    s.infer_bonds(restrict_residues=True)

    sequon = "(N)(?=[A-OQ-Z][ST])"
    residues = s.find(sequon)

    # f = "/Users/oahhk/GIT/glycosylator/test.prot"
    # s.save(f)

    # s = Scaffold.load(f)
    # _s = Scaffold(s.structure)
    # _s.bonds = s.bonds
    # s = _s

    mol = gl.Molecule.from_pdb(
        "/Users/noahhk/GIT/glycosylator/support/examples/man9.pdb"
    )
    mol.infer_bonds(restrict_residues=False)
    mol.reindex()

    mol.root_atom = 1

    import time

    t1 = time.time()
    s.attach(
        mol,
        remove_atoms=("HD22",),
        mol_remove_atoms=("HO1", "O1"),
        at_atom="ND2",
        sequon="N-linked",
    )
    t2 = time.time()
    print("Oh yeah - it works!")
    from datetime import datetime

    s.to_pdb(f"final_scaffold_superduper.{datetime.now()}.pdb")
    print(f"Time: {t2 - t1}")
    pass
