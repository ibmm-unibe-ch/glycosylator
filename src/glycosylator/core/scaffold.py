"""
The `Scaffold` class is used to represent a scaffold structure such as a protein or membrane on to which a modification in form of one or more :class:`glycosylator.core.molecule.Molecule`s can be added.

"""

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

        self._molecule_residues = set()
        self._internal_residues = []
        self._molecule_connections = set()

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

    def exclude_residue(self, residue: Union[int, str, bio.Residue.Residue]):
        """
        Exclude a residue of the scaffold from modification
        by a molecule

        Parameters
        ----------
        residue : int or str or Residue
            The residue to exclude
        """
        res = self.get_residue(residue)
        if res is not None:
            self._excluded_residues.add(res)
        else:
            raise ValueError(f"Residue '{residue}' not found")

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
                self.residues[idx]
                for idx in indices
                if self.residues[idx] not in self._excluded_residues
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

        to_remove = [
            res
            for res in self.residues
            if res not in residues and res not in self._molecule_residues
        ]

        self._internal_residues = self.remove_residues(*to_remove)

    def fill(self):
        """
        Fill the structure with the residues that were removed by the `hollow_out` method.
        """
        self.add_residues(*self._internal_residues, chain=self.chains[0])

    def attach(
        self,
        mol: "Molecule",
        remove_atoms: tuple,
        mol_remove_atoms: tuple,
        sequon: str = None,
        at_atom: str = None,
        chain=None,
        _copy: bool = False,
    ):
        """
        Attach a molecule to the scaffold via the root atoms.

        Parameters
        ----------
        mol : Molecule
            The molecule to attach
        remove_atoms : tuple
            The atoms to remove from the scaffold while stitching the molecule to it. These must be the atom ids (e.g. "HO4")
            and they must be part of the scaffold's root residue.
        mol_remove_atoms : tuple
            The atoms to remove from the molecule while stitching it to the scaffold. These must be the atom ids (e.g. "HO4")
            and they must be part of the molecule's root residue.
        sequon : str
            Instead of providing one or multiple residues at which to attach the molecule, you can also provide a
            sequon, which will be used to find the residues at which to attach the molecule.
        at_atom: str
            In case a sequon is provided, this specifies the atom of the matching residues that shall be used for stitching.
            This must be the atom's id (e.g. "CA").

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

        if _copy:
            scaffold = deepcopy(self)
        else:
            scaffold = self

        if sequon is not None:
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

            res_dict = scaffold.find(sequon)
            for chain_id in res_dict:
                residues = res_dict[chain_id]
                for res in residues:
                    root = res.child_dict.get(at_atom)
                    if root is None:
                        raise ValueError(
                            f"No atom with id '{at_atom}' found in residue {res}"
                        )

                    scaffold.root_atom = root
                    scaffold.attach_residue = scaffold.root_residue
                    scaffold = scaffold.attach(
                        mol,
                        remove_atoms,
                        mol_remove_atoms,
                        chain=_chain,
                        _copy=False,
                    )

            scaffold.root_atom = _root
            scaffold.attach_residue = _attach_residue
            return scaffold

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
        connections = set((i.serial_number, j.serial_number) for i, j in connections)
        scaffold._molecule_connections.update(connections)

        scaffold.add_residues(*_mol.residues, chain=chain)
        scaffold._bonds.extend(_mol.bonds)
        scaffold._AtomGraph.add_edges_from(_mol._AtomGraph.edges)
        scaffold.add_bond(s._anchors[0].serial_number, s._anchors[1].serial_number)

        return scaffold

    def __repr__(self) -> str:
        return f"Scaffold({self.id})"


if __name__ == "__main__":
    import glycosylator as gl

    f1 = "/Users/noahhk/GIT/glycosylator/support/examples/4tvp.prot.pdb"
    s = Scaffold.from_pdb(f1)
    s.reindex()
    s.infer_bonds(restrict_residues=True)

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
    s.to_pdb("final_scaffold_superduper999.pdb")
    print(f"Time: {t2 - t1}")
    pass
