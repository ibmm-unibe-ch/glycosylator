"""
Attaching to residues matching a sequence pattern
-------------------------------------------------

The `Protein` class offers the `find_glycosylation_sites` method which accepts a `regex` pattern to search 
the protein sequence for matching residues. This is useful for attaching to all residues that match a specific sequence pattern,
such as all N-glycosylation sequons. The `find_glycosylation_sites` method returns a dictionary keyed by chain-id containing lists of matching residues.
For convenience, the `Protein`'s `attach` method accepts a `sequon` argument where a regex pattern can directly be provided to do the searching automatically.
For this to work, the sequon must contain a single capturing group that specifies the amino acid to glycosylate (in the example below the `(N)`) and a `link` must be provided
if not a corresponding `{}-glyco` linkage is defined in the default linkages (e.g. here an `ASN-glyco` since `N` is one letter code for Asparagine, which is available by default).

.. code-block:: python

    import glycosylator as gl 

    my_prot = gl.protein("my_protein.pkl")
    my_glycan = gl.glycan("my_glycan.pkl")

    # attach to all N-glycosylation sites

    # (using the default N-linked sequon)
    my_prot.attach(my_glycan, sequon="N-linked")

    # (or by providing a regex pattern manually)
    my_prot.attach(my_glycan, sequon="(N)(?=[^P][ST])")

If you have non-canonical glycosylation that you wish to perform on an amino acid other than ASN, SER, or THR, or you wish to perform a special kind of
connection between the amino acid and the glycan, you can also provide a `link` that contains the necessary instructions.

.. code-block:: python

    # first define a linkage to be used when glycosylating
    # (this is actually the default ASN-glyco linkage)
    link = gl.linkage("ND2", "C1", ["HD22"], ["O1", "HO1"])

    # get all N-glycosylation sites in chain A
    glyco_sites = my_prot.find_glycosylation_sites("(N)(?=[^P][ST])")
    glyco_sites = glyco_sites["A"] # only the sites in chain A

    # attach the glycan to the sites
    my_prot.attach(my_glycan, link=link, residues=glyco_sites)
"""

from typing import Union

from glycosylator.core.scaffold import Scaffold
from glycosylator.core.Glycan import Glycan
import glycosylator.resources as resources
import glycosylator.utils.visual as visual


from buildamol.core import Linkage, molecule, base_classes


import re

__all__ = ["Protein", "protein"]


def protein(input) -> "Protein":
    """
    Create a Protein from a file or a string.

    Parameters
    ----------
    input : str
        The input file or string

    Returns
    -------
    Protein
        The Protein object
    """
    try:
        s = molecule(input)
        return Protein.from_molecule(s)
    except Exception as e:
        raise ValueError(f"Could not create a Protein from input: {input}") from e


class Protein(Scaffold):
    """
    The :class:`Protein` class is used to represent a protein scaffold structure onto which a
    modification in form of one or more :class:`glycosylator.core.Glycan.Glycan`s can be added.
    """

    def __init__(self, structure, model: int = 0) -> None:
        super().__init__(structure, model)
        self.type = "protein"

    @property
    def seq(self) -> str:
        """
        Returns the sequences of the protein
        for each chain in the structure

        Returns
        -------
        str
            The sequence of the structure. Different chains are delimited by a ":" character.
        """
        ":".join(self.get_sequence().values())

    def get_sequence(
        self, chains: Union[str, list] = None, ignore_unknown: bool = True
    ) -> Union[str, dict]:
        """
        Returns the sequences of the scaffold (if it is a protein)
        for each chain in the structure

        Parameters
        ----------
        chains : str or list
            The chain(s) for which to get the sequence. If None, all chains are considered.
        ignore_unknown : bool
            Whether to ignore unknown residues in the sequence. If False, unknown residues are
            included as an 'X' in the sequence.

        Returns
        -------
        str or dict
            If a single chain is provided as a string, the string of the sequence, otherwise
            a dictionary of the sequences of the scaffold for each chain in the structure
        """
        compounds = resources.get_default_compounds()
        seqs = {}
        if chains is None:
            chain_iter = self._base_struct.get_chains()
        elif isinstance(chains, (str, base_classes.Chain)):
            chain_iter = iter((self.get_chain(chains),))
        elif isinstance(chains, (tuple, list)):
            chain_iter = iter([self.get_chain(c) for c in chains])
        else:
            raise ValueError(
                f"Invalid chains argument. Must be a string, list or None, but got: {chains}"
            )
        for chain in chain_iter:
            ids = [residue.resname for residue in chain.get_residues()]
            ids = compounds.translate_ids_3_to_1(ids, ignore_unknown=ignore_unknown)
            seqs[chain] = "".join(ids)

        if isinstance(chains, str):
            return seqs[chains]

        return seqs

    def find_glycosylation_sites(self, sequon: str = "N-linked") -> list:
        """
        Find the residues at which Glycans can be attached to the protein

        Parameters
        ----------
        sequon : str
            A regex pattern sequon to search for, which must contain a single capturing group
            specifying the residue at which the Glycan should be attached. Defaults to the "N-linked"
            sequon which matches Asn-X-Ser/Thr sites. Alternatively, "O-linked" can be used, which matches
            Ser/Thr-X-Ser/Thr sites. If a custom sequon is provided, it must contain a single capturing group
            specifying the residue at which the Glycan should be attached. To search for all registered sequons
            use "all".

        Returns
        -------
        dict
            A dictionary of lists of matching residues in each chain.
        """
        if sequon == "all":
            sequons = resources.SEQUONS
            sequons = [sequons[key] for key in sequons]
            sequons = [sequon.pattern for sequon in sequons]
            sequon = "|".join(sequons)

        if sequon in resources.SEQUONS:
            sequon = resources.SEQUONS[sequon]
        elif isinstance(sequon, str):
            sequon = resources.Sequon("new", sequon)

        seqs = self.get_sequence()
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
        if isinstance(sequon, str):
            sequon = resources.get_sequon(sequon)

        if sequon is not None:
            scaffold = self if inplace else self.copy()
            res_dict = scaffold.find_glycosylation_sites(sequon)
            for chain_id in res_dict:
                residues = res_dict[chain_id]
                scaffold = Scaffold.attach(
                    scaffold,
                    mol,
                    link=link,
                    residues=residues,
                    chain=chain,
                    inplace=True,
                    _topology=_topology,
                )
        else:
            scaffold = Scaffold.attach(
                self,
                mol,
                link=link,
                residues=residues,
                chain=chain,
                inplace=inplace,
                _topology=_topology,
            )
        return scaffold

    glycosylate = attach

    def snfg(
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

    draw2d = snfg

    def py3dmol(
        self,
        style: str = "cartoon",
        color: str = "spectrum",
        glycan_style: str = "stick",
        glycan_color: str = None,
        glycans: bool = True,
    ):
        """
        Visualize the scaffold in a 3D viewer using py3Dmol.

        Parameters
        ----------
        style : str
            The style of the protein representation. Defaults to "cartoon".
        color : str
            The color of the protein representation. Defaults to "spectrum".
        glycan_style : str
            The style of the glycan representation. Defaults to "stick".
        glycan_color : str
            The color of the glycan representation. Defaults to None.
        glycans : bool
            If True, the glycans are drawn.
        """
        viewer = super().py3dmol(style, color)
        if glycans:
            for glycan in self.get_glycans().values():
                v = glycan.py3dmol(style=glycan_style, color=glycan_color)
                viewer.add(v)
        return viewer


if __name__ == "__main__":
    import glycosylator as gl

    s = Protein.from_pdb(
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
    import matplotlib.pyplot as plt

    plt.show()

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
