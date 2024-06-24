"""
The `Membrane` class is used to represent a membrane scaffold structure onto which one or more glycans can be attached.
In contrast to amino acids, lipids have a less well defined atom nomenclature which makes it difficult to provide a set of default
linkages for different lipids. However, as glycans are usually attached via a hydroxyl group, the `Membrane` class can automatically
try and infer a proper linkage to an unknown lipid residue by searching for the right hydroxyl group. Naturally, due to the inferential
nature of this process, there may be cases where a linkage still needs to be manually provided for proper results.
"""

from glycosylator.core.scaffold import Scaffold
from glycosylator.core.Glycan import Glycan
import glycosylator.resources as resources
import glycosylator.utils.visual as visual

import buildamol.structural.neighbors as neighbors

constraints = neighbors.constraints
from buildamol.core import Molecule, linkage, molecule


__all__ = ["Membrane", "membrane"]


def membrane(input) -> "Membrane":
    """
    Create a membrane from an input file or string

    Parameters
    ----------
    input : str
        The input file or string

    Returns
    -------
    Membrane
        A membrane object
    """
    try:
        s = molecule(input)
        return Membrane.from_molecule(s)
    except Exception as e:
        raise ValueError(f"Could not create a Membrane from input: {input}") from e


class Membrane(Scaffold):
    """
    The :class:`Membrane` class is used to represent a membrane scaffold structure onto which a
    modification in form of one or more :class:`glycosylator.core.Glycan.Glycan`s can be added.
    """

    def __init__(self, structure, model: int = 0) -> None:
        super().__init__(structure, model)
        self.type = "membrane"

    def attach(
        self,
        mol: "Glycan",
        link: "Linkage" = None,
        residues: list = None,
        chain=None,
        inplace: bool = True,
        _topology=None,
    ):
        if residues is None:
            if self.attach_residue is None:
                raise ValueError(
                    "No residue to attach glycan to. Provide a list of residues or set an attach_residue"
                )
            residues = [self.attach_residue]

        # first make sure we have linkages for every
        # lipid residue where the glycan should be attached
        if link is None:
            for residue in residues:
                _create_linkage_for_missing_lipid(residue, mol)

        return Scaffold.attach(
            self,
            mol,
            link=link,
            residues=residues,
            chain=chain,
            inplace=inplace,
            _topology=_topology,
        )

    glycosylate = attach

    def py3dmol(
        self,
        style: str = "sphere",
        color: str = "beige",
        glycan_style: str = "stick",
        glycan_color: str = None,
        glycans: bool = True,
    ):
        """
        Visualize the scaffold in a 3D viewer using py3Dmol.

        Parameters
        ----------
        style : str
            The style of the membrane representation. Defaults to "sphere".
        color : str
            The color of the membrane representation. Defaults to "beige".
        glycan_style : str
            The style of the glycan representation. Defaults to "stick".
        glycan_color : str
            The color of the glycan representation. Defaults to None.
        glycans : bool
            Whether to include glycans in the visualization. Defaults to True.
        """
        tmp = Molecule.empty()
        tmp.add_residues(
            *(i for i in self.get_residues() if i not in self.glycan_residues),
            adjust_seqid=False,
            _copy=True,
        )
        viewer = Scaffold.py3dmol(tmp, style, color)
        if glycans:
            for glycan in self.get_glycans().values():
                v = glycan.py3dmol(style=glycan_style, color=glycan_color)
                viewer.add(v)
        return viewer


_lipid_to_glycan_hydroxyl_link_constraints = [
    constraints.has_neighbor_hist({"C": 1, "O": 1, "H": 2}),
    constraints.extended_has_n_neighbors(2, 5),
]
"""
The constraints to help identify the hydroxyl group that is used to attach
glycans
"""


def _create_linkage_for_missing_lipid(lipid, glycan):
    """
    Create a "{}-glyco" linkage for a lipid that does not
    yet have one
    """
    if resources.get_linkage(lipid.resname + "-glyco"):
        return
    _lipid = Molecule.empty()
    _lipid.add_residues(lipid, _copy=True)
    m = _lipid.search_by_constraint(_lipid_to_glycan_hydroxyl_link_constraints)
    if len(m):
        m = m[0]
        lipid_binder = m[1]
        lipid_delete = _lipid.get_hydrogen(lipid_binder)
        if lipid_delete is None:
            raise ValueError(
                f"Could not find a suitable hydroxyl group to attach a glycan to lipid {lipid}."
            )

        glycan_binder = glycan.root_atom or "C1"
        glycan_deletes = None
        if glycan.root_atom:
            O = next(
                iter(
                    glycan.get_neighbors(
                        glycan.root_atom, filter=lambda x: x.element == "O"
                    )
                ),
                None,
            )
            H = None if not O else glycan.get_hydrogen(O)
            if O and H:
                glycan_deletes = (O, H)
            elif O:
                glycan_deletes = (O,)

        if not glycan_deletes:
            glycan_deletes = ("O1", "HO1")

        link = linkage(
            lipid_binder,
            glycan_binder,
            (lipid_delete,),
            glycan_deletes,
            id=lipid.resname + "-glyco",
        )
        resources.add_linkage(link)
    else:
        raise ValueError(
            f"Could not find a suitable hydroxyl group to attach a glycan to lipid {lipid}."
        )


if __name__ == "__main__":
    f = "/Users/noahhk/GIT/glycosylator/docs/source/examples/files/membrane.pkl"
    m = Membrane.load(f)
    g = Glycan.from_iupac(None, "Gal(b1-4)GlcNAc")

    m.attach(g, residues=m.get_residues("CER", by="name", chain="U"), chain="new")
    m.to_pdb("membrane_glycan.pdb")
