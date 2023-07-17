"""
The basic glycan molecule class
"""
from typing import Union
import warnings

import biobuild.core as core
import biobuild.structural as structural
import biobuild.resources as resources

import glycosylator.utils as utils


def glycan(id: str, g: Union[str, list] = None):
    """
    The toplevel function to generate an entire glycan molecule from either an IUPAC/SNFG string, a list of residues from a graph structure, or just get a single residue glycan molecule (e.g. one Glucose, etc.).

    Parameters
    ----------
    id : str
        The id of the molecule to create.
        If no further arguments are passed, this function is used like the `biobuild.molecule` function to obtain a single residue molecule (e.g. `glycan("glucose")`).

    g : str or list
        The glycan string to parse.
        The string must be in IUPAC condensed format - currently, neither extended nor short formats are supported (refer to the `read_snfg` function for more information).
        Alternatively, a list of residues from a graph structure can be passed (refer to the `read_graph` function for more information).

    Returns
    -------
    molecule : Glycan
        The created Glycan molecule.
    """
    if g is None:
        mol = core.molecule(id)
        mol = Glycan(mol)
        return mol

    if isinstance(g, str):
        try:
            return read_snfg(id, g)
        except:
            try:
                mol = bb.molecule(g)
                mol.id = id
                mol = Glycan(mol)
                return mol
            except:
                raise ValueError(
                    f"Failed to interpret input '{g}'. Could not parse as IUPAC/SNFG string and not get a molecule from it. If it is a IUPAC/SNFG string, perhaps one of the residues could not be found in the database? Try building the glycan directly."
                )
    elif isinstance(g, list):
        return read_graph(id, g)


def read_snfg(id: str, snfg: str, _topology=None) -> "Glycan":
    """
    Make a molecule from an IUPAC/SNFG glycan string in condensed format.

    Note
    ----


    Parameters
    ----------
    id : str
        The id of the molecule to create
    g : str
        The glycan string to parse.
        The string must be in IUPAC condensed format - currently, neither extended nor short formats are supported.
    _topology
        A particular topology to use. If None, the default topology is used.

    Returns
    -------
    molecule : Glycan
        The created Glycan molecule.

    Examples
    --------
    To generate a small glycan of the structure:
    ```
    ~ --- NAG                  MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN
    ```
    the IUPAC/SNFG string would be:
    >>> iupac = "Man(a1-6)[Man(a1-3)]b-Man(a1-4)GlcNAc(b1-4)GlcNAc(b1-" # notice the final "b1-" to indicate where the glycan attaches to a scaffold

    Which can be parsed into a molecule with:
    >>> mol = read_snfg("my_glycan", iupac)
    """
    if isinstance(g, str):
        g = utils.IUPACParser().parse(snfg)
    mol = _parse_iupac_graph(id, g, _topology)
    return mol


def read_graph(id: str, g: list, _topology=None) -> "Glycan":
    """
    Build a molecule from a glycan graph.

    Parameters
    ----------
    id : str
        The id of the molecule to create
    g : list
        A list of tuples in the form (parent, child, linkage), where parent and child are strings
        designating the residues and an "@{id}" suffix to distinguish individual residues. The linkage must be
        a valid id of any defined linkage in the provided or default topology, e.g. "14bb" or "16ab" for the default
        CHARMM topology (see example below).
    _topology
        A particular topology to use. If None, the default topology is used.

    Returns
    -------
    molecule : Glycan
        The created Glycan molecule.

    Examples
    --------
    To generate a small glycan of the structure:
    ```
    ~ --- NAG                  MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN
    ```
    We can formulate a graph structure for the glycan above as follows:
    >>> graph = [
    ("NAG@1", "NAG@2", "14bb"),
    ("NAG@2", "BMA@1", "14bb"),
    ("BMA@1", "MAN@1", "13ab"),
    ("BMA@1", "MAN@2", "16ab"),
    ] # notice the @{id} after each residue name 

    The `@{id}` suffix is used to distinguish between residues with the same name, for example the two mannoses (`MAN@1` and `MAN@2`).
    The in the example above the ids reflect the number of residues with the same name, hence `NAG@2` connects to `BMA@1` "the second NAG connecting to the first BMA".
    However, this is not a strict requirement. Any numeric or string value that will mark each residue as a unique node is acceptable - that is, each combination of one particuar residue is identified by a unique "{name}@{id}".
    Hence, also the following graph is valid where the index simply reflects the order of the residues in the molecule: 
    
    >>> graph = [
    ("NAG@1", "NAG@2", "14bb"),
    ("NAG@2", "BMA@3", "14bb"),
    ("BMA@3", "MAN@4", "13ab"),
    ("BMA@3", "MAN@5", "16ab"),
    ]
    or even 
    >>> graph = [
    ("NAG@a", "NAG@b", "14bb"),
    ("NAG@b", "BMA@c", "14bb"),
    ("BMA@c", "MAN@d", "13ab"),
    ("BMA@c", "MAN@e", "16ab"),
    ] # here the ids are simply letters

    we can then create a molecule using:
    >>> mol = read_graph("my_glycan", graph)
    """
    if not isinstance(g, list):
        raise ValueError("g must be either a string or a list")
    if not isinstance(g[0], (tuple, list)) or len(g[0]) != 3:
        raise ValueError("g must be a list of tuples of length 3")
    mol = _parse_iupac_graph(id, g, _topology)
    return mol


class GlycanTree:
    """
    Stores the meta-glycan structure as a tree structure for easier access when plotting
    glycan schematics.
    """

    def __init__(self):
        # REMARK
        # Originally, the whole thing was supposed to be a networkx graph, or a dictionary,
        # but strangely, during removal residues could not be found in the graph using `i in graph`
        # even though iterations over the graph clearly showed them being present `i == some_residue for i in graph.nodes`
        # for instance showed some True values. Since the issues persisted with various other implementations, this one finally seemed to solve
        # the problem, so we stick with it even though it is not particularly elegant.
        self._segments = []
        self._linkages = []

    def add(self, residue1, residue2, linkage):
        self._segments.append((residue1, residue2))
        self._linkages.append(linkage)

    def remove(self, residue1, residue2):
        idx = self._segments.index((residue1, residue2))
        self._segments.pop(idx)
        self._linkages.pop(idx)

    def get_linkage(self, residue1, residue2):
        idx = self._segments.index((residue1, residue2))
        return self._linkages[idx]

    def update(self, other):
        self._segments.extend(other._segments)
        self._linkages.extend(other._linkages)

    def __iter__(self):
        return iter(self._segments)

    def __len__(self):
        return len(self._segments)


class Glycan(core.Molecule):
    """
    A glycan molecule

    Parameters
    ----------
    id : str
        The id of the molecule
    structure: biobuild.Structure or biopython.PDB.Structure or biobuild.Molecule
        The structure of the molecule
    """

    def __init__(
        self,
        structure,
        root_atom: Union[str, int, core.base_classes.Atom] = None,
        model: int = 0,
        chain: str = None,
    ):
        if isinstance(structure, core.Molecule):
            model = model or structure.model.id
            chain = chain or structure._chain.id
            super().__init__(
                structure.structure,
                structure.root_atom,
                model,
                chain,
            )
            self.add_bonds(*structure._bonds)
        else:
            super().__init__(structure, root_atom, model, chain)
        self._glycan_tree = GlycanTree()

    @classmethod
    def from_compound(cls, compound: str):
        new = core.Molecule.from_compound(compound)
        # this is the only method that does not use cls() to create a new instance in the biobuild implementation...
        return cls(new)

    def attach(
        self,
        other: "core.Molecule",
        link: Union[str, "core.Linkage"] = None,
        at_residue: Union[int, "core.base_classes.Residue"] = None,
        other_residue: Union[int, "core.base_classes.Residue"] = None,
        inplace: bool = True,
        other_inplace: bool = False,
        _topology=None,
    ) -> "Glycan":
        out = core.Molecule.attach(
            self,
            other,
            link,
            at_residue,
            other_residue,
            inplace,
            other_inplace,
            _topology,
        )

        # now we need to add a new strip to the glycan tree
        # since we cannot know beforehand which residues were actually connected
        # we need to retrace some steps that are taken inside `attach` because the patcher or stitcher
        # is the only object that actually remembers the residues that it connected
        if isinstance(link, str):
            _topology = _topology or resources.get_default_topology()
            link = _topology.get_patch(link)
        factory = (
            structural.__default_keep_keep_patcher__
            if link.has_IC
            else structural.__default_keep_keep_stitcher__
        )
        at_residue = factory._target_residue
        other_residue = factory._source_residue

        # add a new strip to the glycan tree
        out._glycan_tree.add(at_residue, other_residue, link)

        # if we have additional glycan strips from the attached structure, add these as well...
        if isinstance(other, Glycan):
            out._glycan_tree.update(factory.source._glycan_tree)

        return out

    def remove_residues(
        self, *residues: Union[int, "core.base_classes.Residue"]
    ) -> list:
        residues = core.Molecule.remove_residues(self, *residues)
        _strips_to_remove = set()
        for strip in self._glycan_tree:
            for residue in residues:
                if strip[0] == residue or strip[1] == residue:
                    _strips_to_remove.add(strip)
        for s in _strips_to_remove:
            self._glycan_tree.remove(*s)
        return residues

    def __repr__(self):
        return f"Glycan({self.id})"


def _parse_iupac_graph(id, glycan_segments, _topology=None):
    """
    Make a molecule from a list of glycan segments that were generated by the IUPACParser class

    Parameters
    ----------
    id : str
        The id of the molecule
    glycan_segments : list
        A list of glycan segments

    Returns
    -------
    Molecule
        The molecule
    """
    if not _topology:
        _topology = bb.get_default_topology()

    # Check that all segments have a known patch
    for segment in glycan_segments:
        link = segment[-1]

        if not _topology.has_patch(link):
            # ---------------------------- TO REMOVE LATER ----------------------------
            # STILL SOME DEBUG STUFF HERE
            # ---------------------------- TO REMOVE LATER ----------------------------
            f = core.molecule(segment[0])
            t = core.molecule(segment[1])
            f.show()
            t.show()
            print("MISSING: ", link)
            if input("Next: (y)").lower() == "y":
                continue
            raise ValueError(
                f"No patch/recipe available for linkage: {link}. Try adding a patch or recipe to the topology."
            )

    mol = None
    first_mol = None
    second_mol = None
    at_residue = None
    other_residue = None
    residue_id_mapping = {}
    for i, segment in enumerate(glycan_segments):
        first, second, link = segment

        first_name = first.split("@")[0]
        second_name = second.split("@")[0]

        if first in residue_id_mapping:
            at_residue = residue_id_mapping[first]
            first_mol = mol
        else:
            first_mol = core.molecule(first_name)
            if isinstance(first_mol, list):
                first_mol = first_mol[0]
            # if we did not get the compound from the PDBE compounds,
            # we probably got them from PubChem, in which case we need to autolabel them
            if not bb.has_compound(first_mol.id):
                first_mol.autolabel()
            residue_id_mapping[first] = len(residue_id_mapping) + 1
            at_residue = None

        if second in residue_id_mapping:
            other_residue = residue_id_mapping[second]
            second_mol = mol
        else:
            second_mol = core.molecule(second_name)
            if isinstance(second_mol, list):
                second_mol = second_mol[0]
            if not bb.has_compound(second_mol.id):
                second_mol.autolabel()
            residue_id_mapping[second] = len(residue_id_mapping) + 1
            other_residue = None

        if not mol:
            mol = first_mol
        first_mol.attach(
            second_mol,
            link,
            at_residue=at_residue,
            other_residue=other_residue,
            _topology=_topology,
        )

    mol.id = id
    return mol


__all__ = ["Glycan", "read_snfg", "read_graph", "glycan"]

if __name__ == "__main__":
    import biobuild as bb

    bb.load_sugars()
    glc = Glycan.from_compound("GLC")
    read_snfg = bb.connect(glc, glc, "14bb")
    read_snfg = bb.connect(read_snfg, glc, "16ab")
    read_snfg = bb.connect(read_snfg, glc, "12ab", at_residue_a=-2)
    out = bb.connect(read_snfg, read_snfg, "13ab")
    # out.show()
    _removed = out.remove_residues(4)
    print(out._glycan_tree)
    pass
