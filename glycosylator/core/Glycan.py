"""
Glycans are represented by the `Glycan` class. The `Glycan` is actually just a derivative of a `buildamol.Molecule` with some additional methods and attributes to handle glycan-specific operations. 
The `Glycan` class is used to represent glycan molecules and is used to generate glycan molecules from IUPAC strings, graph structures, or to extend/crop glycan molecules to match a given IUPAC string. 
The `Glycan` class also has methods to draw 2D (SNFG) and 3D representations of the glycan molecule.

Making Glycans from IUPAC strings
---------------------------------
The IUPAC nomenclature has been widely established to describe glycans in textual form because their SMILES tend to be very long and cumbersome to work with.
`Glycans` have a `from_iupac` method that can read an IUPAC string and produce a 3D model for the glycan. Alternatively one may use the top-level `read_iupac` function
to the same effect, or pass an IUPAC string to the `glycan` function. 

.. code-block:: python

    import glycosylator as gl 

    my_glycan = gl.read_iupac("my new glycan", "Man(b1-4)Glc")

The `Glycan` class also has a `to_iupac` method that can convert a glycan molecule to an IUPAC string.

.. code-block:: python

    iupac_string = my_glycan.to_iupac()

Making Glycans from other inputs
--------------------------------
Since the `Glycan` class inherits from the `buildamol.Molecule` it supports many more inputs such as RDKit molecules, SMILES, or PDB files.
You can either use a dedicated classmethod such as `from_pdb` to get a desired glycan molecule, or trust that the top-level `glycan` function
can automatically figure out what kind of input you are providing. This function is very versatile and can be provided with a variety of inputs
which are automatically processed behind the scenes. It is the most convenient way for users to obtain a glycan structure.

.. code-block:: python

    # make a non-standard glucose from SMILES
    smiles = "NCC1OC(O)C(O)C(O)C1O"
    nitrogen_glucose = gl.glycan(smiles)

Modifying individual sugars
---------------------------
If you do not feel like working with SMILES in order to make small chemical changes on individual sugars or even whole glycans,
you can use the flexibility of `BuildAMol`, which Glycosylator is built upon, to help you out. Let's say we we want to make a phospho-glucose
we can do something like this:

.. code-block:: python

    import glycosylator as gl
    import buildamol as bam

    # get a glucose 
    glc = gl.glycan("GLC")

    # now use buildamol directly to modify
    # the molecule
    bam.phosphorylate(glc, at_atom="C6")
    glc.remove_atoms("O6", "HO6")

"""

from typing import Union
import warnings

import glycosylator.utils as utils
import buildamol.structural as structural
import glycosylator.resources as resources
import buildamol.core as core

from buildamol.extensions.bio.glycans.glycan import (
    _segments_to_mol as _parse_iupac_graph,
)

import pandas as pd


def glycan(g: Union[str, list], id: str = None, _topology=None):
    """
    The toplevel function to generate an entire glycan molecule from either an IUPAC/SNFG string, a list of residues from a graph structure, or just get a single residue glycan molecule (e.g. one Glucose, etc.).

    Parameters
    ----------
    g : str or list
        The glycan string to parse.
        The string may be a single sugar residue's name, in which case this function is used like `buildamol.molecule`, or can be an entire glycan structure in IUPAC/SNFG condensed format - currently, neither extended nor short formats are supported (refer to the `read_iupac` function for more information).
        Alternatively, a list of residues from a graph structure can be passed (refer to the `read_graph` function for more information).
    id : str
        The id of the molecule to create.
        If not provided, the id will be the same as the input string.
    _topology
        A particular topology to use. If None, the default topology is used.

    Returns
    -------
    Glycan
        The created Glycan molecule.
    """
    if isinstance(g, str):
        try:
            return read_iupac(id, g, _topology)
        except:
            try:
                mol = core.molecule(g)
                if isinstance(mol, list):
                    mol = mol[0]
                mol.id = id if id is not None else g
                mol = Glycan(mol)
                return mol
            except:
                try:
                    mol = Glycan.from_glycosmos(g, _topology)
                    return mol
                except:
                    raise ValueError(
                        f"Failed to interpret input '{g}'. Could not parse as IUPAC/SNFG string and not get a molecule from it. If it is a IUPAC/SNFG string, perhaps one of the residues could not be found in the database or a linkage is not available in the used topology? Try building the glycan directly."
                    )

    elif isinstance(g, list):
        mol = read_graph(id, g)
        return mol


def read_iupac(id: str, s: str, _topology=None) -> "Glycan":
    """
    Make a molecule from an IUPAC glycan string in condensed format.

    Parameters
    ----------
    id : str
        The id of the molecule to create
    s : str
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
    >>> mol = read_iupac("my_glycan", iupac)
    """
    if isinstance(s, str):
        s = utils.__default_IUPACParser__(s)
    else:
        raise TypeError(
            f"Expected string, got {type(s)} instead. Perhaps you meant to use the `read_graph` function?"
        )
    mol = Glycan(_parse_iupac_graph(id, s, _topology))
    return mol


# Alias
read_snfg = read_iupac


def write_iupac(mol: "Glycan") -> str:
    """
    Write a molecule as an IUPAC string in condensed format.

    Parameters
    ----------
    mol : Glycan
        The molecule to write.

    Returns
    -------
    iupac : str
        The IUPAC/SNFG string in condensed format.
    """
    return mol.to_iupac()


write_snfg = write_iupac


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
        self._connectivity = {}

    @property
    def segments(self):
        return self._segments

    @property
    def linkages(self):
        return self._linkages

    def add(self, residue1, residue2, linkage):
        self._segments.append((residue1, residue2))
        self._linkages.append(linkage)
        if residue1 not in self._connectivity:
            self._connectivity[residue1] = []
        self._connectivity[residue1].append(residue2)

    def remove(self, residue1, residue2):
        idx = self._segments.index((residue1, residue2))
        self._segments.pop(idx)
        self._linkages.pop(idx)

    def get_linkage(self, residue1, residue2):
        idx = self._segments.index((residue1, residue2))
        return self._linkages[idx]

    def has_matching_segment(self, residue1, residue2):
        return (residue1, residue2) in self._segments or any(
            residue1.matches(j) and residue2.matches(l) for j, l in self._segments
        )

    def get_linkage_with_matching(self, residue1, residue2):
        for i, (j, l) in enumerate(self._segments):
            if residue1.matches(j) and residue2.matches(l):
                return self._linkages[i]
        return None

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
    structure: buildamol.Structure or biopython.PDB.Structure or buildamol.Molecule
        The structure of the molecule
    """

    def __init__(
        self,
        structure,
        root_atom: Union[str, int, core.Atom] = None,
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
        self._scaffold = None

    @classmethod
    def from_pdb(cls, filename: str):
        new = core.Molecule.from_pdb(filename)
        new = cls(new)
        try:
            new.infer_glycan_tree()
        except:
            pass
        return new

    @classmethod
    def from_compound(cls, compound: str, by: str = None, root_atom=None):
        new = core.Molecule.from_compound(compound, by, root_atom)
        # this is the only method that does not use cls() to create a new instance in the buildamol implementation...
        return cls(new)

    @classmethod
    def from_iupac(cls, id: str, iupac: str, _topology=None) -> "Glycan":
        """
        Generate a glycan molecule from an IUPAC/SNFG string
        """
        new = read_iupac(id, iupac, _topology)
        try:
            new.infer_glycan_tree()
        except:
            pass
        return new

    @classmethod
    def from_glycosmos(cls, id: str, _topology=None) -> "Glycan":
        """
        Generate a glycan molecule from a GlyCosmos/GlyTouCan ID
        """
        iupac = resources.get_iupac_from_glycosmos(id)
        return cls.from_iupac(id, iupac, _topology)

    from_glytoucan = from_glycosmos
    from_snfg = from_iupac

    def get_glytoucan_id(self) -> str:
        """
        Get the GlyTouCan ID of the glycan molecule

        Returns
        -------
        str
            The GlyTouCan ID (if available)
        """
        return resources.get_glytoucan_id_from_iupac(self.to_iupac())

    get_glycosmos_id = get_glytoucan_id

    def search_glytoucan_ids(self) -> list:
        """
        Find GlyTouCan IDs for glycans that are partial matches of the glycan.

        Returns
        -------
        list
            A list of GlyTouCan IDs (if available)
        """
        return resources.find_glytoucan_ids_from_iupac(self.to_iupac())

    search_glycosmos_ids = search_glytoucan_ids

    def find_glytoucan_ids(self) -> list:
        warnings.warn(
            "This method is deprecated. Use `search_glytoucan_ids` instead.",
            DeprecationWarning,
        )
        return self.search_glytoucan_ids()

    find_glycosmos_ids = find_glytoucan_ids

    def to_iupac(self, add_terminal_conformation: bool = True) -> str:
        """
        Generate an IUPAC/SNFG string from the glycan molecule

        Parameters
        ----------
        add_terminal_conformation : bool
            Whether to add the terminal conformation of the first residue as `(a1-` or `(b1-` to the end string.
        """
        return utils.make_iupac_string(self, add_terminal_conformation)

    to_snfg = to_iupac

    def infer_glycan_tree(self):
        """
        Infer the glycan tree connectivity in case of a glycan molecule that was loaded externally
        """
        if not self.root_atom:
            self.root_atom = 1

        # now be sure to fill the glycan tree
        connections = self.get_residue_connections(triplet=False)
        connections = self._AtomGraph.direct_edges(self.root_atom, edges=connections)
        for bond in connections:
            a, b = bond
            res_a, res_b = a.parent, b.parent

            id = utils.iupac.make_link_id(bond)
            _link = core.Linkage(id=id, description="<auto-generated: labels only>")
            _link.atom1 = a
            _link.atom2 = b
            self._glycan_tree.add(res_a, res_b, _link)

    def hist(self) -> pd.DataFrame:
        """
        Get a histogram of the glycan residues

        Returns
        -------
        pd.DataFrame
            The histogram
        """
        hist = {}
        for res in self.get_residues():
            name = resources.id_to_name(res.resname)
            if name is None:
                name = res.resname
            hist[name] = hist.get(name, 0) + 1

        return pd.DataFrame(
            {"residue": [i for i in hist.keys()], "count": [i for i in hist.values()]}
        )

    def attach(
        self,
        other: "core.Molecule",
        link: Union[str, "core.Linkage"] = None,
        at_residue: Union[int, "core.Residue"] = None,
        other_residue: Union[int, "core.Residue"] = None,
        use_patch: bool = True,
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
            use_patch,
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
            if link.has_IC and use_patch
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

    def extend_to(self, full_iupac: str, inplace: bool = True, _topology=None):
        """
        Extend a partial glycan to a glycan matching the provided IUPAC/SNFG string.
        Note that this can only extend the glycan from the root-onward. It cannot retro-fit
        a root given some leaf glycans.

        Parameters
        ----------
        full_iupac : str
            The full IUPAC/SNFG string which the glycan should be extended to.
        inplace : bool
            Whether to extend the glycan in place or return a new glycan.
        _topology
            A particular topology to use. If None, the default topology is used.
        """
        full = Glycan.from_iupac(None, full_iupac, _topology)
        if self.root_residue is not None:
            full.rename_chain("A", self.root_residue.parent.id)

        obj = self if inplace else self.copy()

        present_segments = [
            (i.resname, j.resname) for i, j in obj._glycan_tree.segments
        ]
        for segment in full._glycan_tree:

            # if the segment is already present in the current glycan, skip it
            _s = (segment[0].resname, segment[1].resname)
            if _s in present_segments:
                present_segments.remove(_s)
                continue

            # if obj._glycan_tree.has_matching_segment(*segment):
            #     continue

            # get the linkage to attach with
            link = full._glycan_tree.get_linkage(*segment)

            # get the residue to attach to from the loaded compounds
            # they should be available since all standard IUPAC-represented glycans should be there
            incoming = Glycan.from_compound(segment[1].resname, by="id")

            # get the residue in the current glycan to attach to
            attach_residue = None
            candidates = (i for i in obj.get_residues() if i.matches(segment[0]))
            for candidate in candidates:
                if link.can_apply(obj, incoming, target_residue=candidate):
                    attach_residue = candidate
                    break

            if attach_residue is None:
                raise ValueError(
                    f"Could not find a residue matching {segment[0]} to attach to in the current glycan"
                )

            # attach the residue
            obj.attach(
                incoming,
                link,
                at_residue=attach_residue,
                other_inplace=True,
                _topology=_topology,
            )

        return obj

    def crop_to(self, full_iupac: str, inplace: bool = True, _topology=None):
        """
        Crop a glycan to a glycan matching the provided IUPAC/SNFG string.
        Note that this can only crop the glycan from the leaves-onward. It cannot retro-fit
        a root given some leaf glycans.

        Parameters
        ----------
        full_iupac : str
            The full IUPAC/SNFG string which the glycan should be cropped to.
        inplace : bool
            Whether to crop the glycan in place or return a new glycan.
        _topology
            A particular topology to use. If None, the default topology is used.
        """
        full = Glycan.from_iupac(None, full_iupac, _topology)

        obj = self if inplace else self.copy()

        present_segments = [
            (i.resname, j.resname) for i, j in obj._glycan_tree.segments
        ]
        for segment in present_segments:
            if not full._glycan_tree.has_matching_segment(*segment):
                obj.remove_residues(
                    *[i for i in obj.get_residues() if i.matches(segment[0])]
                )

        return obj

    def remove_residues(self, *residues: Union[int, "core.Residue"]) -> list:
        residues = core.Molecule.remove_residues(self, *residues)
        _strips_to_remove = set()
        for strip in self._glycan_tree:
            for residue in residues:
                if strip[0] == residue or strip[1] == residue:
                    _strips_to_remove.add(strip)
        for s in _strips_to_remove:
            self._glycan_tree.remove(*s)
        return residues

    def clashes_with_scaffold(
        self,
        clash_threshold: float = 1.0,
        ignore_hydrogens: bool = True,
        coarse_precheck: bool = True,
    ) -> bool:
        """
        Check if the glycan clashes with the scaffold

        Parameters
        ----------
        clash_threshold : float
            The minimum distance to consider a clash
        ignore_hydrogens : bool
            Whether to ignore hydrogens
        coarse_precheck : bool
            Whether to use a coarse pre-check to speed up the process. This may lead to false negatives,
            especially if the scaffold has very large residues (e.g. lipids with long carbon chains).

        Returns
        -------
        bool
            Whether the glycan clashes with the scaffold
        """
        if self._scaffold is None:
            return False
        return (
            len(
                self.find_clashes_with(
                    self._scaffold,
                    clash_threshold=clash_threshold,
                    ignore_hydrogens=ignore_hydrogens,
                    coarse_precheck=coarse_precheck,
                )
            )
            > 0
        )

    def snfg(
        self,
        ax=None,
        axis="y",
        **kwargs,
    ) -> "utils.visual.plt.Axes":
        """
        Draw the SNFG 2D schematic of the glycan

        Parameters
        ----------
        ax : matplotlib.Axes
            The axes to draw on
        axis : str
            The orientation of the glycan y-axis = vertical, x-axis = horizontal.

        Returns
        -------
        matplotlib.Axes
            The axes
        """
        v = utils.visual.GlycanViewer2D(self)
        return v.draw(
            ax=ax,
            axis=axis,
            **kwargs,
        )

    draw2d = snfg

    # def show2d(
    #     self,
    #     axis="y",
    #     **kwargs,
    # ):
    #     """
    #     Draw and show the 2D schematic of the glycan

    #     Parameters
    #     ----------
    #     axis : str
    #         The orientation of the glycan y-axis = vertical, x-axis = horizontal.
    #     node_size : int
    #         The size of the nodes
    #     """
    #     ax = self.draw2d(
    #         axis=axis,
    #         **kwargs,
    #     )
    #     if kwargs.get("svg", False):
    #         return ax
    #     utils.visual.plt.show()

    def __repr__(self):
        return f"Glycan({self.id})"


__all__ = [
    "Glycan",
    "read_iupac",
    "read_snfg",
    "write_snfg",
    "write_iupac",
    "read_graph",
    "glycan",
]

if __name__ == "__main__":
    # g = glycan("G78791QP")
    # g.show2d()

    #     # s = "Gal(b1-3)GlcNAc(b1-3)[Gal(b1-3)GlcNAc(b1-3)[Gal(b1-4)GlcNAc(b1-6)]Gal(b1-4)GlcNAc(b1-6)][Gal(b1-4)GlcNAc(b1-2)]Gal(b1-4)Glc"
    #     # glc = Glycan.from_iupac(None, s)
    #     # glc2 = glc % "14bb" + glc

    import glycosylator as gl

    g = Glycan.from_iupac(None, "Gal(b1-4)GlcNAc")
    g.extend_to("Gal(b1-4)[Gal(b1-4)Glc(b1-6)]GlcNAc")

    g.root = g.get_atom(1)
    graph = g.make_atom_graph()
    edges = g.get_residue_connections(triplet=False, direct_by="root")
    env = gl.optimizers.DistanceRotatron(graph, edges)
    out = gl.optimizers.optimize(g, env)

    g.show()
    # g.show_snfg()

#     print(g.hist())
#     print(g)

# f = "/Users/noahhk/GIT/glycosylator/examples/my_first_glycan.pdb"
# g = Glycan.from_pdb(f)
# g.infer_glycan_tree()
# g.show2d()

# _read_snfg = gl.connect(glc, glc, "14bb")
# _read_snfg = gl.connect(_read_snfg, glc, "16ab")
# _read_snfg = gl.connect(_read_snfg, glc, "12ab", at_residue_a=-2)
# out = gl.connect(_read_snfg, _read_snfg, "13ab")
# out.show()
# _removed = out.remove_residues(4)
# print(out._glycan_tree)
# out.show2d()
# pass
