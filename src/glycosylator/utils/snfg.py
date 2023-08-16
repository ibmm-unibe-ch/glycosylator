"""
Functions to work with the IUPAC glycan nomenclature.
"""
import re
from typing import Any
import networkx as nx
import glycosylator.resources.names as names


def reformat_link(string):
    """
    Reformat a linkage string from the IUPAC format to the Glycosylator (CHARMM force field) format.

    Parameters
    ----------
    string : str
        The linkage string in IUPAC format.

    Returns
    -------
    str
        The linkage string in Glycosylator format.
    """

    type, a, b = re.match("([ab])(\d+)-(\d+)", string).groups()
    return a + b + type * 2


def reverse_format_link(string, pretty: bool = False):
    """
    Reverse format a linkage string from the glycosylator (CHARMM force field) format to IUPAC format.
    If the string cannot be formatted, it is returned as is.

    Parameters
    ----------
    string : str
        The linkage string in Glycosylator format.
    pretty : bool
        If True, returns a prettier formatted version with `α` and `β` instead of `a` and `b`, and an actual arrow instead of a dash.

    Returns
    -------
    str
        The linkage string in IUPAC format if possible, otherwise the original string.
    """
    match = re.match("(\d+)(\d+)([ab]{2})", string)
    if match is None:
        return string
    a, b, type = match.groups()
    if pretty:
        pre = "α" if type[1] == "a" else "β"
        return pre + "(" + a + "→" + b + ")"
    return type[1] + a + "-" + b


class SNFGParser:
    """
    A parser for condensed IUPAC glycan nomenclature strings. This class will generate
    a list of connecting glycan segments from a string from which a Molecule can be built.
    """

    def __init__(self):
        self._string = ""
        self.reset()

    @property
    def _current(self):
        return self._string[-self._idx]

    @property
    def _next(self):
        return self._string[-self._idx - 1]

    @property
    def _can_shift(self):
        return self._idx < len(self._string)

    @property
    def _is_at_bracket(self):
        return self._current in ("(", ")")

    @property
    def _is_at_open_square_bracket(self):
        return self._current == "]"

    @property
    def _is_at_close_square_bracket(self):
        return self._current == "["

    @property
    def _can_store(self):
        return (
            not self._is_still_parsing
            and len(self._latest_residue) >= 1
            and len(self._latest_linkage) >= 3
            and len(self._second_latest_residue) >= 1
        )

    @property
    def _is_still_parsing(self):
        if not self._can_shift:
            return False
        return self._current not in ("(", ")", "[", "]")

    def parse(self, string):
        """
        Parse a string of IUPAC glycan nomenclature into a list of glycan segments.

        Parameters
        ----------
        string : str
            The IUPAC glycan nomenclature string.

        Returns
        -------
        list
            A list of tuples where each segment is a tuple of (residue1, residue2, linkage).
        """
        self._string = self._prep_greek_letters(string)
        self.reset()
        self._parse()
        return self._glycan

    def reset(self):
        """
        Reset the parser.
        """
        self._glycan = []
        self._idx = 1
        self._residue_counts = {}
        self._latest_residue = ""
        self._latest_conformation = ""
        self._latest_linkage = ""
        self._second_latest_residue = ""
        self._second_latest_conformation = ""
        self._latest_residue_before_square_bracket = ""
        self._latest_conformation_before_square_bracket = ""

    def _shift(self):
        self._idx += 1

    def _shift_residue(self):
        self._second_latest_residue = self._latest_residue
        self._latest_residue = ""
        self._second_latest_conformation = self._latest_conformation
        self._latest_conformation = ""

    def _parse(self):
        self._crop_end()
        while self._can_shift:
            if self._can_store:
                self._store()
                continue
            if self._is_at_bracket:
                self._shift_residue()
                self._latest_linkage = self._parse_linkage()
                self._shift()
                continue
            if self._is_at_open_square_bracket:
                self._latest_residue_before_square_bracket = self._fit_residue(
                    self._latest_residue, increase_count=False
                )
                self._latest_residue = self._latest_residue_before_square_bracket
                self._latest_conformation_before_square_bracket = (
                    self._latest_conformation
                )
                self._shift()
                continue
            if self._is_at_close_square_bracket:
                self._latest_residue = self._latest_residue_before_square_bracket
                self._latest_conformation = (
                    self._latest_conformation_before_square_bracket
                )
                self._second_latest_residue = ""
                self._second_latest_conformation = ""
                self._shift()
                continue
            self._latest_residue += self._current
            self._shift()
        self._latest_residue += self._current
        self._store()

    def _store(self):
        second = self._second_latest_residue
        if "@" not in second:
            second = self._fit_residue(second)
            self._second_latest_residue = second

        latest = self._latest_residue
        if not "@" in latest:
            latest = self._fit_residue(latest)
            self._latest_residue = latest

        branch = (second, latest, self._reformat_link(self._latest_linkage))
        self._glycan.append(branch)
        self._latest_linkage = ""

    def _fit_residue(self, r, increase_count=True):
        if "@" in r:
            return r

        r = r[::-1]
        if r in self._residue_counts:
            if increase_count:
                self._residue_counts[r] += 1
            r += "@" + str(self._residue_counts[r])
        else:
            self._residue_counts[r] = 1
            r += "@1"
        return r

    def _parse_linkage(self):
        self._shift()
        linkage = ""
        while self._can_shift and not self._is_at_bracket:
            linkage += self._current
            self._shift()
        self._latest_conformation = linkage[-1]
        linkage = linkage[:-1]
        return linkage

    def _crop_end(self):
        if self._string[-1] == "-":
            while self._next != "(":
                self._shift()
            self._latest_conformation = self._current
            self._shift()
            self._string = self._string[: -self._idx]
            self._idx = 1
        else:
            self._latest_conformation = (
                "a"  # assume alpha conformation if not specified
            )

    def _reformat_link(self, link):
        link = link[::-1].replace("-", "")
        link = link + self._second_latest_conformation + self._latest_conformation
        return link

    def _prep_greek_letters(self, string):
        string = string.replace("α", "a").replace("β", "b")
        return string

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.parse(*args, **kwds)


def parse_snfg(string):
    """
    Parse a string of IUPAC/SNFG glycan nomenclature into a list of glycan segments.

    Parameters
    ----------
    string : str
        The IUPAC glycan nomenclature string.

    Returns
    -------
    list
        A list of tuples where each segment is a tuple of (residue1, residue2, linkage).
    """
    p = SNFGParser()
    return p.parse(string)


parse_iupac = parse_snfg


class SNFGStringMaker:
    """
    This class generates IUPAC/SNFG strings from a list of glycan molecule.
    """

    def write_string(self, glycan: "Glycan") -> str:
        """
        Write the IUPAC string for a glycan molecule.

        Parameters
        ----------
        glycan : Glycan
            The glycan molecule to write the IUPAC string for.

        Returns
        -------
        str
            The IUPAC string for the glycan molecule.
        """
        # setup the class
        self._setup(glycan)

        # generate the string (in reverse order)
        self._string = self._write_string(self.root_residue)

        # fill in all branch-place holders
        for key, value in self._sub_branch_mapping.items():
            self._string = self._string.replace(key, "]" + value + "[")

        return self._string[::-1]

    def _setup(self, glycan):
        """Setup the class for writing the IUPAC string."""
        self.root = glycan.get_root() or glycan.get_atom(1)
        self.root_residue = self.root.get_parent()
        self.graph = nx.Graph(
            glycan._glycan_tree._segments,
        )
        nx.set_edge_attributes(
            self.graph,
            {
                i: j
                for i, j in zip(
                    glycan._glycan_tree._segments, glycan._glycan_tree._linkages
                )
            },
            "linkage",
        )
        self._string = ""
        self._child_mapping = self._make_child_mapping()
        self._sub_branch_mapping = {}
        self.idx = 0

    def _make_child_mapping(self):
        """Make a mapping of each residue to its children."""
        graph = nx.dfs_tree(self.graph, self.root_residue)
        child_mapping = {residue: [] for residue in graph.nodes}
        for node in child_mapping.keys():
            children = list(graph[node].keys())
            child_mapping[node] = children
        return child_mapping

    def _write_string(self, parent, residue=None):
        """Write the IUPAC string for a residue."""
        # start of the string at the root residue
        if parent and not residue:
            string = self._get_name(parent)
            children = self._child_mapping[parent]
            if len(children) > 1:
                cdx = 0
                for child in children:
                    if cdx < len(children) - 1:
                        string += f"<{self.idx}>"
                        self._sub_branch_mapping[f"<{self.idx}>"] = self._write_string(
                            parent, child
                        )
                        self.idx += 1
                        cdx += 1
                    else:
                        string += self._write_string(parent, child)

            return string

        # recursive part for residues that are children of other residues
        name_residue = self._get_name(residue)
        linkage = self._get_linkage(self.graph[parent][residue]["linkage"])
        string = linkage + name_residue

        children = self._child_mapping[residue]
        if len(children) != 0:
            for child in children:
                # if we have a branch, we just set a placeholder and later fill them all in at once...
                if self._should_open_branch(residue, child):
                    string += f"<{self.idx}>"
                    self._sub_branch_mapping[f"<{self.idx}>"] = self._write_string(
                        residue, child
                    )
                    self.idx += 1
                else:
                    string += self._write_string(residue, child)
        return string

    def _should_open_branch(self, parent, residue=None) -> bool:
        """Check if we should open a branch."""
        child_mapping = self._child_mapping[parent]
        verdict = len(child_mapping) > 1
        if residue:
            verdict = verdict and child_mapping.index(residue) != len(child_mapping) - 1
        return verdict

    @staticmethod
    def _get_name(residue):
        return names.id_to_name(residue.resname)[::-1]

    @staticmethod
    def _get_linkage(linkage):
        return f"({reverse_format_link(linkage.id)})"[::-1]


def make_snfg_string(glycan):
    """
    Make an IUPAC/SNFG string from a glycan molecule.

    Parameters
    ----------
    glycan : Glycan
        The glycan molecule to make the IUPAC string for.

    Returns
    -------
    str
        The IUPAC string for the glycan molecule.
    """
    if len(glycan._glycan_tree._segments) == 0:
        raise ValueError(
            "This glycan does not have any entries in its glycan tree. Make sure it was generated using glycosylator! Currently, only glycans generated using glycosylator are supported."
        )
    m = SNFGStringMaker()
    return m.write_string(glycan)


make_iupac_string = make_snfg_string

__all__ = ["SNFGParser", "SNFGStringMaker", "parse_snfg", "parse_iupac", "make_iupac_string", "make_snfg_string"]

if __name__ == "__main__":
    # string2 = "Man(b1-6)[Man(b1-3)]BMA(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    # string = "F(b1-4)[E(a2-3)D(a1-4)]C(a1-6)B(b1-4)A"
    # string3 = "Glc(b1-3)[Gal(a1-4)Man(b1-6)]GlcNAc(a1-4)Man(b1-4)Glc(b1-"

    # p = SNFGParser()
    # g = p.parse(string3)
    # print(g)

    import glycosylator as gl

    s = "GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)[GalNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"

    mol = gl.glycan(s)
    out = SNFGStringMaker().write_string(mol)
    mol2 = gl.glycan(out)
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2)
    mol.draw2d(ax=axs[0])
    mol2.draw2d(ax=axs[1])
    axs[0].set_title("Original")
    axs[1].set_title("Reconstructed")
    plt.show()
    print(out)
    pass
