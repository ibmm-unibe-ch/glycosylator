"""
Functions to work with the IUPAC glycan nomenclature.
"""

import re
from typing import Any
import networkx as nx
import glycosylator.resources.names as names


def make_link_id(bond: "Bond") -> str:
    """
    Make a link id string in Glycosylator format
    from a Bond object.

    Parameters
    ----------
    bond : Bond
        The bond from which to make an id string.

    Returns
    -------
    str
        The linkage id string in Glycosylator format.
    """
    atom1, atom2 = bond
    if atom1.element != "O" or atom2.element != "C":
        raise ValueError(
            f"This linkage does not seem to be a glycosidic bond and cannot be converted to a linkage id string.\nExpected a bond between an oxygen atom and a carbon atom, but got a bond between a {atom1.element} atom and a {atom2.element} atom."
        )
    numeral_atom1 = atom1.id[1:]
    numeral_atom2 = atom2.id[1:]
    conf_res1 = "b" if names.is_beta(atom1.parent.resname) else "a"
    conf_res2 = "b" if names.is_beta(atom2.parent.resname) else "a"

    return numeral_atom2 + numeral_atom1 + conf_res2 + conf_res1


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


def reverse_format_link(string, pretty: bool = False, small: bool = False):
    """
    Reverse format a linkage string from the glycosylator (CHARMM force field) format to IUPAC format.
    If the string cannot be formatted, it is returned as is.

    Parameters
    ----------
    string : str
        The linkage string in Glycosylator format.
    pretty : bool
        If True, returns a prettier formatted version with `α` and `β` instead of `a` and `b`, and an actual arrow instead of a dash.
    small : bool
        If True, returns a smaller version of the string with only the a/b and the number of the second atom.

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
    else:
        pre = type[1]
    if small:
        return pre + b
    if pretty:
        return pre + "(" + a + "→" + b + ")"
    return pre + a + "-" + b


class IUPACParser:
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
        self._glycan = self._adjust_betas(self._glycan)
        return self._glycan

    def reset(self):
        """
        Reset the parser.
        """
        self._glycan = []
        self._idx = 1
        self._residue_counts = {}
        self._past_residue_before_square_bracket = []
        self._past_conformations_before_square_bracket = []
        self._latest_residue = ""
        self._latest_conformation = ""
        self._latest_linkage = ""
        self._second_latest_residue = ""
        self._second_latest_conformation = ""
        # self._latest_residue_before_square_bracket = ""
        # self._latest_conformation_before_square_bracket = ""

    @property
    def _latest_residue_before_square_bracket(self):
        return self._past_residue_before_square_bracket[-1]

    @property
    def _latest_conformation_before_square_bracket(self):
        return self._past_conformations_before_square_bracket[-1]

    def _push_residue_before_square_bracket(self, residue):
        self._past_residue_before_square_bracket.append(residue)

    def _push_conformation_before_square_bracket(self, conformation):
        self._past_conformations_before_square_bracket.append(conformation)

    def _pop_residue_before_square_bracket(self):
        return self._past_residue_before_square_bracket.pop()

    def _pop_conformation_before_square_bracket(self):
        return self._past_conformations_before_square_bracket.pop()

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
                _latest_residue_before_square_bracket = self._fit_residue(
                    self._latest_residue, increase_count=False
                )
                self._latest_residue = _latest_residue_before_square_bracket
                self._push_residue_before_square_bracket(self._latest_residue)
                self._push_conformation_before_square_bracket(self._latest_conformation)
                # self._latest_conformation_before_square_bracket = (
                #     self._latest_conformation
                # )
                self._shift()
                continue
            if self._is_at_close_square_bracket:
                self._latest_residue = self._pop_residue_before_square_bracket()
                self._latest_conformation = (
                    self._pop_conformation_before_square_bracket()
                )
                # self._latest_residue = self._latest_residue_before_square_bracket
                # self._latest_conformation = (
                #     self._latest_conformation_before_square_bracket
                # )
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

    def _adjust_betas(self, segments):
        for i in range(len(segments)):
            a, b, link = segments[i]
            if link[2] == "b" and not a.startswith("b-"):
                a = "b-" + a
            if link[3] == "b" and not b.startswith("b-"):
                b = "b-" + b
            segments[i] = (a, b, link)
        return segments

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.parse(*args, **kwds)


class IUPACStringMaker:
    """
    This class generates IUPAC/SNFG strings from a list of glycan molecule.
    """

    double_beta_pattern = r"(b-(\w+\(b))"
    double_beta_replacement = r"\2"
    terminal_beta_pattern = r"b-(\w+)$"
    terminal_beta_replacement = r"\1(b1-"
    terminal_alpha_pattern = r"([\]\)]\w+)$"
    terminal_alpha_replacement = r"\1(a1-"

    def write_string(
        self, glycan: "Glycan", add_terminal_conformation: bool = True
    ) -> str:
        """
        Write the IUPAC string for a glycan molecule.

        Parameters
        ----------
        glycan : Glycan
            The glycan molecule to write the IUPAC string for.
        add_terminal_conformation : bool
            If True, adds the terminal conformation as `last(a1-` or `last(b1-` to the string.

        Returns
        -------
        str
            The IUPAC string for the glycan molecule.
        """
        # setup the class
        self._setup(glycan)

        # generate the string (in reverse order)
        self._string = self._init_string(self.root_residue)

        # for child in self._child_mapping[self.root_residue]:
        #     self._string += self._continue_string(self.root_residue, child)

        # # fill in all branch-place holders
        for key, value in self._sub_branch_mapping.items():
            for key2, value2 in self._sub_branch_mapping.items():
                if key in value2:
                    # self._string = self._string.replace(key, key2)
                    self._sub_branch_mapping[key2] = value2.replace(
                        key, "]" + value + "["
                    )
                    # self._sub_branch_mapping[key] = self._sub_branch_mapping[key2]
        for key, value in self._sub_branch_mapping.items():
            self._string = self._string.replace(key, "]" + value + "[")
        # for key, value in self._sub_branch_mapping.items():
        #     self._string = self._string.replace(key, "]" + value + "[")
        string = self._string[::-1]
        string = re.sub(self.double_beta_pattern, self.double_beta_replacement, string)
        if add_terminal_conformation:
            string = re.sub(
                self.terminal_alpha_pattern, self.terminal_alpha_replacement, string
            )
            string = re.sub(
                self.terminal_beta_pattern, self.terminal_beta_replacement, string
            )
        return string

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

    def _continue_string(self, parent, child):
        """Continue the string with a new residue."""
        string = self._make_one_string(parent, child)
        children = self._child_mapping[child]
        if children:
            parent = child
            for child in children:
                if self._should_open_branch(parent, child):
                    self.idx += 1
                    key = f"<{self.idx}>"
                    string += key
                    self._sub_branch_mapping[key] = self._continue_string(parent, child)
                else:
                    string += self._continue_string(parent, child)
        return string

    def _make_one_string(self, parent, child):
        """Make the string component for one residue."""
        name_residue = self._get_name(child)
        linkage = self._get_linkage(self.graph[parent][child]["linkage"])
        string = linkage + name_residue
        return string

    def _make_child_mapping(self):
        """Make a mapping of each residue to its children."""
        graph = nx.dfs_tree(self.graph, self.root_residue)
        child_mapping = {residue: [] for residue in graph.nodes}
        for node in child_mapping.keys():
            children = list(graph[node].keys())
            child_mapping[node] = children
        return child_mapping

    def _init_string(self, parent):
        """Start writing the string from the root residue."""
        string = self._get_name(parent)
        children = self._child_mapping[parent]
        for child in children:
            if self._should_open_branch(parent, child):
                key = f"<{self.idx}>"
                string += key
                self._sub_branch_mapping[key] = self._continue_string(parent, child)
                self.idx += 1
            else:
                string += self._continue_string(parent, child)
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

    def __call__(self, glycan):
        return self.write_string(glycan)


__default_IUPACParser__ = IUPACParser()
"""
The default IUPACParser instance.
"""

__default_IUPACStringMaker__ = IUPACStringMaker()
"""
The default IUPACStringMaker instance.
"""


def parse_iupac(string):
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
    return __default_IUPACParser__.parse(string)


parse_snfg = parse_iupac


def make_iupac_string(glycan, add_terminal_conformation: bool = True):
    """
    Make an IUPAC/SNFG string from a glycan molecule.

    Parameters
    ----------
    glycan : Glycan
        The glycan molecule to make the IUPAC string for.
    add_terminal_conformation : bool
        If True, adds the terminal conformation as `last(a1-` or `last(b1-` to the end of the string.

    Returns
    -------
    str
        The IUPAC string for the glycan molecule.
    """
    if len(glycan._glycan_tree._segments) == 0 and len(glycan.residues) > 1:
        raise ValueError(
            "This glycan does not have any entries in its glycan tree. Make sure it was generated using glycosylator! Currently, only glycans generated using glycosylator are supported."
        )
    elif len(glycan.residues) == 1:
        _id = glycan.residues[0].resname
        string = names.id_to_name(_id)
        if add_terminal_conformation:
            prefix = "b" if names.is_beta(_id) else "a"
            string += f"({prefix}1-"
        return string
    return __default_IUPACStringMaker__.write_string(glycan, add_terminal_conformation)


make_snfg_string = make_iupac_string

__all__ = [
    "IUPACParser",
    "IUPACStringMaker",
    "parse_snfg",
    "parse_iupac",
    "make_iupac_string",
    "make_snfg_string",
    "make_link_id",
    "__default_IUPACParser__",
    "__default_IUPACStringMaker__",
]

if __name__ == "__main__":
    string2 = "Man(b1-6)[Man(b1-3)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
    # string = "F(b1-4)[E(a2-3)D(a1-4)]C(a1-6)B(b1-4)A"
    string3 = "Gal(b1-3)GlcNAc(b1-3)[Gal(b1-3)GlcNAc(b1-3)[Gal(b1-4)GlcNAc(b1-6)]Gal(b1-4)GlcNAc(b1-6)][Gal(b1-4)GlcNAc(b1-2)]Gal(b1-4)Glc"

    p = IUPACParser()
    g = p.parse(string2)
    # print(g)
    # g = p.parse(string3)
    # print(g)

    import glycosylator as gl

    s1 = "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)[Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"

    # s = "Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Neu5Gc(a2-8)Neu5Gc(a2-3/6)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc".replace(
    #     "/6", ""
    # )

    mol = gl.read_iupac(s1, s1)
    out = IUPACStringMaker().write_string(mol)
    print(out)
    #  Neu5Gc(a2-8)Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc
    s = "Neu5Gc(a2-8)Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
    # g = p.parse(s)
    # print(g)
    link28 = gl.linkage("O8", "C2", ["HO8"], ["O2", "HO2"], id="28aa")
    gl.add_linkage(link28)
    link23 = gl.linkage("C3", "O2", ["O3", "HO3"], ["HO2"], id="23ba")
    gl.add_linkage(link23)
    mol = gl.read_iupac(s, s)
    out = IUPACStringMaker().write_string(mol)

    s2 = "Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc(a1-4)Glc"
    mol2 = gl.glycan(s2)
    out2 = IUPACStringMaker().write_string(mol2)
    print(out2)

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
