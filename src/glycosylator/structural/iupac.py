"""
Functions to work with the IUPAC glycan nomenclature.
"""
import re
from typing import Any


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


class IUPACParser:
    """
    A parser for IUPAC glycan nomenclature strings. This class will generate
    a list of connecting glycan segments from a string from which a Molecule can be built.
    """

    def __init__(self):
        self._string = ""
        self.reset()

    @property
    def _current(self):
        return self._string[-self._idx]

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
            and len(self._latest_residue) >= 3
            and len(self._latest_linkage) > 3
            and len(self._second_latest_residue) >= 3
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
        self._string = string
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
        self._latest_linkage = ""
        self._second_latest_residue = ""
        self._latest_residue_before_square_bracket = ""

    def _shift(self):
        self._idx += 1

    def _shift_residue(self):
        self._second_latest_residue = self._latest_residue
        self._latest_residue = ""

    def _parse(self):
        self._crop_end()
        while self._can_shift:
            if self._can_store:
                self._store()
                continue
            if self._is_at_bracket:
                self._latest_linkage = self._parse_linkage()
                self._shift_residue()
                self._shift()
                continue
            if self._is_at_open_square_bracket:
                self._latest_residue_before_square_bracket = self._fit_residue(
                    self._latest_residue, increase_count=False
                )
                self._latest_residue = self._latest_residue_before_square_bracket
                self._shift()
                continue
            if self._is_at_close_square_bracket:
                self._latest_residue = self._latest_residue_before_square_bracket
                self._second_latest_residue = ""
                self._shift()
                continue
            self._latest_residue += self._current
            self._shift()
        self._latest_residue += self._current
        self._store()

    def _store(self):
        latest = self._latest_residue
        if not "@" in latest:
            latest = self._fit_residue(latest)
            self._latest_residue = latest

        second = self._second_latest_residue
        if "@" not in second:
            second = self._fit_residue(second)
            self._second_latest_residue = second

        branch = (second, latest, reformat_link(self._latest_linkage[::-1]))
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
        return linkage

    def _crop_end(self):
        if self._string[-1] == "-":
            while self._current != "(":
                self._idx -= 1
        self._string = self._string[: self._idx - 1]
        self._idx = 1

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.parse(*args, **kwds)


if __name__ == "__main__":
    string2 = "Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"

    p = IUPACParser()
    g = p.parse(string2)
    print(g)
