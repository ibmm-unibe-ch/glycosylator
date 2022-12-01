from dataclasses import dataclass
from typing import NamedTuple


class ResiduePath(NamedTuple):
    """The sequence of patches that go from a glycan's root to the residue"""

    residue_name: str
    linking_atom: str
    patches: tuple[str]

    @classmethod
    def from_str(cls, path_string: str):
        path_string = path_string.removeprefix("UNIT ")
        residue_name, *data = path_string.split(" ")
        if data:
            linking_atom, *patches = data
            patches = tuple(patches)
        else:
            linking_atom, patches = "", ()
        return cls(residue_name, linking_atom, patches)

    def __str__(self) -> str:
        patches_str = " ".join(self.patches)
        return f"{self.residue_name} {self.linking_atom} {patches_str}"

    def __len__(self) -> int:
        return len(self.patches)


@dataclass
class GlycanTopology:
    """Class containing all the paths to residues in a glycan"""

    name: str
    paths: list[ResiduePath]
    num_residues: int

    def __str__(self) -> str:
        name_string = "RESI " + self.name
        path_strings = ["UNIT " + str(path) for path in self.paths]
        return "\n".join([name_string, *path_strings])

    def __len__(self) -> int:
        return len(self.paths)

    def __iter__(self):
        return iter(self.paths)

    @classmethod
    def from_patch_strings(cls, name: str, strings: list[str]):
        paths = [ResiduePath.from_str(s) for s in strings]
        paths.sort(key=len)
        return cls(name, paths, len(paths))

    @classmethod
    def from_tuples(cls, name: str, patch_tuples: tuple[str, str, list[str]]):
        # TODO: extract the linking atom
        paths = [
            ResiduePath(name, atom, tuple(patches))
            for name, atom, patches in patch_tuples
        ]
        paths.sort(key=len)
        return cls(name, paths, len(paths))
