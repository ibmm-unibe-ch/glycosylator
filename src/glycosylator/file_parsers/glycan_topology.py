import copy
import sqlite3
from dataclasses import dataclass
from typing import NamedTuple


class ResiduePath(NamedTuple):
    """The sequence of patches that go from a glycan's root to the residue"""

    residue_name: str
    linking_atom: str
    patches: list[str]

    @classmethod
    def from_str(cls, path_string: str):
        path_string = path_string.removeprefix("UNIT ")
        residue_name, *data = path_string.split(" ")
        if data:
            linking_atom, *patches = data
        else:
            linking_atom, patches = "", [""]
        return cls(residue_name, linking_atom, patches)

    def __str__(self) -> str:
        patches_str = " ".join(self.patches)
        return f"{self.residue_name} {self.linking_atom} {patches_str}"


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


def read_glycan_topology(file_path: str) -> dict[str, GlycanTopology]:
    """read glycan.top file"""
    with open(file_path, "r") as f:
        # remove comments
        lines = [line.split("!")[0] for line in f.readlines()]
        # remove whitespace
        lines = [line.strip() for line in lines]
        # remove blank lines
        lines = [line for line in lines if line != ""]
        # rejoin and split by each new glycan entry
        # discard everything before first glycan entry (should be just empty string)
        _, *glycans = "\n".join(lines).split("RESI ")

    topologies = {}
    for glycan in glycans:
        glycan_name, *lines = glycan.splitlines()
        paths = [ResiduePath.from_str(line) for line in lines]
        topology = GlycanTopology(glycan_name, paths, len(paths))
        topologies[topology.name] = topology
    return topologies


def write_glycan_topologies(
    file_path: str, topologies: dict[str, GlycanTopology] | GlycanTopology
):
    # handle case if a single GlycanTopology is provided instead of a dict of GlycanTopologies
    if type(topologies) == GlycanTopology:
        topologies = {topologies.name: topologies}

    output = "\n\n\n".join([str(topology) for topology in topologies.values()])

    with open(file_path, "w") as f:
        f.write(output)


def import_connectivity_topology(filename: str) -> dict[str, GlycanTopology]:
    """Import connectivity topology from sql database
    This function will initialize connect_topology
    Parameters:
        filename: path to database
    """
    try:
        conn = sqlite3.connect(filename)
    except:
        print(f"Error while connecting to the database {filename}")
        return -1
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM glycans")
    glycans = cursor.fetchall()

    topologies = {}
    for glycan in glycans:
        glycan_name, tree = glycan
        paths = [ResiduePath.from_str(line) for line in tree.split("|")]
        topology = GlycanTopology(glycan_name, paths, len(paths))
        topologies[topology.name] = topology
    return topologies


def export_connectivity_topology(filename, connect_topology: dict[str, GlycanTopology]):
    """Export connectivity topology to sql database
    This function will
    """
    try:
        conn = sqlite3.connect(filename)
    except:
        print(f"Error while connecting to the database {filename}")
        return -1
    cursor = conn.cursor()
    tn = "glycans"
    gn = "glycan_name"
    gt = "glycan_tree"
    cursor.execute(f"DROP TABLE IF EXISTS {tn}")
    cursor.execute(f"CREATE TABLE {tn} ({gn} text, {gt} text)")

    for key, value in connect_topology.items():
        glycan = [str(path) for path in value.paths]
        glycan = "|".join(glycan)
        cursor.execute(f"INSERT INTO {tn} VALUES ('{key}', '{glycan}')")

    conn.commit()
    conn.close()


if __name__ == "__main__":
    t = read_glycan_topology(
        "/Users/borisgusev/Bioinformatics/glycosylator/support/topology/man9.top"
    )
