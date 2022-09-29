import networkx as nx
from prody import AtomGroup

ELEMENTS = ["C", "H", "N", "O"]
BOND_LENGTHS = {
    "C-H": 1.20,
    "H-N": 1.20,
    "H-O": 1.20,
    "C-C": 1.7,
    "C-N": 1.7,
    "C-O": 1.7,
}


def build_bond_graph(atom_group: AtomGroup, default_bond_length=1.6):
    bond_graph = nx.Graph()
    for a in atom_group:
        sel = ""
        for e in ELEMENTS:
            key = "-".join(sorted([a.getElement(), e]))
            if key in BOND_LENGTHS:
                if sel:
                    sel += " or "
                sel += f'((element "{e}") and (within "{BOND_LENGTHS[key]}" of serial "{a.getSerial()}"))'

        if not sel:
            sel = f"within {default_bond_length} of serial {a.getSerial()}"

        sel = f"({sel}) and (not serial {a.getSerial()})"

        bonds = [
            (a.getIndex(), neighbour.getIndex()) for neighbour in atom_group.select(sel)
        ]
        bond_graph.add_edges_from(bonds)

    return bond_graph  # bond_graph.edges() ?
