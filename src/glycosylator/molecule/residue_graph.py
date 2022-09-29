import networkx as nx
from prody import AtomGroup, KDTree

BOND_LENGTH = 1.7


def build_residue_graph(atom_group: AtomGroup, root_residue):
    # this is also finding bonds -> redundant as have bond_graph?
    kd = KDTree(atom_group.getCoords())
    kd.search(BOND_LENGTH)
    atom_pairs = kd.getIndices()

    res_n, res_names = atom_group.getResnums(), atom_group.getResnames()

    graph = nx.Graph()

    for a1, a2 in atom_pairs:
        a1_num, a2_num = res_n[a1], res_n[a2]
        a1_name, a2_name = res_names[a1], res_names[a2]
        if a1_num != a2_num:
            graph.add_node(a1_num, resname=a1_name)
            graph.add_node(a2_num, resname=a2_name)
            graph.add_edge(a1_num, a2_num)

    # create directed graph and remove all unnecessary edges
    if graph:
        residue_graph = graph.to_directed()
        for edges in nx.dfs_edges(graph, root_residue):
            edges = list(edges)
            edges.reverse()
            residue_graph.remove_edge(*edges)

    return residue_graph
