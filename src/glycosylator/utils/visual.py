"""
Visualizations for glycosylator
"""
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

import numpy as np
from biobuild.utils.visual import *
import glycosylator.utils.iupac as iupac
import glycosylator.resources.icons as icons

edge_label_defaults = {
    "font_size": 5,
    "font_color": "k",
    "font_weight": "normal",
    "alpha": 1.0,
    "rotate": True,
    "horizontalalignment": "center",
    "clip_on": False,
}
"""
Default edge label drawing parameters
"""


class GlycanViewer2D:
    """
    A 2D scematic viewer for glycan structures
    """

    def __init__(self, glycan):
        self.glycan = glycan
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

    def layout(self, graph) -> dict:
        """
        Make a vertical fork layout

        Parameters
        ----------
        graph : nx.Graph
            The graph to layout

        Returns
        -------
        dict
            A dictionary of node positions
        """
        nodes = [self.root_residue]
        levels = {self.root_residue: 0}
        child_mapping = {self.root_residue: []}
        idx = 0
        graph = nx.dfs_tree(graph, self.root_residue)
        while len(levels) < len(graph):
            parent_node = nodes[idx]
            for node in nx.descendants_at_distance(graph, parent_node, 1):
                levels.setdefault(node, levels[parent_node] + 1)
                nodes.append(node)
                child_mapping.setdefault(parent_node, [])
                child_mapping[parent_node].append(node)
            idx += 1
        parent_mapping = {i: k for k, v in child_mapping.items() for i in v}

        nodes_per_level = {}
        for node, level in levels.items():
            nodes_per_level.setdefault(level, []).append(node)
        nodes_per_level.pop(0)

        _all_nodes = list(graph.nodes)
        n_levels = len(nodes_per_level)
        self._max_nodes = max(map(len, nodes_per_level.values()))
        node_grid = np.zeros((n_levels + 1, self._max_nodes * 2 + 1), dtype=int)
        node_grid[0, self._max_nodes] = _all_nodes.index(self.root_residue) + 1

        for child, parent in parent_mapping.items():
            parent_level = levels[parent]
            parent_idx = _all_nodes.index(parent) + 1
            parent_grid_idx = int(np.argwhere(node_grid[parent_level] == parent_idx))
            child_level = levels[child]
            y = self._y_pos(
                node_grid,
                child_level,
                parent_grid_idx,
            )
            node_grid[child_level, y] = _all_nodes.index(child) + 1

        positions = {
            node: (
                level,
                int(np.argwhere(node_grid[level] == _all_nodes.index(node) + 1)),
            )
            for node, level in levels.items()
        }

        return positions

    @staticmethod
    def _y_pos(node_grid, level, parent_idx):
        row = node_grid[level]
        dist = 0
        while True:
            if row[parent_idx - dist] == 0:
                return parent_idx - dist
            elif row[parent_idx + dist] == 0:
                return parent_idx + dist
            dist += 1

    def draw(
        self,
        ax=None,
        axis="y",
        node_size: int = 20,
        draw_edge_labels: bool = False,
        edge_label_kws: dict = None,
        **kwargs
    ) -> plt.Axes:
        """
        Draw the glycan structure

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on
        axis : str
            The axis to draw on, either "x" or "y".
        node_size : int
            The size of the nodes
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        edge_label_kws : dict
            Keyword arguments to pass to the edge label drawing function.
            The defaults are defined in `edge_label_defaults`.
        kwargs : dict
            Keyword arguments to pass to `nx.draw_networkx_edges`

        Returns
        -------
        matplotlib.axes.Axes
        """
        aspect = 0.8 if axis == "y" else 1  # 0.5 if axis == "x" else 1
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect(aspect)
        pos = self.layout(self.graph)
        if axis == "y":
            pos = {k: (v[1], self._max_nodes + v[0]) for k, v in pos.items()}

        nx.draw_networkx_edges(self.graph, pos=pos, ax=ax, **kwargs)
        ax.axis("off")

        if draw_edge_labels:
            if edge_label_kws:
                edge_label_defaults.update(edge_label_kws)

            edge_labels = nx.get_edge_attributes(self.graph, "linkage")
            edge_labels = {
                k: iupac.reverse_format_link(v.id, pretty=True)
                for k, v in edge_labels.items()
            }
            nx.draw_networkx_edge_labels(
                self.graph,
                pos=pos,
                edge_labels=edge_labels,
                ax=ax,
                **edge_label_defaults,
            )

        # Add the respective image to each node
        for n in self.graph.nodes:
            icon = icons.get_icon(n.resname)
            icon = OffsetImage(icon, zoom=node_size / ax.figure.dpi)
            icon.image.axes = ax
            ab = AnnotationBbox(icon, pos[n], xycoords="data", frameon=False)
            ax.add_artist(ab)
        return ax

    def show(
        self,
        ax=None,
        axis="y",
        node_size: int = 20,
        draw_edge_labels: bool = False,
        edge_label_kws: dict = None,
        **kwargs
    ):
        """
        Draw and show the glycan structure

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on
        axis : str
            The axis to draw on, either "x" or "y".
        node_size : int
            The size of the nodes
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        edge_label_kws : dict
            Keyword arguments to pass to the edge label drawing function.
            The defaults are defined in `edge_label_defaults`.
        kwargs : dict
            Keyword arguments to pass to `nx.draw_networkx_edges`
        """
        self.draw(
            ax=ax,
            axis=axis,
            node_size=node_size,
            draw_edge_labels=draw_edge_labels,
            edge_label_kws=edge_label_kws,
            **kwargs,
        )
        plt.show()


if __name__ == "__main__":
    import glycosylator as gl

    man = gl.glycan("MAN")
    glc = gl.glycan("GLC")
    gulnac = gl.glycan("GulNAc")

    glycan = man % "14bb" + glc + glc + man + glc + glc
    glycan @ -4
    glycan % "16ab"
    glycan += man + glc + gulnac
    glycan % "12bb"
    glycan += man + gulnac + glc + man
    glycan = glycan @ -6 + glycan

    print(glycan)
    viewer = GlycanViewer2D(glycan)
    viewer.draw(axis="y", node_size=20)
    plt.show()
    pass
