"""
Visualizations for glycosylator
"""
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

import numpy as np
import glycosylator.utils.iupac as iupac
import glycosylator.resources.icons as icons
from biobuild.utils.visual import *

import os

try:
    from glycowork.motif.draw import GlycoDraw
except ImportError:
    GlycoDraw = None


use_glycowork = GlycoDraw is not None
"""
Whether to use GlycoWork from glycowork to draw glycan structures
"""

edge_label_defaults = {
    "font_size": 7,
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
        self.graph = self._make_graph(glycan)

        if use_glycowork:
            self.draw = self._glycowork_draw
            self.show = self._glycowork_show
        else:
            self.draw = self._native_draw
            self.show = self._native_show

    def disable_glycowork(self):
        """
        Disable the use of GlycoWork for drawing (if available)
        """
        self.draw = self._native_draw
        self.show = self._native_show

    def enable_glycowork(self):
        """
        Enable the use of GlycoWork for drawing (if available)
        """
        if GlycoDraw is None:
            raise ImportError("GlycoWork is not available")
        self.draw = self._glycowork_draw
        self.show = self._glycowork_show

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
            parent_grid_idx = int(np.argwhere(node_grid[parent_level] == parent_idx)[0])
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
                int(np.argwhere(node_grid[level] == _all_nodes.index(node) + 1)[0]),
            )
            for node, level in levels.items()
        }

        return positions

    @staticmethod
    def _y_pos(node_grid, level, parent_idx):
        row = node_grid[level]
        dist = 0
        _evals = (
            lambda parent_idx, dist: row[parent_idx - dist] == 0,
            lambda parent_idx, dist: row[parent_idx + dist] == 0,
        )
        _new_dists = (
            lambda parent_idx, dist: parent_idx - dist,
            lambda parent_idx, dist: parent_idx + dist,
        )
        _dir = [0, 1]
        while True:
            if _evals[_dir[0]](parent_idx, dist):
                return _new_dists[_dir[0]](parent_idx, dist)
            elif _evals[_dir[1]](parent_idx, dist):
                return _new_dists[_dir[1]](parent_idx, dist)
            _dir = _dir[::-1]
            dist += 1
            if dist > 3:
                np.roll(row, 1)
                dist = 0
                node_grid[level] = row

    def _native_draw(
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
                _edge_kws = dict(edge_label_defaults).update(edge_label_kws)
            else:
                _edge_kws = edge_label_defaults
            edge_labels = nx.get_edge_attributes(self.graph, "linkage")
            edge_labels = {
                k: iupac.reverse_format_link(getattr(v, "id", v), pretty=True)
                for k, v in edge_labels.items()
            }
            nx.draw_networkx_edge_labels(
                self.graph,
                pos=pos,
                edge_labels=edge_labels,
                ax=ax,
                **_edge_kws,
            )

        # Add the respective image to each node
        for n in self.graph.nodes:
            icon = icons.get_icon(n.resname)
            icon = OffsetImage(icon, zoom=node_size / ax.figure.dpi)
            icon.image.axes = ax
            ab = AnnotationBbox(icon, pos[n], xycoords="data", frameon=False)
            ax.add_artist(ab)
        return ax

    def _native_show(
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
            The axis to draw along, either "x" (horizontal) or "y" (vertical).
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
        self._native_draw(
            ax=ax,
            axis=axis,
            node_size=node_size,
            draw_edge_labels=draw_edge_labels,
            edge_label_kws=edge_label_kws,
            **kwargs,
        )
        plt.show()

    def _glycowork_draw(
        self,
        ax=None,
        axis="x",
        compact: bool = False,
        draw_edge_labels: bool = True,
        svg: bool = False,
    ):
        """
        Draw the glycan structure using GlycoWork

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on
        axis : str
            The axis to draw along, either "x" (horizontal) or "y" (vertical).
        compact : bool
            Whether to draw the glycan in a compact form.
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        svg : bool
            Whether to draw the structure as an SVG image. If True,
            this method will return the SVG as created by GlycoDraw.
            If False, it will return a matplotlib.Axes object.
        """
        vertical = axis == "y"
        img = GlycoDraw(
            self.glycan.to_snfg(),
            compact=compact,
            vertical=vertical,
            show_linkage=draw_edge_labels,
        )
        if not svg:
            img = self._glycowork_to_plt(img)
            if ax is None:
                fig, ax = plt.subplots()
            ax.imshow(img)
            ax.axis("off")
            return ax
        return img

    def _glycowork_show(
        self,
        ax=None,
        axis="x",
        compact: bool = False,
        draw_edge_labels: bool = True,
        svg: bool = False,
    ):
        """
        Draw and show the glycan structure using GlycoWork

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on
        axis : str
            The axis to draw along, either "x" (horizontal) or "y" (vertical).
        compact : bool
            Whether to draw the glycan in a compact form.
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        svg : bool
            Whether to draw the structure as an SVG image. If True,
            this method will return the SVG as created by GlycoDraw.
            If False, it will show the matplotlib figure.

        """
        img = self._glycowork_draw(
            ax=ax,
            axis=axis,
            compact=compact,
            draw_edge_labels=draw_edge_labels,
            svg=svg,
        )
        if svg:
            return img
        plt.show()

    @staticmethod
    def _glycowork_to_plt(img):
        """Convert a drawing to a matplotlib image"""
        img.save_png("tmp.glycan.png")
        img = plt.imread("tmp.glycan.png")
        os.remove("tmp.glycan.png")
        return img

    def _make_graph(self, glycan):
        src = glycan._glycan_tree._segments
        label_src = (i.id for i in glycan._glycan_tree._linkages)
        if len(src) == 0:
            if len(glycan.residues) == 1:
                src = [(glycan.residues[0], glycan.residues[0])]
                label_src = ["<NaN>"]
            else:
                src = glycan.make_residue_graph().edges
                label_src = ("<NaN>" for _ in src)
        graph = nx.Graph(src)

        nx.set_edge_attributes(
            graph,
            {
                i: j
                for i, j in zip(
                    src, [iupac.reverse_format_link(i, pretty=True) for i in label_src]
                )
            },
            "linkage",
        )
        return graph


if __name__ == "__main__":
    import glycosylator as gl

    # man = gl.glycan("MAN")
    # glc = gl.glycan("GLC")
    # gulnac = gl.glycan("GulNAc")

    # glycan = man % "14bb" + glc + glc + man + glc + glc
    # glycan @ -4
    # glycan % "16ab"
    # glycan += man + glc + gulnac
    # glycan % "12bb"
    # glycan += man + gulnac + glc + man
    # glycan = glycan @ -6 + glycan

    # viewer = GlycanViewer2D(glycan)
    # viewer.show()
    use_glycowork = False
    mol = gl.glycan(
        # "Gal(a1-2)Man(a1-3)[Gal(a1-2)Man(a1-3)]Glc(a1-4)[Gal(a1-2)Man(a1-3)][Gal(a1-2)Man(a1-3)]GlcNAc(b1-4)Man(b1-4)Glc(b1-",
        "Neu5Gc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Gc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-8)Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    )
    v = GlycanViewer2D(mol)
    v.show()
    plt.show()
    pass
