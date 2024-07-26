"""
Visualizations for glycosylator
"""

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# from matplotlib.offsetbox import OffsetImage, AnnotationBbox
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import numpy as np
import glycosylator.utils.iupac as iupac
import glycosylator.resources.icons as icons
from buildamol.utils.visual import *

import os
import importlib

# try:
#     from glycowork.motif.draw import GlycoDraw
# except ImportError:
#     GlycoDraw = None


USE_GLYCOWORK = False
"""
Whether to use GlycoDraw from GlycoWork to draw glycan structures
"""


def use_glycowork():
    """
    Use GlycoDraw from GlycoWork to draw glycan structures
    """
    global USE_GLYCOWORK
    USE_GLYCOWORK = True


def dont_use_glycowork():
    """
    Don't use GlycoDraw from GlycoWork to draw glycan structures
    """
    global USE_GLYCOWORK
    USE_GLYCOWORK = False


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
        self.canvas = None

        if len(glycan._glycan_tree._segments) == 0:
            glycan.infer_glycan_tree()

        if USE_GLYCOWORK:
            if importlib.util.find_spec("glycowork") is not None:
                global GlycoDraw
                from glycowork.motif.draw import GlycoDraw
            else:
                raise ImportError("GlycoWork is not available")

            self.draw = self._glycowork_draw
            self.show = self._glycowork_show
        else:
            self.draw = self._native_draw
            self.show = self._native_show
            glycan.set_root(self.root)
            self.canvas = Canvas.from_glycan(glycan)

    def disable_glycowork(self):
        """
        Disable the use of GlycoWork for drawing (if available)
        """
        self.draw = self._native_draw
        self.show = self._native_show
        if self.canvas is None:
            self.glycan.set_root(self.root)
            self.canvas = Canvas.from_glycan(self.glycan)

    def enable_glycowork(self):
        """
        Enable the use of GlycoWork for drawing (if available)
        """
        if GlycoDraw is None:
            raise ImportError("GlycoWork is not available")
        self.draw = self._glycowork_draw
        self.show = self._glycowork_show

    def _native_draw(
        self,
        ax=None,
        axis="x",
        draw_edge_labels: bool = True,
        pretty_labels: bool = True,
        small_labels: bool = False,
        edge_label_kws: dict = None,
        edge_kws: dict = None,
        adjust_axes: bool = True,
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
        pretty_labels : bool
            Whether to use the pretty format for the edge labels. Including greek letters and arrows.
        small_labels : bool
            Whether to use the small format for the edge labels, without arrows.
        edge_label_kws : dict
            Keyword arguments to pass to the edge label drawing function.
            The defaults are defined in `edge_label_defaults`.

        Returns
        -------
        matplotlib.axes.Axes
        """
        ax = ax or plt.gca()
        self.canvas.tree_layout()
        if axis == "x":
            self.canvas.rotate90()
        self.canvas.render(
            ax,
            draw_edge_labels=draw_edge_labels,
            pretty_labels=pretty_labels,
            small_labels=small_labels,
            label_kws=edge_label_kws,
            edge_kws=edge_kws,
            adjust_axes=adjust_axes,
        )
        return ax

    def _native_draw_no_layout(
        self,
        ax,
        draw_edge_labels: bool = True,
        pretty_labels: bool = True,
        small_labels: bool = False,
        edge_label_kws: dict = None,
        edge_kws: dict = None,
        adjust_axes: bool = True,
    ) -> plt.Axes:
        """
        Draw the glycan structure

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on
        node_size : int
            The size of the nodes
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        pretty_labels : bool
            Whether to use the pretty format for the edge labels. Including greek letters and arrows.
        small_labels : bool
            Whether to use the small format for the edge labels, without arrows.
        edge_label_kws : dict
            Keyword arguments to pass to the edge label drawing function.
            The defaults are defined in `edge_label_defaults`.

        Returns
        -------
        matplotlib.axes.Axes
        """
        self.canvas.render(
            ax,
            draw_edge_labels=draw_edge_labels,
            pretty_labels=pretty_labels,
            small_labels=small_labels,
            label_kws=edge_label_kws,
            edge_kws=edge_kws,
            adjust_axes=adjust_axes,
        )
        return ax

    def _native_show(
        self,
        ax=None,
        axis="x",
        draw_edge_labels: bool = True,
        pretty_labels: bool = True,
        small_labels: bool = False,
        edge_label_kws: dict = None,
        edge_kws: dict = None,
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
        pretty_labels : bool
            Whether to use the pretty format for the edge labels. Including greek letters and arrows.
        small_labels : bool
            Whether to use the small format for the edge labels, without arrows.
        edge_label_kws : dict
            Keyword arguments to pass to the edge label drawing function.
            The defaults are defined in `edge_label_defaults`.
        """
        self._native_draw(
            ax=ax,
            axis=axis,
            draw_edge_labels=draw_edge_labels,
            pretty_labels=pretty_labels,
            small_labels=small_labels,
            edge_label_kws=edge_label_kws,
            edge_kws=edge_kws,
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
            iupac.make_iupac_string(self.glycan, add_terminal_conformation=False),
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


class ScaffoldViewer2D:
    """
    A 2D schematic for a protein scaffold

    Parameters
    ----------
    scaffold : Scaffold
        The protein scaffold to visualize.
    bar_color : str
        The color to use for the bar representing the protein sequence.
    bar_height : int
        The height of the bar representing the protein sequence.
    height : float
        The total height of each subplot (chain) in the figure.
    """

    def __init__(
        self,
        scaffold,
        bar_color: str = "xkcd:sandy",
        bar_height: int = 6,
        height: float = 20,
        figsize: tuple = None,
    ):
        self.scaffold = scaffold
        self._chain_lengths = {
            chain: sum(1 for i in seq if i != "X")
            for chain, seq in scaffold.get_sequence().items()
            if chain not in scaffold._excluded_chains
        }
        for i in scaffold._excluded_chains:
            self._chain_lengths[i] = 0
        _remove = []
        for chain, length in self._chain_lengths.items():
            if length == 0:
                _remove.append(chain)
        for i in _remove:
            del self._chain_lengths[i]

        self.fig, self.subplots = plt.subplots(
            len(self._chain_lengths), 1, figsize=figsize
        )
        if len(self._chain_lengths) == 1:
            self.subplots = [self.subplots]
        self._subplot_dict = {
            chain: ax for chain, ax in zip(self._chain_lengths.keys(), self.subplots)
        }
        self.bar_color = bar_color
        self.bar_height = bar_height
        self.height = height
        self.ylim = -(self.height - self.bar_height), self.height
        self._setup()

    def _setup(self):
        """
        Setup the bar to represent the protein sequence
        """
        for i, ((chain, base), ax) in enumerate(
            zip(self._chain_lengths.items(), self.subplots)
        ):
            ax.axhspan(
                0,
                self.bar_height,
                color=self.bar_color,
                alpha=0.5,
            )
            ax.yaxis.set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.set_xticks(np.arange(0, base, base // 10))
            ax.set_ylim(*self.ylim)
            ax.set_title("chain: " + str(chain.id))

        self.fig.suptitle(self.scaffold.id if self.scaffold.id else "")
        self.fig.tight_layout()

    def get_figure(self):
        """
        Get the figure of the viewer
        """
        return self.fig

    def save(self, filename: str, dpi: int = 300, **kwargs):
        """
        Save the figure to a file

        Parameters
        ----------
        filename : str
            The filename to save to.
        dpi : int
            The resolution of the image.
        """
        self.fig.savefig(filename, dpi=dpi, **kwargs)

    def title(self, t: str):
        """
        Set the title of the figure

        Parameters
        ----------
        t : str
            The title to set
        """
        self.fig.suptitle(t)

    def draw_glycans(self, scalar: float = 1):
        """
        Draw the glycans

        Parameters
        ----------
        scalar : float
            A scalar to use for artificially upscaling the glycans.
        """
        self._got_flipped = False
        drawn = set()
        scales = {chain: self.height for chain in self.scaffold.chains}
        flipped = {chain: False for chain in self.scaffold.chains}
        for root, glycan in self.scaffold.get_glycans().items():

            chain, rdx = self.scaffold.index(root.parent)
            _chain = self.scaffold.get_chain(chain)
            if _chain in self.scaffold._excluded_chains:
                continue

            ax = self._subplot_dict[_chain]

            v = GlycanViewer2D(glycan)
            v.disable_glycowork()
            v.canvas.tree_layout()

            # seems like the canvas is inverted already
            if v.canvas.ylim[1] == 0 and v.canvas.ylim[0] < 0:
                v.canvas.flip(vertical=True)
                # v.canvas.root_nodes[0].move_to(
                #     v.canvas.root_nodes[0].x, -v.canvas.scalar * 0.3
                # )
                # self._got_flipped = True
            _y = max(1, v.canvas.ylim[1])
            scale = scalar * (self.height / _y)
            v.canvas.scale(scale)
            _y = max(1, v.canvas.ylim[1])
            scales[_chain] = max(scales[_chain], _y)

            self._switch_direction(v, drawn, chain, rdx)
            flipped[_chain] = self._got_flipped

            v._native_draw_no_layout(ax, draw_edge_labels=False, adjust_axes=False)

        for chain, ax in self._subplot_dict.items():
            ytop = scales[chain]
            _got_flipped = flipped[chain]
            ybottom = (
                -(1.2 * self.bar_height)
                if not _got_flipped
                else -(self.bar_height + ytop * 1.05)
            )
            ax.set_xlim(0, self._chain_lengths[chain])
            ax.set_ylim(ybottom, (self.bar_height + ytop * 1.2))
            # ax.set_box_aspect(0.3)

    def draw_glycosylation_sites(
        self,
        n_linked_color="xkcd:sky blue",
        o_linked_color="xkcd:medium pink",
        draw_legend: bool = True,
    ):
        """
        Highlight glycosylation sites

        Parameters
        ----------
        n_linked_color : str
            The color to use for N-linked glycosylation sites.
        o_linked_color : str
            The color to use for O-linked glycosylation sites.
        draw_legend : bool
            Whether to draw a legend for the glycosylation sites.
        """
        n_linked = self.scaffold.find_n_linked_sites()
        o_linked = self.scaffold.find_o_linked_sites()

        for _dict, _color in zip(
            (n_linked, o_linked), (n_linked_color, o_linked_color)
        ):
            for chain, residues in _dict.items():
                ax = self._subplot_dict.get(chain, None)
                if not ax:
                    continue
                for residue in residues:
                    _, idx = self.scaffold.index(residue)
                    rect = Rectangle(
                        (idx - 1, 0),
                        2,
                        self.bar_height,
                        facecolor=_color,
                        alpha=0.5,
                    )
                    ax.add_patch(rect)
        if draw_legend:
            self.subplots[0].legend(
                labels=["N-linked", "O-linked"],
                handles=[
                    Rectangle((0, 0), 1, 1, color=n_linked_color),
                    Rectangle((0, 0), 1, 1, color=o_linked_color),
                ],
                loc="upper left",
                title="Glycosylation Sites",
                frameon=False,
                bbox_to_anchor=(1.0, 1.1),
            )
            self.fig.tight_layout()

    def show(self):
        """
        Show the figure
        """
        self.fig.show()

    def _switch_direction(self, v, drawn, chain, rdx):
        """
        Make sure that close-by glycans are drawn in different directions
        """

        min_diff = min((20, *(rdx - i for _, i, d in drawn if _ == chain)))
        direction = 1
        self._got_flipped = False

        if min_diff < 20:
            closest = next(i for i in drawn if i[0] == chain and rdx - i[1] == min_diff)
            ref_direction = closest[2]
            if ref_direction == 1:
                v.canvas.flip(vertical=True)
                v.canvas.root_nodes[0].move_to(rdx, -v.canvas.scalar * 0.3)
                direction = 0
                self._got_flipped = True
            else:
                v.canvas.root_nodes[0].move_to(
                    rdx, self.bar_height + v.canvas.scalar * 0.3
                )

        else:
            v.canvas.root_nodes[0].move_to(rdx, self.bar_height + v.canvas.scalar * 0.3)

        drawn.add((chain, rdx, direction))

        # for _, i, direction in drawn:
        #     if _ != chain:
        #         continue
        #     if rdx - i < 20 and direction == 1:
        #         v.canvas.flip(vertical=True)
        #         v.canvas.root_nodes[0].move_to(rdx, -v.canvas.scalar * 0.3)
        #         direction = 0
        #         break
        # else:
        #     direction = 1
        #     v.canvas.root_nodes[0].move_to(rdx, self.bar_height + v.canvas.scalar * 0.3)
        # drawn.add((chain, rdx, direction))


# =============================================================================


class Node:
    """
    A node in a tree structure.

    Parameters
    ----------
    id : str
        The id of the node.
    icon : np.ndarray
        The icon to represent the node in SNFG format.

    Attributes
    ----------
    id : str
        The id of the node.
    icon : np.ndarray
        The icon of the node.
    children : list
        The children of the node.
    parent : Node
        The parent of the node.
    canvas : Canvas
        The canvas that the node is placed on.
    x : int
        The x coordinate of the node.
    y : int
        The y coordinate of the node.
    """

    def __init__(self, id, icon=None):
        self.id = id
        self.icon = icon
        self.children = []
        self.parent = None
        self.canvas = None

        self.x = 0
        self.y = 0

    @property
    def depth(self):
        """
        The depth of the node.
        """
        if self.parent is None:
            return 0
        else:
            return self.parent.depth + 1

    @property
    def siblings(self) -> list:
        """
        Sibling nodes that share the same parent.
        """
        if self.parent is None:
            return []
        else:
            return [i for i in self.parent.children if i != self]

    def __eq__(self, other):
        if isinstance(other, Node):
            return self.id == other.id
        elif isinstance(other, str):
            return self.id == other
        else:
            raise TypeError("Node can only be compared with Node or str")

    def __repr__(self):
        return f"Node({self.id})"

    def move_to(self, x: int, y: int):
        """
        Move this node to (x, y) on the canvas.
        """
        if self.canvas is None:
            raise ValueError("Node has no canvas")
        dx = x - self.x
        dy = y - self.y
        self.shift_down(dx, dy)
        # self.shift_up(dx, dy)

    def shift_down(self, dx, dy):
        """
        Shift this node and all children by dx and dy.
        """
        self.x += dx
        self.y += dy
        for child in self.children:
            child.shift_down(dx, dy)

    def shift_up(self, dx, dy):
        """
        Shift this node and all parents by dx and dy.
        This will also affet all siblings and their children.
        """
        if self.parent is not None:
            for child in self.parent.children:
                child.shift_down(dx, dy)
            self.parent.shift_up(dx, dy)

    def place(self, canvas: "Canvas", x: int, y: int):
        """
        Place this node at (x, y) on the canvas.

        Parameters
        ----------
        canvas : Canvas
            The canvas to place the node on.
        x : int
            The x coordinate of the node.
        y : int
            The y coordinate of the node.
        """
        self.x = x
        self.y = y
        self.canvas = canvas

    def has_valid_placement(self) -> bool:
        """
        Check if a node has a valid placement on the canvas.
        Placements are valid if there is no other node at the same position.
        """
        if self.canvas is None:
            raise ValueError("Node has no canvas")
        for _node in self.canvas.all_nodes:
            if self == _node:
                continue
            if _node.x == self.x and _node.y == self.y:
                return False
        return True

    def add_child(self, child: "Node", x: int = None, y: int = None):
        """
        Add a child node to this node.

        Parameters
        ----------
        child : Node
            The child node to add.
        x : int
            The x coordinate of the child node.
        y : int
            The y coordinate of the child node.
        """
        x = self.x if x is None else x
        y = self.y + 1 if y is None else y
        child.place(self.canvas, x, y)
        self.children.append(child)
        child.parent = self
        self.canvas.all_nodes.append(child)

    def remove_child(self, child: "Node"):
        """
        Remove a child node from this node.
        """
        self.children.remove(child)
        self.canvas.all_nodes.remove(child)
        child.parent = None

    def has_descendant(self, node: "Node") -> bool:
        """
        Check if a node is a descendant of this node.
        """
        if node in self.children:
            return True
        else:
            for child in self.children:
                if child.is_descendant(node):
                    return True
        return False

    def get_descendants(self) -> list:
        """
        Get the descendants of this node.
        """
        descendants = []
        for child in self.children:
            descendants.append(child)
            descendants.extend(child.get_descendants())
        return descendants

    def get_ancestors(self) -> list:
        """
        Get the ancestors of this node.
        """
        ancestors = []
        if self.parent is not None:
            ancestors.append(self.parent)
            ancestors.extend(self.parent.get_ancestors())
        return ancestors

    def get_left_siblings(self):
        """
        Get the left siblings of this node.
        """
        if self.parent is None:
            return []
        else:
            return self.parent.children[: self.parent.children.index(self)]

    def get_right_siblings(self):
        """
        Get the right siblings of this node.
        """
        if self.parent is None:
            return []
        else:
            return self.parent.children[self.parent.children.index(self) + 1 :]

    def get_left_cousins(self):
        """
        Get the left cousins of this node.
        """
        if self.parent is None:
            return []
        else:
            return self.parent.get_left_siblings()

    def get_right_cousins(self):
        """
        Get the right cousins of this node.
        """
        if self.parent is None:
            return []
        else:
            return self.parent.get_right_siblings()


class Canvas:
    """
    A canvas of 2D slots for placing nodes.

    Attributes
    ----------
    root_nodes : list
        The root nodes on the canvas.
    all_nodes : list
        The nodes on the canvas (not regarding connectivity or hierarchy).
    """

    def __init__(self):
        self.root_nodes = []
        self.all_nodes = []
        self.data = {}
        self.scalar = 1

    @property
    def leaves(self) -> list:
        """
        The leaves of the canvas.
        """
        return [node for node in self.all_nodes if len(node.children) == 0]

    @property
    def xlim(self) -> tuple:
        """
        The x limits of the canvas.
        """
        return min([node.x for node in self.all_nodes]), max(
            [node.x for node in self.all_nodes]
        )

    @property
    def ylim(self) -> tuple:
        """
        The y limits of the canvas.
        """
        return min([node.y for node in self.all_nodes]), max(
            [node.y for node in self.all_nodes]
        )

    def get_xlim(self):
        """
        Get the x limits of the canvas.
        """
        return self.xlim

    def get_ylim(self):
        """
        Get the y limits of the canvas.
        """
        return self.ylim

    @staticmethod
    def _make_nodes_from_dict(parent_node, child_dict):
        for key, value in child_dict.items():
            node = Node(key)
            parent_node.add_child(node)
            if isinstance(value, dict):
                Canvas._make_nodes_from_dict(node, value)
            elif isinstance(value, (list, tuple, set)):
                for i in value:
                    node.add_child(Node(i))
            else:
                node.add_child(Node(value))

    @classmethod
    def from_glycan(cls, glycan):
        new = cls()
        root = Node(glycan.root_residue, icons.get_icon(glycan.root_residue.resname))
        new.add(root)

        new.data["bonds"] = {}

        nodes = {glycan.root_residue: root}
        for root, children in glycan._glycan_tree._connectivity.items():
            for child in children:
                if child not in nodes:
                    node = Node(child, icons.get_icon(child.resname))
                    nodes[child] = node
                else:
                    node = nodes[child]

                nodes[root].add_child(node)
                link = glycan._glycan_tree.get_linkage(root, child)
                new.data["bonds"][(root, child)] = link
        return new

    @classmethod
    def from_dict(cls, __dict):
        """
        Make a canvas from a dictionary.

        Parameters
        ----------
        __dict : dict
            The dictionary to make the canvas from.

        Returns
        -------
        canvas : Canvas
            The canvas made from the dictionary.
        """
        new = cls()
        for key, value in __dict.items():
            node = Node(key)
            new.add(node)
            if isinstance(value, dict):
                cls._make_nodes_from_dict(node, value)
            elif isinstance(value, (list, tuple, set)):
                for i in value:
                    node.add_child(Node(i))
            else:
                node.add_child(Node(value))
        return new

    def flip(self, horizontal: bool = False, vertical: bool = True):
        """
        Flip the canvas horizontally or vertically.

        Parameters
        ----------
        horizontal : bool
            Whether to flip the canvas horizontally.
        vertical : bool
            Whether to flip the canvas vertically.
        """
        if horizontal:
            for node in self.all_nodes:
                node.x = -node.x
        if vertical:
            for node in self.all_nodes:
                node.y = -node.y

    def rotate90(self, clockwise: bool = True):
        """
        Rotate the canvas 90 degrees.

        Parameters
        ----------
        clockwise : bool
            Whether to rotate the canvas clockwise.
        """
        for node in self.all_nodes:
            node.x, node.y = (node.y, -node.x) if clockwise else (-node.y, node.x)

    def has_node(self, node_or_id):
        """
        Check if a node is on the canvas.

        Parameters
        ----------
        node_or_id : Node or str
            The node or id to check for.
        """
        if isinstance(node_or_id, Node):
            return node_or_id in self.root_nodes
        elif isinstance(node_or_id, str):
            return node_or_id in [node.id for node in self.root_nodes]
        else:
            raise TypeError("Node can only be compared with Node or str")

    def add(self, node: Node, x: int = None, y: int = None):
        """
        Add a node at a specific location on the canvas.

        Parameters
        ----------
        node : Node
            The node to place.
        x : int
            The x coordinate to place the node at.
        y : int
            The y coordinate to place the node at.
        """
        if x is not None and y is not None:
            node.place(self, x, y)
        else:
            node.canvas = self
        self.root_nodes.append(node)
        self.all_nodes.append(node)

    def remove(self, node: Node):
        """
        Remove a node from the canvas.

        Parameters
        ----------
        node : Node
            The node to remove.
        """
        self.all_nodes.remove(node)
        self.root_nodes.remove(node)
        node.canvas = None

    def is_free(self, x: int, y: int) -> bool:
        """
        Check if a slot is free.

        Parameters
        ----------
        x : int
            The x coordinate of the slot.
        y : int
            The y coordinate of the slot.

        Returns
        -------
        is_free : bool
            Whether the slot is free.
        """
        for node in self.all_nodes:
            if node.x == x and node.y == y:
                return False
        return True

    def get_clashing_nodes(self) -> iter:
        """
        Get an iterator of nodes that are clashing with each other.
        """
        found_clashes = False
        while True:
            for node in self.all_nodes:
                for other_node in self.all_nodes:
                    if node == other_node:
                        continue
                    if node.x == other_node.x and node.y == other_node.y:
                        found_clashes = True
                        yield node, other_node

            if not found_clashes:
                break
            else:
                found_clashes = False

    def scale(self, factor: float):
        """
        Scale the canvas by a factor.
        """
        for node in self.all_nodes:
            node.x *= factor
            node.y *= factor
        self.scalar = factor

    def icicle_layout(self):
        """
        Make an icicle layout of the nodes.
        """
        nodes = self.root_nodes.copy()
        while len(nodes) > 0:
            node = nodes.pop()
            if not node.has_valid_placement():
                self._autoplace_node(node)
            nodes.extend(node.children)

        root = self.root_nodes[0]
        root.shift_down(-root.x, -root.y)

    def tree_layout(self):
        """
        Make a tree layout of the nodes.
        """
        for nodes in self.get_clashing_nodes():
            node, other_node = nodes
            self._autoplace_node(node)
            self._autoplace_node(other_node)

        root = self.root_nodes[0]
        root.shift_down(-root.x, -root.y)

    def render(
        self,
        ax,
        draw_edge_labels: bool = True,
        pretty_labels: bool = True,
        small_labels: bool = False,
        label_kws: dict = None,
        edge_kws: dict = None,
        adjust_axes: bool = True,
    ):
        """
        Render the canvas on a matplotlib axis.
        This requires that all nodes have a matplotlib patch associated with them.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to render the canvas on.
        draw_edge_labels : bool
            Whether to draw the id labels of the connecting edges.
        pretty_labels : bool
            Whether to use the pretty format for the edge labels. Including greek letters and arrows.
        small_labels : bool
            Whether to use the small format for the edge labels, without arrows.
        label_kws : dict
            The keyword arguments to pass to the edge label drawing function.
        edge_kws : dict
            The keyword arguments to pass to matplotlib.patches.ConnectionPatch.
        """
        edge_kws = {"color": "black"} if edge_kws is None else edge_kws
        label_kws = (
            dict(
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=10,
                bbox=dict(
                    boxstyle="square",
                    pad=0,
                    edgecolor="none",
                    alpha=1,
                    facecolor="white",
                ),
            )
            if label_kws is None
            else label_kws
        )
        for node in self.all_nodes:
            for child in node.children:
                ax.plot([node.x, child.x], [node.y, child.y], zorder=1, **edge_kws)
                if draw_edge_labels:
                    link = self.data["bonds"][(node.id, child.id)]
                    ax.text(
                        (node.x + child.x) / 2,
                        (node.y + child.y) / 2,
                        iupac.reverse_format_link(
                            link.id, pretty=pretty_labels, small=small_labels
                        ),
                        **label_kws,
                    )
            ax.imshow(
                node.icon,
                extent=(
                    node.x - self.scalar * 0.3,
                    node.x + self.scalar * 0.3,
                    node.y - self.scalar * 0.3,
                    node.y + self.scalar * 0.3,
                ),
                zorder=2,
            )
        if adjust_axes:
            ax.axis("off")
            ax.set_aspect("equal")
            _minx = min(node.x for node in self.all_nodes) - 1
            _maxx = max(node.x for node in self.all_nodes) + 1

            if abs(_minx) < abs(_maxx):
                _minx = -_maxx
            else:
                _maxx = -_minx

            ax.set_xlim(_minx, _maxx)
            ax.set_ylim(
                min(node.y for node in self.all_nodes) - 1,
                max(node.y for node in self.all_nodes) + 1,
            )

        # ax.set_xlim(
        #     min(node.x for node in self.all_nodes) - 1,
        #     max(node.x for node in self.all_nodes) + 1,
        # )
        # ax.set_ylim(
        #     min(node.y for node in self.all_nodes) - 1,
        #     max(node.y for node in self.all_nodes) + 1,
        # )

    def _autoplace_node(self, node):
        """
        Place a node at a free slot.
        """
        for i in range(-1, 2):
            for j in range(0, 1):
                if self.is_free(node.x + i, node.y + j):
                    node.shift_down(i, j)
                    return
        self._autoplace_node(node.parent)


if __name__ == "__main__":
    import glycosylator as gl

    # # man = gl.glycan("MAN")
    # # glc = gl.glycan("GLC")
    # # gulnac = gl.glycan("GulNAc")

    # # glycan = man % "14bb" + glc + glc + man + glc + glc
    # # glycan @ -4
    # # glycan % "16ab"
    # # glycan += man + glc + gulnac
    # # glycan % "12bb"
    # # glycan += man + gulnac + glc + man
    # # glycan = glycan @ -6 + glycan

    # # viewer = GlycanViewer2D(glycan)
    # # viewer.show()
    # USE_GLYCOWORK = False
    # mol = gl.Glycan.from_iupac(
    #     None,
    #     # "Gal(a1-2)Man(a1-3)[Gal(a1-2)Man(a1-3)]Glc(a1-4)[Gal(a1-2)Man(a1-3)][Gal(a1-2)Man(a1-3)]GlcNAc(b1-4)Man(b1-4)Glc(b1-",
    #     "Gal(a1-2)Man(a1-3)[Gal(a1-2)Man(a1-3)]Glc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Gc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-8)Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    # )
    # mol.set_root(1)
    # # c = Canvas.from_glycan(mol)
    # # c.tree_layout()
    # # fig, ax = plt.subplots()
    # # c.render(ax)
    # # fig.show()

    # v = GlycanViewer2D(mol)
    # v.show(axis="x")
    # plt.show()
    # # pass

    s = gl.Protein.load(
        "/Users/noahhk/GIT/glycosylator/docs/source/examples/files/protein_optimized.pkl"
    )
    # s.find_glycans()
    # s.save("/Users/noahhk/GIT/glycosylator/final_scaffold_superduper.2023-07-26 20:51:40.416918.pkl")
    # exit()
    # s = gl.Scaffold.load("/Users/noahhk/GIT/glycosylator/scaffold.3.glycan.pkl")
    # s.exclude_chain("B")

    v = ScaffoldViewer2D(s, figsize=(10, 10))

    v.draw_glycans(1)
    v.draw_glycosylation_sites(draw_legend=True)
    v.fig.tight_layout()
    plt.show()
