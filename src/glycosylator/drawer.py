import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Polygon, Rectangle

from utils import alphanum_sort


class Drawer:
    """Class to display glycans."""

    def __init__(self):
        self.Symbols = {}
        self.Symbols["GLC"] = ["circle", "blue"]
        self.Symbols["MAN"] = ["circle", "green"]
        self.Symbols["BMA"] = ["circle", "green"]
        self.Symbols["GAL"] = ["circle", "yellow"]
        self.Symbols["NAG"] = ["square", "blue"]
        self.Symbols["NGA"] = ["square", "yellow"]
        self.Symbols["FUC"] = ["triangle", "red"]
        self.Symbols["SIA"] = ["rhombus", "purple"]

        self.Colors = {}
        self.Colors["blue"] = np.array([0, 0, 204])
        self.Colors["green"] = np.array([0, 204, 0])
        self.Colors["yellow"] = np.array([255, 255, 102])
        self.Colors["red"] = np.array([204, 0, 0])
        self.Colors["purple"] = np.array([167, 13, 206])
        self.side = 0.25
        self.scaling = 5.0

    def get_shape(self, name, pos, direction=1, axis=0):
        """Returns a Matplotlib patch
        name: name of shape
        pos: vector defining where to draw object
        direction: direction for drawing glycan: 1 positive y; -1 negative y
        axis: axis along which to draw glycan. 0 for x and 1 for y.
        """
        pos = np.array(pos)
        if name in self.Symbols:
            shape, color = self.Symbols[name]
            color = self.Colors[color] / 255.0
            if shape == "circle":
                patch = Circle(pos, self.side)
                return patch, color
            if shape == "square":
                return Rectangle(pos - self.side, 2 * self.side, 2 * self.side), color
            if shape == "triangle":
                return (
                    Polygon(
                        [
                            (pos[0] - self.side, pos[1]),
                            (pos[0] + self.side, pos[1] - self.side),
                            (pos[0] + self.side, pos[1] + self.side),
                        ]
                    ),
                    color,
                )
            if shape == "rhombus":
                if axis:
                    angle = 135 if direction == -1 else -45.0
                else:
                    angle = -135 if direction == -1 else 45.0
                return Rectangle(pos, 2 * self.side, 2 * self.side, angle=angle), color

    def draw_protein_fragment(
        self,
        pos=[0.0, 0.0],
        length=10,
        protein_color=[0.7, 0.7, 0.7],
        sequon_color=[0.5, 0.5, 0.5],
        ax=None,
        axis=0,
    ):
        """
        Parameters:
            pos: position of sequon
            length: length of fragment
            color: color of protein
            axis: axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
        shape = [0, 0]
        shape[axis] = length / self.scaling
        shape[np.mod(axis + 1, 2)] = 2 * self.side
        p_fragment = np.array([0.0, 0.0])
        p_fragment[axis] = -shape[axis] / 2
        p_fragment[np.mod(axis + 1, 2)] = -(self.side * 3.5)
        protein = Rectangle(p_fragment - self.side, shape[0], shape[1])
        patches = []
        colors = []
        patches.append(protein)
        colors.append(protein_color)
        new_pos = np.array([0.0, 0.0])
        new_pos[axis] = pos[axis]
        new_pos[np.mod(axis + 1, 2)] = pos[np.mod(axis + 1, 2)] - (self.side * 3.5)
        patches.append(Rectangle(new_pos - self.side, 2 * self.side, 2 * self.side))
        colors.append(sequon_color)

        p = PatchCollection(patches, zorder=3)
        p.set_edgecolors([0.2, 0.2, 0.2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)

        return ax

    def draw_glycoprotein(
        self,
        length,
        start_resid,
        sequons,
        pos=[0.0, 0.0],
        protein_color=[0.7, 0.7, 0.7],
        sequon_color={},
        ax=None,
        axis=0,
        trees=None,
        names=None,
    ):
        """
        Parameters:
            length: length of the protein
            start_resid: Index of frist residue
            sequons: position of sequons
            pos: starting position of sequence
            protein_color: color of protein (list)
            sequon_color: dictionary of colors for sequons, with sequon id as key and color code as value. Default color [.5, .5, .5]
            axis: axis along which to draw glycan. 0 for x and 1 for y.
            trees: dictionary with sequon as key and list of root and glycan graph as value
            names: dictionary with glycan nodes as keys and resnames as values
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
        shape = [0, 0]
        shape[axis] = length / self.scaling
        shape[np.mod(axis + 1, 2)] = 2 * self.side
        protein = Rectangle(np.array(pos) - self.side, shape[0], shape[1])
        patches = []
        colors = []
        patches.append(protein)
        colors.append(protein_color)
        direction = -1
        sequons = alphanum_sort(sequons)
        for seq in sequons:
            s, c, rr, i = seq.split(",")
            if seq in sequon_color:
                color = sequon_color[seq]
            else:
                color = [0.5, 0.5, 0.5]
            r = rr
            if i:
                r = str(float(r) + ord(i.upper()) - 64)
            r = (float(r) - start_resid) / self.scaling
            new_pos = np.array([0.0, 0.0])
            new_pos[axis] = pos[axis] + r
            new_pos[np.mod(axis + 1, 2)] = pos[np.mod(axis + 1, 2)]
            patches.append(Rectangle(new_pos - self.side, 2 * self.side, 2 * self.side))
            colors.append(color)

            new_pos[np.mod(axis + 1, 2)] = direction * (self.side * 3.5)
            if axis:
                h_align = "right" if direction == -1 else "left"
                v_align = "center"
                rotation = "horizontal"
            else:
                h_align = "center"
                v_align = "top" if direction == -1 else "bottom"
                rotation = "vertical"
            ax.text(
                new_pos[0],
                new_pos[1],
                rr + i,
                verticalalignment=v_align,
                horizontalalignment=h_align,
                rotation=rotation,
            )
            direction *= -1

            if seq in trees:
                r_pos = [0, 0]
                r_pos[axis] = r
                r_pos[np.mod(axis + 1, 2)] = direction * (self.side * 3.5)
                root, tree = trees[seq]
                self.draw_tree(
                    tree, root, names, r_pos, direction=direction, ax=ax, axis=axis
                )

        p = PatchCollection(patches, zorder=3)
        p.set_edgecolors([0.2, 0.2, 0.2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)

        return ax

    def draw_all_trees(self, trees, start_resid, names, ax=None, axis=0):
        """Draws all glycans in a protein
        Parameters:
            trees: dictionary with sequon as key and list of root and glycan graph as value
            start_resid: first residue number of protein sequence
            names: dictionary with glycan nodes as keys and resnames as values
            ax: axis handle
            axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()
        direction = 1
        keys = alphanum_sort(trees.keys())
        for k in keys:
            s, c, r, i = k.split(",")
            if i:
                r += ord(i.upper()) - 64
            r_pos = [0, 0]
            r_pos[axis] = (float(r) - start_resid) / self.scaling
            r_pos[np.mod(axis + 1, 2)] = direction * (self.side * 3.5)
            root, tree = trees[k]
            self.draw_tree(
                tree, root, names, r_pos, direction=direction, ax=ax, axis=axis
            )
            direction *= -1

        plt.axis("equal")
        plt.axis("off")
        fig = plt.gcf()
        fig.tight_layout()
        return ax

    def draw_tree(
        self, tree, root, names, root_pos=[0, 0], direction=1, ax=None, axis=0
    ):
        """Draws a tree representation of a glycan
        Parameters:
            tree: netwrokX graph representation of glycan
            root: name of root node
            names: dictionary with nodes as keys and resnames as values
            root_pos: position of the root node (default [0,0])
            direction: direction for drawing glycan: 1 positive y; -1 negative y
            axis: axis along which to draw glycan. 0 for x and 1 for y.
            ax: axes handle
         Returns:
            ax: axis of plot
        """
        if not ax:
            fig, ax = plt.subplots()

        pos = self.compute_positions(tree, root, root_pos, direction, axis)
        # draw edges
        for edge in tree.edges():
            x1, y1 = pos[edge[0]]
            x2, y2 = pos[edge[1]]
            ax.add_line(mlines.Line2D([x1, x2], [y1, y2], lw=2.0, c=[0.2, 0.2, 0.2]))
        # draw nodes
        patches = []
        colors = []
        for node in tree.nodes():
            patch, color = self.get_shape(names[node], pos[node], direction, axis)
            patches.append(patch)
            colors.append(color)
        p = PatchCollection(patches, zorder=3)
        p.set_edgecolors([0.2, 0.2, 0.2])
        p.set_linewidth(2)
        p.set_facecolors(colors)
        ax.add_collection(p)

        return ax

    def tree_to_text(self, tree, root, names, glyPAC="", visited=[]):
        """Returns UIPAC like string of a glycan.
        []: branch
        (): connectivity
        Example: Man3 == NAG(14bb)NAG(14bb)BMA[(16ab)MAN](13ab)MAN
        Parameters:
            tree: directed graph of a rooted tree
            root: id of root vertex
        Returns:
            glyPAC: string representing the glycan
        """
        if not glyPAC:
            glyPAC = names[root]
        while 1:
            visited.append(root)
            successors = [n for n in tree.neighbors(root) if n not in visited]
            if not successors:
                break
            branches = len(successors)
            for s in successors:
                patch = tree[root][s]["patch"]
                if branches > 1:
                    glyPAC += "[(" + patch + ")" + names[s]
                    glyPAC = self.tree_to_text(tree, s, names, glyPAC, visited)
                    glyPAC += "]"
                    branches -= 1
                else:
                    glyPAC += "(" + patch + ")" + names[s]
            root = s
        return glyPAC

    def compute_positions(self, tree, root, root_pos=[0, 0], direction=1, axis=0):
        r"""Computes positions of vertices for drawing a tree
        1-2--3-8
           \-4--5
              \-6-7
        Parameters:
            tree: directed graph of a rooted tree
            root: id of root vertex
            root_pos: position of the root vertex (default: 0,0)
            direction: -1/1 direction to draw tree
            axis: axis along which to draw glycan. 0 for x and 1 for y.
        Returns:
            pos: dictionary with nodes as keys and positions as values
        """
        nodes = list(nx.dfs_preorder_nodes(tree, root))
        successors = nx.dfs_successors(tree, root)
        nodes.reverse()
        branches = {}
        for n in nodes:
            if n in successors:
                branches[n] = len(successors[n])
                for s in successors[n]:
                    if branches[s] > 1:
                        branches[n] += branches[s] - 1
            else:
                branches[n] = 1

        nodes.reverse()
        pos = {}

        for n in nodes:
            if not n in pos:
                pos[n] = root_pos

            if n in successors:
                nbr_branches = []
                for s in successors[n]:
                    nbr_branches.append(branches[s])
                c = 0
                for s in np.argsort(nbr_branches):
                    x = pos[n]
                    x_new = [0, 0]
                    x_new[axis] = x[axis] + c
                    x_new[np.mod(axis + 1, 2)] = x[np.mod(axis + 1, 2)] + direction
                    pos[successors[n][s]] = x_new
                    c += nbr_branches[s]

        return pos
