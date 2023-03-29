"""
The basic Class for Molecular Graphs
"""

from abc import abstractmethod

import Bio.PDB as bio
import networkx as nx
import numpy as np
from scipy.spatial.transform import Rotation


class BaseGraph(nx.Graph):
    """
    The basic class for molecular graphs
    """

    def __init__(self, id, bonds: list):
        if len(bonds) == 0:
            raise ValueError("Cannot create a graph with no bonds")
        super().__init__(bonds)
        self.id = id
        self._structure = None
        self._neighborhood = None

    @property
    def structure(self):
        """
        Returns the underlying `bio.PDB.Structure` object
        """
        if not self._structure:
            self._structure = self._get_structure()
        return self._structure

    @property
    def residues(self):
        """
        Returns the residues in the molecule
        """
        return list(self.structure.get_residues())

    @property
    def atoms(self):
        """
        Returns the atoms in the molecule
        """
        return list(self.structure.get_atoms())

    @property
    def bonds(self):
        """
        Returns the bonds in the molecule
        """
        return list(self.edges)

    def to_pdb(self, filename: str):
        """
        Save to a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        structure = bio.Structure.Structure(self.id)
        model = bio.Model.Model(0)
        chain = bio.Chain.Chain("A")
        structure.add(model)
        model.add(chain)

        for atom in self._structure:
            chain.add(atom)

        io = bio.PDBIO()
        io.set_structure(structure)
        io.save(filename)

    @abstractmethod
    def get_neighbors(self, node, n: int = 1, mode="upto"):
        """
        Get the neighbors of a node

        Parameters
        ----------
        node
            The target node

        n : int, optional
            The number of edges to separate the node from its neighbors.

        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` edges
            - "exact": get all neighbors exactly `n` edges away

        Returns
        -------
        set
            The neighbors of the node
        """
        raise NotImplementedError

    def get_descendants(self, node_1, node_2):
        """
        Get all descendant nodes that come after a specific edge
        defined in the direction from node1 to node2 (i.e. get all
        nodes that come after node2). This method is directed
        in contrast to the `get_neighbors()` method, which will get all neighboring
        nodes of an anchor node irrespective of direction.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge

        Returns
        -------
        set
            The descendants of the node

        Examples
        --------
        In case of this graph
        ```
        A---B---C---D---E
             \\
              F---H
              |
              G
        ```
        >>> graph.get_descendants("B", "C")
        {"D", "E"}
        >>> graph.get_descendants("B", "F")
        {"H", "G"}
        >>> graph.get_descendants("B", "A")
        set() # because in this direction there are no other nodes
        """
        if node_1 == node_2:
            raise ValueError("Cannot get descendants if only one node is given (no direction)!")

        neighbors = self.get_neighbors(node_2)
        neighbors.remove(node_1)
        if len(neighbors) == 0:
            return neighbors

        _new_neighbors = neighbors.copy()
        _seen = set((node_1, node_2))

        while len(_new_neighbors) > 0:

            neighbor = _new_neighbors.pop()
            descendants = self.get_neighbors(neighbor)

            descendants -= _seen
            _seen.add(neighbor)
            neighbors.add(neighbor)

            if len(descendants) == 0:
                continue
            _new_neighbors.update(descendants)

        return neighbors

    def rotate_around_edge(self, node_1, node_2, angle: float, descendants_only: bool = False):
        """
        Rotate descending nodes around a specific edge by a given angle.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge around which to rotate.
        angle: float
            The angle to rotate by, in radians.
        descendants_only: bool, optional
            Whether to only rotate the descending nodes, by default False, in which case the entire graph
            will be rotated.
        """

        # get the node coordinates as a dictionary
        # # we do this in order to update the node attributes later...
        # node_dict = nx.get_node_attributes(self, 'coord')

        if node_1 not in self.nodes or node_2 not in self.nodes:
            raise ValueError("One or more nodes not in graph!")

        # we need to get a reference node index to normalise the rotated
        # coordinates to the original coordinate system
        indices = list(self.nodes)
        idx_1 = indices.index(node_1)

        # define the axis of rotation as the cross product of the edge's vectors
        edge_vector = node_2.coord - node_1.coord
        edge_vector /= np.linalg.norm(edge_vector)

        # create the rotation matrix
        r = Rotation.from_rotvec(angle * edge_vector)

        # create a numpy array of the node coordinates
        if descendants_only:
            nodes = {i: i.coord for i in self.get_descendants(node_1, node_2)}
            nodes[node_2] = node_2.coord
        else:
            nodes = nx.get_node_attributes(self, 'coord')

        node_coords = np.array(tuple(nodes.values()))

        indices = list(nodes.keys())
        idx_2 = indices.index(node_2)

        # apply the rotation matrix to the node coordinates
        node_coords_rotated = r.apply(node_coords)

        # now adjust for the translational shift around the axis
        _diff = node_coords_rotated[idx_2] - node_coords[idx_2]
        node_coords_rotated -= _diff

        # update the node coordinates in the graph
        for i, node in enumerate(nodes):
            node.coord = node_coords_rotated[i]

        new_coords = {i: i.coord for i in self.nodes}

        # set the node attributes in the graph
        nx.set_node_attributes(self, new_coords, 'coord')

    def _get_structure(self):
        """
        Get the underlying `bio.PDB.Structure` object
        """
        if not hasattr(list(self.nodes)[0], "get_parent"):
            raise ValueError("Nodes are not Biopython entities with linked parents!")
        structure = list(self.nodes)[0].get_parent()
        while not isinstance(structure, bio.Structure.Structure):
            structure = structure.get_parent()
        return structure
