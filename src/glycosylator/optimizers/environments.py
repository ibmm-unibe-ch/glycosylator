from copy import deepcopy
import gymnasium as gym
import numpy as np
from scipy.spatial.transform import Rotation

import glycosylator.utils.visual as visual


class Rotatron(gym.Env):
    """
    The base class for rotation environments

    Parameters
    ----------
    graph : glycosylator.graphs.Graph
        The graph to optimize
    rotateable_edges : list, optional
        A list of edges to rotate around
    """

    def __init__(self, graph, rotateable_edges) -> None:
        super().__init__()
        self.graph = graph
        self.rotateable_edges = [i for i in rotateable_edges if i in graph.edges]

        self._make_coord_array()
        self._make_descendant_masks()
        self._make_edge_ref_masks()

        self._effector_mask = np.ones(len(self._nodes), dtype=bool)
        self.__orig_coords = deepcopy(self._coords)

        self.__renderer__ = None

    @property
    def effector_coords(self):
        """
        Get the effector coordinates
        """
        return self._coords[self._effector_mask]

    @property
    def effector_nodes(self):
        """
        Get the effector nodes
        """
        return np.array(self._nodes)[self._effector_mask]

    @property
    def state(self):
        """
        Get the current state of the environment
        """
        return self._coords

    @state.setter
    def state(self, state):
        """
        Set the current state of the environment
        """
        self._coords = deepcopy(state)

    def set_state(self, state):
        """
        Set the current state of the environment
        """
        self.state = state

    def get_state(self):
        """
        Get the current state of the environment
        """
        return self.state

    def reset(self):
        """
        Reset the environment
        """
        self._coords = deepcopy(self.__orig_coords)
        self.apply_to_graph()
        return self.state

    def render(self, **kwargs):
        """
        Render the environment
        """
        raise NotImplementedError

    def step(self, action):
        """
        Take a step in the environment
        """
        raise NotImplementedError

    def apply_to_graph(self):
        """
        Apply the current state to the graph
        """
        for i, node in enumerate(self._nodes):
            node.coord = deepcopy(self._coords[i])

    def mask_effector_coords(self, mask):
        """
        Mask the effector coordinates for which the reward is calculated.
        By default, all coordinates are used.

        Parameters
        ----------
        mask : np.ndarray
            The mask to apply
        """
        self._effector_mask = mask

    def rotate(self, bond: int, angle: float):
        """
        Rotate a coordinate array around an axis by a specified angle

        Parameters
        ----------
        bond : int
            The bond to rotate around (as sampled from the action space)
        angle : float
            The angle to rotate by in radians
        """

        desc_mask = self._descendant_masks[bond]
        coords = self._coords[desc_mask]

        # for c in coords:
        #     self._v.draw_point(None, c, color="red", opacity=0.3, showlegend=False)

        # Get the axis to rotate around
        edge_masks = self._edge_ref_masks[bond]
        a = self._coords[edge_masks[0]][0]
        b = self._coords[edge_masks[1]][0]
        axis = b - a
        axis = axis / np.linalg.norm(axis)

        # get the reference coordinate for the bond
        ref_coord = b

        # self._v.draw_point("anchor", ref_coord, color="limegreen")

        # translate the coordinates so that the reference coordinate is at the origin
        coords -= ref_coord

        # Get the rotation matrix
        R = Rotation.from_rotvec(angle * axis)

        # Rotate the coordinates
        rotated_coords = R.apply(coords)

        # Translate the coordinates back to the original position
        rotated_coords += ref_coord

        # Update the coordinate array
        self._coords[desc_mask] = rotated_coords

        # for c in self._coords[desc_mask]:
        #     self._v.draw_point(None, c, color="purple", opacity=0.5, showlegend=False)

        return rotated_coords

    def blank(self):
        """
        Return a blank action
        """
        raise NotImplementedError

    def compute_reward(self, coords=None):
        """
        Compute the reward of the given or current coordinate array
        """
        if coords is None:
            coords = self.effector_coords

        # Compute the inter-residue distances
        dists = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

        # Mask the diagonal
        np.fill_diagonal(dists, np.nan)
        dist = np.nanmin(dists)

        # now also include a penalty for any atoms that are too close
        # to each other
        _dists = dists[np.where(~np.isnan(dists))]
        penalty = np.sum(_dists < 0.85)
        penalty *= 100

        reward = dist - penalty
        return reward

    def _make_coord_array(self):
        """
        Make the coordinate array
        """
        self._nodes = list(self.graph.nodes)
        self._coords = np.array([i.coord for i in self._nodes])

    def _make_descendant_masks(self):
        """
        Map the descendant nodes for each rotateable edge
        """
        self._descendant_masks = {
            idx: np.array([j in self.graph.get_descendants(*edge) for j in self._nodes])
            for idx, edge in enumerate(self.rotateable_edges)
        }

    def _make_edge_ref_masks(self):
        """
        Map the coordinates of one end of each rotateable edge from the _coords array
        """
        self._edge_ref_masks = {
            idx: (
                np.array([j == node1 for j in self._nodes]),
                np.array([j == node2 for j in self._nodes]),
            )
            for idx, (node1, node2) in enumerate(self.rotateable_edges)
        }


class SingleBondRotatron(Rotatron):
    """
    Optimize the geometry of a structure by rotating around single bonds successively
    """

    def __init__(self, graph, connections):
        super().__init__(graph, connections)

        self._bond_space = gym.spaces.Discrete(len(self.rotateable_edges))
        self._angle_space = gym.spaces.Box(low=-np.pi, high=np.pi, shape=(1,))
        self.action_space = gym.spaces.Tuple((self._bond_space, self._angle_space))

    def step(self, action):
        """
        Take a step in the environment
        """
        bond, angle = action
        self.rotate(bond, angle)
        reward = self.compute_reward()
        done = False
        info = {}
        return self.state, reward, done, info

    def blank(self):
        """
        Return a blank action
        """
        return (0, 0)

    def render(self, **kwargs):
        """
        Render the environment
        """
        print(f"Current reward: {self.compute_reward()}")


__all__ = ["SingleBondRotatron"]


if __name__ == "__main__":
    from glycosylator import Molecule
    import glycosylator as gl
    from time import time

    mol = Molecule.from_pdb("/Users/noahhk/GIT/glycosylator/support/examples/man9.pdb")
    mol.infer_bonds(restrict_residues=False)

    # for c in mol.get_residue_connections():

    #     angle = np.random.randint(120, 180)
    #     mol.rotate_around_bond(*c, angle)

    # connections = mol.get_residue_connections()

    # g = mol.make_residue_graph(detailed=True)
    # v = visual.MoleculeViewer3D(g)

    # env = SingleBondRotatron(g, connections)
    # env._v = v

    # bond, angle = env.action_space.sample()
    # for i in range(10):
    #     env.step((bond, angle))

    # # v.draw_point("root", mol.get_atom(1).coord)

    # # for c in env.rotateable_edges:
    # #     v.draw_vector(None, c[0].coord, 1.2 * (c[1].coord - c[0].coord), color="red")

    # # t1 = time()
    # # # bond, _ = env.action_space.sample()
    # # # _bond = env.bonds[bond]
    # # # v.draw_vector(
    # # #     None, _bond[0].coord, 1.2 * (_bond[1].coord - _bond[0].coord), color="limegreen"
    # # # )

    # # start_reward = env.compute_reward()
    # # print(f"Starting reward: {start_reward}")

    # # colors = ["yellow", "orange", "red", "magenta", "blue", "cyan", "green"]
    # # angles = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]
    # # for run in range(3):

    # #     global_best = -1000
    # #     global_best_coords = None
    # #     env.reset()

    # #     for i in range(3):
    # #         c, reward, *_ = env.step((3, angles[run]))

    # #         if reward > global_best:
    # #             global_best = reward
    # #             global_best_coords = c

    # #         env.apply_to_graph()
    # #         v.draw_edges(env.graph.edges, color=colors[run])
    # #         v.draw_point(
    # #             None,
    # #             env.rotateable_edges[3][0].coord,
    # #             color="limegreen",
    # #         )

    # #     v.draw_point(
    # #         "root", mol.get_residue(2).coord, color="magenta", showlegend=False
    # #     )

    # #     print(f"=== run {run} ===")
    # #     print(f"Best: {global_best}")

    # # env.set_state(global_best_coords)
    # # env.apply_to_graph()
    # # v.draw_edges(env.graph.edges, color="magenta")

    # # t2 = time()
    # # print(f"Time: {t2 - t1}")

    # v.show()
    # env.reset()

    connections = mol.get_residue_connections()

    env = SingleBondRotatron(mol.make_residue_graph(detailed=True), connections)
    env._v = gl.utils.visual.MoleculeViewer3D(env.graph)

    initial_coords = env._coords.copy()
    nodes = list(env.graph.nodes)

    initial_distances = []
    for node in nodes:
        neighbors = env.graph.get_neighbors(node)
        dists = np.array(
            [np.linalg.norm(node.coord - neighbor.coord) for neighbor in neighbors]
        )
        initial_distances.append(dists)

    v = gl.utils.visual.MoleculeViewer3D(mol.make_residue_graph(detailed=True))
    for run in range(5):

        for i in range(30):
            action = env.action_space.sample()
            descendant_mask = env._descendant_masks[action[0]]
            descendant_coords = env._coords[descendant_mask].copy()
            ancestor_coords = env._coords[~descendant_mask].copy()

            env.step(action)

            assert not np.allclose(descendant_coords, env._coords[descendant_mask])
            assert np.allclose(ancestor_coords, env._coords[~descendant_mask])

            new_distances = []
            for node in nodes:
                dists = []
                neighbors = env.graph.get_neighbors(node)
                for neighbor in neighbors:
                    dists.append(np.linalg.norm(node.coord - neighbor.coord))
                new_distances.append(np.asarray(dists))

            for i in range(len(nodes)):
                assert np.all(np.abs(initial_distances[i] - new_distances[i]) < 1e-3)

        env.apply_to_graph()
        v.draw_edges(env.graph.edges, color="cyan", linewidth=2)

        env.reset()
        assert np.allclose(initial_coords, env.state)

    v.show()
