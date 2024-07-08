from buildamol.optimizers import Rotatron, ConstraintRotatron
import numpy as np
from scipy.spatial.distance import cdist


def _compute_glycan_scaffold_distances(
    each_glycan_nodes: list,
    each_glycan_scaffold_vicinity_nodes: list,
    state: np.ndarray,
):
    _glycan_scaffold_distances = np.zeros(len(each_glycan_nodes))
    _glycan_glycan_distances = np.zeros(
        (len(each_glycan_nodes), len(each_glycan_nodes))
    )
    for i, glycan_nodes in enumerate(each_glycan_nodes):
        if len(glycan_nodes) == 0:
            continue
        glycan_coords = state[glycan_nodes]
        scaffold_coords = state[each_glycan_scaffold_vicinity_nodes[i]]

        d = cdist(glycan_coords, scaffold_coords)
        _glycan_scaffold_distances[i] = 1 / d.mean()

        for j, glycan_nodes_j in enumerate(each_glycan_nodes):
            if i == j:
                _glycan_glycan_distances[i, j] = -1
                continue
            if _glycan_glycan_distances[j, i] > 0:
                _glycan_glycan_distances[i, j] = _glycan_glycan_distances[j, i]
                continue
            glycan_coords_j = state[glycan_nodes_j]
            _glycan_glycan_distances[i, j] = (
                1 / cdist(glycan_coords, glycan_coords_j).min()
            )

    a = np.mean(_glycan_scaffold_distances)
    b = np.apply_along_axis(_along_axis_min, 1, _glycan_glycan_distances).mean()
    return a, b


def _along_axis_min(x):
    mask = x > -1
    if mask.any():
        return x[mask].min()
    return 0


class ScaffoldRotatron(ConstraintRotatron):
    """
    A Rotatron that is specific to glycosylated scaffolds. The heuristic is designed to increase the distance of glycan residues from the scaffold's residues.

    Parameters
    ----------
    base_rotatron: Rotatron
        The Rotatron to use for the base optimization. This needs to be a fully initialized Rotatron.
    """

    def __init__(self, base_rotatron: "Rotatron", *args, **kwargs):
        super().__init__(base_rotatron, self.constraint, **kwargs)
        self.mol = base_rotatron.graph._molecule
        self.graph = base_rotatron.graph
        self.hyperparameters = self.rotatron.hyperparameters
        self._identify_nodes()
        self._compute_glycan_scaffold_distances = (
            self._normal_compute_glycan_scaffold_distances
        )
        self.n_edges = len(self.rotatron.rotatable_edges)
        self.rotatable_edges = self.rotatron.rotatable_edges

    def _identify_nodes(self):
        """
        Identify which nodes belong to the scaffold and which belong to the glycans.
        """
        _glycan_nodes = []
        _each_glycan_nodes = []
        for glycan in self.mol.get_glycans().values():
            _each_glycan_nodes.append([])
            for residue in glycan.get_residues():
                _glycan_nodes.append(residue)
                _glycan_nodes.extend(residue.child_list)
                _each_glycan_nodes[-1].extend(residue.child_list)

        self.scaffold_nodes = []
        self.glycan_nodes = []
        for node in self.graph.nodes:
            if node in _glycan_nodes:
                self.glycan_nodes.append(self.rotatron.node_dict[node])
            else:
                self.scaffold_nodes.append(self.rotatron.node_dict[node])

        _each_glycan_scaffold_vicinity_nodes = []
        scaff_coords = self.rotatron.state[self.scaffold_nodes]
        for i, glycan in enumerate(self.mol.get_glycans().values()):

            x = np.array(
                [
                    self.rotatron.node_dict[node]
                    for node in _each_glycan_nodes[i]
                    if node in self.rotatron.node_dict
                ]
            )
            if len(x) == 0:
                _each_glycan_scaffold_vicinity_nodes.append([])
                _each_glycan_nodes[i].clear()
                continue
            glycan_mean = glycan.root_residue.coord
            _each_glycan_scaffold_vicinity_nodes.append(
                np.where(cdist([glycan_mean], scaff_coords) < 8.0)[1]
            )
            _each_glycan_nodes[i] = x

        self.scaffold_nodes = np.array(self.scaffold_nodes)
        self.glycan_nodes = np.array(self.glycan_nodes)
        self.each_glycan_nodes = [i for i in _each_glycan_nodes if len(i) > 0]
        self.each_glycan_scaffold_vicinity_nodes = [
            i for i in _each_glycan_scaffold_vicinity_nodes if len(i) > 0
        ]
        pass

    def constraint(self, rotatron, state, **kwargs):
        scaffold_dist, glycan_dist = self._compute_glycan_scaffold_distances(state)
        return 5 * scaffold_dist + 10 * glycan_dist

    def _normal_compute_glycan_scaffold_distances(self, state: np.ndarray):
        """
        Compute the distances of glycans to the scaffold...
        """
        return _compute_glycan_scaffold_distances(
            self.each_glycan_nodes,
            self.each_glycan_scaffold_vicinity_nodes,
            state,
        )


if __name__ == "__main__":
    import glycosylator as gl

    rot = gl.utils.load_pickle(
        "/Users/noahhk/GIT/glycosylator/test.scaf.optimize.distrot-full.pkl"
    )

    rot = ScaffoldRotatron(rot, numba=False)
    e = rot.eval(rot.rotatron.state)
    print("allright")
    final = gl.optimizers.optimize(rot.mol.copy(), rot)
    final.to_pdb("scaf.opt-newScafRot.pdb")
    final.show(residue_graph=True)
    pass
