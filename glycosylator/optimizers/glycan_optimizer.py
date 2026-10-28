"""
Glycan-wise optimization of glycosylated scaffolds.

Instead of optimizing all torsions of all glycans at once in one large environment, each glycan
is optimized individually against a static environment (the scaffold and all other glycans in their current conformation).
The glycans are visited in order of how badly they clash, and the procedure can be repeated over multiple sweeps.
Because the current conformation of each glycan is always part of the search population and new conformations are only accepted
if they improve the evaluation, a glycan can never end up worse than it started.

The evaluation of a conformation is a soft clash penalty against the environment and within the glycan itself,
plus a (weak) reward for glycans that protrude away from the scaffold surface.

.. code-block:: python

    import glycosylator as gl

    protein = gl.Protein.from_pdb("my_protein.pdb")
    # ... attach glycans
    optimized = gl.optimize_glycans(protein)
"""

import warnings

import numpy as np
from scipy.optimize import differential_evolution
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import cKDTree

__all__ = ["optimize_glycans"]


def optimize_glycans(
    scaffold,
    glycans: list = None,
    include_n_ancestors: int = 2,
    sweeps: int = 2,
    coarse_first: bool = True,
    popsize: int = 24,
    maxiter: int = 40,
    clash_distance: float = 2.8,
    protrusion_weight: float = 0.05,
    seed: int = None,
    inplace: bool = False,
    verbose: bool = False,
):
    """
    Optimize the conformations of the glycans on a scaffold one glycan at a time.

    Parameters
    ----------
    scaffold : Scaffold
        The glycosylated scaffold.
    glycans : list
        The glycans to optimize. These can be Glycan objects or the scaffold atoms to which the glycans are attached.
        If None, all glycans are optimized. Glycans that are not optimized are part of the static environment.
    include_n_ancestors : int
        The number of torsions upstream of the glycan root connection (i.e. in the side chain of the anchor residue)
        that may also be rotated.
    sweeps : int
        The number of times to iterate over all glycans.
    coarse_first : bool
        If True, each glycan is first only optimized around its root connection (and anchor ancestors) and the
        glycan-internal torsions are only optimized if the glycan is still clashing afterwards.
    popsize : int
        The population size of the differential evolution for each glycan.
    maxiter : int
        The maximal number of generations of the differential evolution for each glycan.
    clash_distance : float
        Heavy atoms closer than this distance (in Angstrom) are penalized.
    protrusion_weight : float
        The weight of the reward for moving glycan atoms away from the scaffold surface (along the outward normal at the anchor).
        Set to 0 to only resolve clashes.
    seed : int
        A random seed for reproducibility.
    inplace : bool
        If True, the scaffold is modified in place, otherwise a copy is optimized.
    verbose : bool
        If True, print the evaluation of each glycan before and after optimization.

    Returns
    -------
    Scaffold
        The scaffold with optimized glycans.
    """
    if not inplace:
        scaffold = scaffold.copy()

    rng = np.random.default_rng(seed)

    all_glycans = scaffold.get_glycans()
    if glycans is None:
        to_optimize = list(all_glycans.items())
    else:
        to_optimize = []
        for g in glycans:
            for root, glycan in all_glycans.items():
                if g is glycan or g is root:
                    to_optimize.append((root, glycan))
                    break
            else:
                raise ValueError(f"Glycan {g} is not attached to the scaffold")

    units = [
        _GlycanUnit(scaffold, root, glycan, include_n_ancestors)
        for root, glycan in to_optimize
    ]
    units = [u for u in units if u.n_edges > 0]
    if len(units) == 0:
        return scaffold

    # every heavy atom that is not part of a unit is static
    moving = set().union(*(u.atoms for u in units))
    static_coords = np.array(
        [a.coord for a in scaffold.get_atoms() if a.element != "H" and a not in moving]
    )

    for unit in units:
        unit.setup_scaffold_normal(static_coords)

    def environment(unit):
        """The static environment of one unit: all static atoms + all other units in their current conformation"""
        others = [u.heavy_coords for u in units if u is not unit]
        return cKDTree(np.concatenate([static_coords, *others]))

    for sweep in range(sweeps):
        costs = np.array([u.clash_cost(environment(u), clash_distance) for u in units])
        order = np.argsort(costs)[::-1]
        for udx in order:
            unit = units[udx]
            if sweep > 0 and costs[udx] == 0:
                continue
            tree = environment(unit)

            if coarse_first:
                stages = [unit.root_edges, None]
            else:
                stages = [None]

            for edges in stages:
                if edges is not None and len(edges) == 0:
                    continue
                before = unit.evaluate(
                    np.zeros((1, unit.n_edges)), tree, clash_distance, protrusion_weight
                )[0]
                if (
                    edges is None
                    and coarse_first
                    and unit.clash_cost(tree, clash_distance) == 0
                ):
                    break
                angles = unit.search(
                    tree,
                    edges,
                    clash_distance,
                    protrusion_weight,
                    popsize,
                    maxiter,
                    rng,
                )
                after = unit.evaluate(
                    angles[None], tree, clash_distance, protrusion_weight
                )[0]
                if after < before:
                    unit.apply(angles)
                if verbose:
                    stage = "root" if edges is not None else "full"
                    print(
                        f"sweep {sweep} glycan {unit.glycan.id} ({stage}): {before:.3f} -> {min(before, after):.3f}"
                    )

    for unit in units:
        unit.write_back()

    return scaffold


class _GlycanUnit:
    """
    The local system of one glycan: the glycan atoms, the atoms of its anchor residue,
    the rotatable edges and the precomputed masks to rotate and evaluate conformations in batches.
    """

    def __init__(self, scaffold, root, glycan, include_n_ancestors):
        self.glycan = glycan
        self.root = root
        anchor = root.parent

        atoms = list(glycan.get_atoms()) + [
            a for a in anchor.get_atoms() if a not in glycan.get_atoms()
        ]
        self.atom_list = atoms
        self.atoms = set(atoms)
        index = {a: i for i, a in enumerate(atoms)}

        internal = glycan.get_residue_connections(direct_by="root")
        incoming = scaffold.get_glycan_root_connections(
            glycan=glycan, include_scaffold_ancestors=include_n_ancestors
        )

        graph = scaffold._AtomGraph
        edges, masks, root_edges = [], [], []
        for edge in list(incoming) + list(internal):
            a, b = edge
            if a not in index or b not in index:
                continue
            descendants = graph.get_descendants(a, b)
            if not descendants or not descendants <= self.atoms:
                warnings.warn(
                    f"Skipping edge {a.full_id} - {b.full_id} because it does not separate the glycan from the scaffold"
                )
                continue
            if (a, b) in edges:
                continue
            mask = np.zeros(len(atoms), dtype=bool)
            mask[[index[d] for d in descendants]] = True
            if edge in incoming:
                root_edges.append(len(edges))
            edges.append((a, b))
            masks.append(mask)

        # rotations are applied from the scaffold outwards
        order = np.argsort([-m.sum() for m in masks], kind="stable")
        self.edges = [edges[i] for i in order]
        self.masks = np.array([masks[i] for i in order]).reshape(len(edges), len(atoms))
        self.root_edges = [int(np.where(order == i)[0][0]) for i in root_edges]
        self.n_edges = len(self.edges)
        self.axes = np.array(
            [[index[a], index[b]] for a, b in self.edges], dtype=int
        ).reshape(-1, 2)

        self.coords = np.array([a.coord for a in atoms], dtype=np.float64)
        self.heavy = np.array([a.element != "H" for a in atoms])
        # static atoms of the anchor residue (e.g. the backbone) are bonded to their neighbours
        # in the environment and would only contribute a constant clash
        self.query = self.heavy & self.masks.any(axis=0)

        # intra-glycan pairs that can actually move relative to each other and that are
        # at least 4 bonds apart (closer pairs have their distances fixed by the geometry)
        heavy_idx = np.where(self.heavy)[0]
        n = len(atoms)
        adjacency = np.zeros((n, n), dtype=bool)
        for a, b in graph.edges(atoms):
            if a in index and b in index:
                adjacency[index[a], index[b]] = adjacency[index[b], index[a]] = True
        topological = shortest_path(adjacency, unweighted=True, directed=False)
        i, j = np.triu_indices(len(heavy_idx), k=1)
        i, j = heavy_idx[i], heavy_idx[j]
        relative = (self.masks[:, i] != self.masks[:, j]).any(axis=0)
        self.pairs = np.stack([i, j], axis=1)[relative & (topological[i, j] > 3)]

    def setup_scaffold_normal(self, static_coords):
        """The outward normal at the anchor, based on the static atoms around it"""
        anchor = self.root.coord
        local = static_coords[np.linalg.norm(static_coords - anchor, axis=1) < 12.0]
        if len(local) == 0:
            self.normal = np.zeros(3)
        else:
            normal = anchor - local.mean(axis=0)
            self.normal = normal / (np.linalg.norm(normal) or 1.0)
        self.origin = anchor.copy()
        self.glycan_heavy = np.array(
            [a in self.glycan.get_atoms() and a.element != "H" for a in self.atom_list]
        )

    @property
    def heavy_coords(self):
        return self.coords[self.heavy]

    def rotate(self, angles):
        """
        Rotate the local system for a batch of torsion angles.

        Parameters
        ----------
        angles : np.ndarray
            (n_candidates, n_edges) angles in radians

        Returns
        -------
        np.ndarray
            (n_candidates, n_atoms, 3) coordinates
        """
        X = np.repeat(self.coords[None], len(angles), axis=0)
        for k in range(self.n_edges):
            theta = angles[:, k]
            if np.all(np.abs(theta) < 1e-6):
                continue
            a, b = self.axes[k]
            mask = self.masks[k]
            ref = X[:, a][:, None]
            axis = X[:, b] - X[:, a]
            axis /= np.linalg.norm(axis, axis=1, keepdims=True)
            axis = axis[:, None]
            v = X[:, mask] - ref
            cos = np.cos(theta)[:, None, None]
            sin = np.sin(theta)[:, None, None]
            # Rodrigues' rotation formula
            X[:, mask] = (
                v * cos
                + np.cross(axis, v) * sin
                + axis * np.sum(axis * v, axis=2, keepdims=True) * (1 - cos)
                + ref
            )
        return X

    def _clash(self, X, tree, clash_distance):
        """Soft clash penalty for a batch of coordinates (n_candidates, n_atoms, 3)"""
        n = len(X)
        H = X[:, self.query]
        d, _ = tree.query(
            H.reshape(-1, 3), k=4, distance_upper_bound=clash_distance, workers=-1
        )
        d = np.where(np.isfinite(d), clash_distance - d, 0.0).reshape(n, -1)
        cost = np.square(d).sum(axis=1)

        if len(self.pairs):
            D = np.linalg.norm(X[:, self.pairs[:, 0]] - X[:, self.pairs[:, 1]], axis=2)
            cost += np.square(np.clip(clash_distance - D, 0, None)).sum(axis=1)
        return cost

    def evaluate(self, angles, tree, clash_distance, protrusion_weight):
        X = self.rotate(angles)
        cost = self._clash(X, tree, clash_distance)
        if protrusion_weight:
            protrusion = ((X[:, self.glycan_heavy] - self.origin) @ self.normal).mean(
                axis=1
            )
            cost -= protrusion_weight * protrusion
        return cost

    def clash_cost(self, tree, clash_distance):
        return float(self._clash(self.coords[None], tree, clash_distance)[0])

    def search(
        self, tree, edges, clash_distance, protrusion_weight, popsize, maxiter, rng
    ):
        """
        Search for better torsion angles using differential evolution.
        If edges are given, only those are optimized.
        """
        edges = np.arange(self.n_edges) if edges is None else np.asarray(edges)
        n = len(edges)

        def f(x):
            # scipy passes (n_dims, n_candidates) in vectorized mode
            angles = np.zeros((x.shape[1], self.n_edges))
            angles[:, edges] = x.T
            return self.evaluate(angles, tree, clash_distance, protrusion_weight)

        init = rng.uniform(-np.pi, np.pi, (max(popsize, 5), n))
        init[0] = 0.0  # the current conformation is always part of the population
        result = differential_evolution(
            f,
            [(-np.pi, np.pi)] * n,
            init=init,
            maxiter=maxiter,
            tol=1e-8,
            polish=False,
            vectorized=True,
            updating="deferred",
            seed=int(rng.integers(2**31)),
        )
        angles = np.zeros(self.n_edges)
        angles[edges] = result.x
        return angles

    def apply(self, angles):
        """Make the given torsions the new current conformation"""
        self.coords = self.rotate(angles[None])[0]

    def write_back(self):
        for atom, coord in zip(self.atom_list, self.coords):
            atom.coord = coord
