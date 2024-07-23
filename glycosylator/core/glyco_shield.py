"""
The `GlycoShield` class offers a convenient interface to perform a quick
simulation of the shielding effect of glycans on a scaffold molecule, and can be
used to obtain numeric data or visualizations of the scaffold at different glycan conformations.
"""

from copy import deepcopy

from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
import pandas as pd

import glycosylator.core as core
import buildamol.structural as structural
import glycosylator.utils as utils

__all__ = ["quickshield", "GlycoShield"]


def quickshield(
    scaffold,
    repeats: int = 1,
    angle_step: int = 50,
    save_conformations_to: str = None,
    verbose: bool = False,
) -> "GlycoShield":
    """
    A convenience function to perform a quick simulation of the glycan shielding on a scaffold.
    This function will set up a GlycoShield and perform a standard simulation. Use a manual procedure
    to retain more control.

    Parameters
    ----------
    scaffold
        The scaffold molecule with attached glycans
    repeats: int
        The number of times to repeat the edge sampling and simulation process for each glycan.
    angle_step: int
        The step size in degrees to rotate each sampled rotatable edge during the simulation.
        The smaller this is the more fine-grained the simulation will be, but the longer it will take.
    save_conformations_to: str
        If a path for a directory is provided, the conformations will be saved to that path as PDB files (one file per glycan, containing models for each accepted conformation).
    verbose: bool
        If True, a progress bar will be shown during the simulation.

    Returns
    -------
    GlycoShield
        The GlycoShield object that stores the analyses results such as exposure dataframe (`df` attribute)
        or any available visualizations.
    """
    shield = GlycoShield(scaffold)
    shield.simulate(
        repeats=repeats,
        angle_step=angle_step,
        save_conformations_to=save_conformations_to,
        verbose=verbose,
    )
    return shield


class ExposureData:
    """
    A class to house the data of per-residue exposure values for a scaffold molecule.
    """

    def __init__(self):
        self.df = None
        self.scaffold = None

    def plot(self, backend: str = "matplotlib", cmap: str = "viridis"):
        """
        Plot the exposure for each residue using Matplotlib or Plotly

        Parameters
        ----------
        backend: str
            The plotting backend to use. Either "matplotlib" or "plotly".
        cmap: str
            The colormap to use when coloring the bars. This must be an acceptable
            colormap name for either matplotlib or plotly (depending on the backend).

        Returns
        -------
        fig
            The figure, either a Matplotlib figure or a Plotly figure.
        """
        if self.df is None:
            raise ValueError("No exposure data available. Run simulate() first.")

        if backend == "matplotlib":
            return self.matplotlib(cmap)

        elif backend == "plotly":
            return self.plotly(cmap)

        else:
            raise ValueError(f"Unknown plotting backend: {backend}")

    def plotly(self, cmap: str = "viridis"):
        """
        Plot the exposure for each residue using Plotly

        Parameters
        ----------
        cmap: str
            The colormap to use when coloring the bars. This must be an acceptable
            colormap name for Plotly.

        Returns
        -------
        fig
            The Plotly figure.
        """
        import plotly.express as px
        import plotly.graph_objs as go

        fig = px.bar(
            self.df,
            x="serial",
            y="exposure",
            facet_row="chain",
            color="exposure",
            color_continuous_scale=cmap,
        )

        return fig

    def matplotlib(self, cmap: str = "viridis"):
        """
        Plot the exposure for each residue using Matplotlib

        Parameters
        ----------
        cmap: str
            The colormap to use when coloring the bars. This must be an acceptable
            colormap name for Matplotlib.

        Returns
        -------
        fig
            The Matplotlib figure.
        """
        global plt
        import matplotlib.pyplot as plt
        import seaborn as sns

        _, get_color = _prepare_colormap(self.df, cmap)

        sns.set_style("ticks")
        fig, axes = plt.subplots(len(self.df["chain"].unique()), 1)
        for i, chain in enumerate(sorted(self.df["chain"].unique())):
            ax = axes[i]
            chain_df = self.df[self.df["chain"] == chain]
            ax.bar(
                chain_df["serial"],
                chain_df["exposure"],
                color=[get_color(e, hex=True) for e in chain_df["exposure"].values],
            )
            ax.set_title(f"Chain {chain}")
            ax.set_xticks(chain_df["serial"][::10])

        fig.suptitle(f"Exposure of Residues in {self.scaffold.id}")
        fig.supylabel("Exposure")
        fig.supxlabel("Residue Serial")
        sns.despine()
        plt.tight_layout()
        return fig

    def show(self, backend: str = "matplotlib", cmap: str = "viridis"):
        """
        Show the exposure plot using Matplotlib or Plotly

        Parameters
        ----------
        backend: str
            The plotting backend to use. Either "matplotlib" or "plotly".
        cmap: str
            The colormap to use when coloring the bars. This must be an acceptable
            colormap name for either matplotlib or plotly (depending on the backend).
        """
        fig = self.plot(backend, cmap)
        if backend == "matplotlib":
            plt.show(fig)
        elif backend == "plotly":
            fig.show()

    def save(self, path: str, backend: str = "matplotlib", cmap: str = "viridis"):
        """
        Save the exposure plot to a file using Matplotlib or Plotly

        Parameters
        ----------
        path: str
            The path to save the plot to.
        backend: str
            The plotting backend to use. Either "matplotlib" or "plotly".
        cmap: str
            The colormap to use when coloring the bars. This must be an acceptable
            colormap name for either matplotlib or plotly (depending on the backend).
        """
        fig = self.plot(backend, cmap)
        if backend == "matplotlib":
            fig.savefig(path)
        elif backend == "plotly":
            fig.write_html(path)


class GlycoShield:
    """
    The GlycoShield class can be used to perform a quick simulation of the shielding effect of glycans on a scaffold molecule.

    Attributes
    ----------
    df: pandas.DataFrame
        The dataframe with exposure values for each scaffold residue
    conformation_viewer: Py3dmolViewer
        A Py3Dmolviewer with all accepted conformations overlaid (if conformations were visualized during simulation)
    conformations: Scaffold
        A copy of the scaffold with additional models; one for each accepted conformation (if conformations were recorded during simulation)

    Parameters
    ----------
    scaffold
        A scaffold molecule with glycans to analyze
    exposure_metric: callable
        The metric function to use for evaluating "exposure". This function will receive three arguments (in order): `self` (the GlycoShield object itself), `scaffold_coords` (the coordinates of the scaffold residues), and `glycan_coords` (the coordinates of the glycan residues).
        The function must return a 1D array of exposure values for each scaffold residue. By default `line_of_sight_exposzure` is used.
    concatenation_function: callable
        The function to use to concatenate the exposure values from each iteration of the simulation into the final dataframe. The default is `np.add` which will sum the exposure values from each iteration.
        This function should take two 1D arrays (the total exposure and the exposure from the current iteration) and return a new 1D array (the updated total exposure).
    """

    def __init__(
        self,
        scaffold: core.Scaffold,
        exposure_metric: callable = None,
        concatenation_function: callable = None,
    ):
        if scaffold.count_glycans() == 0:
            raise ValueError("The scaffold must have at least one glycan attached")
        self.scaffold = (
            scaffold.copy()
        )  # the copy is a precaution since the simulation will modify the scaffold

        self.exposure = ExposureData()
        self.exposure.scaffold = self.scaffold

        self.conformation_viewer = None
        self.conformations = None

        self._scaffold_residues = list(self.scaffold.scaffold_residues)
        self._scaffold_coords = np.array([i.coord for i in self._scaffold_residues])

        glycan_coords = [
            [i.coord for i in glycan.get_residues()] for glycan in scaffold.glycans
        ]
        glycan_coords = np.concatenate(glycan_coords).reshape(-1, 3)
        glycan_index_masks = []
        start = 0
        for glycan in scaffold.glycans:
            end = start + glycan.count_residues()
            mask = np.zeros(len(glycan_coords), dtype=bool)
            mask[start:end] = True
            glycan_index_masks.append(mask)
            start = end

        self._glycan_coords = glycan_coords
        self._glycan_index_masks = glycan_index_masks

        if exposure_metric is None:
            exposure_metric = line_of_sight_exposure
            setup_line_of_sight_exposure(self)
        if concatenation_function is None:
            concatenation_function = lambda x, y: x + y

        self.exposure_metric = exposure_metric
        self.concatenation_function = concatenation_function

    @property
    def df(self):
        return self.exposure.df

    def simulate(
        self,
        repeats: int = 3,
        edge_samples: int = 3,
        angle_step: int = 30,
        coarse_precheck: bool = True,
        visualize_conformations: bool = False,
        save_conformations_to: str = None,
        capture_full_scaffold: bool = False,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """
        Perform a quick simulation of the shielding effect of glycans on a scaffold molecule.

        Parameters
        ----------
        repeats: int
            The number of times to repeat the edge sampling and simulation process for each glycan. If edge_samples is "all" this will be ignored.
        edge_samples: int or str
            The number of rotatable edges within the glycan to sample. Either give a number
            of edges to sample from. Alternatively, provide "all" to use all available rotatable edges in the glycan (e.g. even bonds to protruding OH-groups),
            or "connections" to use all inter-residue connections (i.e. only bonds connecting different sugars).
        angle_step: int
            The step size in degrees to rotate each sampled rotatable edge during the simulation.
            The smaller this is the more fine-grained the simulation will be, but the longer it will take.
        coarse_precheck: bool
            When computing whether a specific conformation clashes with the scaffold (and should therefore be ignored)
            perform a coarse-grained precomputation to speed things up. This will greatly enhance the performance but may
            miss clashes if individual scaffold residues are particularly large (e.g. membrane lipids). If you find the results
            to be inaccurate, set this to False to perform a more accurate (albeit much slower) check.
        visualize_conformations: bool
            If True, each accepted conformation will be visualized in Py3DMol and a combined view (total overlay of all conformations) will be returned at the end together with the dataframe.
            This will make the computation slower.
        save_conformations_to: str
            If a path for a directory is provided, the conformations will be saved to that path as PDB files (one file per glycan, containing models for each accepted conformation).
        capture_full_scaffold: bool
            If True the conformations that are saved are the full scaffold. Otherwise only the glycan in question is saved.
        verbose: bool
            If True, a progress bar will be shown during the simulation.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the "exposure" of each residue in the scaffold at each glycan conformation.
            The more "exposed" a residue is the further away it is from the shielding effect of glycans.
            Exposure is measured as the minimal euclidian distance of each residue to any glycan residue.
        Py3DmolViewer
            The visualization of the glycan conformations (if visualize_conformations is True). Otherwise only the dataframe is returned!
        """
        if edge_samples == "all":
            edge_sampler = (
                lambda glycan, glycan_graph: glycan_graph.find_rotatable_edges()
            )
            repeats = 1
        elif edge_samples == "connections":
            edge_sampler = lambda glycan, glycan_graph: glycan.get_residue_connections()
        elif isinstance(edge_samples, int):
            edge_sampler = lambda glycan, glycan_graph: (
                glycan_graph.sample_edges(
                    glycan.get_residue_connections(), m=edge_samples
                )
                if glycan.count_residues() > 1
                else []
            )
        else:
            raise ValueError(
                f"edge_samples must be 'all', 'connections', or an integer, got {edge_samples}"
            )

        angle_range = range(0, 360 + angle_step, angle_step)
        root_angle_range = range(0, 360 + min(angle_step, 30), min(angle_step, 30))

        # we need to convert to radians since the _rotate_around_bond function expects radians
        angle_range = np.deg2rad(angle_range)
        root_angle_range = np.deg2rad(root_angle_range)

        exposures = np.zeros(len(self._scaffold_residues))

        if visualize_conformations:
            v = self.scaffold.py3dmol()
            self.conformation_viewer = v
        else:
            v = _dummy()

        metric = self.exposure_metric
        concatenate = self.concatenation_function

        record_conformations = save_conformations_to is not None
        if record_conformations:
            save_conformations_to = Path(save_conformations_to)
            self.conformations = []

            def save_model(model):
                model = model.copy()
                self.conformations.append(model)

        else:

            def save_model(model):
                pass

        total = len(self.scaffold.glycans) * repeats * len(root_angle_range)
        if verbose:
            _bar = utils.progress_bar
        else:
            _bar = utils.DummyBar
        with _bar(total) as bar:

            for gdx, (root, glycan) in enumerate(self.scaffold.get_glycans().items()):
                glycan_graph = glycan._AtomGraph

                # a backup so we can reset the glycan to its original conformation
                glycan_atoms = list(glycan.get_atoms())
                current_coords = [deepcopy(i.coord) for i in glycan_atoms]

                if record_conformations:
                    self.conformations = []

                for angle in root_angle_range:
                    # first rotate around the root connection
                    self.scaffold._rotate_around_bond(
                        root, glycan.root_atom, angle, descendants_only=True
                    )
                    if glycan.clashes_with_scaffold(coarse_precheck=coarse_precheck):
                        bar.update(repeats)
                        continue

                    glycan_coords = np.array([i.coord for i in glycan.get_residues()])
                    self._glycan_coords[self._glycan_index_masks[gdx]] = glycan_coords

                    exposure = metric(self, self._scaffold_coords, self._glycan_coords)
                    exposures = concatenate(exposures, exposure)

                    v.add(glycan, style={"stick": {"color": "lightblue"}})

                    if capture_full_scaffold:
                        save_model(self.scaffold._model)
                    else:
                        save_model(glycan._model)

                    for _ in range(repeats):
                        sampled_edges = edge_sampler(glycan, glycan_graph)
                        sampled_edges = glycan_graph.direct_edges(
                            glycan.root_atom, sampled_edges
                        )

                        # then rotate around the sampled edges
                        for edge in sampled_edges:
                            for angle in angle_range:

                                glycan._rotate_around_bond(
                                    *edge, angle, descendants_only=True
                                )
                                if (
                                    glycan.count_clashes()
                                    or glycan.clashes_with_scaffold(
                                        coarse_precheck=coarse_precheck
                                    )
                                ):
                                    continue

                                glycan_coords = np.array(
                                    [i.coord for i in glycan.get_residues()]
                                )
                                self._glycan_coords[self._glycan_index_masks[gdx]] = (
                                    glycan_coords
                                )

                                exposure = metric(
                                    self, self._scaffold_coords, self._glycan_coords
                                )
                                exposures = concatenate(exposures, exposure)

                                v.add(glycan, style={"stick": {"color": "lightblue"}})
                                if capture_full_scaffold:
                                    save_model(self.scaffold._model)
                                else:
                                    save_model(glycan._model)
                        bar.update(1)

                    # now reset the glycan to its original conformation
                    for idx, atom in enumerate(glycan_atoms):
                        atom.coord = current_coords[idx]

                if record_conformations:
                    _id = glycan.id.split("/")[-1]
                    if capture_full_scaffold:
                        confs = self.scaffold.copy()
                    else:
                        confs = glycan.copy()

                    for i in self.conformations[1:]:
                        confs.add_model(i)

                    confs.to_pdb(save_conformations_to / f"{_id}_conformations.pdb")

        exposures /= len(self.scaffold.glycans) * repeats * len(angle_range)

        _df = {
            "chain": [i.parent.id for i in self._scaffold_residues],
            "resname": [i.resname for i in self._scaffold_residues],
            "serial": [i.serial_number for i in self._scaffold_residues],
            "exposure": exposures,
        }
        self.exposure.df = pd.DataFrame(_df)

        if visualize_conformations:
            return self.df, v
        return self.df

    def py3dmol(
        self,
        add_to_conformation_view: bool = True,
        color: str = "white",
        glycan_color: str = None,
        colormap: str = "viridis",
        opacity: float = 0.8,
    ):
        """
        Visualize the exposure in Py3DMol

        Parameters
        ----------
        add_to_conformation_view: bool
            If conformations were already recorded during `simulate`
            add the exposure onto that same view. Otherwise a new view is made.
        color: str
            The color of the scaffold residues
        glycan_color: str
            The color of the glycan residues
        colormap: str
            The name of the colormap to use for visualizing the exposure values.
            This can be any matplotlib-recognized colormap.
        opacity: float
            The opacity for the highlighted exposure values (will be shown in "sphere" style on top of the normal scaffold visualization).

        Returns
        -------
        Py3DmolViewer
            The viewer containing the visualization
        """
        if self.df is None:
            raise ValueError("No exposure data available. Run simulate() first.")

        get_color, _ = _prepare_colormap(self.df, colormap)

        if self.conformation_viewer is not None and add_to_conformation_view:
            v = self.conformation_viewer
        else:
            v = self.scaffold.py3dmol(color=color, glycan_color=glycan_color)

        for residue in self._scaffold_residues:
            pdb = utils.pdb.make_atoms_table(residue)
            v.add(
                pdb, style={"sphere": {"color": get_color(residue), "opacity": opacity}}
            )
        return v

    def plotly(self, glycan_color: str = "black", colormap: str = "viridis"):
        """
        Visualize the exposure in Plotly

        Parameters
        ----------
        glycan_color: str
            The color of the glycan residues
        colormap: str
            The name of the colormap to use for visualizing the exposure values.
            This can be any matplotlib-recognized colormap.

        Returns
        -------
        MoleculeViewer3D
            The viewer containing the visualization
        """
        v = utils.visual.MoleculeViewer3D()
        get_color, _ = _prepare_colormap(self.df, colormap)

        for residue in self._scaffold_residues:
            v.draw_point(
                id=residue.get_id(),
                coord=residue.coord,
                color=get_color(residue, hex=True),
                showlegend=False,
            )

        for glycan in self.scaffold.glycans:
            v += glycan.plotly(line_color=glycan_color, atoms=False)

        return v

    def draw(self, *args, **kwargs):
        if utils.visual.DEFAULT_BACKEND == "py3dmol":
            return self.py3dmol(*args, **kwargs)
        elif utils.visual.DEFAULT_BACKEND == "plotly":
            return self.plotly(*args, **kwargs)
        else:
            raise RuntimeError(
                "The GlycoShield only supports Py3DMol and Plotly default backends. Use a dedicated method directly or change the backend."
            )

    def show(self, *args, **kwargs):
        self.draw(*args, **kwargs).show()


def _prepare_colormap(df, cmap_name):
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib as mpl

    data = df["exposure"]

    norm = mcolors.Normalize(vmin=data.min(), vmax=data.max())
    cmap = mpl.colormaps.get_cmap(cmap_name)

    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    def getter_from_residue(residue, hex=False):
        value = df.query(
            f"chain == '{residue.parent.id}' and serial == {residue.serial_number}"
        )["exposure"].values[0]
        r, g, b, _ = mapper.to_rgba(value)
        if hex:
            return mpl.colors.to_hex((r, g, b))
        return f"rgb({r*255}, {g*255}, {b*255})"

    def getter_from_value(value, hex=False):
        r, g, b, _ = mapper.to_rgba(value)
        if hex:
            return mpl.colors.to_hex((r, g, b))
        return f"rgb({r*255}, {g*255}, {b*255})"

    return getter_from_residue, getter_from_value


class _dummy:
    def add(*args, **kwargs):
        pass


def _line_of_sight_vector(scaf_coords, idx, length=20) -> np.ndarray:
    res_coord = scaf_coords[idx]
    close_by_mask = np.linalg.norm(scaf_coords - res_coord, axis=1) < 15
    if not np.any(close_by_mask):
        return np.zeros(3)

    close_by_coords = scaf_coords[close_by_mask]

    vec = structural.plane_of_points(close_by_coords)
    if np.linalg.norm(close_by_coords - res_coord - vec) < np.linalg.norm(
        close_by_coords - res_coord + vec
    ):
        vec = -vec

    vec *= length
    return vec


def _line_to_linspace(start, vec, n_points):
    return np.linspace(start, start + vec, n_points)


def _distance_to_line(points, line_linspace):
    dists = cdist(points, line_linspace)
    dists = dists.min()
    return dists


def line_of_sight_exposure(self, scaffold_coords, glycan_coords):
    """
    Line of sight exposure measures the distances of glycan residue coordinates to
    perpendicular cylinders that protrude from each scaffold residue's surface upward (the residue's "line of sight").
    The exposure is measured as the minimal distance of each line of sight to any glycan residue.
    """
    straight_line_coverage = np.zeros(len(scaffold_coords))
    for rdx in range(len(scaffold_coords)):
        linspace = self._line_of_sight_linspaces[rdx]
        distance = _distance_to_line(glycan_coords, linspace)
        straight_line_coverage[rdx] = distance

    return straight_line_coverage


def setup_line_of_sight_exposure(glycoshield, length: int = 20):
    """
    Perpare the Glycoshield object for line of sight exposure calculations.
    This is necessary since the function will draw on some precomputed values for
    enhanced performance.
    (This function is called automatically if a GlycoShield is setup using default settings)
    """
    glycoshield._line_of_sight_linspaces = []
    for rdx in range(len(glycoshield._scaffold_coords)):
        vec = _line_of_sight_vector(glycoshield._scaffold_coords, rdx, length=length)
        linspace = _line_to_linspace(
            glycoshield._scaffold_coords[rdx], vec, length // 3
        )
        glycoshield._line_of_sight_linspaces.append(linspace)
    glycoshield._line_of_sight_linspaces = np.array(
        glycoshield._line_of_sight_linspaces
    )


if __name__ == "__main__":
    prot = utils.load_pickle(
        "/Users/noahhk/GIT/glycosylator/__projects__/solf3/solF_plus_3G_rsr017_coot_30_man5.pdb_glycosylated_optimized.pkl"
    )
    shield = GlycoShield(prot)
    from time import time

    print("shielding")
    t1 = time()
    df = shield.simulate(edge_samples=3, repeats=1)
    print("done, took " + str(time() - t1) + " seconds")
    fig = shield.plot_exposure()
    fig.savefig("exposure.png")
