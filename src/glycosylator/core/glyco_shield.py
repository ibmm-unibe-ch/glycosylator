"""
The `GlycoShield` class offers a convenient interface to perform a quick
simulation of the shielding effect of glycans on a scaffold molecule, and can be
used to obtain numeric data or visualizations of the scaffold at different glycan conformations.
"""

import numpy as np
from scipy.spatial.distance import cdist
import pandas as pd

import glycosylator.core as core
import glycosylator.utils as utils

__all__ = ["quickshield", "GlycoShield"]


def quickshield(scaffold) -> "GlycoShield":
    """
    A convenience function to perform a quick simulation of the glycan shielding on a scaffold.
    This function will set up a GlycoShield and perform a standard simulation. Use a manual procedure
    to retain more control.

    Parameters
    ----------
    scaffold
        The scaffold molecule with attached glycans

    Returns
    -------
    GlycoShield
        The GlycoShield object that stores the analyses results such as exposure dataframe (`df` attribute)
        or any available visualizations.
    """
    shield = GlycoShield(scaffold)
    shield.simulate()
    return shield


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
    """

    def __init__(self, scaffold: core.Scaffold):
        if scaffold.count_glycans() == 0:
            raise ValueError("The scaffold must have at least one glycan attached")
        self.scaffold = scaffold
        self.df = None
        self.conformation_viewer = None
        self.conformations = None
        self._scaffold_residues = list(self.scaffold.scaffold_residues)
        self._scaffold_coords = np.array([i.coord for i in self._scaffold_residues])

        glycan_coords = [
            [i.coord for i in glycan.get_residues()] for glycan in scaffold.glycans
        ]
        glycan_coords = np.array(glycan_coords).reshape(-1, 3)
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

    def simulate(
        self,
        repeats: int = 3,
        edge_samples: int = 3,
        angle_step: int = 30,
        coarse_precheck: bool = True,
        visualize_conformations: bool = False,
        record_conformations: bool = False,
    ) -> pd.DataFrame:
        """
        Perform a quick simulation of the shielding effect of glycans on a scaffold molecule.

        Parameters
        ----------
        repeats: int
            The number of times to repeat the simulation process for each glycan.
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
        record_conformations: bool
            If True, each accepted conformation is saved and returned at the end. This will make the computation considerably slower!
            Also since each model will be a copy of the entire scaffold this will consume a lot of memory! If you do not absolutely require these models
            do *not* set this to true!

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
        elif edge_samples == "connections":
            edge_sampler = lambda glycan, glycan_graph: glycan.get_residue_connections()
        elif isinstance(edge_samples, int):
            edge_sampler = lambda glycan, glycan_graph: glycan_graph.sample_edges(
                glycan.get_residue_connections(), m=edge_samples
            )
        else:
            raise ValueError(
                f"edge_samples must be 'all', 'connections', or an integer, got {edge_samples}"
            )

        angle_range = range(0, 360, angle_step)

        exposures = np.zeros(len(self._scaffold_residues))

        if visualize_conformations:
            v = self.scaffold.py3dmol()
            self.conformation_viewer = v
        else:
            v = _dummy()

        if record_conformations:
            self.conformations = self.scaffold.copy()

            def save_model(model):
                model = model.copy()
                self.conformations._base_struct.child_list.append(model)
                self.conformations._base_struct.child_dict[model.get_id()] = model

        else:

            def save_model(model):
                pass

        for gdx, glycan in enumerate(self.scaffold.glycans):
            glycan_graph = glycan._AtomGraph

            for _ in range(repeats):
                sampled_edges = edge_sampler(glycan, glycan_graph)
                sampled_edges = glycan_graph.direct_edges(
                    glycan.root_atom, [tuple(i) for i in sampled_edges]
                )

                for edge in sampled_edges:
                    for angle in angle_range:

                        glycan.rotate_descendants(*edge, angle)
                        if glycan.count_clashes() or glycan.clashes_with_scaffold(
                            coarse_precheck=coarse_precheck
                        ):
                            continue

                        glycan_coords = np.array(
                            [i.coord for i in glycan.get_residues()]
                        )
                        self._glycan_coords[self._glycan_index_masks[gdx]] = (
                            glycan_coords
                        )

                        distances = cdist(self._scaffold_coords, self._glycan_coords)
                        exposure = distances.min(axis=1)
                        exposures += exposure

                        v.add(glycan, style={"stick": {"color": "lightblue"}})
                        save_model(self.scaffold._model)

        exposures /= len(self.scaffold.glycans) * repeats * len(angle_range)

        _df = {
            "chain": [i.parent.id for i in self._scaffold_residues],
            "resname": [i.resname for i in self._scaffold_residues],
            "serial": [i.serial_number for i in self._scaffold_residues],
            "exposure": exposures,
        }
        self.df = pd.DataFrame(_df)

        to_return = (self.df,)

        if visualize_conformations:
            to_return = to_return + (v,)
        if record_conformations:
            to_return = to_return + (self.conformations,)

        if len(to_return) == 1:
            to_return = to_return[0]

        return to_return

    def plot_exposure(self, backend: str = "matplotlib", cmap: str = "viridis"):
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

        elif backend == "plotly":
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
            v += glycan.draw(line_color=glycan_color, atoms=False)

        return v

    # alias
    draw = plotly


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


if __name__ == "__main__":
    prot = utils.load_pickle(
        "/Users/noahhk/GIT/glycosylator/__projects__/solf3/solF_plus_3G_rsr017_coot_30_man5.pdb_glycosylated_optimized.pkl"
    )
    shield = GlycoShield(prot)
    from time import time

    print("shielding")
    t1 = time()
    df = shield.simulate(record_conformations=True, repeats=1)
    print("done, took " + str(time() - t1) + " seconds")
    fig = shield.plot_exposure()
    fig.savefig("exposure.png")
