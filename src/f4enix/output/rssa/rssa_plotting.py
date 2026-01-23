from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import polars as pl
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from f4enix.egroups import GROUP_STRUCTURES
from f4enix.output.rssa.rssa_helpers import NEUTRON_INDICATOR
from f4enix.output.rssa.rssa_reader import FileParameters


@dataclass
class PlotParameters:
    title: str = ""
    xlabel: str = "X [cm]"
    ylabel: str = "Y [cm]"
    legend_label: str = ""
    legend_orientation: Literal["vertical", "horizontal"] = "horizontal"
    number_of_colors: int = 10
    norm: Literal["linear", "log"] = "log"
    vmin: float | None = None
    vmax: float | None = None


@dataclass
class SpectraInfo:
    normalized_counts: np.ndarray
    energy_bins: np.ndarray
    label: str = ""


class PlottingFunctions(ABC):
    tracks: pl.LazyFrame
    rssa_parameters: FileParameters

    def set_particle(self, particle_type: Literal["n", "p"]):
        """Set the particle type to filter the tracks."""
        if particle_type == "n":
            self.tracks = self.tracks.filter(pl.col("b") == NEUTRON_INDICATOR)
        elif particle_type == "p":
            self.tracks = self.tracks.filter(pl.col("b") != NEUTRON_INDICATOR)
        return self

    def set_surface_ids(self, surface_ids: list[int]):
        """Set the surface IDs to filter the tracks."""
        valid_surface_ids = [s.id for s in self.rssa_parameters.surfaces]
        if not all(sid in valid_surface_ids for sid in surface_ids):
            raise ValueError(
                f"Some surface IDs are not valid. Valid IDs are: {valid_surface_ids}"
            )
        self.tracks = self.tracks.filter(pl.col("c").is_in(surface_ids))
        return self

    def set_z_limits(self, vmin: float, vmax: float):
        """Set the z limits for the plot."""
        self.tracks = self.tracks.filter(
            pl.col("z").is_between(vmin, vmax, closed="both")
        )
        return self

    def set_perimeter_limits(self, vmin: float, vmax: float):
        """Set the limits for the perimeter positions."""
        if "perimeter_pos" not in self.tracks.collect_schema().names():
            self.calculate_perimeter_positions()
        self.tracks = self.tracks.filter(
            pl.col("perimeter_pos").is_between(vmin, vmax, closed="both")
        )
        return self

    def calculate_perimeter_positions(self) -> None:
        """
        It adds a new column to the tracks DataFrame called 'perimeter_pos' that takes
        the X and Y coordinates and calculates the position as a perimeter coordinate.
        The perimeter position is calculated as theta * r, where theta is the angle in
        radians and r is the average radius of all the points.
        """
        radius = (pl.col("x").pow(2) + pl.col("y").pow(2)).sqrt()
        thetas = pl.arctan2(pl.col("y"), pl.col("x"))
        perimeter_pos = (thetas * radius).alias("perimeter_pos")
        self.tracks = self.tracks.with_columns(perimeter_pos)

    @abstractmethod
    def get_plot(self) -> tuple[Figure, Axes]: ...

    def save_figure(self, out_path: Path | str):
        """Save the figure to the specified path."""
        fig, _ax = self.get_plot()

        out_path = Path(out_path)
        if out_path.suffix != ".png":
            out_path = Path(out_path).with_suffix(".png")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        return self

    def show(self):
        """Show the plot."""
        fig, _ax = self.get_plot()
        fig.show()
        return self


class RSSAPlot(PlottingFunctions):
    def __init__(
        self,
        tracks: pl.DataFrame,
        rssa_parameters: FileParameters,
        x_col: str = "x",
        y_col: str = "y",
    ):
        self.tracks = tracks.lazy()
        self.rssa_parameters = rssa_parameters
        self.x_col = x_col
        self.y_col = y_col
        self._x_bins: pl.Series | None = None
        self._y_bins: pl.Series | None = None
        self.raster: np.ndarray | None = None
        self.plot_parameters = PlotParameters()

    @property
    def x_bins(self) -> pl.Series:
        """Returns the x bins for the plot."""
        if self._x_bins is None:
            raise ValueError("X bins are not set. Call set_bins() or calculate_bins().")
        return self._x_bins

    @property
    def y_bins(self) -> pl.Series:
        """Returns the y bins for the plot."""
        if self._y_bins is None:
            raise ValueError("Y bins are not set. Call set_bins() or calculate_bins().")
        return self._y_bins

    def set_plot_parameters(self, plot_parameters: PlotParameters) -> "RSSAPlot":
        """Set the plot parameters for the plot."""
        self.plot_parameters = plot_parameters
        return self

    def set_bins(self, x_bins: Sequence[float], y_bins: Sequence[float]) -> "RSSAPlot":
        """Set the x and y bins for the plot."""
        self._x_bins = pl.Series("x_bins", x_bins).sort()
        self._y_bins = pl.Series("y_bins", y_bins).sort()
        return self

    def calculate_bins(self, bin_width: float = 10.0) -> "RSSAPlot":
        """Automatically calculate the bins for the x and y coordinates of the plot by
        giving a bin width in cm. Instead use `set_bins()` to apply custom bins."""
        collected_tracks = self.tracks.collect()
        if collected_tracks.is_empty():
            raise ValueError("The tracks DataFrame is empty at this point.")
        x_min = collected_tracks[self.x_col].min()
        x_max = collected_tracks[self.x_col].max()
        y_min = collected_tracks[self.y_col].min()
        y_max = collected_tracks[self.y_col].max()
        return self.set_bins(
            np.arange(x_min, x_max + bin_width, bin_width),  # type: ignore
            np.arange(y_min, y_max + bin_width, bin_width),  # type: ignore
        )

    def get_particle_current(self, source_intensity: float) -> "RSSAPlot":
        """Calculate the particle current from the tracks. It automatically divides the
        weight by the nps value."""
        self.apply_source_intensity(source_intensity)
        self.divide_by_nps()
        raster = self._get_2d_grid_of_weights(agg_func="sum")
        self.raster = raster / calculate_areas(self.x_bins, self.y_bins)
        return self

    def get_particle_current_errors(self) -> "RSSAPlot":
        """Calculate the particle current errors from the tracks as the square root of
        the number of tracks in each bin.
        """
        raster = self._get_2d_grid_of_weights(agg_func="count")
        self.raster = 1 / (raster**0.5)
        return self

    def apply_source_intensity(self, source_intensity: float) -> "RSSAPlot":
        """Apply a source intensity to the weights of the tracks."""
        self.tracks = self.tracks.with_columns(
            (pl.col("wgt") * source_intensity).alias("wgt")
        )
        return self

    def divide_by_nps(self) -> "RSSAPlot":
        """Divide the weights by the number of histories. The nps value used is the one
        read in the header abs(np1)."""
        self.tracks = self.tracks.with_columns(
            (pl.col("wgt") / abs(self.rssa_parameters.np1)).alias("wgt")
        )
        return self

    def get_ratio_to(self, other: "RSSAPlot") -> "RSSAPlot":
        """Calculate the ratio of the current to another RSSAPlot instance in %
        difference. Calculated as (other - self) / self * 100. The other RSSAPlot
        instance should have undergone the same processing as the current one."""
        if self.raster is None or other.raster is None:
            raise ValueError(
                "Raster data is not calculated. Call get_particle_current()"
                " or other first for both `self` and `other`."
            )
        self.raster = (other.raster - self.raster) / self.raster * 100
        return self

    def get_plot(self) -> tuple[Figure, Axes]:
        """Returns the Matplotlib figure and axes of the plot for further manual
        customization."""
        if self.raster is None:
            raise ValueError(
                "Raster data is not calculated. Call get_particle_current()"
                " or other first."
            )
        norm = (
            colors.LogNorm if self.plot_parameters.norm == "log" else colors.Normalize
        )

        fig, ax = plt.subplots()
        im = ax.pcolormesh(
            self.x_bins.to_numpy(),
            self.y_bins.to_numpy(),
            self.raster,
            cmap=plt.get_cmap("jet", self.plot_parameters.number_of_colors),
            norm=norm(self.plot_parameters.vmin, self.plot_parameters.vmax),
        )
        fig.colorbar(
            im,
            ax=ax,
            label=self.plot_parameters.legend_label,
            orientation=self.plot_parameters.legend_orientation,
        )
        ax.set_aspect("equal")
        ax.set_title(self.plot_parameters.title)
        ax.set_xlabel(self.plot_parameters.xlabel)
        ax.set_ylabel(self.plot_parameters.ylabel)
        return fig, ax

    def _get_2d_grid_of_weights(
        self,
        agg_func: Literal["sum", "count"] = "sum",
    ) -> np.ndarray:
        # Remove points outside of the bins
        filtered_df = (
            self.tracks.select([self.x_col, self.y_col, "wgt"])
            .filter(
                pl.col(self.x_col).is_between(
                    self.x_bins[0],
                    self.x_bins[-1],
                    closed="left",
                ),
                pl.col(self.y_col).is_between(
                    self.y_bins[0],
                    self.y_bins[-1],
                    closed="left",
                ),
            )
            .collect()
        )

        # Decide if the weights are summed or the number of tracks are counted
        agg_expression = (
            pl.col("wgt").sum() if agg_func == "sum" else pl.col("wgt").count()
        )
        grid = (
            filtered_df.lazy()
            # Find the bin indices for each point
            .with_columns(
                (
                    self.x_bins.search_sorted(filtered_df[self.x_col], side="right") - 1
                ).alias("bin_x"),
                (
                    self.y_bins.search_sorted(filtered_df[self.y_col], side="right") - 1
                ).alias("bin_y"),
            )
            # Group by the bin indices and sum the weights
            .group_by(["bin_x", "bin_y"])
            .agg(agg_expression)
            .collect()
        )
        raster = _get_raster(grid, self.x_bins, self.y_bins)
        return raster

    # Override the return types of PlottingFunctions methods
    # necessary due to the lack of typing.Self in Python <3.11
    def set_particle(self, particle_type: Literal["n", "p"]) -> "RSSAPlot":
        return super().set_particle(particle_type)

    def set_surface_ids(self, surface_ids: list[int]) -> "RSSAPlot":
        return super().set_surface_ids(surface_ids)

    def set_z_limits(self, vmin: float, vmax: float) -> "RSSAPlot":
        return super().set_z_limits(vmin, vmax)

    def set_perimeter_limits(self, vmin: float, vmax: float) -> "RSSAPlot":
        return super().set_perimeter_limits(vmin, vmax)

    def save_figure(self, out_path: Path | str) -> "RSSAPlot":
        return super().save_figure(out_path)

    def show(self) -> "RSSAPlot":
        return super().show()


def _get_raster(
    grid: pl.DataFrame,
    x_bins: pl.Series,
    y_bins: pl.Series,
) -> np.ndarray:
    # The dimensions of our grid are determined by the number of bins
    num_y_bins = len(y_bins) - 1
    num_x_bins = len(x_bins) - 1
    raster = np.zeros((num_y_bins, num_x_bins), dtype=np.float64)

    # Extract columns to NumPy and use advanced indexing to fill the raster
    bin_x = grid.get_column("bin_x").to_numpy()
    bin_y = grid.get_column("bin_y").to_numpy()
    wgt = grid.get_column("wgt").to_numpy()

    # The bin indices from Polars directly correspond to the raster indices
    raster[bin_y, bin_x] = wgt

    return raster


def calculate_areas(
    x_bins: Sequence[float] | pl.Series, y_bins: Sequence[float] | pl.Series
) -> np.ndarray:
    """Calculate the areas of the bins in a 2D histogram.

    Parameters
    ----------
    x_bins: Sequence[float] | pl.Series
        The x-axis bin edges.
    y_bins: Sequence[float] | pl.Series
        The y-axis bin edges.

    Returns
    -------
    np.ndarray
        A 2D array with the areas of each bin.
    """
    x_edges = np.array(x_bins)
    y_edges = np.array(y_bins)
    dx = np.diff(x_edges)
    dy = np.diff(y_edges)
    return np.outer(dy, dx)


@dataclass
class SpectraPlotParameters:
    title: str = "Particle energy spectra"
    xlabel: str = "Energy [eV]"
    ylabel: str = "Normalized counts per unit of lethargy"
    label: str = ""
    min_energy: float | None = None
    max_energy: float | None = None


class RSSASpectraPlot(PlottingFunctions):
    def __init__(self, tracks: pl.DataFrame, rssa_parameters: FileParameters):
        self.tracks = tracks.lazy()
        self.rssa_parameters = rssa_parameters
        self.energy_bins: np.ndarray = GROUP_STRUCTURES["VITAMIN-J-175"]
        self.plot_parameters = SpectraPlotParameters(
            title="Energy spectra",
            xlabel="Energy [eV]",
            ylabel="Normalized counts per unit of lethargy",
        )

    def set_plot_parameters(
        self, plot_parameters: SpectraPlotParameters
    ) -> "RSSASpectraPlot":
        self.plot_parameters = plot_parameters
        return self

    def set_energy_bins(self, energy_bins: Sequence[float]) -> "RSSASpectraPlot":
        """Set the energy bins for the spectra plot."""
        self.energy_bins = np.asarray(energy_bins)
        return self

    def get_spectra_info(self) -> SpectraInfo:
        # Get the energies and weights
        data = self.tracks.select(
            (pl.col("erg") * 1e6).alias("erg"),  # Convert energy from MeV to eV,
            pl.col("wgt"),
        ).collect()

        # Group the data by energy bins
        weighted_counts, _ = np.histogram(
            data["erg"],
            bins=self.energy_bins,
            weights=data["wgt"],  # each contribution is weighted by the particle weight
        )

        # Calculate the counts per lethargy
        lethargies = np.log(self.energy_bins[1:] / self.energy_bins[:-1])
        normalized_counts = (
            weighted_counts
            / lethargies  # The Y axis is in units of energy per lethargy
            / np.sum(weighted_counts)  # Normalize to the total number of counts to 1
        )

        return SpectraInfo(
            normalized_counts=normalized_counts,
            energy_bins=self.energy_bins,
            label=self.plot_parameters.label,
        )

    def get_plot(self) -> tuple[Figure, Axes]:
        spectra_info = self.get_spectra_info()

        # Create the plot
        fig, ax = plt.subplots()
        ax.step(
            spectra_info.energy_bins[:-1],
            spectra_info.normalized_counts,
            label=spectra_info.label,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        if spectra_info.label:
            ax.legend()
        ax.set_xlabel(self.plot_parameters.xlabel)
        ax.set_ylabel(self.plot_parameters.ylabel)
        ax.set_title(self.plot_parameters.title)
        ax.set_xlim(self.plot_parameters.min_energy, self.plot_parameters.max_energy)
        ax.grid(True)
        return fig, ax

    def get_combined_plot_with_other_spectras(
        self, *other_spectra_info: SpectraInfo
    ) -> tuple[Figure, Axes]:
        """Combine the current spectra plot with other spectra for comparison."""
        current_spectra_info = self.get_spectra_info()

        fig, ax = plt.subplots()
        ax.step(
            current_spectra_info.energy_bins[:-1],
            current_spectra_info.normalized_counts,
            label=current_spectra_info.label,
        )

        for spectra_info in other_spectra_info:
            ax.step(
                spectra_info.energy_bins[:-1],
                spectra_info.normalized_counts,
                label=spectra_info.label,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()
        ax.set_xlabel(self.plot_parameters.xlabel)
        ax.set_ylabel(self.plot_parameters.ylabel)
        ax.set_title(self.plot_parameters.title)
        ax.set_xlim(self.plot_parameters.min_energy, self.plot_parameters.max_energy)
        ax.grid(True)

        return fig, ax

    # Override the return types of PlottingFunctions methods
    # necessary due to the lack of typing.Self in Python <3.11
    def set_particle(self, particle_type: Literal["n", "p"]) -> "RSSASpectraPlot":
        return super().set_particle(particle_type)

    def set_surface_ids(self, surface_ids: list[int]) -> "RSSASpectraPlot":
        return super().set_surface_ids(surface_ids)

    def set_z_limits(self, vmin: float, vmax: float) -> "RSSASpectraPlot":
        return super().set_z_limits(vmin, vmax)

    def set_perimeter_limits(self, vmin: float, vmax: float) -> "RSSASpectraPlot":
        return super().set_perimeter_limits(vmin, vmax)

    def save_figure(self, out_path: Path | str) -> "RSSASpectraPlot":
        return super().save_figure(out_path)

    def show(self) -> "RSSASpectraPlot":
        return super().show()
