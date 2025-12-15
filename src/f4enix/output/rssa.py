"""This module is related to the parsing of D1S-UNED meshinfo files."""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Literal

import numpy as np
import polars as pl
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

BYTE = np.byte
CHAR = np.char
INT = np.int32
FLOAT = np.float64
LONG = np.int64
NUMBER_OF_EXPECTED_VALUES = 11  # Number of values recorded for each particle
NEUTRON_INDICATOR = 8  # Neutron indicator in the packed variable


@dataclass
class _SurfaceParameters:
    id: int
    info: int
    type: int
    num_parameters: int
    parameters: list[int]


@dataclass
class _FileParameters:
    np1: int  # Number of histories of the simulation, given as a negative number
    nrss: int  # Number of tracks recorded
    nrcd: int  # Number of values recorded for each particle, it should be 11
    njsw: int  # Number of surfaces in JASW
    niss: int  # Number of different histories that reached the SSW surfaces
    niwr: int  # Number of cells in RSSA file
    mipts: int  # Source particle type
    kjaq: int  # Flag for macrobodies surfaces
    surfaces: list[_SurfaceParameters]  # List with surface ids that appear in the file


class RSSA:
    def __init__(self, path: Path | str):
        """Representation of a RSSA file.

        Parameters
        ----------
        path: Path
            Path to the RSSA file.

        Attributes
        ----------
        path: Path
            Path to the RSSA file.
        parameters: _FileParameters
            Parameters extracted from the RSSA file header.

            np1   # Number of histories of the simulation, given as a negative number
            nrss  # Number of tracks recorded
            nrcd  # Number of values recorded for each particle, it should be 11
            njsw  # Number of surfaces in JASW
            niss  # Number of different histories that reached the SSW surfaces
            niwr  # Number of cells in RSSA file
            mipts  # Source particle type
            kjaq  # Flag for macrobodies surfaces
            surfaces
        tracks: pl.DataFrame
            DataFrame containing the tracks recorded in the RSSA file.

            Each row of the table has 11 values
            0 a,  # History number of the particle, negative if uncollided
            1 b,  # Packed variable, the sign is the sign of the third direction cosine
                  # starts with 8 = neutron, 16 = photon
            2 wgt,
            3 erg,
            4 tme,
            5 x,
            6 y,
            7 z,
            8 u,  # Particle direction cosine with X-axis
            9 v,  # Particle direction cosine with Y-axis, to calculate w (Z-axis) use
                  # the sign from b
            10 c  # Surface id

        Examples
        --------
        >>> from f4enix.output.rssa import RSSA
        ... my_rssa = RSSA('small_cyl.w')
        ... print(my_rssa)
        RSSA file small_cyl.w was recorded using the following surfaces:
          Surface ID: 1, type: 1
        The total number of tracks recorded is 72083.
        Neutrons: 72083 photons: 0.
        The simulation that produced this RSSA run 100000 histories.
        The amount of independent histories that reached the RSSA surfaces was 70797.
        """
        self.path = Path(path)
        with open(path, "rb") as infile:
            self.parameters = _parse_header(infile)
            self.tracks = _parse_tracks(infile)

        # Modify the value of "b" for fast filtering of neutrons and photons
        self.tracks = self.tracks.with_columns(
            (pl.col("b").abs() / (10 ** pl.col("b").abs().log10().floor()))
            .cast(int)
            .alias("b")
        )

    def __repr__(self) -> str:
        return self.get_summary()

    def __str__(self) -> str:
        return self.get_summary()

    def get_summary(self) -> str:
        """Returns a summary of the RSSA file."""
        summary = f"RSSA file {self.path.name} was recorded using the following"
        summary += " surfaces:\n"
        for surface in self.parameters.surfaces:
            summary += f"  Surface ID: {surface.id}, type: {surface.type}\n"

        summary += f"The total number of tracks recorded is {self.parameters.nrss}.\n"
        summary += f"Neutrons: {self.neutron_tracks.shape[0]}"
        summary += f" photons: {self.photon_tracks.shape[0]}, "

        summary += "The simulation that produced this RSSA run "
        summary += f"{abs(self.parameters.np1)} histories.\n"
        summary += "The amount of independent histories that reached the RSSA surfaces "
        summary += f"was {self.parameters.niss}.\n"
        return summary

    @property
    def neutron_tracks(self) -> pl.DataFrame:
        """Returns the neutron tracks from the RSSA file."""
        return self.tracks.filter(pl.col("b") == NEUTRON_INDICATOR)

    @property
    def photon_tracks(self) -> pl.DataFrame:
        """Returns the photon tracks from the RSSA file."""
        return self.tracks.filter(pl.col("b") != NEUTRON_INDICATOR)

    @property
    def x(self) -> pl.Series:
        """Returns the x coordinates of the tracks."""
        return self.tracks["x"]

    @property
    def y(self) -> pl.Series:
        """Returns the y coordinates of the tracks."""
        return self.tracks["y"]

    @property
    def z(self) -> pl.Series:
        """Returns the z coordinates of the tracks."""
        return self.tracks["z"]

    @property
    def energies(self) -> pl.Series:
        """Returns the energies of the tracks."""
        return self.tracks["erg"]

    @property
    def wgt(self) -> pl.Series:
        """Returns the weights of the tracks."""
        return self.tracks["wgt"].abs()

    @property
    def histories(self) -> pl.Series:
        """Returns the history numbers of the tracks."""
        return self.tracks["a"]

    def plot_plane(self) -> "RSSAPlot":
        """Returns an instance of RSSAPlot to plot the RSSA data assumin an XY plane."""
        return RSSAPlot(self)

    def plot_cyl(self) -> "RSSAPlot":
        """Returns an instance of RSSAPlotCyl to plot the RSSA data asuming a
        cylindrical geometry with an axis following the Z-coordinate axis."""
        return RSSAPlot(self, x_col="perimeter_pos", y_col="z")


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


class RSSAPlot:
    def __init__(self, rssa: RSSA, x_col: str = "x", y_col: str = "y"):
        self.tracks = rssa.tracks.lazy()
        self.rssa_parameters = rssa.parameters
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

    def set_particle(self, particle_type: Literal["n", "p"]) -> "RSSAPlot":
        """Set the particle type to filter the tracks."""
        if particle_type == "n":
            self.tracks = self.tracks.filter(pl.col("b") == NEUTRON_INDICATOR)
        elif particle_type == "p":
            self.tracks = self.tracks.filter(pl.col("b") != NEUTRON_INDICATOR)
        return self

    def set_surface_ids(self, surface_ids: list[int]) -> "RSSAPlot":
        """Set the surface IDs to filter the tracks."""
        valid_surface_ids = [s.id for s in self.rssa_parameters.surfaces]
        if not all(sid in valid_surface_ids for sid in surface_ids):
            raise ValueError(
                f"Some surface IDs are not valid. Valid IDs are: {valid_surface_ids}"
            )
        self.tracks = self.tracks.filter(pl.col("c").is_in(surface_ids))
        return self

    def set_z_limits(self, vmin: float, vmax: float) -> "RSSAPlot":
        """Set the z limits for the plot."""
        self.tracks = self.tracks.filter(
            pl.col("z").is_between(vmin, vmax, closed="both")
        )
        return self

    def set_perimeter_limits(self, vmin: float, vmax: float) -> "RSSAPlot":
        """Set the limits for the perimeter positions."""
        if "perimeter_pos" not in self.tracks.collect_schema().names():
            self.calculate_perimeter_positions()
        self.tracks = self.tracks.filter(
            pl.col("perimeter_pos").is_between(vmin, vmax, closed="both")
        )
        return self

    def calculate_perimeter_positions(self) -> "RSSAPlot":
        radius = (pl.col("x").pow(2) + pl.col("y").pow(2)).sqrt().mean()
        thetas = pl.arctan2(pl.col("y"), pl.col("x"))
        perimeter_pos = (thetas * radius).alias("perimeter_pos")
        self.tracks = self.tracks.with_columns(perimeter_pos)
        return self

    def set_bins(self, x_bins: Sequence[float], y_bins: Sequence[float]) -> "RSSAPlot":
        """Set the x and y bins for the plot."""
        self._x_bins = pl.Series("x_bins", x_bins).sort()
        self._y_bins = pl.Series("y_bins", y_bins).sort()
        return self

    def calculate_bins(self, bin_width: float = 10.0) -> "RSSAPlot":
        """Automatically calculate the bins for the x and y coordinates of the plot by
        giving a bin width in cm. Instead use `set_bins()` to apply custom bins."""
        self._ensure_that_columns_are_set()
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

    def _ensure_that_columns_are_set(self) -> None:
        if (
            self.x_col == "perimeter_pos"
            and self.x_col not in self.tracks.collect_schema().names()
        ):
            self.calculate_perimeter_positions()

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

    def set_plot_parameters(self, plot_parameters: PlotParameters) -> "RSSAPlot":
        """Set the plot parameters for the plot."""
        self.plot_parameters = plot_parameters
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

    def save_figure(self, out_path: Path | str) -> "RSSAPlot":
        """Save the figure to the specified path."""
        _fig, _ax = self.get_plot()

        out_path = Path(out_path)
        if out_path.suffix != ".png":
            out_path = Path(out_path).with_suffix(".png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        return self

    def show(self) -> "RSSAPlot":
        """Show the plot."""
        _fig, _ax = self.get_plot()
        plt.show()
        return self

    def _get_2d_grid_of_weights(
        self,
        agg_func: Literal["sum", "count"] = "sum",
    ) -> np.ndarray:
        self._ensure_that_columns_are_set()

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


def _parse_header(infile: BinaryIO) -> _FileParameters:
    first_record = _read_fortran_record(infile)
    # The first line of the file with information like the code version, date and title
    formatted_record_id = first_record.tobytes().decode("UTF-8")
    if "d1suned" in formatted_record_id:
        _last_dump = np.frombuffer(first_record[-4:], INT)
    elif "SF_00001" in formatted_record_id:
        _header = _read_fortran_record(infile)  # code version and other info
    else:
        raise NotImplementedError(
            f"The code that generated this RSSA file has not been implemented"
            f" in this parser, see the code here: {formatted_record_id}..."
        )

    second_record = _read_fortran_record(infile)
    np1 = np.frombuffer(second_record, LONG, 1, 0)[0]
    nrss = np.frombuffer(second_record, LONG, 1, 8)[0]
    nrcd = np.frombuffer(second_record, INT, 1, 16)[0]
    njsw = np.frombuffer(second_record, INT, 1, 20)[0]
    niss = np.frombuffer(second_record, LONG, 1, 24)[0]
    if abs(nrcd) != NUMBER_OF_EXPECTED_VALUES:
        raise NotImplementedError(
            "The amount of values recorded for each particle should be 11 instead of"
            f" {nrcd}..."
        )

    if np1 < 0:
        third_record = _read_fortran_record(infile)
        niwr, mipts, kjaq = np.frombuffer(third_record, INT, 3)
    else:
        raise NotImplementedError("The np1 value is not negative...")

    surfaces = []
    for _ in range(njsw):
        data = _read_fortran_record(infile)
        surf_id = np.frombuffer(data, INT, 1, 0)[0]
        surf_info = np.frombuffer(data, INT, 1, 4)[0] if kjaq == 1 else -1
        surf_type = np.frombuffer(data, INT, 1, 8)[0]
        num_parameters = np.frombuffer(data, INT, 1, 12)[0]
        parameters = np.frombuffer(data, INT, offset=16).tolist()
        surfaces.append(
            _SurfaceParameters(
                id=surf_id,
                info=surf_info,
                type=surf_type,
                num_parameters=num_parameters,
                parameters=parameters,
            )
        )

    # we read any extra records as determined by njsw+niwr...
    # no known case of their actual utility
    for _j in range(njsw, njsw + niwr):
        _read_fortran_record(infile)
        raise NotImplementedError(
            "njsw + niwr values are bigger than njsw, behavior not explained"
        )

    # Summary record
    _data = _read_fortran_record(infile)
    # Summary record not processed, its information does not interest us for now

    return _FileParameters(
        np1=np1,  # Number of histories of the simulation, given as a negative number
        nrss=nrss,  # Number of tracks recorded
        nrcd=nrcd,  # Number of values recorded for each particle, it should be 11
        njsw=njsw,  # Number of surfaces in JASW
        niss=niss,  # Number of different histories that reached the SSW surfaces
        niwr=niwr,  # Number of cells in RSSA file
        mipts=mipts,  # Source particle type
        kjaq=kjaq,  # Flag for macrobodies surfaces
        surfaces=surfaces,
    )


def _parse_tracks(file: BinaryIO) -> pl.DataFrame:
    # Read the whole remaining of the file at once, store all the bytes as a 1D np array
    data = np.fromfile(file, BYTE)

    # Reshape the array so each index holds the information of a single particle
    # we can do this because we know that the particle records have always the same
    # length, 96 bytes
    data = data.reshape(-1, 96)

    # Remove the first and last 4 bytes, these are two integers that tell the record is
    # 88 bytes long
    data = data[:, 4:-4]

    # Convert the array into a 1D array of float numbers instead of simply bytes
    data = np.frombuffer(data.flatten(), FLOAT)

    # Reshape the array so each index holds the information of a single particle
    # all the data is already converted from bytes to floats
    data = data.reshape(-1, 11)

    return pl.DataFrame(
        data,
        schema={
            "a": int,
            "b": int,
            "wgt": float,
            "erg": float,
            "tme": float,
            "x": float,
            "y": float,
            "z": float,
            "u": float,
            "v": float,
            "c": int,
        },
    )


def _read_fortran_record(infile: BinaryIO):
    count_1 = np.fromfile(infile, INT, 1)[0]
    data = np.fromfile(infile, np.byte, count_1)
    count_2 = np.fromfile(infile, INT, 1)[0]
    if count_1 != count_2:
        raise ValueError(
            "The integers that go before and after the Fortran record are not equal..."
        )
    return data
