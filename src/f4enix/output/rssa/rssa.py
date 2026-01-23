"""This module is related to the parsing of D1S-UNED meshinfo files."""

from collections.abc import Sequence
from pathlib import Path

import polars as pl

from f4enix.output.rssa.rssa_helpers import NEUTRON_INDICATOR
from f4enix.output.rssa.rssa_plotting import RSSAPlot, RSSASpectraPlot
from f4enix.output.rssa.rssa_reader import parse_header, parse_tracks


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
            self.parameters = parse_header(infile)
            self.tracks = parse_tracks(infile)

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

    def get_energy_spectra(self, energy_bins: Sequence[float]) -> pl.DataFrame:
        raise NotImplementedError()

    def plot_plane(self) -> RSSAPlot:
        """Returns an instance of RSSAPlot to plot data assuming an XY plane."""
        return RSSAPlot(self.tracks, self.parameters)

    def plot_cyl(self) -> RSSAPlot:
        """Returns an instance of RSSAPlotCyl to plot data asuming a cylindrical
        geometry with an axis following the Z-coordinate axis."""
        rssa_plot = RSSAPlot(
            self.tracks, self.parameters, x_col="perimeter_pos", y_col="z"
        )
        rssa_plot.calculate_perimeter_positions()
        return rssa_plot

    def plot_spectra(self) -> RSSASpectraPlot:
        return RSSASpectraPlot(self.tracks, self.parameters)
