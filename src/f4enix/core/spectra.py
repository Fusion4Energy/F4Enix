import numpy as np
from f4enix.core.constants import PathLike
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import logging

LINESTYLES = ["--", "-.", ":"] * 20


class SpectraParsingError(Exception):
    """Custom error for handling parsing issues in Spectra-related operations."""

    pass


class Spectra:
    def __init__(
        self, ebins: np.ndarray, spectra_values: np.ndarray, name: str | None = None
    ) -> None:
        """Class representing a neutron spectra. Values are expcted as normalized to 1.

        Parameters
        ----------
        ebins : np.ndarray
            N energy bins including first and last
        spectra_values : np.ndarray
            N-1 spectra values, one for each energy bin (in eV)
        name : str, None
            Name of the spectra. "Unnamed_Spectra" is set if None.

        Attributes
        ----------
        ebins : np.ndarray
            N energy bins including first and last
        spectra_values : np.ndarray
            N-1 spectra values, one for each energy bin
        name : str
            Name of the spectra.

        Raises
        ------
        ValueError
            If the spectra values are not normalized to 1
        ValueError
            If the spectra values length does not match the energy bins length minus 1
        """
        if name is None:
            name = "Unnamed_Spectra"

        # Spectra must be normalized to 1
        if not np.isclose(np.sum(spectra_values), 1.0, atol=5e-3):
            logging.warning(
                f"Spectra values must be normalized to 1. Current sum: {np.sum(spectra_values)}"
                " Normalizing automatically."
            )
            spectra_values = spectra_values / np.sum(spectra_values)

        if not len(ebins) - 1 == len(spectra_values):
            raise ValueError(
                f"Spectra values length {len(spectra_values)} does not match "
                f"energy bins length {len(ebins)-1}."
            )
        self.ebins = ebins  # bin edges including lowest and highest
        self.spectra_values = spectra_values  # values per bin
        self.name = name

    def to_fispact_fluxes(self, outdir: PathLike, format: str | None = None) -> None:
        """Save the spectra to a file in fispact "fluxes" format.

        Parameters
        ----------
        outdir : PathLike
            Path to the output directory.
        format : str, None
            format for numbers. If None, default '.3e' is used.
        """
        if format is None:
            format = ".3e"
        with open(Path(outdir, self.name), "w") as f:
            for flux in np.flip(self.spectra_values):
                f.write(f"{flux:{format}}\n")
            # Add dummy wall load
            f.write(f"{1:{format}}\n")
            # add flux name
            f.write(f"{self.name}\n")

    def print_collapse_inp(self, outdir: PathLike) -> Path:
        """Print a fispact collapse input file to perform XS collapse with this spectra.

        Parameters
        ----------
        outdir : PathLike
            Path to the output directory.

        Returns
        -------
        Path
            Path to the generated input file.
        """

        text = f"""MONITOR 1
CLOBBER
GETXS 1 {len(self.ebins) - 1}
FISPACT
* COLLAPSE 
END
* END OF RUN
"""
        outfile = Path(outdir, "collapse.i")
        with open(outfile, "w") as f:
            f.write(text)
        return outfile

    def print_condensed_inp(self, outdir: PathLike) -> Path:
        """Print a fispact condense input file to perform XS condensing with this
        spectra.

        Parameters
        ----------
        outdir : PathLike
            Path to the output directory.
        Returns
        -------
        Path
            Path to the generated input file.
        """
        text = f"""MONITOR 1
CLOBBER
SPEK
GETDECAY 1
FISPACT
* CONDENSE 
END
* END OF RUN
"""
        outfile = Path(outdir, "condense.i")
        with open(outfile, "w") as f:
            f.write(text)
        return outfile

    def to_arb_flux_file(self, outdir: PathLike) -> Path:
        """print the spectra to a fispact "arb_flux" format to be used in a spectra
        conversion

        Parameters
        ----------
        outdir : Path
            _description_

        Returns
        -------
        Path
            _description_
        """
        outfile = Path(outdir, f"arb_flux_{self.name}.txt")
        with open(outfile, "w") as f:
            # first write the energies in formatted columns of 6
            for i, e in enumerate(np.flip(self.ebins)):
                f.write(f"{e:.6e} ")
                if (i + 1) % 6 == 0:
                    f.write("\n")
            # if there is no new line at the end, add it
            if (i + 1) % 6 != 0:
                f.write("\n")
            # then write the flux values in formatted columns of 6
            for i, flux in enumerate(np.flip(self.spectra_values)):
                f.write(f"{flux:.6e} ")
                if (i + 1) % 6 == 0:
                    f.write("\n")
            # if there is no new line at the end, add it
            if (i + 1) % 6 != 0:
                f.write("\n")
            f.write(f"1\n")  # must have a wall load, likely not used
            f.write(f"{self.name}\n")

        return outfile

    @classmethod
    def from_fispact(cls, ebins: np.ndarray, file: PathLike) -> "Spectra":
        """Read a "fluxes" fispact file and get a Spectra object. Wall load must always
        be present (even if it will be ignored by F4Enix). A name also needs to be
        provided. This reduces possible wrong binning errors by the user.

        Parameters
        ----------
        ebins : np.ndarray
            ebins have to be externally provided. It is recommended to use the ones
            available at f4enix.core.egroups
        file : PathLike
            Path to the fispact fluxes file to read.

        Returns
        -------
        Spectra
            Spectra object read from the file.
        """
        # Read and flatten all tokens (handles multi-column and single-column)
        lines = []
        with open(file, "r") as f:
            for line in f:
                lines.append(line)

        # the last non empty line should be the title
        removal_counter = 0
        for line in reversed(lines):
            removal_counter += 1
            if line.strip():
                if len(line) > 1:
                    title_line = line
                    break

        # replace eventual spaces in the title
        title_line = title_line.strip().replace(" ", "_")

        # now go for the tokens
        tokens = []
        for line in lines[:-removal_counter]:
            tokens.extend(line.strip().split())

        n_bins = len(ebins) - 1

        if not len(tokens) == n_bins + 1:
            raise SpectraParsingError(
                f"File {file} does not contain the expected number of tokens "
                f"for the provided energy bins. Possible mismatch."
            )

        flux_values = np.flip(np.array(tokens[:n_bins], dtype=float))

        return Spectra(ebins=ebins, spectra_values=flux_values, name=title_line)

    def get_by_lethargy(self) -> np.ndarray:
        """Get converted spectra values by unit lethargy."""
        # Energies for lethargy computation
        flux = self.spectra_values / np.log(self.ebins[1:] / self.ebins[:-1])

        return flux

    def plot(
        self, lethargy: bool = False, add_spectra: list["Spectra"] | None = None
    ) -> tuple[Figure, Axes]:
        """Plot the spectra.

        Parameters
        ----------
        lethargy : bool, optional
            If True, plot the spectra by unit lethargy. Default is False.
        add_spectra : list[Spectra] | None, optional
            Additional spectra to plot for comparison.
            If None, only the current spectra is plotted. Default is None.

        Returns
        -------
        tuple[Figure, Axes]
            The matplotlib Figure and Axes objects containing the plot.
        """

        if lethargy:
            ylabel = "Spectra by unit lethargy"
        else:
            ylabel = "Spectra"

        fig, ax = plt.subplots()

        to_plot = [self]
        if add_spectra is not None:
            to_plot = [self] + add_spectra

        for i, spec in enumerate(to_plot):
            if lethargy:
                values = spec.get_by_lethargy()
            else:
                values = spec.spectra_values
            if len(to_plot) > 1:
                linestyle = LINESTYLES[i]
            else:
                linestyle = "-"
            ax.step(
                spec.ebins[1:] * 1e-6,
                values,
                label=spec.name,
                where="pre",
                linestyle=linestyle,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Energy (MeV)")
        ax.set_ylabel(ylabel)
        ax.grid()
        if add_spectra is not None:
            ax.set_title("Spectra comparison")
        else:
            ax.set_title(f"Spectrum: {self.name}")
        ax.legend()

        return fig, ax
