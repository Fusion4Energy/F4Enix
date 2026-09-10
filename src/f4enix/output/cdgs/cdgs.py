from enum import Enum
from io import TextIOWrapper
import logging
from matplotlib.pyplot import grid
import numpy as np
import pandas as pd
from pathlib import Path
import pyvista as pv
import actigamma as ag
from f4enix.output.cdgs import InterpolationKernel, RegularMeshDefinition
import f4enix
from datetime import datetime
from typing import Iterable

ACTIVITY_TAG = "_specific_activity [Bq/m3]"
INTENSITY_TAG = "_intensity [ph/s]"
ATOM_DENSITY_TAG = "_specific_atom_density [atom/m3]"


class CDGS_ENERGY_TYPE(Enum):
    LINE = "line"
    BINS = "bins"


class CDGS_MESH_TYPE(Enum):
    RECTANGULAR = "rec"
    # CYLINDRICAL = 'cyl'
    # UNSTRUCTURED = 'unstr'


class CDGS:
    def __init__(
        self,
        mesh: pv.StructuredGrid,
        energy_type: CDGS_ENERGY_TYPE = CDGS_ENERGY_TYPE.LINE,
        particle: str = "gamma",
        e_bins: np.ndarray | None = None,
    ):
        """CDGS object. Handles all operations with respect to .cdgs and .vtk format.
        Computes lines and emission probabilities for isotopes.

        Parameters
        ----------
        mesh : pv.StructuredGrid
            The mesh containing the activity and atom density data for isotopes.
        energy_type : CDGS_ENERGY_TYPE, optional
            The type of energy representation, by default CDGS_ENERGY_TYPE.LINE
        particle : str, optional
            The type of particle, by default "gamma". Other supported one is
            "neutron". Be sure that the radioisotopes you are using have the
            corresponding particle emission data in the actigamma database.
        e_bins : np.ndarray, optional
            The energy bins for the CDGS object, by default None.
            Mandatory if energy_type is CDGS_ENERGY_TYPE.BINS. Should be a 1D array of
            bin edges in eV.
        """

        # TODO: this may be extended to other type of geometries
        self.mesh: pv.StructuredGrid = mesh
        self.mesh_type: CDGS_MESH_TYPE = CDGS_MESH_TYPE.RECTANGULAR
        self.voxel_vol: float = np.min(np.abs(mesh.compute_cell_sizes()["Volume"]))

        self.energy_type: CDGS_ENERGY_TYPE = energy_type

        if particle not in ["gamma", "neutron"]:
            raise ValueError(
                f"Particle type {particle} not supported. Supported types are 'gamma' and 'neutron'."
            )
        self._particle: str = particle

        self._translation: np.ndarray = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self._cooling_time: float = 0.0

        isotopes = []
        for name in self.mesh.array_names:
            if ACTIVITY_TAG in name:
                isotopes.append(name.replace(ACTIVITY_TAG, ""))

        self.isotopes: dict[str, dict[str, np.ndarray]] = {}
        db = ag.Decay2012Database()
        self._db = db

        if energy_type == CDGS_ENERGY_TYPE.LINE:
            # Compute the lines for each isotope
            for isotope in isotopes:
                lines = db.getenergies(isotope, spectype=particle)  # eV
                intensities = db.getintensities(isotope, spectype=particle)
                self.isotopes[isotope] = {
                    "lines": np.array(lines),
                    "intensities": np.array(intensities),
                }

            # merge all lines and intensities for all isotopes
            to_concat = []
            for isotope in isotopes:
                df = pd.DataFrame()
                df["lines"] = self.isotopes[isotope]["lines"]
                df["intensities"] = self.isotopes[isotope]["intensities"]
                df["isotope"] = isotope
                to_concat.append(df)
            all_df = pd.concat(to_concat, ignore_index=True).sort_values(by="lines")
            self.all_df = all_df

        elif energy_type == CDGS_ENERGY_TYPE.BINS:
            e_bins = np.array(e_bins)  # be sure e_bins is a numpy array
            if e_bins is None:
                raise ValueError(
                    "Energy bins must be provided when energy_type is CDGS_ENERGY_TYPE.BINS."
                )
            self.e_bins = e_bins
            grid = ag.EnergyGrid(bounds=e_bins)
            self._lc = ag.LineAggregator(db, grid)

            # create an aggegated inventory for each isotope first
            for isotope in isotopes:
                inv = ag.UnstablesInventory(
                    # 1 Bq
                    data=[(db.getzai(isotope), 1)]
                )
                intensities, _ = self._lc(inv, spectype=self._particle)
                self.isotopes[isotope] = {"intensities": np.array(intensities)}

        # add the intensities to the mesh cell data for each isotope
        for isotope, vals in self.isotopes.items():
            values = (
                self.mesh.cell_data[f"{isotope}{ACTIVITY_TAG}"]  # Bq/m3
                * self.voxel_vol  # m3
                * np.sum(self.isotopes[isotope]["intensities"])  # photons/decay
            )  # photons/s
            self.mesh.cell_data[f"{isotope}{INTENSITY_TAG}"] = values

        # Store original mesh at shutdown
        self._mesh_zero_cooling = mesh.copy()

    def particle(self) -> str:
        return self._particle

    def __repr__(self) -> str:
        return f"CDGS(mesh={self.mesh}, energy_type={self.energy_type}, isotopes={list(self.isotopes.keys())}, particle={self._particle}, cooling_time={self.cooling_time})"

    def __str__(self) -> str:
        return f"CDGS object with mesh of dimensions {self.mesh.dimensions} and energy type {self.energy_type.value}"

    # @property
    # def translation(self) -> np.ndarray:
    #     """Get the translation matrix of the CDGS object.

    #     Returns
    #     -------
    #     np.ndarray
    #         The translation matrix of the CDGS object.
    #     """
    #     return self._translation

    # @translation.setter
    # def translation(self, value: np.ndarray) -> None:
    #     """Set the translation matrix of the CDGS object.

    #     Parameters
    #     ----------
    #     value : np.ndarray
    #         The new translation matrix to set.
    #     """
    #     if not isinstance(value, np.ndarray):
    #         raise TypeError("Translation must be a numpy ndarray.")
    #     if value.shape != (3, 3):
    #         raise ValueError("Translation matrix must be of shape (3, 3).")
    #     self._translation = value

    def to_vtk(self, outfile: str | Path) -> None:
        """Write the CDGS object to a VTK file.

        Parameters
        ----------
        outfile : str | Path
            Path to the output VTK file.
        """
        mesh = self.mesh.copy()
        # mesh.translate(self.translation[0] * self.translation[1] * self.translation[2])
        mesh.save(outfile)

    @property
    def cooling_time(self) -> float:
        """Get the cooling time of the CDGS object.

        Returns
        -------
        float
            The cooling time of the CDGS object.
        """
        return self._cooling_time

    @cooling_time.setter
    def cooling_time(self, value: float) -> None:
        """Set the cooling time of the CDGS object in seconds.
        This will update the mesh activities based on the cooling time.

        Parameters
        ----------
        value : float
            The new cooling time to set.
        """
        if not isinstance(value, (int, float)):
            raise TypeError("Cooling time must be a number.")
        if value < 0:
            raise ValueError("Cooling time cannot be negative.")
        self._cooling_time = float(value)

        if value == 0:
            self.mesh = self._mesh_zero_cooling.copy()
        else:
            self._update_cooling(value)

    def _update_cooling(self, cooling_time: float) -> None:
        """Update the mesh activities based on the cooling time."""
        mesh = self._mesh_zero_cooling.copy()
        for isotope in self.isotopes.keys():
            hl = self._db.gethalflife(isotope)  # seconds
            lambda_decay = np.log(2) / hl  # decay constant

            # Update the atom density and activity based on the cooling time
            for tag in [ACTIVITY_TAG, ATOM_DENSITY_TAG, INTENSITY_TAG]:
                mesh.cell_data[f"{isotope}{tag}"] = mesh.cell_data[
                    f"{isotope}{tag}"
                ] * np.exp(-lambda_decay * cooling_time)
        self.mesh = mesh

    def get_total_intensity(self, isotope: str) -> float:
        """Get the total intensity (photons/s) of the specified isotope in the mesh
        or of all isotopes together if "all" is specified.
        The total intensity is calculated as the sum of the specific activities in the
        mesh multiplied by the voxel volume and the sum of the gamma line
        intensities for the isotope.

        Parameters
        ----------
        isotope : str
            The name of the isotope. If "all" is specified, combine all isotopes.

        Returns
        -------
        float
            The total intensity (photons/s) of the specified isotope in the mesh.
        """
        if isotope == "all":
            return np.sum(
                [
                    np.sum(self.mesh.cell_data[f"{iso}{ACTIVITY_TAG}"])
                    * self.voxel_vol
                    * np.sum(self.isotopes[iso]["intensities"])
                    for iso in self.isotopes
                ]
            )
        if isotope not in self.isotopes:
            raise ValueError(f"Isotope {isotope} not found in the CDGS object.")
        return (
            np.sum(self.mesh.cell_data[f"{isotope}{ACTIVITY_TAG}"])
            * self.voxel_vol
            * np.sum(self.isotopes[isotope]["intensities"])
        )

    @classmethod
    def from_cloud_point(
        cls,
        cloud_point: str | Path | pd.DataFrame,
        isotopes: dict[str, str],
        interpolation_kernel: InterpolationKernel,
        mesh_definition: RegularMeshDefinition,
        col_names: dict[str, str] | None = None,
        e_bins: np.ndarray | None = None,
        particle: str = "gamma",
    ) -> "CDGS":
        """Read a csv (typically produced from fluent) that contains a cloud
        point of data and convert it into a CDGS object.

        Parameters
        ----------
        cloud_point : str | Path | pd.DataFrame
            Path to the CSV file containing the cloud point data or a DataFrame directly.
        isotopes : dict[str, str]
            Dictionary mapping isotope names to their column in the CSV file. Isotopes
            names should be like "Co60".
        interpolation_kernel : InterpolationKernel
            An instance of an InterpolationKernel subclass that defines how to distribute
            the activity from the cloud point to the regular mesh.
        mesh_definition : RegularMeshDefinition
            An instance of a RegularMeshDefinition subclass that defines the mesh structure.
        col_names : dict[str, str] | None, optional
            Dictionary mapping column names in the CSV to desired names in the CDGS
            object. If None, default names will be used. Default dictionary is:
            {
                "x": "x-coordinate",
                "y": "y-coordinate",
                "z": "z-coordinate",
                "vol": "cell-volume",
            }
        e_bins : np.ndarray, optional
            The energy bins for the CDGS object, by default None. If provided, the
            emission lines will be computed in the energy bins. Otherwise,
            pure emission lines will be used.
        particle : str, optional
            The type of particle, by default "gamma". Other supported one is
            "neutron". Be sure that the radioisotopes you are using have the
            corresponding particle emission data in the actigamma database.

        Returns
        -------
        CDGS
            A CDGS object populated with data from the CSV file.
        """
        if isinstance(cloud_point, pd.DataFrame):
            df = cloud_point.copy()
        else:
            df = pd.read_csv(cloud_point)
        df.columns = (
            df.columns.str.strip()
        )  # remove any leading/trailing whitespace from column names
        if col_names is None:
            col_names = {
                "x": "x-coordinate",
                "y": "y-coordinate",
                "z": "z-coordinate",
                "vol": "cell-volume",
            }
        coords = df[[col_names["x"], col_names["y"], col_names["z"]]].to_numpy()

        # Check if there are duplicated points, if so, drop them and warn the user.
        # pandas.duplicated() hashes rows instead of sorting them like
        # np.unique(axis=0) would, which is much cheaper for millions of points
        dup_cols = [col_names["x"], col_names["y"], col_names["z"]]
        if df.duplicated(subset=dup_cols).any():
            print(
                "Warning: Duplicated points found in the cloud point data. "
                "These will be dropped."
            )
            df = df.drop_duplicates(subset=dup_cols)
            coords = df[[col_names["x"], col_names["y"], col_names["z"]]].to_numpy()

        UM_vols = df[col_names["vol"]].to_numpy()

        # build the regular mesh
        mesh = mesh_definition._build_grid(coords)
        vol = np.min(np.abs(mesh.compute_cell_sizes()["Volume"]))
        # strip the pyvista_ndarray subclass: its __array_finalize__/__array_wrap__
        # hooks fire on every slice and add up over millions of kernel calls
        centroids_voxels = np.asarray(mesh.cell_centers().points)

        # Query the KDTree to find neighboring points
        neighbours_idx = interpolation_kernel.get_neighbours(
            source_coords=coords,
            dest_coords=centroids_voxels,
            voxel_size=vol ** (1 / 3),
        )

        # add the different isotopes activities
        db = ag.Decay2012Database()

        for isotope, col in isotopes.items():
            # initialize the activity array for each isotope and isotope counter
            activity_array = np.zeros(len(centroids_voxels))
            atoms_array = np.zeros(len(centroids_voxels))
            # extract once as a numpy array: df[col].iloc[i] rebuilds a Series
            # on every call and dominates runtime over millions of rows
            col_values = df[col].to_numpy()
            atoms_per_point = col_values * UM_vols  # atoms

            if getattr(interpolation_kernel, "supports_batch", False):
                # activity_from_atoms is a linear scaling of atoms, so it
                # broadcasts fine over the whole array in a single call
                activities_per_point = ag.activity_from_atoms(
                    db, isotope, atoms_per_point
                )  # Bq
                interpolation_kernel._kernel_batch(
                    coords,
                    neighbours_idx,
                    centroids_voxels,
                    activity_array,
                    activities_per_point,
                )
                interpolation_kernel._kernel_batch(
                    coords,
                    neighbours_idx,
                    centroids_voxels,
                    atoms_array,
                    atoms_per_point,
                )
            else:
                # Cycle on all UM centroids and distribute the activity
                for i in range(len(coords)):
                    source_coord = coords[i]
                    atoms = atoms_per_point[i]
                    activity = ag.activity_from_atoms(db, isotope, atoms)  # Bq
                    interpolation_kernel._kernel(
                        source_coord,
                        neighbours_idx[i],
                        centroids_voxels,
                        activity_array,
                        activity,
                    )
                    interpolation_kernel._kernel(
                        source_coord,
                        neighbours_idx[i],
                        centroids_voxels,
                        atoms_array,
                        atoms,
                    )
            mesh.cell_data[f"{isotope}{ACTIVITY_TAG}"] = activity_array / vol  # Bq/m3
            mesh.cell_data[f"{isotope}{ATOM_DENSITY_TAG}"] = (
                atoms_array / vol
            )  # atoms/m3

        if e_bins is not None:
            e_type = CDGS_ENERGY_TYPE.BINS
        else:
            e_type = CDGS_ENERGY_TYPE.LINE

        return cls(mesh, energy_type=e_type, particle=particle, e_bins=e_bins)

    def to_cdgs(self, outfile: str | Path, isotope: str) -> None:
        # Assume only one mesh for now
        num_mesh = 1
        tot_intensity = self.get_total_intensity(isotope)
        mesh_intensity = tot_intensity
        header = f"""
num_meshes {num_mesh}
global_source {tot_intensity:.5e}
mesh_id {1}
CDGS produced by F4Enix {f4enix.__version__} {datetime.now()}
cooling_time {self.cooling_time:.5e}
total_source {mesh_intensity:.5e}
energy_type {self.energy_type.value}
"""
        with open(outfile, "w") as f:
            f.write(header)

            # Energy boundaries
            self._write_energy_boundaries(f, isotope)

            # Mesh boundaries
            self._write_mesh_boundaries(f)

            # Values
            f.write("source_data\n")
            self._write_values(f, isotope)
            f.write("end_source_data\n")

        logging.info(f"CDGS file written to {outfile}")

    def _write_energy_boundaries(self, f: TextIOWrapper, isotope: str) -> None:
        """Write the energy boundaries section for the specified isotope to the file."""
        # Energy boundaries
        if self.energy_type == CDGS_ENERGY_TYPE.LINE:
            if isotope == "all":
                energies = self.all_df["lines"] * 1e-6  # convert eV to MeV
            else:
                energies = self.isotopes[isotope]["lines"] * 1e-6  # convert eV to MeV
            f.write(f"energy_boundaries {len(energies)}\n")
            f.write(_floats_to_multiline_string(energies))

        elif self.energy_type == CDGS_ENERGY_TYPE.BINS:
            f.write(f"energy_boundaries {len(self.e_bins)}\n")
            f.write(
                _floats_to_multiline_string(self.e_bins * 1e-6)
            )  # convert eV to MeV

        else:
            raise NotImplementedError(
                f"Energy type {self.energy_type} not implemented yet."
            )

    def _write_mesh_boundaries(self, f: TextIOWrapper) -> None:
        """Write the mesh type, boundaries, and translation matrix to the file."""
        # write the translation matrix
        f.write(f"mesh_type {self.mesh_type.value}\n")

        if self.mesh_type == CDGS_MESH_TYPE.RECTANGULAR:
            nx, ny, nz = self.mesh.dimensions
            f.write(f"mesh_boundaries {nx} {ny} {nz}\n")
            # Bin limits (grid edges) per direction
            edges1 = sorted(np.unique(self.mesh.points[:, 0]))
            edges2 = sorted(np.unique(self.mesh.points[:, 1]))
            edges3 = sorted(np.unique(self.mesh.points[:, 2]))
        else:
            raise NotImplementedError(
                f"Mesh type {self.mesh_type} not implemented yet."
            )

        # self._write_translation_matrix(f)

        f.write(_floats_to_multiline_string(edges1))
        f.write(_floats_to_multiline_string(edges2))
        f.write(_floats_to_multiline_string(edges3))

    def _write_values(self, f: TextIOWrapper, isotope: str) -> None:
        """Write the per-voxel source data values for the given isotope to the file."""
        if isotope == "all":
            tot_intensity = np.zeros(self.mesh.n_cells)
            for iso, _ in self.isotopes.items():
                tot_intensity += (
                    self.mesh.cell_data[f"{iso}{INTENSITY_TAG}"]  # ph/s
                )
        else:
            tot_intensity = self.mesh.cell_data[f"{isotope}{INTENSITY_TAG}"]  # ph/s

        for i, val in enumerate(tot_intensity):
            if val > 0:
                # printing element id, element emission intensity, element volume,
                # number of active cells under mesh element (assumed to be 1)
                f.write(f"{i} {val:.5e} {self.voxel_vol:.5e} 1\n")
                # Cell ID, volume fraction, intensity (ph/s).
                f.write(f"0 1.0 {val:.5e}\n")  # All homogenous

                probabilities = self._get_probabilities(isotope, i)
                intensities = probabilities * val  # ph/s
                f.write(_floats_to_multiline_string(intensities))
                # List of gamma source uncertainties is set to 0 for all values
                f.write(_floats_to_multiline_string([0] * len(intensities)))

    def _get_probabilities(self, isotope: str, idx: int) -> np.ndarray:
        """Return normalised per-line emission probabilities for a voxel and isotope."""
        if isotope == "all":
            if self.energy_type == CDGS_ENERGY_TYPE.LINE:
                df = self.all_df.copy()
                df["intensity"] = 0.0  # Initialize the intensity column
                # assign intensities to the dataframe using the isotope column
                for iso in self.isotopes.keys():
                    df.loc[df["isotope"] == iso, "intensity"] = self.mesh.cell_data[
                        f"{iso}{INTENSITY_TAG}"
                    ][idx]  # ph/s
                df["probability"] = df["intensity"] * df["intensities"]
                df["probability"] /= df["probability"].sum()  # Normalize to sum to 1
                return df["probability"].to_numpy()
            elif self.energy_type == CDGS_ENERGY_TYPE.BINS:
                data = []
                for iso in self.isotopes.keys():
                    data.append(
                        (
                            self._db.getzai(iso),
                            self.mesh.cell_data[f"{iso}{ACTIVITY_TAG}"][idx]
                            * self.voxel_vol,
                        )
                    )
                inv = ag.UnstablesInventory(data)
                hist, _ = self._lc(inv, spectype=self._particle)
                return hist / hist.sum()  # Normalize to sum to 1
            else:
                raise NotImplementedError(
                    f"Energy type {self.energy_type} not implemented yet."
                )
        else:
            prob = self.isotopes[isotope]["intensities"]
            return prob / prob.sum()  # Normalize to sum to 1

    # def _write_translation_matrix(self, f: TextIOWrapper) -> None:
    #     """Write the three rows of the translation matrix to the file."""
    #     f.write(_floats_to_multiline_string(self.translation[0, :]))
    #     f.write(_floats_to_multiline_string(self.translation[1, :]))
    #     f.write(_floats_to_multiline_string(self.translation[2, :]))


def _floats_to_multiline_string(
    values: Iterable[float], max_line_length: int = 100
) -> str:
    """Format an iterable of floats as a space-separated, line-wrapped string."""
    lines = []
    current_line = []

    for value in values:
        formatted_value = f"{value:.7e}"
        # Check if adding the next value exceeds the line length limit
        if (
            sum(len(v) for v in current_line)
            + len(formatted_value)
            + len(current_line)
            - 1
            < max_line_length
        ):
            current_line.append(formatted_value)
        else:
            # Join the current line and add it to the lines list
            lines.append(" ".join(current_line))
            current_line = [formatted_value]

    # Add the last line if it's not empty
    if current_line:
        lines.append(" ".join(current_line))

    return "\n".join(lines) + "\n"
