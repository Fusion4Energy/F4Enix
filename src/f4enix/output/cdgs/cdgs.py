from enum import Enum
from io import TextIOWrapper
import logging
import numbers
import numpy as np
import pandas as pd
from pathlib import Path
import pyvista as pv
import actigamma as ag
from matplotlib import pyplot as plt
from f4enix.output.cdgs import (
    InterpolationKernel,
    KNearestKernel,
    MeshByBinSize,
    RegularMeshDefinition,
)
import f4enix
from datetime import datetime
from typing import Iterable

ACTIVITY_TAG = "_specific_activity [Bq/m3]"
INTENSITY_TAG = "_intensity [ph/s]"
ATOM_DENSITY_TAG = "_specific_atom_density [atom/m3]"

# conversion factors to m for the supported input length units. Internally,
# CDGS works only in SI units
LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3}
# MCNP works in cm
M_TO_CM = 1e2


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

        All quantities are handled internally in SI units (m, s, Bq). Energies
        are in eV, following the actigamma convention.

        Parameters
        ----------
        mesh : pv.StructuredGrid
            The mesh containing the activity and atom density data for isotopes.
            Mesh coordinates must be in m, activities in Bq/m3 and atom
            densities in atom/m3.
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
        # copy so that the intensity arrays added below do not leak into the
        # caller's mesh
        self.mesh: pv.StructuredGrid = mesh.copy()
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
        if len(isotopes) == 0:
            raise ValueError(
                f"No isotope activity arrays (named '<isotope>{ACTIVITY_TAG}') "
                "found in the mesh cell data."
            )

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
            if e_bins is None:
                raise ValueError(
                    "Energy bins must be provided when energy_type is CDGS_ENERGY_TYPE.BINS."
                )
            e_bins = np.array(e_bins)  # be sure e_bins is a numpy array
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

                # lines falling outside e_bins are dropped by the aggregator,
                # so warn the user about the emission that is not accounted for
                full = np.sum(db.getintensities(isotope, spectype=particle))
                binned = np.sum(intensities)
                if full > 0 and binned < full * (1 - 1e-6):
                    logging.warning(
                        f"{isotope}: {(1 - binned / full) * 100:.2f}% of the "
                        f"{particle} emission ({full - binned:.4e} per decay) "
                        "falls outside the provided energy bins and will be "
                        "ignored."
                    )

        # add the intensities to the mesh cell data for each isotope
        for isotope, vals in self.isotopes.items():
            values = (
                self.mesh.cell_data[f"{isotope}{ACTIVITY_TAG}"]  # Bq/m3
                * self.voxel_vol  # m3
                * np.sum(self.isotopes[isotope]["intensities"])  # photons/decay
            )  # photons/s
            self.mesh.cell_data[f"{isotope}{INTENSITY_TAG}"] = values

        # Store original mesh at shutdown
        self._mesh_zero_cooling = self.mesh.copy()

    @property
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
        # numbers.Real also covers numpy scalars (e.g. np.int64)
        if not isinstance(value, numbers.Real) or isinstance(value, bool):
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
        length_unit: str = "m",
    ) -> "CDGS":
        """Read a csv (typically produced from fluent) that contains a cloud
        point of data and convert it into a CDGS object.

        Parameters
        ----------
        cloud_point : str | Path | pd.DataFrame
            Path to the CSV file containing the cloud point data or a DataFrame directly.
        isotopes : dict[str, str]
            Dictionary mapping isotope names to their column in the CSV file. Isotopes
            names should be like "Co60". Each column must contain the atom density
            of the isotope, in atoms per cubic length_unit (e.g. atom/cm3 if
            length_unit is "cm").
        interpolation_kernel : InterpolationKernel
            An instance of an InterpolationKernel subclass that defines how to distribute
            the activity from the cloud point to the regular mesh. Any length
            parameter of the kernel (e.g. radius) must be given in m.
        mesh_definition : RegularMeshDefinition
            An instance of a RegularMeshDefinition subclass that defines the mesh
            structure. Any length parameter of the definition (e.g. bin_size,
            bbox) must be given in m.
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
            The energy bin edges in eV for the CDGS object, by default None. If
            provided, the emission lines will be computed in the energy bins.
            Otherwise, pure emission lines will be used.
        particle : str, optional
            The type of particle, by default "gamma". Other supported one is
            "neutron". Be sure that the radioisotopes you are using have the
            corresponding particle emission data in the actigamma database.
        length_unit : str, optional
            Length unit used in the cloud point data, one of "m", "cm" or "mm",
            by default "m". Coordinates are expected in length_unit, cell
            volumes in length_unit^3 and atom densities in atom/length_unit^3.
            Everything is converted to SI (m) internally.

        Returns
        -------
        CDGS
            A CDGS object populated with data from the CSV file.
        """
        df, coords, UM_vols, to_m = _read_cloud_point(
            cloud_point, col_names, length_unit
        )

        # extract once as numpy arrays: df[col].iloc[i] rebuilds a Series
        # on every call and dominates runtime over millions of rows
        atoms_per_point = {
            isotope: df[col].to_numpy() / to_m**3 * UM_vols  # atom/m3 * m3 = atoms
            for isotope, col in isotopes.items()
        }
        mesh, atoms, vol = _interpolate_cloud_point(
            coords, atoms_per_point, interpolation_kernel, mesh_definition
        )

        # add the different isotopes activities
        db = ag.Decay2012Database()
        for isotope, atoms_array in atoms.items():
            # activity_from_atoms is a linear scaling of atoms, so interpolating
            # the atoms once and converting afterwards is equivalent to
            # interpolating the activity separately
            activity_array = ag.activity_from_atoms(db, isotope, atoms_array)  # Bq
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
        """Write the CDGS object to a .cdgs file readable by MCNP.

        Parameters
        ----------
        outfile : str | Path
            Path to the output .cdgs file.
        isotope : str
            The isotope to write. If "all" is specified, combine all isotopes.
        """
        MultiCDGS([self]).to_cdgs(outfile, isotope)

    def _write_mesh_block(
        self, f: TextIOWrapper, isotope: str, mesh_id: int, mesh_intensity: float
    ) -> None:
        """Write the self-contained block describing this mesh to the file."""
        header = f"""mesh_id {mesh_id}
CDGS produced by F4Enix {f4enix.__version__} {datetime.now()}
cooling_time {self.cooling_time:.5e}
total_source {mesh_intensity:.5e}
energy_type {self.energy_type.value}
"""
        f.write(header)

        # Energy boundaries
        self._write_energy_boundaries(f, isotope)

        # Mesh boundaries
        self._write_mesh_boundaries(f)

        # Values
        f.write("source_data\n")
        self._write_values(f, isotope)
        f.write("end_source_data\n")

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
            # Bin limits (grid edges) per direction, converted from m to cm
            # as the file is read by MCNP
            edges1 = np.unique(self.mesh.points[:, 0]) * M_TO_CM
            edges2 = np.unique(self.mesh.points[:, 1]) * M_TO_CM
            edges3 = np.unique(self.mesh.points[:, 2]) * M_TO_CM
        else:
            raise NotImplementedError(
                f"Mesh type {self.mesh_type} not implemented yet."
            )

        # self._write_translation_matrix(f)

        f.write(_floats_to_multiline_string(edges1))
        f.write(_floats_to_multiline_string(edges2))
        f.write(_floats_to_multiline_string(edges3))

    def _write_values(
        self, f: TextIOWrapper, isotope: str, chunk_size: int = 10000
    ) -> None:
        """Write the per-voxel source data values for the given isotope to the file."""
        if isotope == "all":
            tot_intensity = np.zeros(self.mesh.n_cells)
            for iso, _ in self.isotopes.items():
                tot_intensity += (
                    self.mesh.cell_data[f"{iso}{INTENSITY_TAG}"]  # ph/s
                )
        else:
            tot_intensity = self.mesh.cell_data[f"{isotope}{INTENSITY_TAG}"]  # ph/s

        tot_intensity = np.asarray(tot_intensity)
        voxel_vol_cm3 = self.voxel_vol * M_TO_CM**3  # MCNP works in cm
        active_idx = np.flatnonzero(tot_intensity > 0)
        # probabilities are computed for a chunk of voxels at a time: a single
        # (n_voxels, n_lines) matrix for the whole mesh could exhaust memory
        for start in range(0, len(active_idx), chunk_size):
            chunk = active_idx[start : start + chunk_size]
            vals = tot_intensity[chunk]
            probabilities = self._get_probabilities(isotope, chunk)
            intensities = probabilities * vals[:, None]  # ph/s
            zeros = _floats_to_multiline_string([0] * intensities.shape[1])

            for i, val, voxel_intensities in zip(chunk, vals, intensities):
                # printing element id, element emission intensity, element volume,
                # number of active cells under mesh element (assumed to be 1)
                f.write(f"{i} {val:.5e} {voxel_vol_cm3:.5e} 1\n")
                # Cell ID, volume fraction, intensity (ph/s).
                f.write(f"0 1.0 {val:.5e}\n")  # All homogenous

                f.write(_floats_to_multiline_string(voxel_intensities))
                # List of gamma source uncertainties is set to 0 for all values
                f.write(zeros)

    def _get_probabilities(
        self, isotope: str, idx: int | np.ndarray
    ) -> np.ndarray:
        """Return normalised per-line (or per-bin) emission probabilities for the
        given isotope in the voxel(s) idx.

        If idx is an int the result is a 1D array of shape (n_lines,), if it is
        an array of voxel indices the result has shape (len(idx), n_lines).
        """
        if isotope != "all":
            # spectrum shape of a single isotope does not depend on the voxel
            prob = self.isotopes[isotope]["intensities"]
            prob = prob / prob.sum()  # Normalize to sum to 1
            return np.broadcast_to(prob, (*np.shape(idx), len(prob)))

        if self.energy_type == CDGS_ENERGY_TYPE.LINE:
            # emission rate of each line is activity * line intensity (per decay)
            # of the isotope it belongs to
            iso_order = {iso: i for i, iso in enumerate(self.isotopes)}
            line_iso = self.all_df["isotope"].map(iso_order).to_numpy()
            per_decay = self.all_df["intensities"].to_numpy()
            weights = self._get_activities(idx)[..., line_iso] * per_decay
        elif self.energy_type == CDGS_ENERGY_TYPE.BINS:
            # self.isotopes[iso]["intensities"] is the binned spectrum for 1 Bq
            # and the aggregation is linear in the activity, so the spectrum of
            # the mixture is the activity-weighted sum of the isotope spectra
            per_decay = np.stack(
                [vals["intensities"] for vals in self.isotopes.values()]
            )  # (n_isotopes, n_bins)
            weights = self._get_activities(idx) @ per_decay
        else:
            raise NotImplementedError(
                f"Energy type {self.energy_type} not implemented yet."
            )
        # Normalize to sum to 1. The voxel volume would cancel out here, so
        # activities (Bq/m3) are used directly
        return weights / weights.sum(axis=-1, keepdims=True)

    def _get_activities(self, idx: int | np.ndarray) -> np.ndarray:
        """Return the specific activities (Bq/m3) of all isotopes in the voxel(s)
        idx, with shape (n_isotopes,) or (len(idx), n_isotopes)."""
        return np.stack(
            [
                np.asarray(self.mesh.cell_data[f"{iso}{ACTIVITY_TAG}"])[idx]
                for iso in self.isotopes
            ],
            axis=-1,
        )

    # def _write_translation_matrix(self, f: TextIOWrapper) -> None:
    #     """Write the three rows of the translation matrix to the file."""
    #     f.write(_floats_to_multiline_string(self.translation[0, :]))
    #     f.write(_floats_to_multiline_string(self.translation[1, :]))
    #     f.write(_floats_to_multiline_string(self.translation[2, :]))


class MultiCDGS:
    def __init__(self, cdgs_list: Iterable[CDGS], check_overlap: bool = True):
        """Collection of CDGS objects to be written as a single multi-mesh .cdgs
        file. Each CDGS object is written as an independent mesh block and the
        global source is the sum of the total intensities of all meshes.

        Each mesh block carries its own energy boundaries, so particle,
        energy type and energy bins can differ across meshes. All CDGS objects
        must share the same cooling time.

        Parameters
        ----------
        cdgs_list : Iterable[CDGS]
            The CDGS objects to combine, one per mesh.
        check_overlap : bool, optional
            Warn if the bounding boxes of the meshes overlap, by default True.
            Overlapping meshes are fine only if they were built from disjoint
            parts of the cloud point, otherwise the source in the overlap is
            counted twice.
        """
        self.cdgs_list: list[CDGS] = list(cdgs_list)
        if len(self.cdgs_list) == 0:
            raise ValueError("At least one CDGS object must be provided.")
        self._check_consistency()
        if check_overlap:
            self._check_overlap()

    def __len__(self) -> int:
        return len(self.cdgs_list)

    def __getitem__(self, idx: int) -> CDGS:
        return self.cdgs_list[idx]

    def __repr__(self) -> str:
        return f"MultiCDGS(n_meshes={len(self)}, isotopes={self.isotopes}, cooling_time={self.cooling_time})"

    def __str__(self) -> str:
        return f"MultiCDGS object with {len(self)} meshes"

    def _check_consistency(self) -> None:
        """Check that all the CDGS objects can be written in the same file."""
        ref = self.cdgs_list[0]
        for i, cdgs in enumerate(self.cdgs_list[1:], start=1):
            if cdgs.cooling_time != ref.cooling_time:
                raise ValueError(
                    f"Mesh {i} has cooling time {cdgs.cooling_time} s while "
                    f"mesh 0 has {ref.cooling_time} s."
                )

    def _check_overlap(self) -> None:
        """Warn if the bounding boxes of the meshes overlap."""
        bounds = [np.array(cdgs.mesh.bounds).reshape(3, 2) for cdgs in self.cdgs_list]
        for i in range(len(bounds)):
            for j in range(i + 1, len(bounds)):
                low = np.maximum(bounds[i][:, 0], bounds[j][:, 0])
                high = np.minimum(bounds[i][:, 1], bounds[j][:, 1])
                if np.all(high > low):
                    logging.warning(
                        f"Meshes {i} and {j} overlap: unless they were built "
                        "from disjoint parts of the cloud point, the source in "
                        "the overlapping region will be counted twice."
                    )

    @property
    def isotopes(self) -> list[str]:
        """All the isotopes found in at least one of the meshes."""
        return list(dict.fromkeys(iso for c in self.cdgs_list for iso in c.isotopes))

    @property
    def cooling_time(self) -> float:
        """Get the cooling time (s), shared by all the meshes."""
        return self.cdgs_list[0].cooling_time

    @cooling_time.setter
    def cooling_time(self, value: float) -> None:
        """Set the cooling time (s) of all the meshes."""
        for cdgs in self.cdgs_list:
            cdgs.cooling_time = value

    def get_total_intensity(self, isotope: str) -> float:
        """Get the total intensity (photons/s) of the specified isotope summed
        over all the meshes, or of all isotopes together if "all" is specified.

        Parameters
        ----------
        isotope : str
            The name of the isotope. If "all" is specified, combine all isotopes.

        Returns
        -------
        float
            The total intensity (photons/s) over all the meshes.
        """
        if isotope != "all" and isotope not in self.isotopes:
            raise ValueError(f"Isotope {isotope} not found in any of the meshes.")
        return float(np.sum(self._get_mesh_intensities(isotope)))

    def _get_mesh_intensities(self, isotope: str) -> list[float]:
        """Total intensity (photons/s) of each mesh, 0 if the isotope is not
        present in the mesh."""
        return [
            cdgs.get_total_intensity(isotope)
            if isotope == "all" or isotope in cdgs.isotopes
            else 0.0
            for cdgs in self.cdgs_list
        ]

    def to_cdgs(self, outfile: str | Path, isotope: str) -> None:
        """Write all the meshes to a single .cdgs file readable by MCNP.
        Meshes with no source for the requested isotope are not written.

        Parameters
        ----------
        outfile : str | Path
            Path to the output .cdgs file.
        isotope : str
            The isotope to write. If "all" is specified, combine all isotopes.
        """
        if isotope != "all" and isotope not in self.isotopes:
            raise ValueError(f"Isotope {isotope} not found in any of the meshes.")

        to_write = []
        for i, (cdgs, intensity) in enumerate(
            zip(self.cdgs_list, self._get_mesh_intensities(isotope))
        ):
            if intensity > 0:
                to_write.append((cdgs, intensity))
            else:
                logging.warning(
                    f"Mesh {i} has no {isotope} source and will not be written."
                )
        if len(to_write) == 0:
            raise ValueError(f"No source found for {isotope} in any of the meshes.")

        global_source = sum(intensity for _, intensity in to_write)
        header = f"""
num_meshes {len(to_write)}
global_source {global_source:.5e}
"""
        with open(outfile, "w") as f:
            f.write(header)
            for mesh_id, (cdgs, intensity) in enumerate(to_write, start=1):
                cdgs._write_mesh_block(f, isotope, mesh_id, intensity)

        logging.info(f"CDGS file written to {outfile}")

    @classmethod
    def from_cloud_point(
        cls,
        cloud_point: str | Path | pd.DataFrame,
        isotopes: dict[str, str],
        col_names: dict[str, str] | None = None,
        e_bins: np.ndarray | None = None,
        particle: str = "gamma",
        length_unit: str = "m",
        split_size: float | None = None,
        min_fill: float = 0.5,
        min_gain: float = 0.1,
        min_points: int = 1000,
        reference_column: str | None = None,
        parameters: dict | None = None,
        recommend_kwargs: dict | None = None,
    ) -> "MultiCDGS":
        """Split a cloud point into spatially compact chunks and build a CDGS
        object for each of them, with the interpolation parameters recommended
        by get_recommended_parameters for that chunk.

        Sparse cloud points (e.g. the fluid in a piping system) fill a tiny
        fraction of their bounding box, so a single mesh would be mostly empty.
        The space is divided in coarse cubes of side split_size and the box
        enclosing the occupied cubes is recursively bisected with an axis
        aligned plane, choosing the plane that minimises the total volume of the
        two resulting boxes. A box is not split further if at least min_fill of
        its coarse cubes are occupied, if the best split reduces its volume by
        less than min_gain, or if a split would leave fewer than min_points
        cloud points on one side.

        Each cloud point belongs to exactly one chunk, so the meshes can overlap
        (e.g. in the padding added around each chunk) without counting any
        source twice.

        Parameters
        ----------
        cloud_point : str | Path | pd.DataFrame
            Path to the CSV file containing the cloud point data or a DataFrame
            directly. See CDGS.from_cloud_point.
        isotopes : dict[str, str]
            Dictionary mapping isotope names to their column in the cloud
            point. See CDGS.from_cloud_point.
        col_names : dict[str, str] | None, optional
            Coordinates and volume column names, see CDGS.from_cloud_point.
        e_bins : np.ndarray, optional
            The energy bin edges in eV, see CDGS.from_cloud_point.
        particle : str, optional
            The type of particle, by default "gamma".
        length_unit : str, optional
            Length unit of the cloud point data, one of "m", "cm" or "mm", by
            default "m".
        split_size : float | None, optional
            Side in m of the coarse cubes used to detect the empty space. By
            default, 10 times the characteristic cell size of the cloud point
            (see get_recommended_parameters). Smaller values follow the
            geometry more closely, producing more and tighter boxes.
        min_fill : float, optional
            Fraction of occupied coarse cubes above which a box is not split
            further, by default 0.5.
        min_gain : float, optional
            Minimum fractional reduction of the box volume for a split to be
            accepted, by default 0.1.
        min_points : int, optional
            Minimum number of cloud points in each chunk, by default 1000.
        reference_column : str | None, optional
            Column used to recommend the interpolation parameters of each chunk.
            By default, the column of the first isotope.
        parameters : dict | None, optional
            {"bin_size": bin size in m, "k": number of neighbours} to be used for
            all the chunks instead of the recommended ones.
        recommend_kwargs : dict | None, optional
            Additional keyword arguments for get_recommended_parameters. If
            slice_axis is not given, the section of each chunk is taken normal
            to its shortest side, where it contains the most points.

        Returns
        -------
        MultiCDGS
            One CDGS object per chunk.
        """
        # avoid a circular import, sensitivity depends on this module
        from f4enix.output.cdgs.sensitivity import (
            AXES,
            _characteristic_cell_size,
            get_recommended_parameters,
        )

        df, coords, vols, _ = _read_cloud_point(cloud_point, col_names, length_unit)
        if split_size is None:
            split_size = 10 * _characteristic_cell_size(vols)
        if reference_column is None:
            reference_column = next(iter(isotopes.values()))
        if recommend_kwargs is None:
            recommend_kwargs = {}

        chunks = _split_cloud_point(coords, split_size, min_fill, min_gain, min_points)
        logging.info(f"Cloud point split in {len(chunks)} chunks.")

        cdgs_list = []
        for i, idx in enumerate(chunks):
            df_chunk = df.iloc[idx]
            if parameters is not None:
                bin_size, k = parameters["bin_size"], parameters["k"]
            else:
                kwargs = dict(recommend_kwargs)
                if "slice_axis" not in kwargs:
                    extent = np.ptp(coords[idx], axis=0)
                    kwargs["slice_axis"] = list(AXES)[int(np.argmin(extent))]
                rec, metrics, fig = get_recommended_parameters(
                    df_chunk,
                    reference_column,
                    col_names=col_names,
                    length_unit=length_unit,
                    **kwargs,
                )
                plt.close(fig)
                if rec is None:
                    # fall back on the candidate that misplaces the least source
                    best = metrics.loc[metrics["misplaced_fraction"].idxmin()]
                    rec = {"bin_size": best["bin_size [m]"], "k": best["k"]}
                    logging.warning(
                        f"Chunk {i}: no candidate below the misplaced fraction "
                        "threshold, using the one with the lowest "
                        f"({best['misplaced_fraction']:.2%})."
                    )
                bin_size, k = float(rec["bin_size"]), int(rec["k"])

            logging.info(f"Chunk {i}: {len(idx)} points, bin size {bin_size} m, k={k}.")
            cdgs_list.append(
                CDGS.from_cloud_point(
                    df_chunk,
                    isotopes,
                    interpolation_kernel=KNearestKernel(k=k),
                    mesh_definition=MeshByBinSize((bin_size,) * 3),
                    col_names=col_names,
                    e_bins=e_bins,
                    particle=particle,
                    length_unit=length_unit,
                )
            )

        return cls(cdgs_list, check_overlap=False)

    def to_vtk(self, outfile: str | Path) -> None:
        """Write all the meshes to a VTK multiblock file (.vtm).

        Parameters
        ----------
        outfile : str | Path
            Path to the output .vtm file.
        """
        multiblock = pv.MultiBlock(
            {f"mesh_{i + 1}": cdgs.mesh.copy() for i, cdgs in enumerate(self.cdgs_list)}
        )
        multiblock.save(outfile)


def _split_cloud_point(
    coords: np.ndarray,
    split_size: float,
    min_fill: float,
    min_gain: float,
    min_points: int,
) -> list[np.ndarray]:
    """Split a cloud point into spatially compact chunks, see
    MultiCDGS.from_cloud_point.

    The recursive bisection works on the occupied coarse cubes rather than on
    the points, which keeps it fast for millions of points.

    Returns
    -------
    list[np.ndarray]
        The indices of the points of each chunk.
    """
    cube_idx = np.floor((coords - coords.min(axis=0)) / split_size).astype(np.int64)
    cubes, point_cube = np.unique(cube_idx, axis=0, return_inverse=True)
    point_cube = point_cube.ravel()
    counts = np.bincount(point_cube, minlength=len(cubes))

    leaves = []
    stack = [np.arange(len(cubes))]
    while stack:
        node = stack.pop()
        sub = cubes[node]
        volume = np.prod(sub.max(axis=0) - sub.min(axis=0) + 1)
        best = None
        if len(node) / volume < min_fill:
            for axis in range(3):
                order = np.argsort(sub[:, axis], kind="stable")
                s = sub[order]
                # split only between cubes with different coordinates, leaving
                # at least min_points on each side
                cut = np.flatnonzero(np.diff(s[:, axis]) > 0) + 1
                left_points = np.cumsum(counts[node[order]])[cut - 1]
                right_points = counts[node].sum() - left_points
                cut = cut[(left_points >= min_points) & (right_points >= min_points)]
                if len(cut) == 0:
                    continue
                # volumes of the boxes enclosing the cubes on each side
                left_lo = np.minimum.accumulate(s, axis=0)[cut - 1]
                left_hi = np.maximum.accumulate(s, axis=0)[cut - 1]
                right_lo = np.minimum.accumulate(s[::-1], axis=0)[::-1][cut]
                right_hi = np.maximum.accumulate(s[::-1], axis=0)[::-1][cut]
                total = np.prod(left_hi - left_lo + 1, axis=1) + np.prod(
                    right_hi - right_lo + 1, axis=1
                )
                j = np.argmin(total)
                if best is None or total[j] < best[0]:
                    best = (total[j], node[order[: cut[j]]], node[order[cut[j] :]])
        if best is not None and best[0] <= (1 - min_gain) * volume:
            stack.extend(best[1:])
        else:
            leaves.append(node)

    # map the cubes of each leaf back to the points
    leaf_of_cube = np.empty(len(cubes), dtype=np.int64)
    for i, leaf in enumerate(leaves):
        leaf_of_cube[leaf] = i
    point_leaf = leaf_of_cube[point_cube]
    order = np.argsort(point_leaf, kind="stable")
    bounds = np.cumsum(np.bincount(point_leaf, minlength=len(leaves)))[:-1]
    return np.split(order, bounds)

def _read_cloud_point(
    cloud_point: str | Path | pd.DataFrame,
    col_names: dict[str, str] | None,
    length_unit: str,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, float]:
    """Read a cloud point, drop duplicated points and convert coordinates and
    volumes to SI.

    Returns
    -------
    tuple[pd.DataFrame, np.ndarray, np.ndarray, float]
        The (deduplicated) dataframe with the original values, the point
        coordinates in m, the point volumes in m3 and the conversion factor
        from length_unit to m.
    """
    if length_unit not in LENGTH_UNITS:
        raise ValueError(
            f"Length unit {length_unit} not supported. Supported units are "
            f"{list(LENGTH_UNITS)}."
        )
    to_m = LENGTH_UNITS[length_unit]

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

    # convert everything to SI
    coords = df[dup_cols].to_numpy() * to_m  # m
    vols = df[col_names["vol"]].to_numpy() * to_m**3  # m3
    return df, coords, vols, to_m


def _interpolate_cloud_point(
    coords: np.ndarray,
    quantities: dict[str, np.ndarray],
    interpolation_kernel: InterpolationKernel,
    mesh_definition: RegularMeshDefinition,
) -> tuple[pv.StructuredGrid, dict[str, np.ndarray], float]:
    """Build the regular mesh and distribute extensive per-point quantities
    (e.g. atoms, volumes) on its voxels.

    Parameters
    ----------
    coords : np.ndarray
        (n_points, 3) coordinates of the cloud point in m.
    quantities : dict[str, np.ndarray]
        Extensive quantities to distribute, one value per point.
    interpolation_kernel : InterpolationKernel
        Kernel used to distribute the quantities.
    mesh_definition : RegularMeshDefinition
        Definition of the regular mesh.

    Returns
    -------
    tuple[pv.StructuredGrid, dict[str, np.ndarray], float]
        The mesh, the per-voxel totals of each quantity and the voxel volume
        in m3.
    """
    mesh = mesh_definition._build_grid(coords, interpolation_kernel)
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

    results = {}
    for name, per_point in quantities.items():
        interpolation = np.zeros(len(centroids_voxels))
        if getattr(interpolation_kernel, "supports_batch", False):
            interpolation_kernel._kernel_batch(
                coords,
                neighbours_idx,
                centroids_voxels,
                interpolation,
                per_point,
            )
        else:
            # Cycle on all UM centroids and distribute the quantity
            for i in range(len(coords)):
                interpolation_kernel._kernel(
                    coords[i],
                    neighbours_idx[i],
                    centroids_voxels,
                    interpolation,
                    per_point[i],
                )
        results[name] = interpolation
    return mesh, results, vol



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
