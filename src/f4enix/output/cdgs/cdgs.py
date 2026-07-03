from enum import Enum
import numpy as np
import pandas as pd
from pathlib import Path
import pyvista as pv
from sklearn.neighbors import KDTree
import actigamma as ag
from f4enix.output.cdgs import InterpolationKernel, RegularMeshDefinition


class CDGS_ENERGY_TYPE(Enum):
    LINE = "line"
    BINS = "bins"


class CDGS_MESH_TYPE(Enum):
    RECTANGULAR = "rec"
    # CYLINDRICAL = 'cyl'
    # UNSTRUCTURED = 'unstr'


class CDGS:
    def __init__(self, mesh: pv.StructuredGrid, translation: np.ndarray | None = None):
        # self._title: str = ""
        # self._cooling_time: float = 0
        # self._mesh_source: float = 0
        # self._energy_type: CDGS_ENERGY_TYPE = None
        # self._energy_boundaries: list[float] = []
        # self._mesh_type: CDGS_MESH_TYPE = None
        # self._mesh_boundaries: np.ndarray = np.array([])
        # if translation is None:
        #     self._translation: np.ndarray = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # else:
        #     self._translation: np.ndarray = np.array(translation)
        self.mesh: pv.StructuredGrid = mesh

    def to_vtk(self, outfile: str | Path) -> None:
        """Write the CDGS object to a VTK file.

        Parameters
        ----------
        outfile : str | Path
            Path to the output VTK file.
        """
        self.mesh.save(outfile)

    @classmethod
    def from_cloud_point(
        cls,
        csv_file: str | Path,
        isotopes: dict[str, str],
        interpolation_kernel: InterpolationKernel,
        mesh_definition: RegularMeshDefinition,
        col_names: dict[str, str] | None = None,
    ) -> "CDGS":
        """Read a csv (typically produced from fluent) that contains a cloud
        point of data and convert it into a CDGS object.

        Parameters
        ----------
        csv_file : str | Path
            Path to the CSV file containing the cloud point data.
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

        Returns
        -------
        CDGS
            A CDGS object populated with data from the CSV file.
        """
        df = pd.read_csv(csv_file)
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

        # Check if there are duplicated points, if so, drop them and warn the user
        if len(coords) != len(np.unique(coords, axis=0)):
            print(
                "Warning: Duplicated points found in the cloud point data. "
                "These will be dropped."
            )
            df = df.drop_duplicates(
                subset=[col_names["x"], col_names["y"], col_names["z"]]
            )
            coords = df[[col_names["x"], col_names["y"], col_names["z"]]].to_numpy()

        UM_vols = df[col_names["vol"]].to_numpy()

        # build the regular mesh
        mesh = mesh_definition._build_grid(coords)
        vol = np.min(np.abs(mesh.compute_cell_sizes()["Volume"]))
        centroids_voxels = mesh.cell_centers().points

        # Query the KDTree to find neighboring points
        neighbours_idx = interpolation_kernel.get_neighbours(
            source_coords=coords,
            dest_coords=centroids_voxels,
            voxel_size=vol ** (1 / 3),
        )

        # add the different isotopes activities
        db = ag.Decay2012Database()

        for isotope, col in isotopes.items():
            # lines = db.getenergies("Co60", spectype="gamma")  # eV
            # intensities = db.getintensities(isotope, spectype="gamma")

            # initialize the activity array for each isotope and isotope counter
            activity_array = np.zeros(len(centroids_voxels))
            atoms_array = np.zeros(len(centroids_voxels))

            # Cycle on all UM centroids and distribute the activity
            for i in range(len(coords)):
                source_coord = coords[i]
                atoms = df[col].iloc[i] * UM_vols[i]  # atoms
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
            mesh.cell_data[isotope] = activity_array / vol  # Bq/m3
            mesh.cell_data[f"{isotope}_atoms"] = atoms_array / vol  # atoms/m3

        return cls(mesh)

    def to_cdgs(outfile: str | Path) -> None:
        pass
