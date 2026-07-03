import pyvista as pv
import numpy as np
from abc import ABC, abstractmethod
from sklearn.neighbors import KDTree


# --- Mesh definitions ---
class RegularMeshDefinition(ABC):
    @abstractmethod
    def _build_grid(
        self,
        coords: np.ndarray,
    ) -> pv.StructuredGrid:
        pass

    def _get_regular_mesh(
        self,
        coords: np.ndarray,
        steps: tuple[float, float, float],
        bbox: dict | None = None,
    ) -> pv.StructuredGrid:
        if bbox is not None:
            x_min, y_min, z_min = bbox["x_min"], bbox["y_min"], bbox["z_min"]
            x_max, y_max, z_max = bbox["x_max"], bbox["y_max"], bbox["z_max"]
        else:
            # bounding box will be min and max of the coordinates + half the step
            x_min, y_min, z_min = coords.min(axis=0) - np.array(steps) / 2
            x_max, y_max, z_max = coords.max(axis=0) + np.array(steps) / 2

        # build the grid
        xrng = np.linspace(x_min, x_max, int((x_max - x_min) / steps[0]) + 1)
        yrng = np.linspace(y_min, y_max, int((y_max - y_min) / steps[1]) + 1)
        zrng = np.linspace(z_min, z_max, int((z_max - z_min) / steps[2]) + 1)
        x, y, z = np.meshgrid(xrng, yrng, zrng, indexing="ij")
        grid = pv.StructuredGrid(x, y, z)
        return grid


class MeshByAverageDistance(RegularMeshDefinition):
    def __init__(self, factor: float, bbox: dict | None = None):
        """Build the structured mesh with constant voxels. Dimension will be based on
        average distance between cloud points.

        Parameters
        ----------
        factor : float
            Factor to scale the average distance between cloud points to
            determine the voxel size.
        bbox : dict, optional
            Bounding box for the mesh. If None, the bounding box will be determined
            from the coordinates of the cloud points. The dictionary should have the
            following keys: "x_min", "x_max", "y_min", "y_max", "z_min", "z_max".
        """
        self.factor = factor
        self.bbox = bbox

    def _build_grid(
        self,
        coords: np.ndarray,
    ) -> pv.StructuredGrid:
        tree = KDTree(coords)
        distances, _ = tree.query(coords, k=2)
        distances = distances[:, 1]
        step = distances.mean() * self.factor
        steps = (step, step, step)
        return self._get_regular_mesh(coords, steps, bbox=self.bbox)


class MeshByBinSize(RegularMeshDefinition):
    def __init__(self, bin_size: tuple[float, float, float], bbox: dict | None = None):
        """Build the structured mesh with constant voxels. Dimension will be based on
        the provided bin size.

        Parameters
        ----------
        bin_size : tuple[float, float, float]
            Size of the bins in each dimension (x, y, z).
        bbox : dict, optional
            Bounding box for the mesh. If None, the bounding box will be determined
            from the coordinates of the cloud points. The dictionary should have the
            following keys: "x_min", "x_max", "y_min", "y_max", "z_min", "z_max".
        """
        self.bin_size = bin_size
        self.bbox = bbox

    def _build_grid(
        self,
        coords: np.ndarray,
    ) -> pv.StructuredGrid:
        return self._get_regular_mesh(coords, self.bin_size, bbox=self.bbox)


class MeshByNumberOfVoxels(RegularMeshDefinition):
    def __init__(self, n_voxels: tuple[int, int, int], bbox: dict):
        """Build the structured mesh with constant voxels. Dimension will be based on
        the provided number of voxels.

        Parameters
        ----------
        n_voxels : tuple[int, int, int]
            Number of voxels in each dimension (x, y, z).
        bbox : dict
            Bounding box for the mesh. The dictionary should have the
            following keys: "x_min", "x_max", "y_min", "y_max", "z_min", "z_max".
        """
        self.n_voxels = n_voxels
        self.bbox = bbox

    def _build_grid(
        self,
        coords: np.ndarray,
    ) -> pv.StructuredGrid:
        x_min, y_min, z_min = (
            self.bbox["x_min"],
            self.bbox["y_min"],
            self.bbox["z_min"],
        )
        x_max, y_max, z_max = (
            self.bbox["x_max"],
            self.bbox["y_max"],
            self.bbox["z_max"],
        )

        steps = (
            (x_max - x_min) / self.n_voxels[0],
            (y_max - y_min) / self.n_voxels[1],
            (z_max - z_min) / self.n_voxels[2],
        )
        return self._get_regular_mesh(coords, steps, bbox=self.bbox)
