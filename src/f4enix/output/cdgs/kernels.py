import numpy as np
from abc import ABC, abstractmethod
import logging
from sklearn.neighbors import KDTree


# --- Interpolation kernels ---
class InterpolationKernel(ABC):
    @abstractmethod
    def _kernel(
        self,
        source_coord: np.ndarray,
        dest_idx: np.ndarray,
        dest_coords: np.ndarray,
        interpolation: np.ndarray,
        activity: float,
    ) -> None:
        pass

    @abstractmethod
    def get_neighbours(
        self, source_coords: np.ndarray, dest_coords: np.ndarray, **kwargs
    ) -> np.ndarray:
        """Given a KDTree of the source coordinates and

        Parameters
        ----------
        source_coords : np.ndarray
            cloud point coordinates of the source activity
        dest_coords : np.ndarray
            coordinates of the destination points

        Returns
        -------
        np.ndarray
            matrix of indices of the neighbors for each source coordinate in the cloud point
        """
        pass


class SphereKernel(InterpolationKernel):
    def __init__(self, radius: float | None = None, n_voxels: int | None = None):
        self.radius = radius
        self.n_voxels = n_voxels
        if self.radius is None and self.n_voxels is None:
            raise ValueError("Either radius or n_voxels must be provided")
        if self.radius is not None and self.n_voxels is not None:
            raise ValueError("Only one of radius or n_voxels should be provided.")

    def get_neighbours(
        self,
        source_coords: np.ndarray,
        dest_coords: np.ndarray,
        voxel_size: float | None = None,
    ) -> np.ndarray:
        if self.radius is not None:
            search_radius = self.radius
        elif self.n_voxels is not None and voxel_size is not None:
            search_radius = self.n_voxels * voxel_size
        else:
            raise ValueError("Insufficient information to determine search radius.")

        tree = KDTree(dest_coords)
        neighbours_idx = tree.query_radius(source_coords, r=search_radius)
        # check if there are empty queries and log a warning
        for i, idx in enumerate(neighbours_idx):
            if len(idx) == 0:
                logging.warning(
                    f"No neighbours found for source coordinate {source_coords[i]} within radius {search_radius}."
                )

        return neighbours_idx

    def _kernel(
        self,
        source_coord: np.ndarray,
        dest_idx: np.ndarray,
        dest_coords: np.ndarray,
        interpolation: np.ndarray,
        activity: float,
    ) -> None:
        """simply distribute the EM force components proportionally to the inverse of the
        distance"""
        interpolation[dest_idx] += activity / len(
            dest_idx
        )  # uniform distribution of activity


class DistanceKernel(SphereKernel):
    def _kernel(
        self,
        source_coord: np.ndarray,
        dest_idx: np.ndarray,
        dest_coords: np.ndarray,
        interpolation: np.ndarray,
        activity: float,
    ) -> None:
        """simply distribute the EM force components proportionally to the inverse of the
        distance"""
        distances = np.linalg.norm(dest_coords[dest_idx, :] - source_coord, axis=1)
        max_distance = np.max(distances)
        weights = 1 / (distances / max_distance)  # inverse distance weights
        weights /= np.sum(weights)  # normalize the weights so that they sum to 1
        interpolation[dest_idx] += weights * activity
