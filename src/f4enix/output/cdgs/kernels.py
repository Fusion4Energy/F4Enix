import numpy as np
from abc import ABC, abstractmethod
import logging
from sklearn.neighbors import KDTree


# --- Interpolation kernels ---
class InterpolationKernel(ABC):
    # subclasses that can distribute all source points in one vectorized call
    # (see KNearestKernel._kernel_batch) should override this to True
    supports_batch: bool = False

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
    def _kernel_batch(
        self,
        source_coords: np.ndarray,
        neighbours_idx: np.ndarray,
        dest_coords: np.ndarray,
        interpolation: np.ndarray,
        activities: np.ndarray,
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

    def _kernel_batch(self, *args, **kwargs) -> None:
        raise NotImplementedError("Batch kernel not implemented for this kernel.")


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
        # manual sqrt(sum of squares) + ndarray methods (.sum()/.max()) avoid the
        # extra argument validation of np.linalg.norm/np.sum/np.max, which is
        # pure overhead when called millions of times on tiny arrays
        diff = dest_coords[dest_idx, :] - source_coord
        distances = np.sqrt((diff * diff).sum(axis=1))
        if (distances == 0).any():
            # Coincident point: assign full weight to it, zero to others
            weights = np.zeros_like(distances)
        else:
            max_distance = distances.max()
            weights = max_distance / distances  # inverse distance weights
            weights /= weights.sum()  # normalize the weights so that they sum to 1
        interpolation[dest_idx] += weights * activity


class KNearestKernel(InterpolationKernel):
    """Distribute activity among a fixed number of destination points, the k
    nearest to the source coordinate, weighted by the inverse of their distance.

    Unlike SphereKernel/DistanceKernel, which search all points within a radius
    (returning a variable neighbour count per source point), KNearestKernel
    always returns exactly k neighbours. This makes get_neighbours return a
    dense (n_points, k) array instead of a ragged one, which is a prerequisite
    for vectorizing the interpolation across all source points at once.
    """

    supports_batch = True

    def __init__(self, k: int):
        if k < 1:
            raise ValueError("k must be a positive integer")
        self.k = k

    def get_neighbours(
        self, source_coords: np.ndarray, dest_coords: np.ndarray, **kwargs
    ) -> np.ndarray:
        if self.k > len(dest_coords):
            raise ValueError(
                f"k ({self.k}) cannot exceed the number of destination points "
                f"({len(dest_coords)})."
            )
        tree = KDTree(dest_coords)
        _, neighbours_idx = tree.query(source_coords, k=self.k)
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
        diff = dest_coords[dest_idx, :] - source_coord
        distances = np.sqrt((diff * diff).sum(axis=1))
        zero_mask = distances == 0
        if zero_mask.any():
            # Coincident point(s): assign full weight to them, zero to others
            weights = zero_mask.astype(distances.dtype)
            weights /= weights.sum()
        else:
            max_distance = distances.max()
            weights = max_distance / distances  # inverse distance weights
            weights /= weights.sum()  # normalize the weights so that they sum to 1
        interpolation[dest_idx] += weights * activity

    def _kernel_batch(
        self,
        source_coords: np.ndarray,
        neighbours_idx: np.ndarray,
        dest_coords: np.ndarray,
        interpolation: np.ndarray,
        activities: np.ndarray,
    ) -> None:
        """Vectorized equivalent of _kernel, distributing all source points at
        once. Requires neighbours_idx to be the dense (n_points, k) array
        returned by get_neighbours (ragged arrays cannot be used here)."""
        dest_pts = dest_coords[neighbours_idx]  # (n_points, k, 3)
        diff = dest_pts - source_coords[:, None, :]
        distances = np.sqrt((diff * diff).sum(axis=2))  # (n_points, k)

        zero_mask = distances == 0
        any_zero = zero_mask.any(axis=1, keepdims=True)

        # inverse distance weights, guarding the coincident-point rows against
        # division by zero (they are discarded by np.where below anyway)
        max_distance = distances.max(axis=1, keepdims=True)
        safe_distances = np.where(zero_mask, 1.0, distances)
        inv_weights = max_distance / safe_distances
        inv_weights /= inv_weights.sum(axis=1, keepdims=True)

        # coincident point(s): full weight to them, zero to others
        coincident_weights = zero_mask.astype(distances.dtype)
        row_sums = coincident_weights.sum(axis=1, keepdims=True)
        coincident_weights = np.divide(
            coincident_weights,
            row_sums,
            out=np.zeros_like(coincident_weights),
            where=row_sums != 0,
        )

        weights = np.where(any_zero, coincident_weights, inv_weights)
        contributions = weights * activities[:, None]
        # scatter-add: several source points can share a destination voxel, so
        # a plain fancy-index assignment would overwrite instead of accumulate
        np.add.at(interpolation, neighbours_idx.ravel(), contributions.ravel())
