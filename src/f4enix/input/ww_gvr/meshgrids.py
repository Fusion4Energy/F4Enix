"""
This file contains functions to create and operate pyvista meshes of cartesian and
cylindrical grids.

These meshes follow the convention of k, j, i vectors.
    j
    --- ---
    | 3 | 4 |
    --- ---
    | 1 | 2 |
    --- ---    i

The values are given in an array of shape (k, j, i).
"""

# flake8: noqa: PLR2004
from typing import Optional

import numpy as np
import pyvista as pv
from numpy.typing import NDArray


def create_cartesian_grid(
    vector_i: NDArray,
    vector_j: NDArray,
    vector_k: NDArray,
    origin: Optional[NDArray] = None,
) -> pv.StructuredGrid:
    """
    Create a cartesian grid with the given vectors and origin.

    Parameters
    ----------
    vector_i : NDArray
        Vector defining the i direction intervals.
    vector_j : NDArray
        Vector defining the j direction intervals.
    vector_k : NDArray
        Vector defining the k direction intervals.
    origin : NDArray | None
        Origin vector of the grid, by default None (0,0,0).
    """
    # meshgrid should receive float arrays to avoid a warning
    vector_i = np.asarray(vector_i, "float32")
    vector_j = np.asarray(vector_j, "float32")
    vector_k = np.asarray(vector_k, "float32")

    x_coordinates, y_coordinates, z_coordinates = np.meshgrid(
        vector_i, vector_j, vector_k, indexing="ij"
    )
    grid = pv.StructuredGrid(x_coordinates, y_coordinates, z_coordinates)

    if origin is not None:
        grid.translate(origin, inplace=True)

    return grid


def create_cylindrical_grid(
    vector_i: NDArray,
    vector_j: NDArray,
    vector_k_revolutions: NDArray,
    origin=None,
    axis=None,
    vec=None,
) -> pv.StructuredGrid:
    """
    Create a cylindrical grid with the given vectors and origin.

    In the case of a single theta interval or two theta intervals, the theta vector is
    artificially extended to 20 intervals as they couldnt be visually represented.

    Parameters
    ----------
    vector_i : NDArray
        Vector defining the radial direction intervals.
    vector_j : NDArray
        Vector defining height direction intervals.
    vector_k_revolutions : NDArray
        Vector defining k direction intervals in revolution units.
    origin : NDArray | None
        Vector pointing to the bottom center of the cylinder, by default None (0,0,0).
    axis : NDArray | None
        Vector defining the axis of the cylinder, by default None (Z-axis).
    vec : NDArray | None
        Vector defining the direction of the radius, by default None (X-axis).
    """
    vector_i = _correct_radius_vector(vector_i)
    vector_k_revolutions = _correct_theta_vector(vector_k_revolutions)

    # Extend the theta vector if less than 3 intervals
    if len(vector_k_revolutions) <= 3:
        vector_k_revolutions = _extend_theta_intervals(vector_k_revolutions)

    # Create the grid on the Z axis
    grid = _create_cylindrical_grid_with_z_axis(
        vector_i, vector_j, vector_k_revolutions
    )

    # Roto translate the cylinder
    if vec is not None:
        _rotate_to_vec(grid, vec)
    if axis is not None:
        _rotate_to_axis(grid, axis)
    if origin is not None:
        grid.translate(origin, inplace=True)

    return grid


def _extend_theta_intervals(vector_k: NDArray, new_theta_ints: int = 20):
    if len(vector_k) == 2:
        extended_vector_k = np.linspace(vector_k[0], vector_k[1], new_theta_ints + 1)

    elif len(vector_k) == 3:
        new_theta_ints = new_theta_ints // 2
        vector_k_0 = np.linspace(vector_k[0], vector_k[1], new_theta_ints + 1)
        vector_k_1 = np.linspace(vector_k[1], vector_k[2], new_theta_ints + 1)[1:]
        extended_vector_k = np.concatenate((vector_k_0, vector_k_1))

    else:
        raise ValueError("Only 1 or 2 theta defined intervals can be extended...")

    return extended_vector_k


def _create_cylindrical_grid_with_z_axis(
    vector_i: NDArray, vector_j: NDArray, vector_k_revolutions: NDArray
) -> pv.StructuredGrid:
    vector_k_radians = np.array(vector_k_revolutions) * 2 * np.pi

    radii, thetas, z = np.array(
        np.meshgrid(vector_i, vector_k_radians, vector_j, indexing="ij")
    )
    x = np.cos(thetas) * radii
    y = np.sin(thetas) * radii
    points = np.c_[x.ravel(order="f"), y.ravel(order="f"), z.ravel(order="f")]

    grid = pv.StructuredGrid()
    grid.points = points
    grid.dimensions = len(vector_i), len(vector_k_revolutions), len(vector_j)

    return grid


def _rotate_to_vec(cylinder: pv.StructuredGrid, vec: NDArray):
    """Rotate a Z axis cylinder by its axis to match the vec."""
    z_axis = np.array([0.0, 0.0, 1.0])
    default_vec = np.array([1.0, 0.0, 0.0])

    # Calculate the angle to rotate
    angle = _angle_between(default_vec, vec)

    # Apply the rotation
    cylinder.rotate_vector(vector=z_axis, angle=angle, inplace=True)


def _rotate_to_axis(cylinder: pv.StructuredGrid, axis: NDArray):
    default_axis = np.array([0.0, 0.0, 1.0])
    normal_to_rotation = np.cross(default_axis, axis)

    # Calculate the angle to rotate
    angle = _angle_between(default_axis, axis)

    # Apply the rotation
    cylinder.rotate_vector(vector=normal_to_rotation, angle=angle, inplace=True)


def _correct_radius_vector(vector_i: NDArray):
    """A radius that start at zero generates degenerate cells at the center."""
    if np.isclose(vector_i[0], 0.0):
        vector_i[0] = vector_i[1] * 0.01
    return vector_i


def _correct_theta_vector(vector_k_revolutions: NDArray) -> NDArray:
    """With floats sometimes the first and last values are not exactly 0 and 1."""
    # If the cylinder starts at exactly 0 revolutions
    if np.isclose(vector_k_revolutions[0], 0.0):
        vector_k_revolutions[0] = 0.0

    # If the cylinder ends at exactly 1 revolution
    if np.isclose(vector_k_revolutions[-1], 1.0):
        vector_k_revolutions[-1] = 1.0

    return vector_k_revolutions


def _angle_between(vec_1: NDArray, vec_2: NDArray) -> float:
    """Returns the angle in degrees between vectors 'v1' and 'v2'"""
    v1_normalized = _normalize_vector(vec_1)
    v2_normalized = _normalize_vector(vec_2)
    radians = np.arccos(np.clip(np.dot(v1_normalized, v2_normalized), -1.0, 1.0))
    return np.rad2deg(radians)


def _normalize_vector(vector: NDArray) -> NDArray:
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)
