"""
These meshes follow the convention of k, j, i vectors.
j
 --- ---
| 3 | 4 |
 --- ---
| 1 | 2 |
 --- ---    i	 
"""
from typing import Optional
import pyvista as pv
import numpy as np
from numpy.typing import NDArray


def create_cartesian_grid(
    vector_i: NDArray,
    vector_j: NDArray,
    vector_k: NDArray,
    origin: Optional[NDArray] = None,
) -> pv.StructuredGrid:
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
    vector_i = correct_radius_vector(vector_i)
    vector_k_revolutions = correct_theta_vector(vector_k_revolutions)

    # Extend the theta vector if less than 3 intervals
    if len(vector_k_revolutions) <= 3:
        vector_k_revolutions = extend_theta_intervals(vector_k_revolutions)

    # Create the grid on the Z axis
    grid = create_cylindrical_grid_with_z_axis(vector_i, vector_j, vector_k_revolutions)

    # Roto translate the cylinder
    if vec is not None:
        rotate_to_vec(grid, vec)
    if axis is not None:
        rotate_to_axis(grid, axis)
    if origin is not None:
        grid.translate(origin, inplace=True)

    return grid


def extend_theta_intervals(vector_k: NDArray, new_theta_ints: int = 20):
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


def create_cylindrical_grid_with_z_axis(
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


def rotate_to_vec(cylinder: pv.StructuredGrid, vec: NDArray):
    """Rotate a Z axis cylinder by its axis to match the vec."""
    z_axis = np.array([0.0, 0.0, 1.0])
    default_vec = np.array([1.0, 0.0, 0.0])

    # Calculate the angle to rotate
    angle = angle_between(default_vec, vec)

    # Apply the rotation
    cylinder.rotate_vector(vector=z_axis, angle=angle, inplace=True)


def rotate_to_axis(cylinder: pv.StructuredGrid, axis: NDArray):
    default_axis = np.array([0.0, 0.0, 1.0])
    normal_to_rotation = np.cross(default_axis, axis)

    # Calculate the angle to rotate
    angle = angle_between(default_axis, axis)

    # Apply the rotation
    cylinder.rotate_vector(vector=normal_to_rotation, angle=angle, inplace=True)


def correct_radius_vector(vector_i: NDArray):
    """A radius that start at zero generates degenerate cells at the center."""
    if np.isclose(vector_i[0], 0.0):
        vector_i[0] = vector_i[1] * 0.01
    return vector_i


def correct_theta_vector(vector_k_revolutions: NDArray) -> NDArray:
    """With floats sometimes the first and last values are not exactly 0 and 1."""
    # If the cylinder starts at exactly 0 revolutions
    if np.isclose(vector_k_revolutions[0], 0.0):
        vector_k_revolutions[0] = 0.0

    # If the cylinder ends at exactly 1 revolution
    if np.isclose(vector_k_revolutions[-1], 1.0):
        vector_k_revolutions[-1] = 1.0

    return vector_k_revolutions


def angle_between(vec_1: NDArray, vec_2: NDArray) -> float:
    """Returns the angle in degrees between vectors 'v1' and 'v2'"""
    v1_normalized = normalize_vector(vec_1)
    v2_normalized = normalize_vector(vec_2)
    radians = np.arccos(np.clip(np.dot(v1_normalized, v2_normalized), -1.0, 1.0))
    return np.rad2deg(radians)


def normalize_vector(vector: NDArray) -> NDArray:
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)
