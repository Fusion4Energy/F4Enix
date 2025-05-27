"""
Includes some numba functions for improved performance of the ratios calculation.
"""

import numpy as np
from numba import njit, prange


@njit(cache=True)
def calculate_max_ratio_array(array: np.ndarray) -> np.ndarray:
    """
    Calculate the maximum ratio of each value to its neighbours in a 3D array.

    Parameters
    ----------
    array : np.ndarray
        The input 3D array.

    Returns
    -------
    np.ndarray
        A 3D array of the same shape as the input, where each value is the maximum
        ratio of the corresponding value in the input array to its neighbours.
    """
    ratios = np.ones_like(array)
    for k in prange(array.shape[0]):
        for j in prange(array.shape[1]):
            for i in prange(array.shape[2]):
                value = array[k, j, i]
                neighbour_values = _find_neighbour_values(array, k, j, i)
                ratios[k, j, i] = _calculate_max_ratio_from_neighbours(
                    value, neighbour_values
                )
    return ratios


@njit(cache=True)
def _find_neighbour_values(array: np.ndarray, k: int, j: int, i: int) -> np.ndarray:
    neighbour_values = np.zeros(6, dtype=float)
    neighbouring_indices = [
        [k, j, i + 1],
        [k, j, i - 1],
        [k, j + 1, i],
        [k, j - 1, i],
        [k + 1, j, i],
        [k - 1, j, i],
    ]
    for i_index, index in enumerate(neighbouring_indices):
        index_is_valid = True

        for index_i, dimension_i in zip(index, array.shape):
            if index_i < 0 or index_i >= dimension_i:
                index_is_valid = False
                break

        if index_is_valid:
            neighbour_values[i_index] = array[index[0], index[1], index[2]]

    return neighbour_values


@njit(cache=True)
def _calculate_max_ratio_from_neighbours(
    value: float, neighbour_values: np.ndarray
) -> float:
    ratios = np.divide(value, neighbour_values)

    # Invert the values that are less than 1
    ratios = np.where(ratios < 1, 1 / ratios, ratios)

    # Remove the infinities
    ratios = ratios[np.isfinite(ratios)]

    if ratios.size > 0:
        return np.max(ratios)
    else:
        return 1
