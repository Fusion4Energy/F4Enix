from __future__ import annotations
from typing import IO, Dict, Tuple, Union

import numpy as np
import pyvista as pv

from f4enix.input.ww_gvr.models import Pathlike


def get_volume_sampling_error_from_cuv(
    file_path: Pathlike, voxel_sampling_points: int
) -> Dict[int, Dict[int, float]]:
    """
    Read the first mesh tally from a CUV file and return the volume sampling error
    """
    with open(file_path) as infile:
        line = infile.readline()
        while line:
            if "INDEX" in line:
                return _read_cuv_data(infile, voxel_sampling_points)

            line = infile.readline()
        return {}


def _read_cuv_data(
    infile: IO, voxel_sampling_points: int
) -> Dict[int, Dict[int, float]]:
    voxel_cell_errors = {}

    while True:
        try:
            voxel_id, cell_errors = _read_voxel(infile, voxel_sampling_points)
            voxel_cell_errors[voxel_id] = cell_errors
        except IndexError:
            break

    return voxel_cell_errors


def _read_voxel(infile: IO, voxel_sampling_points: int) -> Tuple[int, Dict[int, float]]:
    first_line_words = infile.readline().split()
    # 1-based index on the file, we want 0-based
    voxel_id = int(first_line_words[0]) - 1
    amount_of_cells = int(first_line_words[5])

    cell_errors = {}
    for _ in range(amount_of_cells):
        words = infile.readline().split()
        cell_id = int(words[0])
        cell_error = calculate_volume_sampling_error(
            float(words[1]), voxel_sampling_points
        )
        cell_errors[cell_id] = cell_error

    # no need to read the CUV values and relative errors
    for _ in range(2 * amount_of_cells):
        infile.readline()

    return voxel_id, cell_errors


def calculate_volume_sampling_error(
    partial_volume: float, voxel_sampling_points: int
) -> float:
    return (
        (partial_volume - partial_volume**2) / (voxel_sampling_points)
    ) ** 0.5 / partial_volume


pyvista_grid = Union[
    pv.StructuredGrid, pv.UnstructuredGrid, pv.PolyData, pv.RectilinearGrid
]


def add_sampling_error_to_vtk(
    grid: pyvista_grid, cuv_file_path: Pathlike, voxel_sampling_points: int
) -> pyvista_grid:
    """
    Add a new array to the grid with the volume sampling error. In the case of a voxel
    with more than one cell, the error is the biggest one.
    """
    voxel_cell_errors = get_volume_sampling_error_from_cuv(
        cuv_file_path, voxel_sampling_points
    )

    vector_errors = np.zeros(grid.n_cells)
    for i in range(grid.n_cells):
        cell_error = max(voxel_cell_errors.get(i, {}).values(), default=0)
        vector_errors[i] = cell_error

    grid["Volume Sampling Error"] = vector_errors
    return grid
