import numpy as np
import pytest
import pyvista as pv
from numpy.testing import assert_array_almost_equal, assert_array_equal

from f4enix.input.ww_gvr.meshgrids import (
    _correct_theta_vector,
    _create_cylindrical_grid_with_z_axis,
    _extend_theta_intervals,
    create_cartesian_grid,
    create_cylindrical_grid,
)

# ruff: noqa: PLR2004


def test_create_cartesian_grid():
    grid = create_cartesian_grid(
        vector_i=np.array([0, 1, 2]),
        vector_j=np.array([10, 20, 30]),
        vector_k=np.array([5, 6, 7]),
    )
    assert 8 == grid.n_cells

    cell_0 = grid.extract_cells(0)
    cell_1 = grid.extract_cells(1)
    cell_7 = grid.extract_cells(7)

    assert isinstance(cell_0, pv.UnstructuredGrid)
    assert isinstance(cell_1, pv.UnstructuredGrid)
    assert isinstance(cell_7, pv.UnstructuredGrid)

    assert_array_almost_equal([0.5, 15, 5.5], cell_0.center)
    assert_array_almost_equal([1.5, 15, 5.5], cell_1.center)
    assert_array_almost_equal([1.5, 25, 6.5], cell_7.center)


def test_create_cartesian_grid_translated():
    grid = create_cartesian_grid(
        vector_i=np.array([0, 1, 2]),
        vector_j=np.array([10, 20, 30]),
        vector_k=np.array([5, 6, 7]),
        origin=np.array([100, 100, 100]),
    )
    assert 8 == grid.n_cells

    cell_0 = grid.extract_cells(0)
    cell_1 = grid.extract_cells(1)
    cell_7 = grid.extract_cells(7)

    assert isinstance(cell_0, pv.UnstructuredGrid)
    assert isinstance(cell_1, pv.UnstructuredGrid)
    assert isinstance(cell_7, pv.UnstructuredGrid)

    assert_array_almost_equal([100.5, 115, 105.5], cell_0.center)
    assert_array_almost_equal([101.5, 115, 105.5], cell_1.center)
    assert_array_almost_equal([101.5, 125, 106.5], cell_7.center)


def test_create_cylindrical_grid_z_axis():
    grid = _create_cylindrical_grid_with_z_axis(
        vector_i=np.array([0, 1, 2]),
        vector_j=np.array([10, 20, 30]),
        vector_k_revolutions=np.array([0, 0.25, 0.5]),
    )

    assert 8 == grid.n_cells

    cell_0 = grid.extract_cells(0)
    cell_1 = grid.extract_cells(1)
    cell_2 = grid.extract_cells(2)
    cell_7 = grid.extract_cells(7)

    assert isinstance(cell_0, pv.UnstructuredGrid)
    assert isinstance(cell_1, pv.UnstructuredGrid)
    assert isinstance(cell_2, pv.UnstructuredGrid)
    assert isinstance(cell_7, pv.UnstructuredGrid)

    assert_array_almost_equal([0.5, 0.5, 15], cell_0.center)
    assert_array_almost_equal([1.0, 1, 15], cell_1.center)
    assert_array_almost_equal([-0.5, 0.5, 15], cell_2.center)
    assert_array_almost_equal([-1.0, 1, 25], cell_7.center)


def test_extend_theta_intervals_1_initial_int():
    extended_vector = _extend_theta_intervals(np.array([0.0, 1]), 4)
    assert_array_almost_equal([0, 0.25, 0.5, 0.75, 1], extended_vector)


def test_extend_theta_intervals_2_initial_ints():
    extended_vector = _extend_theta_intervals(np.array([0, 0.5, 1]), 4)
    assert_array_almost_equal([0, 0.25, 0.5, 0.75, 1], extended_vector)


def test_extend_theta_intervals_3_initial_ints():
    with pytest.raises(ValueError):
        _extend_theta_intervals(np.array([0, 0.2, 0.4, 0.6]))


def test_correct_theta_vector():
    wrong_vector = np.array([-0.0000000001, 0.5, 0.9999999999995])
    corrected_vector = _correct_theta_vector(wrong_vector)
    assert_array_equal([0, 0.5, 1], corrected_vector)


def test_create_cylindrical_grid_rotated_axis():
    vector_i = np.array([0, 10])
    vector_j = np.array([0, 20])
    vector_k = np.array([1e-25, 1 + 1e-25])
    origin = np.array([100, 100, 100])
    axis = np.array([1, 0, 0])

    grid = create_cylindrical_grid(
        vector_i, vector_j, vector_k, origin=origin, axis=axis
    )

    assert_array_almost_equal([110.0, 100.0, 100.0], np.array(grid.center))
    assert 20 == grid.n_cells


def test_create_cylindrical_grid_rotated_vec():
    vector_i = np.array([0, 10])
    vector_j = np.array([0, 10])
    vector_k = np.array([0, 0.25, 0.5, 0.75, 1])
    origin = np.array([100.0, 100, 100])
    vec = np.array([0.0, 1, 0])

    non_rotated_grid = create_cylindrical_grid(
        vector_i, vector_j, vector_k, origin=origin
    )
    cell_0 = non_rotated_grid.extract_cells(0)
    cell_1 = non_rotated_grid.extract_cells(1)
    assert isinstance(cell_0, pv.UnstructuredGrid)
    assert isinstance(cell_1, pv.UnstructuredGrid)
    assert_array_almost_equal([105.0, 105, 105], cell_0.center)
    assert_array_almost_equal([95.0, 105, 105], cell_1.center)

    # Same grid as before but rotated 90 degrees around the Z axis
    rotated_grid = create_cylindrical_grid(
        vector_i, vector_j, vector_k, origin=origin, vec=vec
    )
    cell_0 = rotated_grid.extract_cells(0)
    cell_1 = rotated_grid.extract_cells(1)
    assert isinstance(cell_0, pv.UnstructuredGrid)
    assert isinstance(cell_1, pv.UnstructuredGrid)
    assert_array_almost_equal([95.0, 105, 105], cell_0.center)
    assert_array_almost_equal([95.0, 95, 105], cell_1.center)
