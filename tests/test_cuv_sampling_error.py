from pathlib import Path

import pytest
import pyvista as pv

from f4enix.output.cuv_sampling_error import (
    add_sampling_error_to_vtk,
    calculate_volume_sampling_error,
    get_volume_sampling_error_from_cuv,
)
from f4enix.output.meshtal import Meshtal

PATH_TO_CUV_FILE = (
    Path(__file__).parent / "resources" / "meshtal" / "test_CuV" / "cuvmsh"
)
ASSUMED_VOXEL_SAMPLING_POINTS = 1000


def test_get_volume_samplig_error_from_cuv():
    voxel_cell_errors = get_volume_sampling_error_from_cuv(
        file_path=PATH_TO_CUV_FILE, voxel_sampling_points=ASSUMED_VOXEL_SAMPLING_POINTS
    )

    assert voxel_cell_errors[1][1] == pytest.approx(0.0)
    assert voxel_cell_errors[6][1] == pytest.approx(0.01918304084569366)
    assert voxel_cell_errors[6][10] == pytest.approx(0.05212937865502626)


@pytest.mark.parametrize(
    ("partial_volume", "voxel_sampling_points", "expected_result"),
    [
        (1.0, 1000, 0.0),
        (0.9, 1000, 0.010540925533894595),
        (0.5, 1000, 0.03162277660168379),
        (0.5, 10, 0.31622776601683794),
    ],
)
def test_calculate_volume_sampling_error(
    partial_volume, voxel_sampling_points, expected_result
):
    result = calculate_volume_sampling_error(partial_volume, voxel_sampling_points)
    assert result == pytest.approx(expected_result)


def test_add_sampling_error_to_vtk(tmpdir):
    # Read a meshtal and create a VTK file
    meshtal = Meshtal(PATH_TO_CUV_FILE, filetype="CUV")
    meshtal.readMesh(norm="ctot")
    meshtal.mesh[44].write(tmpdir)
    grid = pv.read(tmpdir / "Tally_44_vtk.vtr")
    assert type(grid) == pv.RectilinearGrid

    grid_with_errors = add_sampling_error_to_vtk(
        grid=grid, cuv_file_path=PATH_TO_CUV_FILE, voxel_sampling_points=1000
    )

    assert grid_with_errors["Volume Sampling Error"][0] == pytest.approx(0.0)  # type: ignore
    assert grid_with_errors["Volume Sampling Error"][5] == pytest.approx(  # type: ignore
        0.05121707142933009
    )
