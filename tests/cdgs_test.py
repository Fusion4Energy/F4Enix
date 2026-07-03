from importlib.resources import as_file, files
import numpy as np
import pyvista as pv

import tests.resources.cdgs as res
from f4enix.output.cdgs import (
    CDGS,
    SphereKernel,
    DistanceKernel,
    MeshByAverageDistance,
    MeshByNumberOfVoxels,
    MeshByBinSize,
)


class TestKernel:
    def test_get_neighbours_sphere_kernel(self):
        kernel = SphereKernel(radius=0.75)
        source_coords = np.array([[0, 0, 0], [0.1, 0.1, 0.1], [1, 1, 1]])
        dest_coords = np.array([[0, 0, 0], [0.4, 0.4, 0.4], [0.9, 0.9, 0.9]])
        neighbours_idx = kernel.get_neighbours(source_coords, dest_coords)
        assert len(neighbours_idx) == 3
        assert list(neighbours_idx[0]) == [0, 1]
        assert list(neighbours_idx[1]) == [0, 1]
        assert list(neighbours_idx[2]) == [2]


class TestMeshDefinition:
    def test_all_definitions_produce_same_mesh(self):
        # Regular grid of 12 points with unit spacing.
        # Every point's nearest neighbour is at distance 1.0, so the mean
        # nearest-neighbour distance is exactly 1.0 and the auto bounding box
        # expands by half a step on each side:
        #   x in [-0.5, 2.5]  →  3 bins
        #   y in [-0.5, 1.5]  →  2 bins
        #   z in [-0.5, 1.5]  →  2 bins  →  12 cells total
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 1.0, 1.0],
            ]
        )
        # Explicit bbox that matches the auto bbox used by the other two definitions
        bbox = {
            "x_min": -0.5,
            "x_max": 2.5,
            "y_min": -0.5,
            "y_max": 1.5,
            "z_min": -0.5,
            "z_max": 1.5,
        }

        grid_avg = MeshByAverageDistance(factor=1.0)._build_grid(coords)
        grid_bin = MeshByBinSize(bin_size=(1.0, 1.0, 1.0))._build_grid(coords)
        grid_nvox = MeshByNumberOfVoxels(n_voxels=(3, 2, 2), bbox=bbox)._build_grid(
            coords
        )

        assert grid_avg.n_cells == grid_bin.n_cells == grid_nvox.n_cells == 12
        assert np.allclose(grid_avg.points, grid_bin.points)
        assert np.allclose(grid_avg.points, grid_nvox.points)

    def test_mesh_by_avg_distance(self):
        coords = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        mesh_def = MeshByAverageDistance(factor=1.0)
        grid = mesh_def._build_grid(coords)
        assert isinstance(grid, pv.StructuredGrid)
        assert grid.n_cells > 0


class TestCDGS:
    def test_from_cloud_point(self):
        with as_file(files(res).joinpath("test_activity.csv")) as file:
            cdgs = CDGS.from_cloud_point(
                file,
                {"N16": "n16"},
                interpolation_kernel=SphereKernel(n_voxels=10),
                mesh_definition=MeshByAverageDistance(factor=10),
                col_names={"x": "x", "y": "y", "z": "z", "vol": "cell-volume"},
            )
