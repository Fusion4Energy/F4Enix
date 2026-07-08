from importlib.resources import as_file, files
import numpy as np
import pyvista as pv
import io
import pytest

import tests.resources.cdgs as res
from f4enix.output.cdgs import (
    CDGS,
    SphereKernel,
    MeshByAverageDistance,
    MeshByNumberOfVoxels,
    MeshByBinSize,
    CDGS_ENERGY_TYPE,
)
from f4enix.output.cdgs.cdgs import (
    ACTIVITY_TAG,
    ATOM_DENSITY_TAG,
    _floats_to_multiline_string,
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
    @pytest.fixture
    def cdgs(self) -> CDGS:
        with as_file(files(res).joinpath("test_activity.csv")) as file:
            return CDGS.from_cloud_point(
                file,
                {"N16": "n16", "O19": "o19"},
                interpolation_kernel=SphereKernel(n_voxels=10),
                mesh_definition=MeshByAverageDistance(factor=10),
                col_names={"x": "x", "y": "y", "z": "z", "vol": "cell-volume"},
            )

    @pytest.fixture
    def easy_cdgs(self):
        # Build a regular 3x3-cell mesh (4x4x2 points -> 9 cells)
        x = np.arange(4, dtype=float)
        y = np.arange(4, dtype=float)
        z = np.array([0.0, 1.0], dtype=float)
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        mesh = pv.StructuredGrid(xx, yy, zz)

        # Minimal CDGS-like scalar fields on cells
        mesh.cell_data[f"N16{ACTIVITY_TAG}"] = np.ones(mesh.n_cells, dtype=float)
        mesh.cell_data[f"N16{ATOM_DENSITY_TAG}"] = np.full(mesh.n_cells, 1e-6)

        mesh.cell_data[f"Co60{ACTIVITY_TAG}"] = np.ones(mesh.n_cells, dtype=float) * 2
        mesh.cell_data[f"Co60{ATOM_DENSITY_TAG}"] = np.full(mesh.n_cells, 1e-6)

        return CDGS(mesh)

    def test_update_cooling(self, cdgs: CDGS):
        # Test that the cooling time updates the mesh correctly
        initial_mesh = cdgs.mesh.copy()
        cdgs.cooling_time = 10.0  # Set a cooling time
        updated_mesh = cdgs.mesh
        assert not np.allclose(
            initial_mesh.cell_data[f"N16{ACTIVITY_TAG}"],
            updated_mesh.cell_data[f"N16{ACTIVITY_TAG}"],
        )

        # test correct revert if set to zero
        cdgs.cooling_time = 0
        reverted_mesh = cdgs.mesh
        assert np.allclose(
            initial_mesh.cell_data[f"N16{ACTIVITY_TAG}"],
            reverted_mesh.cell_data[f"N16{ACTIVITY_TAG}"],
        )

    def test_write_energy_boundaries(self, easy_cdgs: CDGS):

        cdgs = easy_cdgs

        with io.StringIO() as f:
            cdgs._write_energy_boundaries(f, "N16")
            text = f.getvalue()
            lines = text.splitlines()
            assert len(lines) == 3
            n_boundaries1 = int(lines[0].split()[-1])

        with io.StringIO() as f:
            cdgs._write_energy_boundaries(f, "all")
            text = f.getvalue()
            lines = text.splitlines()
            assert len(lines) > 3
            n_boundaries2 = int(lines[0].split()[-1])

        assert n_boundaries2 > n_boundaries1

    def test_write_mesh_boundaries(self, easy_cdgs: CDGS):
        cdgs = easy_cdgs

        with io.StringIO() as f:
            cdgs._write_mesh_boundaries(f)
            text = f.getvalue()
            lines = text.splitlines()
            assert len(lines) == 5

    def test_write_values(self, easy_cdgs: CDGS):
        cdgs = easy_cdgs

        for isotope, vals in cdgs.isotopes.items():
            for energy in vals["lines"]:
                assert energy in cdgs.all_df["lines"].to_list()

        with io.StringIO() as f:
            cdgs._write_values(f, "N16")
            text1 = f.getvalue()
            lines1 = text1.splitlines()

        with io.StringIO() as f:
            cdgs._write_values(f, "all")
            text2 = f.getvalue()
            lines2 = text2.splitlines()

        # higher intensity
        val1 = lines1[0].split()[1]
        val2 = lines2[0].split()[1]
        assert float(val2) > float(val1)

        # longer file
        assert len(cdgs.all_df["lines"]) == len(cdgs.isotopes["N16"]["lines"]) + len(
            cdgs.isotopes["Co60"]["lines"]
        )
        assert len(text2) > len(text1)

    def test_get_probabilities(self, easy_cdgs: CDGS):
        cdgs = easy_cdgs

        probabilities1 = cdgs._get_probabilities("N16", 0)
        probabilities2 = cdgs._get_probabilities("all", 0)

        assert len(probabilities2) > len(probabilities1)

    @pytest.mark.parametrize(
        ["particle", "e_bins", "isotopes"],
        [
            ("gamma", None, {"N16": "n16", "O19": "o19"}),
            ("neutron", None, {"N17": "n17"}),
            ("gamma", [1, 1e4, 1e6], {"N16": "n16", "O19": "o19"}),
        ],
    )
    def test_exports(self, tmp_path, particle, e_bins, isotopes):
        with as_file(files(res).joinpath("test_activity.csv")) as file:
            cdgs = CDGS.from_cloud_point(
                file,
                isotopes,
                interpolation_kernel=SphereKernel(n_voxels=10),
                mesh_definition=MeshByAverageDistance(factor=10),
                col_names={"x": "x", "y": "y", "z": "z", "vol": "cell-volume"},
                particle=particle,
                e_bins=e_bins,
            )
        cdgs.cooling_time = 100
        cdgs.to_cdgs(tmp_path.joinpath("all.cdgs"), "all")
        cdgs.to_cdgs(tmp_path.joinpath("isotope.cdgs"), list(cdgs.isotopes.keys())[0])
        cdgs.to_vtk(tmp_path.joinpath("cdgs.vtk"))


def test_floats_to_multiline_string():
    floats = np.array([1.0, 2.0, 3.0, 4, 7.5, 1])
    result = _floats_to_multiline_string(floats, max_line_length=70)
    expected = "1.0000000e+00 2.0000000e+00 3.0000000e+00 4.0000000e+00 7.5000000e+00\n1.0000000e+00\n"
    assert result == expected
