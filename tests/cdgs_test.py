from importlib.resources import as_file, files
import logging
import numpy as np
import pandas as pd
import pyvista as pv
import actigamma as ag
import io
import pytest

import tests.resources.cdgs as res
from f4enix.output.cdgs import (
    CDGS,
    CDGS_ENERGY_TYPE,
    SphereKernel,
    DistanceKernel,
    KNearestKernel,
    MeshByAverageDistance,
    MeshByNumberOfVoxels,
    MeshByBinSize,
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

        interp1 = np.zeros(dest_coords.shape[0])
        kernel._kernel(
            source_coord=source_coords[0],
            dest_idx=neighbours_idx[0],
            dest_coords=dest_coords,
            interpolation=interp1,
            activity=1.0,
        )

        interp2 = np.zeros(dest_coords.shape[0])
        kernel = DistanceKernel(radius=0.75)
        kernel._kernel(
            source_coord=source_coords[0],
            dest_idx=neighbours_idx[0],
            dest_coords=dest_coords,
            interpolation=interp2,
            activity=1.0,
        )

        assert np.not_equal(
            interp1, interp2
        ).any()  # The two interpolations should differ

    def test_get_neighbours_k_nearest_kernel(self):
        kernel = KNearestKernel(k=2)
        source_coords = np.array([[0, 0, 0], [0.1, 0.1, 0.1], [1, 1, 1]])
        dest_coords = np.array([[0, 0, 0], [0.4, 0.4, 0.4], [0.9, 0.9, 0.9]])
        neighbours_idx = kernel.get_neighbours(source_coords, dest_coords)

        # always a dense (n_points, k) array, unlike the ragged SphereKernel result
        assert neighbours_idx.shape == (3, 2)
        assert list(neighbours_idx[0]) == [0, 1]
        assert list(neighbours_idx[1]) == [0, 1]
        assert list(neighbours_idx[2]) == [2, 1]

    def test_k_nearest_kernel_invalid_k(self):
        with pytest.raises(ValueError):
            KNearestKernel(k=0)

        dest_coords = np.array([[0, 0, 0], [0.4, 0.4, 0.4]])
        with pytest.raises(ValueError):
            KNearestKernel(k=3).get_neighbours(dest_coords, dest_coords)

    def test_k_nearest_kernel_weights(self):
        kernel = KNearestKernel(k=2)
        dest_coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [3.0, 0, 0]])

        # coincident point: full weight goes to the collocated destination
        interp = np.zeros(dest_coords.shape[0])
        kernel._kernel(
            source_coord=np.array([0.0, 0, 0]),
            dest_idx=np.array([0, 1]),
            dest_coords=dest_coords,
            interpolation=interp,
            activity=2.0,
        )
        assert np.allclose(interp, [2.0, 0.0, 0.0])

        # non-coincident: inverse-distance weights, summing to the full activity
        interp = np.zeros(dest_coords.shape[0])
        kernel._kernel(
            source_coord=np.array([0.5, 0, 0]),
            dest_idx=np.array([0, 1]),
            dest_coords=dest_coords,
            interpolation=interp,
            activity=2.0,
        )
        assert np.isclose(interp.sum(), 2.0)
        assert np.isclose(interp[0], interp[1])  # equidistant -> equal weights

    def test_distance_kernel_coincident_point(self):
        kernel = DistanceKernel(radius=1.0)
        dest_coords = np.array([[0.0, 0, 0], [0.5, 0, 0]])
        interp = np.zeros(dest_coords.shape[0])
        kernel._kernel(
            source_coord=np.array([0.0, 0, 0]),
            dest_idx=np.array([0, 1]),
            dest_coords=dest_coords,
            interpolation=interp,
            activity=10.0,
        )
        # full weight to the collocated destination, activity is conserved
        assert np.allclose(interp, [10.0, 0.0])

    def test_k_nearest_kernel_batch_matches_loop(self):
        kernel = KNearestKernel(k=2)
        dest_coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [3.0, 0, 0], [0.5, 1.0, 0]])
        # includes a source point coincident with a destination point (index 0)
        source_coords = np.array([[0.0, 0, 0], [0.5, 0, 0], [2.0, 0, 0]])
        activities = np.array([2.0, 3.0, 1.5])
        neighbours_idx = kernel.get_neighbours(source_coords, dest_coords)

        interp_loop = np.zeros(dest_coords.shape[0])
        for i in range(len(source_coords)):
            kernel._kernel(
                source_coord=source_coords[i],
                dest_idx=neighbours_idx[i],
                dest_coords=dest_coords,
                interpolation=interp_loop,
                activity=activities[i],
            )

        interp_batch = np.zeros(dest_coords.shape[0])
        kernel._kernel_batch(
            source_coords=source_coords,
            neighbours_idx=neighbours_idx,
            dest_coords=dest_coords,
            interpolation=interp_batch,
            activities=activities,
        )

        assert np.allclose(interp_loop, interp_batch)


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

    @staticmethod
    def _simple_mesh() -> pv.StructuredGrid:
        x = np.arange(3, dtype=float)
        xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
        mesh = pv.StructuredGrid(xx, yy, zz)
        mesh.cell_data[f"Co60{ACTIVITY_TAG}"] = np.ones(mesh.n_cells)
        mesh.cell_data[f"Co60{ATOM_DENSITY_TAG}"] = np.ones(mesh.n_cells)
        return mesh

    def test_init_does_not_modify_input_mesh(self):
        mesh = self._simple_mesh()
        names = list(mesh.array_names)
        CDGS(mesh)
        assert list(mesh.array_names) == names

    def test_init_no_isotopes(self):
        x = np.arange(3, dtype=float)
        mesh = pv.StructuredGrid(*np.meshgrid(x, x, x, indexing="ij"))
        with pytest.raises(ValueError):
            CDGS(mesh)

    def test_init_bins_without_e_bins(self):
        with pytest.raises(ValueError):
            CDGS(self._simple_mesh(), energy_type=CDGS_ENERGY_TYPE.BINS)

    def test_bins_warn_lost_emission(self, caplog):
        # Co60 main lines are at 1.17 and 1.33 MeV, outside a 0-1 MeV grid
        with caplog.at_level(logging.WARNING):
            CDGS(
                self._simple_mesh(),
                energy_type=CDGS_ENERGY_TYPE.BINS,
                e_bins=np.array([0, 1e6]),
            )
        assert "Co60" in caplog.text
        assert "outside the provided energy bins" in caplog.text

        # no warning when all the lines are inside the bins
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            CDGS(
                self._simple_mesh(),
                energy_type=CDGS_ENERGY_TYPE.BINS,
                e_bins=np.array([0, 1e7]),
            )
        assert "outside the provided energy bins" not in caplog.text

    def test_cooling_time_numpy_scalar(self):
        cdgs = CDGS(self._simple_mesh())
        cdgs.cooling_time = np.int64(10)
        assert cdgs.cooling_time == 10.0
        with pytest.raises(TypeError):
            cdgs.cooling_time = True

    def test_activity_consistent_with_atoms(self, cdgs: CDGS):
        db = ag.Decay2012Database()
        for isotope in ["N16", "O19"]:
            lam = np.log(2) / db.gethalflife(isotope)
            assert np.allclose(
                cdgs.mesh.cell_data[f"{isotope}{ACTIVITY_TAG}"],
                cdgs.mesh.cell_data[f"{isotope}{ATOM_DENSITY_TAG}"] * lam,
            )

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

    def test_get_probabilities_all_lines_weighting(self, easy_cdgs: CDGS):
        cdgs = easy_cdgs
        # each isotope's share of the photons is activity * photons per decay
        # (N16 activity is 1, Co60 activity is 2 in every voxel)
        prob = cdgs._get_probabilities("all", 0)
        assert np.isclose(prob.sum(), 1.0)
        yields = {
            "N16": 1 * cdgs.isotopes["N16"]["intensities"].sum(),
            "Co60": 2 * cdgs.isotopes["Co60"]["intensities"].sum(),
        }
        tot = sum(yields.values())
        for iso, y in yields.items():
            mask = (cdgs.all_df["isotope"] == iso).to_numpy()
            assert np.isclose(prob[mask].sum(), y / tot)

    @pytest.mark.parametrize("isotope", ["all", "N16"])
    def test_get_probabilities_vectorized(self, easy_cdgs: CDGS, isotope):
        cdgs = easy_cdgs
        cdgs.mesh.cell_data[f"N16{ACTIVITY_TAG}"] = np.arange(1, cdgs.mesh.n_cells + 1)
        idx = np.array([0, 3, 5])
        batch = cdgs._get_probabilities(isotope, idx)
        for row, i in zip(batch, idx):
            assert np.allclose(row, cdgs._get_probabilities(isotope, i))

    def test_get_probabilities_all_bins(self):
        x = np.arange(2, dtype=float)
        mesh = pv.StructuredGrid(*np.meshgrid(x, x, x, indexing="ij"))
        activities = {"Co60": 3.0, "Cs137": 1.0}
        for iso, act in activities.items():
            mesh.cell_data[f"{iso}{ACTIVITY_TAG}"] = np.full(mesh.n_cells, act)
            mesh.cell_data[f"{iso}{ATOM_DENSITY_TAG}"] = np.ones(mesh.n_cells)
        e_bins = np.linspace(0, 2e6, 21)
        cdgs = CDGS(mesh, energy_type=CDGS_ENERGY_TYPE.BINS, e_bins=e_bins)

        # reference: aggregate the mixed inventory directly with actigamma
        db = ag.Decay2012Database()
        inv = ag.UnstablesInventory(
            data=[(db.getzai(iso), act) for iso, act in activities.items()]
        )
        hist, _ = ag.LineAggregator(db, ag.EnergyGrid(bounds=e_bins))(
            inv, spectype="gamma"
        )
        assert np.allclose(cdgs._get_probabilities("all", 0), hist / np.sum(hist))

    def test_mesh_boundaries_and_volume_in_cm(self, easy_cdgs: CDGS):
        # easy_cdgs mesh is defined in m: x, y edges 0-3 m, z edges 0-1 m,
        # voxels of 1 m3
        with io.StringIO() as f:
            easy_cdgs._write_mesh_boundaries(f)
            lines = f.getvalue().splitlines()
        assert np.allclose(np.array(lines[2].split(), float), [0, 100, 200, 300])
        assert np.allclose(np.array(lines[4].split(), float), [0, 100])

        with io.StringIO() as f:
            easy_cdgs._write_values(f, "N16")
            first = f.getvalue().splitlines()[0].split()
        assert np.isclose(float(first[2]), 1e6)  # cm3

    @pytest.mark.parametrize("length_unit, to_m", [("cm", 1e-2), ("mm", 1e-3)])
    def test_from_cloud_point_length_unit(self, length_unit, to_m):
        with as_file(files(res).joinpath("test_activity.csv")) as file:
            df = pd.read_csv(file)
        df.columns = df.columns.str.strip()
        # same data expressed in a different length unit
        df_scaled = df.copy()
        df_scaled[["x", "y", "z"]] /= to_m
        df_scaled["cell-volume"] /= to_m**3
        df_scaled["n16"] *= to_m**3

        kwargs = dict(
            isotopes={"N16": "n16"},
            interpolation_kernel=KNearestKernel(k=4),
            mesh_definition=MeshByAverageDistance(factor=10),
            col_names={"x": "x", "y": "y", "z": "z", "vol": "cell-volume"},
        )
        ref = CDGS.from_cloud_point(df, **kwargs)
        scaled = CDGS.from_cloud_point(df_scaled, length_unit=length_unit, **kwargs)

        assert np.allclose(ref.mesh.points, scaled.mesh.points)
        assert np.isclose(ref.voxel_vol, scaled.voxel_vol)
        for tag in [ACTIVITY_TAG, ATOM_DENSITY_TAG]:
            assert np.allclose(
                ref.mesh.cell_data[f"N16{tag}"], scaled.mesh.cell_data[f"N16{tag}"]
            )
        assert np.isclose(
            ref.get_total_intensity("all"), scaled.get_total_intensity("all")
        )

    def test_from_cloud_point_invalid_length_unit(self):
        df = pd.DataFrame({"x": [0.0], "y": [0.0], "z": [0.0], "v": [1.0], "a": [1.0]})
        with pytest.raises(ValueError):
            CDGS.from_cloud_point(
                df,
                {"N16": "a"},
                interpolation_kernel=KNearestKernel(k=1),
                mesh_definition=MeshByBinSize((1.0, 1.0, 1.0)),
                col_names={"x": "x", "y": "y", "z": "z", "vol": "v"},
                length_unit="inch",
            )

    def test_write_values_chunking(self, easy_cdgs: CDGS):
        texts = []
        for chunk_size in [2, 10000]:
            with io.StringIO() as f:
                easy_cdgs._write_values(f, "all", chunk_size=chunk_size)
                texts.append(f.getvalue())
        assert texts[0] == texts[1]

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

    def test_from_cloud_point_k_nearest_kernel(self):
        with as_file(files(res).joinpath("test_activity.csv")) as file:
            cdgs = CDGS.from_cloud_point(
                file,
                {"N16": "n16"},
                interpolation_kernel=KNearestKernel(k=4),
                mesh_definition=MeshByAverageDistance(factor=10),
                col_names={"x": "x", "y": "y", "z": "z", "vol": "cell-volume"},
            )
        assert np.sum(cdgs.mesh.cell_data[f"N16{ACTIVITY_TAG}"]) > 0


def test_floats_to_multiline_string():
    floats = np.array([1.0, 2.0, 3.0, 4, 7.5, 1])
    result = _floats_to_multiline_string(floats, max_line_length=70)
    expected = "1.0000000e+00 2.0000000e+00 3.0000000e+00 4.0000000e+00 7.5000000e+00\n1.0000000e+00\n"
    assert result == expected
