import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt

from f4enix.output.cdgs import (
    KNearestKernel,
    MeshByBinSize,
    get_recommended_parameters,
)
from f4enix.output.cdgs.cdgs import _interpolate_cloud_point
from f4enix.output.cdgs.sensitivity import (
    _characteristic_cell_size,
    _full_grid_edges,
    _layer_index,
    _section,
    _slab_layer,
)

COL_NAMES = {"x": "x", "y": "y", "z": "z", "vol": "vol"}


def _lattice(n: int = 12, spacing: float = 1.0, jitter: float = 0.0) -> pd.DataFrame:
    """Cloud point of the centres of a regular lattice of cubic cells, with a
    density increasing along x. Optionally, centres are randomly displaced."""
    c = (np.arange(n) + 0.5) * spacing
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    coords = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    rng = np.random.default_rng(42)
    coords += rng.uniform(-jitter, jitter, coords.shape) * spacing
    return pd.DataFrame(
        {
            "x": coords[:, 0],
            "y": coords[:, 1],
            "z": coords[:, 2],
            "vol": np.full(len(coords), spacing**3),
            "dens": 1 + coords[:, 0],
        }
    )


class TestGetRecommendedParameters:
    def test_recommendation(self):
        rec, metrics, fig = get_recommended_parameters(
            _lattice(jitter=0.2), "dens", col_names=COL_NAMES
        )
        plt.close(fig)

        # 3 inferred bin sizes x 3 default voxel layers
        assert len(metrics) == 9
        assert sorted(metrics["bin_size [m]"].unique()) == [1.0, 2.0, 4.0]
        assert list(metrics["k"].unique()) == [1, 8, 27]
        assert metrics["recommended"].sum() == 1
        best = metrics[metrics["recommended"]].iloc[0]
        assert rec == {
            "bin_size": best["bin_size [m]"],
            "voxel_layers": best["voxel_layers"],
            "k": best["k"],
        }
        assert best["misplaced_fraction"] <= 0.01

        # for a given bin size, more layers always smear more. Note that they do
        # not necessarily misplace less: the inverse distance weights do not
        # spread each point uniformly, so overlapping contributions can
        # overfill voxels that a single layer would fill exactly
        for _, group in metrics.groupby("bin_size [m]"):
            assert group["reach [m]"].is_monotonic_increasing

    def test_aligned_lattice(self):
        # voxels coinciding with the cells: n=1 has no artefact and the
        # smallest reach, so it is recommended
        rec, metrics, fig = get_recommended_parameters(
            _lattice(), "dens", col_names=COL_NAMES, bin_sizes=[1.0]
        )
        plt.close(fig)
        assert metrics.loc[0, "misplaced_fraction"] == pytest.approx(0)
        assert rec == {"bin_size": 1.0, "voxel_layers": 1, "k": 1}

    def test_no_acceptable_candidate(self):
        rec, metrics, fig = get_recommended_parameters(
            _lattice(jitter=0.2),
            "dens",
            col_names=COL_NAMES,
            voxel_layers=[1],
            bin_sizes=[0.9],
            max_misplaced=0.0,
        )
        plt.close(fig)
        assert rec is None
        assert not metrics["recommended"].any()

    def test_length_unit(self):
        # same lattice in cm: bin sizes are always returned in m
        df = _lattice(spacing=100.0)
        df["dens"] = 1.0
        rec, metrics, fig = get_recommended_parameters(
            df, "dens", col_names=COL_NAMES, length_unit="cm"
        )
        plt.close(fig)
        assert sorted(metrics["bin_size [m]"].unique()) == [1.0, 2.0, 4.0]
        assert rec["bin_size"] == 1.0

    def test_figure(self):
        _, _, fig = get_recommended_parameters(
            _lattice(), "dens", col_names=COL_NAMES, voxel_layers=[1, 2],
            slice_axis="z",
        )
        # 3 bin sizes x (reference + 2 x (density, fill)) panels + colorbars
        assert len(fig.axes) >= 3 * 5
        plt.close(fig)

    @pytest.mark.parametrize(
        "kwargs", [{"slice_axis": "w"}, {"voxel_layers": [0]}]
    )
    def test_invalid_inputs(self, kwargs):
        with pytest.raises(ValueError):
            get_recommended_parameters(
                _lattice(), "dens", col_names=COL_NAMES, **kwargs
            )


@pytest.mark.parametrize("n", [1, 2, 3])
@pytest.mark.parametrize("axis", [0, 2])
def test_slab_layer_matches_full_interpolation(n, axis):
    df = _lattice(n=20, jitter=0.3)
    coords = df[["x", "y", "z"]].to_numpy()
    vols = df["vol"].to_numpy()
    density = df["dens"].to_numpy()
    b = 0.9

    edges = _full_grid_edges(coords, b, KNearestKernel(k=n**3))
    layer = _layer_index(edges[axis], np.median(coords[:, axis]))
    dens, fill, _ = _slab_layer(coords, density, vols, edges, axis, layer, n)

    # reference: interpolation of the whole dataset
    mesh, res, vox = _interpolate_cloud_point(
        coords,
        {"q": density * vols, "vol": vols},
        KNearestKernel(k=n**3),
        MeshByBinSize((b, b, b)),
    )
    mesh_edges = [np.unique(mesh.points[:, i]) for i in range(3)]
    for e_full, e_mesh in zip(edges, mesh_edges):
        assert np.allclose(e_full, e_mesh)
    position = (edges[axis][layer] + edges[axis][layer + 1]) / 2
    assert np.allclose(dens, _section(res["q"] / vox, mesh_edges, axis, position))
    assert np.allclose(fill, _section(res["vol"] / vox, mesh_edges, axis, position))


def test_characteristic_cell_size():
    # many small cells and a few big ones holding most of the volume: the
    # volume weighted median follows the big cells
    vols = np.array([1e-3] * 100 + [1.0] * 10)
    assert _characteristic_cell_size(vols) == pytest.approx(1.0)


def test_section():
    edges = [np.arange(5.0), np.arange(4.0), np.arange(3.0)]  # 4 x 3 x 2 voxels
    i, j, k = np.meshgrid(np.arange(4), np.arange(3), np.arange(2), indexing="ij")
    # cell data with x varying fastest, as in pyvista structured grids
    values = (100 * i + 10 * j + k).transpose(2, 1, 0).ravel()

    sec = _section(values, edges, axis=1, position=1.5)  # j = 1
    assert sec.shape == (4, 2)
    assert np.array_equal(sec, 100 * i[:, 1, :] + 10 + k[:, 1, :])

    sec = _section(values, edges, axis=2, position=0.2)  # k = 0
    assert np.array_equal(sec, 100 * i[:, :, 0] + 10 * j[:, :, 0])
