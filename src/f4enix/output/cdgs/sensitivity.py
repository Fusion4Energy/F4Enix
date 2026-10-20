"""
Selection of the interpolation parameters used to build CDGS sources from a
cloud point.

"""

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or -
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at:
    https://eupl.eu/1.2/en/
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from f4enix.output.cdgs.cdgs import _interpolate_cloud_point, _read_cloud_point
from f4enix.output.cdgs.kernels import InterpolationKernel, KNearestKernel
from f4enix.output.cdgs.mesh_definitions import MeshByBinSize, _padding_voxels

AXES = {"x": 0, "y": 1, "z": 2}
# bin sizes tested by default, as multiples of the characteristic cell size
CELL_SIZE_FACTORS = (1, 2, 4)


def get_recommended_parameters(
    cloud_point: str | Path | pd.DataFrame,
    column: str,
    col_names: dict[str, str] | None = None,
    length_unit: str = "m",
    voxel_layers: Sequence[int] = (1, 2, 3),
    bin_sizes: Sequence[float] | None = None,
    slice_axis: str = "y",
    slice_position: float | None = None,
    max_misplaced: float = 0.01,
    vmin: float | None = None,
    vmax: float | None = None,
) -> tuple[dict | None, pd.DataFrame, Figure]:
    """Recommend the bin size of a MeshByBinSize mesh and the k of a
    KNearestKernel to interpolate a cloud point into a CDGS source, comparing
    the candidates on a section of the dataset.

    Each cloud point is spread over the n x n x n cube of its nearest voxels
    (k = n^3), n being the number of voxel layers. Spreading over more layers
    removes the aliasing artefacts that appear when voxels and CFD cells have
    similar sizes, at the cost of blurring the source.

    The candidate bin sizes are by default 1, 2 and 4 times the characteristic
    cell size of the cloud point, i.e. the volume weighted median of the cell
    edge (cell volume^(1/3)).

    To keep it fast, each candidate is evaluated only on the voxel layer
    crossed by the section, using the cloud points of a slab around it that is
    thick enough (n + 2 voxel layers on each side) for the layer to be
    identical to the one obtained interpolating the whole dataset.

    Besides the density column, the cell volumes are distributed with the same
    kernel, giving the fluid "fill" of each voxel (distributed volume / voxel
    volume). A perfect interpolation gives a fill of 1 inside the fluid,
    between 0 and 1 at its boundaries and never above 1. Any excess over 1 is
    source moved to a voxel that cannot physically hold it, which allows to
    quantify the interpolation artefacts without knowing the geometry. For
    each candidate, on the section layer:

    - misplaced_fraction: sum of max(fill - 1, 0) over the total fill, i.e. a
      lower bound of the fraction of the source that has been misplaced.
    - reach [m]: median distance between a cloud point and the farthest voxel
      centroid it contributes to, i.e. how much the source is smeared.

    The recommended candidate is the one with the smallest reach among those
    with misplaced_fraction <= max_misplaced.

    Parameters
    ----------
    cloud_point : str | Path | pd.DataFrame
        Path to the CSV file containing the cloud point data or a DataFrame
        directly. See CDGS.from_cloud_point.
    column : str
        Column of the cloud point containing the density to be interpolated
        (e.g. an isotope atom density), in units of 1/length_unit^3.
    col_names : dict[str, str] | None, optional
        Coordinates and volume column names, see CDGS.from_cloud_point.
    length_unit : str, optional
        Length unit of the cloud point data, one of "m", "cm" or "mm", by
        default "m".
    voxel_layers : Sequence[int], optional
        Numbers of voxel layers n to be tested, by default (1, 2, 3), i.e.
        k = 1, 8, 27.
    bin_sizes : Sequence[float] | None, optional
        Bin sizes to be tested in m. By default, they are inferred from the
        cloud point cell sizes.
    slice_axis : str, optional
        Axis normal to the section, one of "x", "y" or "z", by default "y".
    slice_position : float | None, optional
        Position of the section along slice_axis in m. By default, the median
        coordinate of the cloud point.
    max_misplaced : float, optional
        Maximum misplaced_fraction for a candidate to be recommended, by
        default 0.01. Above about 1%, a speckle is still visible in the
        interpolated density even if the fill looks smooth.
    vmin, vmax : float | None, optional
        Color scale limits for the density, in 1/m3. By default, the 1st and
        99th percentiles of the cloud point in the section.

    Returns
    -------
    tuple[dict | None, pd.DataFrame, Figure]
        The recommended parameters as {"bin_size": bin size in m,
        "voxel_layers": n, "k": n^3} (None if no candidate is acceptable), the
        metrics of all candidates and the figure comparing them on the
        section.

    Examples
    --------
    >>> rec, metrics, fig = get_recommended_parameters(df, "n16")
    >>> cdgs = CDGS.from_cloud_point(
    ...     df,
    ...     {"N16": "n16"},
    ...     interpolation_kernel=KNearestKernel(k=rec["k"]),
    ...     mesh_definition=MeshByBinSize((rec["bin_size"],) * 3),
    ... )
    """
    if slice_axis not in AXES:
        raise ValueError(f"slice_axis must be one of {list(AXES)}")
    if any(n < 1 for n in voxel_layers):
        raise ValueError("voxel_layers must be positive integers")
    axis = AXES[slice_axis]

    df, coords, vols, to_m = _read_cloud_point(cloud_point, col_names, length_unit)
    density = df[column].to_numpy() / to_m**3  # 1/m3
    cell_size = _characteristic_cell_size(vols)
    if bin_sizes is None:
        bin_sizes = [float(f"{f * cell_size:.2g}") for f in CELL_SIZE_FACTORS]
    if slice_position is None:
        slice_position = float(np.median(coords[:, axis]))

    rows = []
    sections = {}
    references = {}
    for b in bin_sizes:
        for n in voxel_layers:
            # the padding of the grid depends on the kernel, but it is made of
            # whole voxels, so the layer crossed by the section is the same
            edges = _full_grid_edges(coords, b, KNearestKernel(k=n**3))
            layer = _layer_index(edges[axis], slice_position)
            references[b] = (coords[:, axis] >= edges[axis][layer]) & (
                coords[:, axis] < edges[axis][layer + 1]
            )
            n_voxels_full = int(np.prod([len(e) - 1 for e in edges]))
            dens, fill, reach = _slab_layer(
                coords, density, vols, edges, axis, layer, n
            )
            rows.append(
                {
                    "bin_size [m]": b,
                    "voxel_layers": n,
                    "k": n**3,
                    "misplaced_fraction": np.sum(np.clip(fill - 1, 0, None))
                    / np.sum(fill),
                    "reach [m]": reach,
                    "n_voxels": n_voxels_full,
                }
            )
            sections[(b, n)] = (dens, fill, edges)

    metrics = pd.DataFrame(rows)
    acceptable = metrics["misplaced_fraction"] <= max_misplaced
    metrics["recommended"] = False
    recommended = None
    if acceptable.any():
        best = metrics.loc[acceptable, "reach [m]"].idxmin()
        metrics.loc[best, "recommended"] = True
        recommended = {
            "bin_size": float(metrics.loc[best, "bin_size [m]"]),
            "voxel_layers": int(metrics.loc[best, "voxel_layers"]),
            "k": int(metrics.loc[best, "k"]),
        }

    fig = _plot(
        coords, density, column, axis, slice_position, cell_size, metrics,
        sections, references, voxel_layers, vmin, vmax,
    )
    return recommended, metrics, fig


def _characteristic_cell_size(vols: np.ndarray) -> float:
    """Volume weighted median of the cell edge (volume^(1/3)). Weighting by
    volume avoids that the many small cells usually found close to walls
    dominate the estimate."""
    sizes = np.cbrt(vols)
    order = np.argsort(sizes)
    cumulative = np.cumsum(vols[order])
    return float(sizes[order][np.searchsorted(cumulative, cumulative[-1] / 2)])


def _full_grid_edges(
    coords: np.ndarray, bin_size: float, kernel: InterpolationKernel
) -> list[np.ndarray]:
    """Edges of the grid that MeshByBinSize((bin_size,) * 3) would build on the
    whole cloud point for the given kernel (see
    RegularMeshDefinition._get_regular_mesh)."""
    steps = np.full(3, bin_size)
    n_pad = _padding_voxels(steps, kernel)
    lo = coords.min(axis=0) - steps * (0.5 + n_pad)
    hi = coords.max(axis=0) + steps * (0.5 + n_pad)
    return [
        np.linspace(lo[i], hi[i], round((hi[i] - lo[i]) / bin_size) + 1)
        for i in range(3)
    ]


def _layer_index(edges: np.ndarray, position: float) -> int:
    """Index of the voxel layer containing position."""
    return int(np.clip(np.searchsorted(edges, position) - 1, 0, len(edges) - 2))


def _slab_layer(
    coords: np.ndarray,
    density: np.ndarray,
    vols: np.ndarray,
    edges: list[np.ndarray],
    axis: int,
    layer: int,
    n: int,
    n_reach_samples: int = 20000,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Interpolate the cloud point with KNearestKernel(k=n^3) on a slab of the
    full grid around the given layer and return the density and fill of that
    layer, together with the reach.

    With k = n^3 a point only contributes to voxels within about n/2 layers
    from it, so the points farther than n + 2 layers cannot reach the target
    layer, and the truncation of the slab cannot affect it.
    """
    margin = n + 2
    first = max(layer - margin, 0)
    last = min(layer + margin, len(edges[axis]) - 2)
    slab_edges = list(edges)
    slab_edges[axis] = edges[axis][first : last + 2]
    in_slab = (coords[:, axis] >= slab_edges[axis][0]) & (
        coords[:, axis] <= slab_edges[axis][-1]
    )
    sub_coords = coords[in_slab]

    # same voxels as the full grid: bbox on its edges and its actual spacing
    bbox = {}
    for i, name in enumerate("xyz"):
        bbox[f"{name}_min"] = slab_edges[i][0]
        bbox[f"{name}_max"] = slab_edges[i][-1]
    steps = tuple(float(e[1] - e[0]) for e in edges)
    kernel = KNearestKernel(k=min(n**3, int(np.prod([len(e) - 1 for e in slab_edges]))))
    mesh, res, vox = _interpolate_cloud_point(
        sub_coords,
        {"q": density[in_slab] * vols[in_slab], "vol": vols[in_slab]},
        kernel,
        MeshByBinSize(steps, bbox=bbox),
    )
    mesh_edges = [np.unique(mesh.points[:, i]) for i in range(3)]
    position = (edges[axis][layer] + edges[axis][layer + 1]) / 2
    dens = _section(res["q"] / vox, mesh_edges, axis, position)
    fill = _section(res["vol"] / vox, mesh_edges, axis, position)

    # reach estimated on the points of the target layer
    centroids = np.asarray(mesh.cell_centers().points)
    in_layer = (sub_coords[:, axis] >= edges[axis][layer]) & (
        sub_coords[:, axis] < edges[axis][layer + 1]
    )
    sample = sub_coords[in_layer]
    if len(sample) > n_reach_samples:
        rng = np.random.default_rng(0)
        sample = sample[rng.choice(len(sample), n_reach_samples, replace=False)]
    if len(sample) == 0:
        return dens, fill, np.nan
    neighbours = kernel.get_neighbours(sample, centroids)
    diff = centroids[neighbours] - sample[:, None, :]
    reach = float(np.median(np.sqrt((diff * diff).sum(axis=2)).max(axis=1)))
    return dens, fill, reach


def _section(
    values: np.ndarray, edges: list[np.ndarray], axis: int, position: float
) -> np.ndarray:
    """Extract the 2D layer of voxel values crossed by the plane
    {axis = position}, with the remaining axes in increasing order."""
    nx, ny, nz = [len(e) - 1 for e in edges]
    # pyvista structured grid cell data are ordered with x varying fastest
    grid = values.reshape(nz, ny, nx).transpose(2, 1, 0)  # [i, j, k]
    return np.take(grid, _layer_index(edges[axis], position), axis=axis)


def _plot(
    coords: np.ndarray,
    density: np.ndarray,
    column: str,
    axis: int,
    slice_position: float,
    cell_size: float,
    metrics: pd.DataFrame,
    sections: dict,
    references: dict,
    voxel_layers: Sequence[int],
    vmin: float | None,
    vmax: float | None,
) -> Figure:
    """One row per bin size: the cloud point in the section layer, then the
    density and fill for each number of voxel layers."""
    ax_u, ax_v = [i for i in range(3) if i != axis]
    names = "xyz"
    bin_sizes = list(references)
    if vmin is None or vmax is None:
        ref = density[np.logical_or.reduce(list(references.values()))]
        ref = ref if len(ref) > 0 else density
        vmin = np.percentile(ref, 1) if vmin is None else vmin
        vmax = np.percentile(ref, 99) if vmax is None else vmax
    dens_kw = dict(cmap="jet", vmin=vmin, vmax=vmax)

    n_cols = 1 + 2 * len(voxel_layers)
    fig, axs = plt.subplots(
        len(bin_sizes), n_cols, figsize=(3.2 * n_cols, 3.4 * len(bin_sizes)),
        squeeze=False, layout="constrained",
    )
    for row, b in zip(axs, bin_sizes):
        mask = references[b]
        sc = row[0].scatter(
            coords[mask, ax_u], coords[mask, ax_v], c=density[mask], s=1, **dens_kw
        )
        row[0].set_title(f"Cloud point in layer\nbin {b * 1e3:.3g} mm", fontsize=9)
        for j, n in enumerate(voxel_layers):
            dens, fill, edges = sections[(b, n)]
            extent = [edges[ax_u][0], edges[ax_u][-1], edges[ax_v][0], edges[ax_v][-1]]
            row[1 + 2 * j].imshow(
                np.ma.masked_equal(dens, 0).T, origin="lower", extent=extent,
                **dens_kw,
            )
            im_fill = row[2 + 2 * j].imshow(
                np.ma.masked_equal(fill, 0).T, origin="lower", extent=extent,
                cmap="RdBu_r", vmin=0, vmax=2,
            )
            m = metrics[(metrics["bin_size [m]"] == b) & (metrics["voxel_layers"] == n)]
            m = m.iloc[0]
            title = (
                f"{'RECOMMENDED ' if m['recommended'] else ''}"
                f"bin {b * 1e3:.3g} mm, n={n} (k={n**3})\n"
                f"misplaced {m['misplaced_fraction']:.1%}, "
                f"reach {m['reach [m]'] * 1e3:.2f} mm"
            )
            row[1 + 2 * j].set_title(
                title, fontsize=9, color="green" if m["recommended"] else "black"
            )
            row[2 + 2 * j].set_title(f"Fill, n={n}", fontsize=9)
    for ax in axs.flat:
        ax.set_aspect("equal")
        ax.set_xlabel(f"{names[ax_u]} [m]")
        ax.set_ylabel(f"{names[ax_v]} [m]")
    fig.colorbar(sc, ax=axs[:, 0].tolist(), label=f"{column} [1/m3]", shrink=0.6)
    fig.colorbar(
        im_fill, ax=axs[:, -1].tolist(),
        label="Fill [-]: 1 full of fluid, <1 partially filled or missing, "
        ">1 overfilled (artefact)",
        shrink=0.6,
    )
    fig.suptitle(
        f"Section at {names[axis]} = {slice_position:.4g} m, characteristic cell "
        f"size {cell_size * 1e3:.3g} mm"
    )
    return fig
