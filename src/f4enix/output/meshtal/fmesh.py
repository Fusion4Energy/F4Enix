from __future__ import annotations

import logging
import re
import warnings
from collections.abc import Iterator
from copy import deepcopy
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import migjorn
import numpy as np
import pandas as pd

from f4enix.core.constants import PathLike

if TYPE_CHECKING:
    import aci_tallies
    import pyvista as pv

COLUMN_LABELS = ("Value - Total", "Error - Total")
PARTICLES = ("neutron", "photon", "electron", "proton", "deuteron")
ALLOWED_OUTPUT_FORMATS = ["point_cloud", "ip_fluent", "csv", "vtk"]
_VTK_EXTENSIONS = {
    "StructuredGrid": ".vts",
    "UnstructuredGrid": ".vtu",
    "RectilinearGrid": ".vtr",
}
# the digits MCNP prints
_VALUE_FORMAT = "%.5e"
# types of the discrete bins of aci-tallies axes
_LABEL_TYPES = {
    "nuclide": "nuc",
    "parent": "nuc",
    "cell": "cel",
    "decay_time": "dtme",
    "line": "line",
}
# MCNP's default time bin, when there is no TMESH
_NO_TMESH = (0.0, 1e36)


@dataclass(frozen=True, eq=False)
class Bins:
    """The energy or time bins of a mesh (or the user bins in their place).

    Attributes
    ----------
    values : np.ndarray
        boundaries of the bins ("erg", "tme"), or the discrete values.
    type : str
        "erg", "tme", or of discrete values "dtme", "nuc", "cel", "line".
    """

    values: np.ndarray
    type: str

    @classmethod
    def from_aci(cls, axis: aci_tallies.Axis, type: str) -> Bins:
        """The bins of an aci-tallies axis: `type` ("erg" or "tme"), unless discrete."""
        if axis.labels is not None:
            return cls(np.asarray(axis.labels, dtype=float), _LABEL_TYPES[axis.kind])
        if axis.bounds is None:
            raise ValueError(f"Unsupported {type} bins: {axis!r}")
        if type == "tme" and tuple(axis.bounds) == _NO_TMESH:
            return cls(np.array([0.0, 1e20]), type)
        return cls(np.asarray(axis.bounds, dtype=float), type)

    @property
    def binbound(self) -> bool:
        """True if the values are bin boundaries"""
        return self.type in ("erg", "tme")

    @property
    def totalbin(self) -> bool:
        """True if there is a Total bin, after the others"""
        return self.type != "dtme" and len(self) > (2 if self.binbound else 1)

    @property
    def n(self) -> int:
        """Number of bins, Total included"""
        return len(self) - self.binbound + self.totalbin

    def select(self, index: int | None) -> tuple[Bins, range]:
        """Bin `index` (None: all of them), and the indices selected."""
        if index is None or (self.n == 1 and index in (0, -1)):
            return self, range(self.n)
        if not -self.n <= index < self.n:
            raise IndexError(f"bin {index} out of {self.n}")
        index %= self.n
        if self.totalbin and index == self.n - 1:
            # the Total of discrete bins is a bin of -1
            values = np.array([self[0], self[-1]] if self.binbound else [-1.0])
        else:
            values = self[index : index + 1 + self.binbound]
        return Bins(values, self.type), range(index, index + 1)

    def suffix(self, letter: str, index: int) -> str:
        """Suffix of the arrays of bin `index`"""
        if self.n == 1:
            return ""
        if self.totalbin and index == self.n - 1:
            return f"_{letter}Tot"
        return f"_{letter}{index + 1:03d}"

    # as an array of its values
    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, key):
        return self.values[key]

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        return np.array(self.values, dtype=dtype, copy=copy)


class Cylinder(NamedTuple):
    """Where a cylindrical mesh is"""

    origin: np.ndarray
    axis: np.ndarray
    vec: np.ndarray | None


@dataclass(eq=False, repr=False)
class Fmesh:
    """A mesh of any kind: FMESH, CuV, CDGS.

    Attributes
    ----------
    grid : pv.DataSet
        PyVista grid of the mesh. Its cell arrays are the values and errors of the
        mesh, the only copy of them.
    tally : str
        mesh label, typically the tally number.
    geom : str
        "rec" or "cyl".
    x1bin, x2bin, x3bin : np.ndarray
        bounds along X, Y, Z or R, Z, theta (in revolutions).
    ebin, tbin : Bins
        energy and time bins, or the user bins in their place.
    trsf : Cylinder | None
        origin, axis and vec of a cylindrical mesh.
    type : str | None
        type of mesh (meshtal, CDGS, CUV).
    particle, comments : str | None
        transported particle and comments of the mesh.
    strength, cooling_time : float | None
        of a CDGS mesh.
    binlabels : tuple[str, str]
        names of the value and error arrays of the grid (without bin suffixes).
    theta_cells : np.ndarray | None
        number of cells that draw each theta bin of a cylindrical mesh.

    Notes
    -----
    Everything uses the arrays of `grid`: modify them, or assign `data`, to change the
    values. Total bins are not recomputed, nor are errors propagated. The arrays must
    keep their names and the grid its cells, and the cells that draw a same voxel of a
    cylindrical mesh must have the same value.
    """

    grid: pv.DataSet
    tally: str
    geom: str
    x1bin: np.ndarray
    x2bin: np.ndarray
    x3bin: np.ndarray
    ebin: Bins
    tbin: Bins
    trsf: Cylinder | None = None
    type: str | None = None
    particle: str | None = None
    comments: str | None = None
    strength: float | None = None
    cooling_time: float | None = None
    binlabels: tuple[str, str] = COLUMN_LABELS
    theta_cells: np.ndarray | None = None

    @classmethod
    def from_aci(cls, tally: aci_tallies.meshtal.MeshTally, **attributes) -> Fmesh:
        """Build an Fmesh from a mesh tally read by aci-tallies, which is not kept.

        Parameters
        ----------
        tally : aci_tallies.meshtal.MeshTally
            the mesh tally.
        **attributes
            attributes of the Fmesh that are not the ones of the tally.
        """
        geometry = tally.geometry
        kwargs = {
            "particle": next((p for p in PARTICLES if p in tally.particle), ""),
            "comments": tally.comment,
        }
        if geometry.kind == "cylindrical":
            vec = None if geometry.vec is None else np.asarray(geometry.vec, float)
            kwargs["trsf"] = Cylinder(
                np.asarray(geometry.origin, float),
                np.asarray(geometry.axis, float),
                vec,
            )
            # aci-tallies draws a coarse theta bin with several cells
            ray = tally.isel(energy=0, time=0, i=0, j=0)
            kwargs["theta_cells"] = np.array(
                [ray.isel(k=k).to_pyvista().n_cells for k in range(tally.shape[4])]
            )
        fmesh = cls(
            tally.to_pyvista(energy=0, time=0, fields=["result"]),
            str(tally.id),
            "cyl" if geometry.kind == "cylindrical" else "rec",
            *(np.asarray(b, dtype=float) for b in geometry.bounds),
            Bins.from_aci(tally.energy, "erg"),
            Bins.from_aci(tally.time, "tme"),
            **(kwargs | attributes),
        )
        fmesh.grid.clear_cell_data()
        values = (tally.result, tally.rel_error)
        for (it, ie, _, i), name in fmesh._arrays():
            # (energy, time, x1, x2, x3) in aci-tallies
            fmesh.grid.cell_data[name] = fmesh._to_cells(values[i][ie, it].T)
        fmesh._check_cells()
        return fmesh

    @property
    def nx1(self) -> int:
        """number of intervals on X/R axis"""
        return len(self.x1bin) - 1

    @property
    def nx2(self) -> int:
        """number of intervals on Y/Z axis"""
        return len(self.x2bin) - 1

    @property
    def nx3(self) -> int:
        """number of intervals on Z/T axis"""
        return len(self.x3bin) - 1

    @property
    def ne(self) -> int:
        """number of energy bins, Total included"""
        return self.ebin.n

    @property
    def nt(self) -> int:
        """number of time bins, Total included"""
        return self.tbin.n

    @property
    def shape(self) -> tuple[int, ...]:
        """shape of `data`: (time, energy, x3, x2, x1, value or error)"""
        return (self.nt, self.ne, self.nx3, self.nx2, self.nx1, 2)

    def array_name(self, time: int, energy: int, error: bool = False) -> str:
        """Name of the array of the grid with the values (or errors) of a bin."""
        suffixes = self.tbin.suffix("t", time) + self.ebin.suffix("e", energy)
        return self.binlabels[error] + suffixes

    def _arrays(self) -> Iterator[tuple[tuple, str]]:
        """Index in `data`, and name, of each array of the grid"""
        for it, ie, i in np.ndindex(self.nt, self.ne, 2):
            yield (it, ie, ..., i), self.array_name(it, ie, bool(i))

    def _to_cells(self, voxels: np.ndarray) -> np.ndarray:
        """Cell array of the grid from an array of the voxels, (x3, x2, x1). A copy."""
        if self.theta_cells is not None:
            # the cells go along r, then theta, then z
            voxels = np.repeat(voxels, self.theta_cells, axis=0).transpose(1, 0, 2)
        return voxels.flatten()

    def _to_voxels(self, name: str) -> np.ndarray:
        """Cell array `name` of the grid as an array of the voxels, (x3, x2, x1)."""
        cells = self._cells(name)
        if self.theta_cells is None:
            return cells.reshape(self.nx3, self.nx2, self.nx1)
        first = np.cumsum(self.theta_cells) - self.theta_cells
        voxels = cells.reshape(self.nx2, -1, self.nx1)[:, first].transpose(1, 0, 2)
        if not np.array_equal(self._to_cells(voxels), cells, equal_nan=True):
            raise ValueError(
                f"array {name!r} of mesh {self.tally} has different values in the cells "
                "of a same voxel"
            )
        return voxels

    def _cells(self, name: str) -> np.ndarray:
        if name not in self.grid.cell_data:
            raise ValueError(
                f"the grid of mesh {self.tally} has no cell array {name!r}: it was "
                "renamed or removed"
            )
        return np.asarray(self.grid.cell_data[name])

    def _check_cells(self) -> None:
        n_theta = self.nx3 if self.theta_cells is None else self.theta_cells.sum()
        n_cells = n_theta * self.nx2 * self.nx1
        if self.grid.n_cells != n_cells:
            raise ValueError(
                f"the grid of mesh {self.tally} has {self.grid.n_cells} cells instead "
                f"of {n_cells}: it is not a grid of this mesh anymore"
            )

    @property
    def data(self) -> np.ndarray:
        """Values and relative errors, shape (time, energy, x3, x2, x1, 2).

        It is a read-only copy of the arrays of `grid`, made each time. Assign a new
        array to change them all.

        Examples
        --------
        >>> data = fmesh.data.copy()
        >>> data[..., 0] *= 2  # values, errors are data[..., 1]
        >>> fmesh.data = data
        """
        self._check_cells()
        data = np.empty(self.shape)
        for index, name in self._arrays():
            data[index] = self._to_voxels(name)
        data.flags.writeable = False
        return data

    @data.setter
    def data(self, data: np.ndarray) -> None:
        data = np.asarray(data, dtype=float)
        if data.shape != self.shape:
            raise ValueError(f"data has shape {data.shape}, it must be {self.shape}")
        self._check_cells()
        for index, name in self._arrays():
            self.grid.cell_data[name] = self._to_cells(data[index])

    def select_bins(
        self,
        ebin: int | None = None,
        tbin: int | None = None,
        binlabels: tuple[str, str] | None = None,
    ) -> Fmesh:
        """A new Fmesh with only energy bin `ebin` and time bin `tbin` (None: all of
        them). Indices count the Total bin, last: -1 is the Total, when there is one.

        Parameters
        ----------
        ebin : int | None, optional
            energy bin, by default None
        tbin : int | None, optional
            time bin, by default None
        binlabels : tuple[str, str] | None, optional
            names of the value and error arrays of the grid, by default None

        Returns
        -------
        Fmesh
            the selected bins of this mesh, with its own grid.
        """
        self._check_cells()
        new_ebin, energies = self.ebin.select(ebin)
        new_tbin, times = self.tbin.select(tbin)
        # copy everything but the arrays of the grid
        grid = self.grid.copy(deep=False)
        grid.clear_data()
        new = deepcopy(replace(self, grid=grid))
        new.ebin, new.tbin = new_ebin, new_tbin
        new.binlabels = binlabels or COLUMN_LABELS
        for (it, ie, _, i), name in new._arrays():
            old_name = self.array_name(times[it], energies[ie], bool(i))
            new.grid.cell_data[name] = self._cells(old_name).copy()
        return new

    def get_etbin_data(self, ebin: int | None = None, tbin: int | None = None) -> Fmesh:
        """Deprecated, use `select_bins`."""
        warnings.warn(
            "get_etbin_data is deprecated, use select_bins",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.select_bins(ebin=ebin, tbin=tbin)

    def print_info(self) -> dict:
        """Mesh information in a dictionary format.

        Returns
        -------
        dict
            fmesh infos
        """
        names = ("tally", "type", "comments", "particle", "x1bin", "x2bin", "x3bin")
        info = {name: getattr(self, name) for name in (*names, "ebin", "tbin")}
        return info | {"geometry": self.geom}

    def convert2tally(self) -> tuple[int, pd.DataFrame, str]:
        """In case the mesh is 1D (only one spatial binning is higher than 1), convert
        the tally data into a pandas DataFrame.

        Returns
        -------
        tuple[int, pd.DataFrame, str]
            A tuple containing the tally number, a pandas DataFrame with the
            converted data, and a comment string.
        """
        axes = {"Cor A": self.x1bin, "Cor B": self.x2bin, "Cor C": self.x3bin}
        axes = {name: bins for name, bins in axes.items() if len(bins) > 2}
        if len(axes) > 1:
            raise RuntimeError(
                "convert2tally can only be used for 1D meshes, this one has shape "
                f"{self.shape[4:1:-1]}"
            )
        if not axes and len(self.ebin) > 2:
            col, bins = "Energy", self.ebin
            data = self.data[0, :-1, 0, 0, 0]
        elif self.ne > 1 or self.nt > 1:
            raise RuntimeError(
                "convert2tally can only be used for 1D meshes without energy or time "
                "bins"
            )
        else:
            col, bins = next(iter(axes.items()), ("Cor C", self.x3bin))
            data = self.data.reshape(-1, 2)
        df = pd.DataFrame({col: bins[1:], "Value": data[:, 0], "Error": data[:, 1]})
        return int(self.tally), df, str(self.comments)

    def write(
        self,
        outpath: PathLike,
        list_array_names: list[str] | None = None,
        out_format: str = "vtk",
        outfile: str | None = None,
    ) -> None:
        """Export the mesh to a file. vtk, csv, fluent (txt) and point cloud
        (txt) formats can be selected.

        Parameters
        ----------
        outpath : os.PathLike | str
            path to the output folder.
        list_array_names : list[str], optional
            arrays to be exported. The default is None, meaning that all the
            available arrays will be used.
        out_format : str, optional
            output format. The allowed ones are ['point_cloud', 'ip_fluent',
            'csv', 'vtk']. Default is .vtk
        outfile : str, optional
            name of the output file. If specified, overrides the default one.
            Do not include the extension of the file here. Default is None.

        Raises
        ------
        KeyError
            raises KeyError if the output format is not allowed.
        """
        if out_format not in ALLOWED_OUTPUT_FORMATS:
            raise KeyError(
                f"Invalid format, these are the ones allowed: {ALLOWED_OUTPUT_FORMATS}"
            )
        path = str(Path(outpath, outfile or f"Tally_{self.tally}_{out_format}"))
        if out_format == "vtk":
            path += _VTK_EXTENSIONS.get(type(self.grid).__name__, ".vtk")
            self.grid.save(path)
            logging.info(f"{path} created successfully!")
            return

        names = list_array_names
        if names is None:
            names = list(self.grid.array_names)
        if all(name in self.grid.cell_data for name in names):
            points = self.grid.cell_centers().points
        elif all(name in self.grid.point_data for name in names):
            points = self.grid.points
        else:
            raise ValueError(
                "mix between cell and point data is only supported for vtk"
            )

        if out_format == "csv":
            path += ".csv"
            table = np.column_stack([points, *(self.grid[name] for name in names)])
            fmt = ["%.3f", " %.3f", " %.3f"] + [f" {_VALUE_FORMAT}"] * len(names)
            with open(path, "w", newline="") as stream:
                np.savetxt(stream, table, fmt, delimiter=",", newline="\r\n")
            logging.info(f"{path} created successfully!")
            return

        for name in names:
            file = path + "_" + _clean_path(name) + ".txt"
            with open(file, "w") as stream:
                if out_format == "point_cloud":
                    stream.write("x, y, z, value\n")
                    table = np.column_stack([points, self.grid[name]])
                    fmt = ["%.3f"] * 3 + [_VALUE_FORMAT]
                    np.savetxt(stream, table, fmt, delimiter=",")
                else:
                    # fluent: the size, then each coordinate and the values in brackets
                    stream.write(f"3\n{points.shape[1]}\n{len(points)}\n1\nuds-0\n")
                    columns = [(column, "%.3f") for column in points.T]
                    for column, fmt in [*columns, (self.grid[name], _VALUE_FORMAT)]:
                        stream.write("(")
                        np.savetxt(stream, column, fmt)
                        stream.write(")\n")
            logging.info(f"{file} created successfully!")

    def sameMesh(self, other_mesh: Fmesh) -> bool:
        """Check if two meshes are the same: same geometry, in the same position, and
        same number of energy and time bins.

        Parameters
        ----------
        other_mesh : Fmesh
            The other mesh to compare with.

        Returns
        -------
        bool
            True if the meshes are the same, False otherwise.
        """
        return (
            self.geom == other_mesh.geom
            and self.shape == other_mesh.shape
            and self.grid.n_points == other_mesh.grid.n_points
            and np.allclose(self.grid.points, other_mesh.grid.points, rtol=1e-9)
        )

    def apply_transformation(self, tr: migjorn.Transform) -> None:
        """Apply a transformation to the mesh object.

        The mesh is assumed to be defined in the auxiliary coordinate system of
        the transformation and is moved to the main one, as MCNP does.

        Parameters
        ----------
        tr : migjorn.Transform
            transformation card to be applied

        Raises
        ------
        migjorn.ReadError
            If the transformation card is not one MCNP would accept
        """
        import pyvista as pv

        matrix = np.eye(4)
        matrix[:3, :3] = tr.affine.rotation
        matrix[:3, 3] = tr.affine.translation
        if isinstance(self.grid, pv.RectilinearGrid):
            # to avoid a bug regarding rotation
            self.grid = self.grid.cast_to_structured_grid()
        self.grid = self.grid.transform(matrix, inplace=False)

    def rescale_dose_map(
        self, scaling_factors_df: pd.DataFrame, cooling_time_col: str
    ) -> pv.DataSet:
        """Rescale the dose map in the mesh using scaling factors for a chosen cooling time.

        Parameters
        ----------
        scaling_factors_df : pd.DataFrame
            DataFrame, index=daughter, columns=cooling times (as string), contains
            the corresponding scaling factors for each daughter and cooling time
        cooling_time_col : str
            str, column name in scaling_factors_df to use for scaling

        Returns
        -------
        pv.DataSet
            new PyVista grid with rescaled arrays and a new total array
        """
        # the first bin is all the other nuclides
        daughters = {int(b): i + 1 for i, b in enumerate(self.ebin) if i > 0}
        return rescale_dose_map_vtk(
            grid=self.grid,
            scaling_factors_df=scaling_factors_df,
            cooling_time_col=cooling_time_col,
            daughter_to_array=daughters,
            arr_prefix=f"{self.binlabels[0]}_e",
        )


def rescale_dose_map_vtk(
    grid: pv.DataSet,
    scaling_factors_df: pd.DataFrame,
    cooling_time_col: str,
    daughter_to_array: dict[int, int],  # mapping from daughter to array number
    arr_prefix: str,  # prefix for array names
) -> pv.DataSet:
    """Rescale the dose map in the mesh using scaling factors for a chosen cooling time.

    Parameters
    ----------
    scaling_factors_df : pd.DataFrame
        DataFrame, index=daughter, columns=cooling times (as string), contains
        the corresponding scaling factors for each daughter and cooling time
    cooling_time_col : str
        str, column name in scaling_factors_df to use for scaling
    daughter_to_array : dict[int, int]
        mapping from daughter to array number
    arr_prefix : str
        prefix for array names

    Returns
    -------
    pv.DataSet
        new PyVista grid with rescaled arrays and a new total array
    """
    new_grid = deepcopy(grid)
    rescaled_arrays = []
    for daughter, arr_num in daughter_to_array.items():
        arr_name = f"{arr_prefix}{str(arr_num).zfill(3)}"
        if arr_name not in new_grid.cell_data:
            raise KeyError(f"Array {arr_name} not found in grid.cell_data")
        scale = scaling_factors_df.loc[daughter, cooling_time_col]
        new_grid.cell_data[arr_name] = grid.cell_data.get(arr_name) * scale
        rescaled_arrays.append(new_grid.cell_data.get(arr_name))
    total_array = np.sum(np.stack(rescaled_arrays, axis=0), axis=0)
    new_grid.cell_data[f"{arr_prefix}Tot"] = total_array
    return new_grid


def _clean_path(name: str) -> str:
    """Remove special characters from the path to make it a valid file name."""
    return re.sub(r"[\[\]\\\/]", "", name)
