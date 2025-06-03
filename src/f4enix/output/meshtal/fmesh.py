import csv
import logging
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv
from numjuggler import parser
from tqdm import tqdm

from f4enix.constants import PathLike
from f4enix.output.meshtal.aux_meshtal_functions import (
    ExtraBin,
    _rectilinear_grid,
    _structured_grid,
)

ALLOWED_OUTPUT_FORMATS = ["point_cloud", "ip_fluent", "csv", "vtk"]


class MeshData:
    """class storing mesh and data information (value + error) from MC simulation"""

    def __init__(
        self,
        geom: str,
        x1bin: np.ndarray,
        x2bin: np.ndarray,
        x3bin: np.ndarray,
        ebin: ExtraBin | None = None,
        tbin: ExtraBin | None = None,
        data: np.ndarray | None = None,
    ):
        if geom not in ("rec", "cyl"):
            raise ValueError("bad geometry label")
        self._geom = geom
        self._x1bin = x1bin
        self._x2bin = x2bin
        self._x3bin = x3bin

        self._nx1 = len(x1bin) - 1
        self._nx2 = len(x2bin) - 1
        self._nx3 = len(x3bin) - 1

        if ebin is None:
            self._ebin = ExtraBin(np.array((0, 100)), "erg")
        else:
            self._ebin = ebin
        self._ne = len(self._ebin) - 1 if self._ebin.binbound else len(self._ebin)
        if self._ebin.totalbin:
            self._ne += 1

        if tbin is None:
            self._tbin = ExtraBin(np.array((0, 1e20)), "tme")
        else:
            self._tbin = tbin
        self._nt = len(self._tbin) - 1 if self._tbin.binbound else len(self._tbin)
        if self._tbin.totalbin:
            self._nt += 1

        shape = (self._nt, self._ne, self._nx3, self._nx2, self._nx1, 2)
        self._data = np.ndarray(shape)

        if data is not None:
            if self._data.shape == data.shape:
                self._data = data[:, :, :, :, :, :]
            else:
                raise TypeError("Data dimension doesn't match expected dimension")

    @property
    def geom(self) -> str:
        """
        Mesh geometry ('rec', 'cyl')
        """
        return self._geom

    @property
    def x1bin(self) -> np.ndarray:
        """
        X/R bins data (rec or cyl geometry)
        """
        return self._x1bin

    @property
    def x2bin(self) -> np.ndarray:
        """
        Y/Z bins data (rec or cyl geometry)
        """
        return self._x2bin

    @property
    def x3bin(self) -> np.ndarray:
        """
        Z/T bins data (rec or cyl geometry)
        """
        return self._x3bin

    @property
    def ebin(self) -> ExtraBin:
        """
        mcnp energy bin (or user bin assingned to mcnp energy bin)
        """
        return self._ebin

    @property
    def tbin(self) -> ExtraBin:
        """
        mcnp timebin (or user bin assingned to mcnp time bin)
        """
        return self._tbin

    @property
    def nx1(self) -> int:
        """
        number of intervals on X/R axis
        """
        return self._nx1

    @property
    def nx2(self) -> int:
        """
        number of intervals on Y/Z axis
        """
        return self._nx2

    @property
    def nx3(self) -> int:
        """
        number of intervals on Z/T axis
        """
        return self._nx3

    @property
    def ne(self) -> int:
        """
        number of energy intervals (or energy values if discrete) on energy grid (or defined user bin)
        """
        return self._ne

    @property
    def nt(self) -> int:
        """
        number of time intervals (or time values if discrete) on time grid (or defined user bin)
        """
        return self._nt

    @property
    def data(self) -> np.ndarray:
        """
        array with mesh values and statistical error
        """
        return self._data

    def get_etbin_data(
        self, ebin: int | None = None, tbin: int | None = None
    ) -> "MeshData":
        """Return a sub mesh with only data of selected energy and time bins"""

        if ebin is None:
            ebin = 0
            new_earray = self._ebin
        else:
            if ebin < 0:
                if -ebin < self._ne + 1:
                    ebin += self._ne
                else:
                    raise TypeError("get_etbin_data: energy bin out of range")

            if ebin < self._ne:
                if self._ebin.totalbin and ebin == self._ne - 1:
                    if self._ebin.binbound:
                        new_earray = ExtraBin(
                            np.array((self._ebin[0], self._ebin[-1])),
                            self._ebin.type,
                        )
                    else:
                        new_earray = ExtraBin(np.array((-1,)), self._ebin.type)
                else:
                    if self._ebin.binbound:
                        new_earray = ExtraBin(
                            self._ebin[ebin : ebin + 2], self._ebin.type
                        )
                    else:
                        new_earray = ExtraBin(
                            self._ebin[ebin : ebin + 1], self._ebin.type
                        )
            else:
                raise TypeError("get_etbin_data: energy bin out of range")

        if tbin is None:
            tbin = 0
            new_tarray = self._tbin
        else:
            if tbin < 0:
                if -tbin < self._nt + 1:
                    tbin += self._nt
                else:
                    raise TypeError("get_etbin_data: time bin out of range")

            if tbin < self._nt:
                if self._tbin.totalbin and tbin == self._nt - 1:
                    if self._tbin.binbound:
                        new_tarray = ExtraBin(
                            np.array((self._tbin[0], self._tbin[-1])),
                            self._tbin.type,
                        )
                    else:
                        new_tarray = ExtraBin(np.array((-1,)), self._tbin.type)
                else:
                    if self._tbin.binbound:
                        new_tarray = ExtraBin(
                            self._tbin[tbin : tbin + 2], self._tbin.type
                        )
                    else:
                        new_tarray = ExtraBin(
                            self._tbin[tbin : tbin + 1], self._tbin.type
                        )
            else:
                raise TypeError("get_etbin_data: time bin out of range")

        bin_mesh = MeshData(
            self._geom, self._x1bin, self._x2bin, self._x3bin, new_earray, new_tarray
        )

        if bin_mesh.ne == 1 and bin_mesh.nt == 1:
            bin_mesh.data[0, 0, :, :, :, :] = self._data[tbin, ebin, :, :, :, :]
        elif bin_mesh.ne == 1:
            bin_mesh.data[:, 0, :, :, :, :] = self._data[:, ebin, :, :, :, :]
        elif bin_mesh.nt == 1:
            bin_mesh.data[0, :, :, :, :, :] = self._data[tbin, :, :, :, :, :]
        else:
            bin_mesh.data[:, :, :, :, :, :] = self._data[:, :, :, :, :, :]

        return bin_mesh


class Fmesh(MeshData):
    def __init__(
        self,
        mesh: MeshData,
        meshLabel: str,
        trsf: None | tuple = None,
        binlabels: tuple[str] | None = None,
    ):
        super().__init__(
            mesh.geom,
            mesh.x1bin,
            mesh.x2bin,
            mesh.x3bin,
            mesh.ebin,
            mesh.tbin,
            mesh.data,
        )
        """""Class storing all kind of mesh Fmesh, CuV, CDGS.

        Parameters
        ----------
        mesh : MeshData
            MeshData object containing the mesh geometry, bins, and data.
        meshLabel : str
            Label for the mesh, typically the tally number.
        trsf : tuple | None
            Mesh transformation, if any. It can be a tuple of transformation parameters
            (Origin vector, rotation matrix) or None if no transformation is applied.
        binlabels : tuple[str] | None, optional
            Labels for the bins in the mesh. If None, default labels are used.
            Default is None.

        
        Attributes
        ----------
        grid : pv.DataSet
            PyVista grid object wrapping the mesh data.

        """

        self._trsf = trsf
        self._tally = meshLabel
        self._type = None
        self._comments = None
        self._particle = None
        self.grid: pv.DataSet = self._create_grid(binlabels=binlabels)
        self._strength = None
        self._cooling_time = None

    @property
    def trsf(self):
        """
        Mesh transformation
        """
        return self._trsf

    @property
    def tally(self) -> str:
        """
        Mesh label
        """
        return self._tally

    @property
    def type(self) -> str | None:
        """
        Mesh type (meshtal, CDGS, CUV, ...)
        """
        return self._type

    @type.setter
    def type(self, value: str | None):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"type should be string type not {type(value)}")
        self._type = value

    @property
    def comments(self) -> str | None:
        """
        Mesh comments
        """
        return self._comments

    @comments.setter
    def comments(self, value: str | None):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"comments should be string type not {type(value)}")
        self._comments = value

    @property
    def particle(self) -> str | None:
        """
        Transported particle
        """
        return self._particle

    @particle.setter
    def particle(self, value: str | None):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"particle should be string type not {type(value)}")
        self._particle = value

    @property
    def strength(self) -> float | None:
        """
        Mesh strength (e.g. for CDGS mesh)
        """
        return self._strength

    @strength.setter
    def strength(self, value: float | None):
        if value is not None and not isinstance(value, (float, int)):
            raise TypeError(f"strength should be float type or None not {type(value)}")
        self._strength = value

    @property
    def cooling_time(self) -> float | None:
        """
        Mesh cooling time (e.g. for CDGS mesh)
        """
        return self._cooling_time

    @cooling_time.setter
    def cooling_time(self, value: float | None):
        if value is not None and not isinstance(value, (float, int)):
            raise TypeError(
                f"cooling_time should be float type or None not {type(value)}"
            )
        self._cooling_time = value

    def print_info(self) -> dict:
        """Print mesh information in a dictionary format.

        Returns
        -------
        dict
            fmesh infos
        """
        info = {
            "tally": self.tally,
            "type": self.type,
            "comments": self.comments,
            "particle": self.particle,
            "x1bin": self.x1bin,
            "x2bin": self.x1bin,
            "x3bin": self.x1bin,
            "ebin": self.ebin,
            "tbin": self.tbin,
            "geometry": self.geom,
        }
        return info

    def convert2tally(self) -> tuple[int, pd.DataFrame, str]:
        """In case the mesh is 1D (only one spatial binning is higher than 1), convert
        the tally data into a pandas DataFrame.

        Returns
        -------
        tuple[int, pd.DataFrame, str]
            A tuple containing the tally number, a pandas DataFrame with the
            converted data, and a comment string.
        """
        val = None
        err = None
        # Check if the mesh is 1D
        if self.nx1 == 1 and self.nx2 == 1 and len(self.ebin) > 2:  # 1D Energy mesh
            if (self.geom == "cyl" and self.nx3 == 2) or self.nx3 == 1:
                col = "Energy"
                idx_vals = self.ebin[1:]
                val = self.data[0, : self.ne - 1, 0, 0, 0, 0]
                err = self.data[0, : self.ne - 1, 0, 0, 0, 1]
        elif self.nx1 == self.nx2 == 1:  # z axis
            col = "Cor C"
            idx = self.nx3
            idx_vals = self.x3bin[1:]
        elif self.nx1 == self.nx3 == 1:  # y axis
            col = "Cor B"
            idx = self.nx2
            idx_vals = self.x2bin[1:]
        elif self.nx2 == self.nx3 == 1:  # x axis
            col = "Cor A"
            idx = self.nx1
            idx_vals = self.x1bin[1:]
        elif self.nx1 == 1 and self.nx3 == 2 and self.geom == "cyl":  # y axis
            col = "Cor B"
            idx = self.nx2
            idx_vals = self.x2bin[1:]
        elif self.nx2 == 1 and self.nx3 == 2 and self.geom == "cyl":  # x axis
            col = "Cor A"
            idx = self.nx1
            idx_vals = self.x1bin[1:]
        else:
            raise RuntimeError(
                (
                    "convert2tally can only be used for 1D meshes",
                    f"Axis lengths x:{len(self.x1bin)}, y:{len(self.x2bin)}, z:{len(self.x3bin)}",
                )
            )
        if val is None:
            val = self.grid["Value - Total"][:idx]
            err = self.grid["Error - Total"][:idx]
        data = {
            col: idx_vals,
            "Value": val,
            "Error": err,
        }
        df = pd.DataFrame(data)

        return int(self.tally), df, str(self.comments)

    def _create_grid(self, binlabels=None) -> pv.PolyData:
        if self._geom == "rec" and self._trsf is None:
            if binlabels:
                gd = _rectilinear_grid(self, labels=binlabels)
            else:
                gd = _rectilinear_grid(self)
        elif binlabels:
            gd = _structured_grid(self, trsf=self._trsf, labels=binlabels)
        else:
            gd = _structured_grid(self, trsf=self._trsf)

        pv_grid = pv.wrap(gd)
        assert isinstance(pv_grid, pv.DataSet), "Could not generate PyVista object..."
        return pv_grid

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
        if list_array_names is None:
            list_array_names = list(self.grid.array_names)

        if outfile is None:
            file_name = f"Tally_{self.tally}_{out_format}"
        else:
            file_name = outfile

        # TODO either all cells or all point data are supported if not a vtk
        if out_format != "vtk":
            len_data = len(self.grid.cell_data)
            len_point = len(self.grid.point_data)
            if len_data > 0 and len_point == 0:
                f_points = self.grid.cell_centers().points
            elif len_point > 0:
                f_points = self.grid.points
            else:
                raise ValueError(
                    "mix between cell and point data is only supported for vtk"
                )

        filepath = os.path.join(outpath, file_name)
        mesh_type = str(type(self.grid)).split(".")[-1][:-2]

        if out_format == "vtk":
            self._write_vtk(filepath, mesh_type)
        elif out_format == "csv":
            self._write_csv(filepath, f_points, list_array_names)
        elif out_format == "point_cloud":
            self._write_point_cloud(filepath, list_array_names, f_points)
        elif out_format == "ip_fluent":
            self._write_fluent(filepath, list_array_names, f_points)
        else:
            raise KeyError(
                f"Invalid format, these are the ones allowed: {ALLOWED_OUTPUT_FORMATS}"
            )

    def _write_vtk(self, filepath: PathLike, mesh_type: str) -> None:
        if mesh_type == "StructuredGrid":
            ext = ".vts"
        elif mesh_type == "UnstructuredGrid":
            ext = ".vtu"
        elif mesh_type == "RectilinearGrid":
            ext = ".vtr"
        else:
            ext = ".vtk"

        self.grid.save(str(filepath) + ext)

    def _write_csv(
        self, filepath: PathLike, f_points: np.ndarray, list_array_names: list[str]
    ) -> None:
        new_name = str(filepath) + ".csv"

        with open(new_name, "w", newline="") as outfile:
            writer = csv.writer(outfile)

            # # TODO This may create some issues...
            # values_type = self.get_array_type(list_array_names[0])
            # if values_type == "cells":  # Take points or centers
            #     f_points = self.centers
            # else:  # Points
            #     f_points = self.points

            for i in tqdm(range(len(f_points)), unit=" Points", desc="Writing"):
                csv_points = [
                    f"{f_points[i][0]:.3f}",
                    f" {f_points[i][1]:.3f}",
                    f" {f_points[i][2]:.3f}",
                ]
                for array_name in list_array_names:
                    csv_points.append(f" {self.grid[array_name][i]:.3f}")
                writer.writerow(csv_points)

            logging.info(f"{new_name} created successfully!")

    def _write_point_cloud(
        self, filepath: PathLike, list_array_names: list[str], f_points: np.ndarray
    ) -> None:
        for array_name in list_array_names:
            values = self.grid[array_name]
            new_name = str(filepath) + f"_{_clean_path(array_name)}.txt"
            with open(new_name, "w") as outfile:
                outfile.write("x, y, z, value\n")
                # TODO this can probably be optmized using
                # outfile.writeline()
                for i in tqdm(range(len(f_points)), unit=" Points", desc="Writing"):
                    outfile.write(f"{f_points[i][0]:.3f},")
                    outfile.write(f"{f_points[i][1]:.3f},")
                    outfile.write(f"{f_points[i][2]:.3f},")
                    outfile.write(f"{values[i]:.3f}\n")
            logging.info(f"{new_name} created successfully!")

    def _write_fluent(
        self, filepath: PathLike, list_array_names: list[str], f_points: np.ndarray
    ) -> None:
        for array_name in list_array_names:
            values = self.grid[array_name]
            new_name = str(filepath) + f"_{_clean_path(array_name)}.txt"
            with open(new_name, "w") as outfile:
                guion1 = "3"
                n_coord = f_points.shape[1]  # self.n_coordinates
                n_values = str(len(f_points))
                guion2 = "1"
                uds = "uds-0"
                beginning = f"{guion1}\n{n_coord}\n{n_values}\n{guion2}\n{uds}\n"
                outfile.write(beginning)
                outfile.write("(")
                for i in tqdm(range(len(f_points)), unit=" x points", desc="Writing x"):
                    outfile.write(f"{f_points[i][0]:.3f}\n")

                outfile.write(")\n")
                outfile.write("(")

                for i in tqdm(range(len(f_points)), unit=" y points", desc="Writing y"):
                    outfile.write(f"{f_points[i][1]:.3f}\n")
                outfile.write(")\n")
                outfile.write("(")
                for i in tqdm(range(len(f_points)), unit=" z points", desc="Writing z"):
                    outfile.write(f"{f_points[i][2]:.3f}\n")
                outfile.write(")\n")
                outfile.write("(")
                for i in tqdm(
                    range(len(f_points)), unit=" values", desc="Writing values"
                ):
                    outfile.write(f"{values[i]:.3f}\n")
                outfile.write(")\n")
            logging.info(f"{new_name} created successfully!")

    def sameMesh(self, other_mesh: MeshData) -> bool:
        """Check if two meshes are the same based on their geometry and bins.

        Parameters
        ----------
        other_mesh : MeshData
            The other mesh to compare with.

        Returns
        -------
        bool
            True if the meshes are the same, False otherwise.
        """
        if self.geom != other_mesh.geom:
            return False
        if self.data.shape != other_mesh.data.shape:
            return False

        mesh1_xbins = (self.x1bin, self.x2bin, self.x3bin)
        mesh2_xbins = (other_mesh.x1bin, other_mesh.x2bin, other_mesh.x3bin)

        for bin1, bin2 in zip(mesh1_xbins, mesh2_xbins):
            if np.any(abs(bin1 - bin2) > 1e-12):
                return False

        return True

    def _read_from_vtk(self, vtk_file: os.PathLike) -> None:
        # This is mostly used for quicker testing
        grid = pv.read(vtk_file)
        self.grid = grid

    def apply_transformation(self, tr: parser.Card) -> None:
        """Apply a transformation to the mesh object

        Parameters
        ----------
        tr : parser.Card
            transformation card to be applied to the fmesh

        Raises
        ------
        ValueError
            If a non-transformation card is passed
        ValueError
            If the transformation card has not 4 (translation) or 13 (affine transformation) values
        """
        if tr.ctype != 5:
            raise ValueError("Numjuggler card is not a transformation")
        if len(tr.values) not in [4, 13]:
            raise ValueError(
                "Numjuggler transformation card has not 4 (translation) or 13 (rototranslation) values"
            )
        transf_values = []

        if len(tr.values) == 13:
            for k, val in enumerate(tr.values[1:]):
                if k < 3:
                    transf_values.append(val[0])
                elif tr.unit == "*":
                    transf_values.append(np.cos(np.radians(val[0])))
                else:
                    transf_values.append(val[0])
            transf_matrix_dcm = np.array(
                [
                    transf_values[3:6],
                    transf_values[6:9],
                    transf_values[9:],
                ]
            )
            # Compute the transpose of the matrix
            transposed_matrix = np.transpose(transf_matrix_dcm)

            # Compute the inverse of the transposed matrix
            try:
                inverted_transposed_matrix = np.linalg.inv(transposed_matrix)
            except np.linalg.LinAlgError:
                print("The transformation matrix is not invertible.")
            transform_matrix = np.array(
                [
                    [
                        inverted_transposed_matrix[0][0],
                        inverted_transposed_matrix[0][1],
                        inverted_transposed_matrix[0][2],
                        transf_values[0],
                    ],
                    [
                        inverted_transposed_matrix[1][0],
                        inverted_transposed_matrix[1][1],
                        inverted_transposed_matrix[1][2],
                        transf_values[1],
                    ],
                    [
                        inverted_transposed_matrix[2][0],
                        inverted_transposed_matrix[2][1],
                        inverted_transposed_matrix[2][2],
                        transf_values[2],
                    ],
                    [0, 0, 0, 1],
                ]
            )
        else:
            transform_matrix = np.array(
                [
                    [1, 0, 0, tr.values[1][0]],
                    [0, 1, 0, tr.values[2][0]],
                    [0, 0, 1, tr.values[3][0]],
                    [0, 0, 0, 1],
                ]
            )
        self.grid = self.grid.transform(transform_matrix, inplace=False)


def _clean_path(name: str) -> str:
    """Remove special characters from the path to make it a valid file name."""
    pat_wrong_file_char = re.compile(r"[\[\]\\\/]")
    fixed_name = pat_wrong_file_char.sub("", name)
    return fixed_name
