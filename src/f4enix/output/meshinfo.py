"""This module is related to the parsing of D1S-UNED meshinfo files.

"""

from __future__ import annotations

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
import json
import os
import re
from typing import Optional
from enum import Enum
import numpy as np
from pathlib import Path
from typing import Optional
import pandas as pd


# PARAMETERS
INDEX_KEY_TIME = "Time"
INDEX_KEY_VOXEL = "Voxel"
INDEX_KEY_CELL = "Cell"
INDEX_KEY_MATERIAL = "Material"
COLUMN_KEY_MASS_GRAMS = "Mass [g]"


class CoordinateType(Enum):
    """Only two acceptable coordinate types exist"""

    CART = "cart"
    CYL = "cyl"


class DataMass:

    def __init__(self, dataframe: pd.DataFrame) -> None:
        """
        The dataframe contains the mass information present in the meshinfo
        file. It can be used either to know the masses of the desired
        Voxel-Cell combinations or to relate the material ids with the MCNP
        cell ids and vice-versa.

        Parameters
        ----------
        dataframe : pd.DataFrame
            mass dataframe

        Attributes
        ----------
        dataframe : pd.DataFrame
            mass dataframe

        Examples
        --------
        Example of mass dataframe

        >>> datamass.df
                                Mass [g]
        Voxel Material Cell
        1     0        327608     0.000000
                    327611     0.000000
                        ...
            110      324893  4977.059766
                    324894  1332.655347
                    324896  8458.281897

        """
        self.df: pd.DataFrame = dataframe.sort_index()
        self._check_dataframe_format()

    def _check_dataframe_format(self) -> None:
        try:
            assert self.df.index.names == [
                INDEX_KEY_VOXEL,
                INDEX_KEY_MATERIAL,
                INDEX_KEY_CELL,
            ]
            assert self.df.columns.tolist() == [COLUMN_KEY_MASS_GRAMS]
        except AssertionError as error:
            raise ValueError(
                "The format of the dataframe is not correct. "
                "It has to be a MultiIndex Dataframe with index: "
                f"['{INDEX_KEY_VOXEL}',"
                f" '{INDEX_KEY_MATERIAL}',"
                f" '{INDEX_KEY_CELL}']"
                f" and a column: '{COLUMN_KEY_MASS_GRAMS}'"
            ) from error

    @property
    def materials(self) -> np.ndarray:
        """Returns the unique values of Material index"""
        material_index = self.df.index.names.index(INDEX_KEY_MATERIAL)
        return self.df.index.levels[material_index].values

    def get_filtered_dataframe(
        self,
        voxels: Optional[list] = None,
        materials: Optional[list] = None,
        cells: Optional[list] = None,
    ) -> pd.DataFrame:
        """Returns a copy of the dataframe, if arguments are given, only the
        given values will appear.

        Parameters
        ----------
        voxels : Optional[list], optional
            voxels to filter for, by default None
        materials : Optional[list], optional
            materials to filter for, by default None
        cells : Optional[list], optional
            cells to filter for, by default None

        Returns
        -------
        pd.DataFrame
            filtered dataframe
        """
        mask = np.full(len(self.df), True)
        for key, filter_values in zip(
            [INDEX_KEY_VOXEL, INDEX_KEY_MATERIAL, INDEX_KEY_CELL],
            [voxels, materials, cells],
        ):
            if filter_values is None:
                continue
            mask *= self.df.index.isin(filter_values, level=key)
        return self.df.iloc[mask]

    def get_cells_from_materials(self, materials: list[int] | None) -> pd.DataFrame:
        """
        Returns an array with all the cell ids that are filled with the
        materials given as argument.

        Parameters
        ----------
        materials : list[int] | None
            materials to be considered for the filter

        Returns
        -------
        pd.DataFrame
            data related to the requested materials
        """
        filtered_dataframe = self.get_filtered_dataframe(materials=materials)
        return filtered_dataframe.index.unique(level=INDEX_KEY_CELL).values

    def save_dataframe_to_hdf5(self, results_path: os.PathLike) -> None:
        """Saves the dataframe for later loading with the load class method.

        Parameters
        ----------
        results_path : os.PathLike
            path to the output folder
        """
        outfile = os.path.join(results_path, "data_mass.hdf5")
        self.df.to_hdf(outfile, key="df")

    @classmethod
    def load(cls, path_to_folder: os.PathLike) -> DataMass:
        """Returns an instance of the class that was saved with the
        :py:meth:`save_dataframe_to_hdf5` method.

        Parameters
        ----------
        path_to_folder : os.PathLike
            path to the folder where the hdf5 file was saved.

        Returns
        -------
        DataMass
            Loaded DataMass object
        """
        infile = os.path.join(path_to_folder, "data_mass.hdf5")
        dataframe = pd.read_hdf(infile)
        return DataMass(dataframe)


class MeshInfo:
    def __init__(
        self,
        coordinates: CoordinateType,
        data_mass: DataMass,
        vector_i: list[float],
        vector_j: list[float],
        vector_k: list[float],
    ) -> None:
        """This class represents the information read in the meshinfo file
        produced by D1SUNED when using a CuV tally. This class represents a
        cartesian mesh or serves as base class for a cylindrical mesh
        (cylindrical meshes have more parameters).

        Parameters
        ----------
        coordinates : CoordinateType
            either cartesian or cylindrical
        data_mass : DataMass
            mass data associated to the mesh
        vector_i : list[float]
            i vector
        vector_j : list[float]
            j vector
        vector_k : list[float]
            k vector

        Attributes
        ----------
        coordinates : CoordinateType
            either cartesian or cylindrical
        data_mass : DataMass
            mass data associated to the mesh
        vector_i : np.array
            i vector
        vector_j : np.array
            j vector
        vector_k : np.array
            k vector
        """

        self.coordinates = coordinates
        self.data_mass = data_mass
        self.vector_i = np.array(vector_i)
        self.vector_j = np.array(vector_j)
        self.vector_k = np.array(vector_k)
        # self.transformation_matrix = None

    def save(self, path_to_results: os.PathLike) -> None:
        """To save the class all the parameters are saved in a JSON file except
        the data_mass which is saved as a HDF5 dataframe.

        Parameters
        ----------
        path_to_results : os.PathLike
            path to the output folder
        """
        json_data = {
            "coordinates": self.coordinates.value,
            "vector_i": self.vector_i.tolist(),
            "vector_j": self.vector_j.tolist(),
            "vector_k": self.vector_k.tolist(),
            # "transformation_matrix": self.transformation_matrix.tolist(),
        }
        outfile = os.path.join(path_to_results, "meshinfo_parameters.json")
        with open(outfile, "w") as infile:
            json.dump(json_data, infile)
        if self.data_mass is not None:
            self.data_mass.save_dataframe_to_hdf5(path_to_results)

    @classmethod
    def load(cls, path_to_folder: os.PathLike) -> MeshInfo | MeshInfoCyl:
        """Given the path to the folder where the class was saved, it reads both
        the JSON file and the data_mass dataframe to recreate the instance of
        the class.

        Parameters
        ----------
        path_to_folder : os.PathLike
            folder where the files to be loaded are stored

        Returns
        -------
        MeshInfo | MeshInfoCyl
            Loaded mesh info
        """
        data_mass = DataMass.load(path_to_folder)
        with open(path_to_folder / "meshinfo_parameters.json", "r") as infile:
            json_data = json.load(infile)
        coordinates = CoordinateType(json_data["coordinates"])
        if coordinates == CoordinateType.CART:
            return MeshInfo(
                coordinates=coordinates,
                data_mass=data_mass,
                vector_i=json_data["vector_i"],
                vector_j=json_data["vector_j"],
                vector_k=json_data["vector_k"],
            )
        return MeshInfoCyl(
            coordinates=coordinates,
            data_mass=data_mass,
            vector_i=json_data["vector_i"],
            vector_j=json_data["vector_j"],
            vector_k=json_data["vector_k"],
            origin=json_data["origin"],
            axis=json_data["axis"],
            vec=json_data["vec"],
        )

    def __eq__(self, __value: MeshInfo | MeshInfoCyl) -> bool:
        if type(__value) not in [MeshInfo, MeshInfoCyl]:
            return False

        try:
            assert self.coordinates == __value.coordinates
            assert (self.vector_i == __value.vector_i).all()
            assert (self.vector_j == __value.vector_j).all()
            assert (self.vector_k == __value.vector_k).all()
            assert self.data_mass.df.equals(__value.data_mass.df)
            if type(self) == MeshInfoCyl:
                assert (self.origin == __value.origin).all()
                assert (self.axis == __value.axis).all()
                assert (self.vec == __value.vec).all()

        except AssertionError:
            return False

        return True


class MeshInfoCyl(MeshInfo):

    def __init__(
        self,
        coordinates: CoordinateType,
        data_mass: Optional[DataMass],
        vector_i: list[float],
        vector_j: list[float],
        vector_k: list[float],
        origin: list[float],
        axis: list[float],
        vec: list[float],
    ) -> None:
        """Same as MeshInfo but adding cylindrical coordinates specific
        parameters like origin, axis and vec.

        Parameters
        ----------
        coordinates : CoordinateType
            either cartesian or cylindrical
        data_mass : DataMass
            mass data associated to the mesh
        vector_i : list[float]
            i vector
        vector_j : list[float]
            j vector
        vector_k : list[float]
            k vector
        origin : list[float]
            x, y, z origin of the mesh
        axis : list[float]
            TBD
        vec : list[float]
            TBD

        Attributes
        ----------
        coordinates : CoordinateType
            either cartesian or cylindrical
        data_mass : DataMass
            mass data associated to the mesh
        vector_i : np.array
            i vector
        vector_j : np.array
            j vector
        vector_k : np.array
            k vector
        origin : np.array
            x, y, z origin of the mesh
        axis : np.array
            TBD
        vec : np.array
            TBD
        """
        super().__init__(
            coordinates=coordinates,
            data_mass=data_mass,
            vector_i=vector_i,
            vector_j=vector_j,
            vector_k=vector_k,
        )
        self.origin = np.array(origin)
        self.axis = np.array(axis)
        self.vec = np.array(vec)

    def save(self, path_to_results: os.PathLike) -> None:
        """To save the class all the parameters are saved in a JSON file except
        the data_mass which is saved as a HDF5 dataframe.

        Parameters
        ----------
        path_to_results : os.PathLike
            _description_
        """
        json_data = {
            "coordinates": self.coordinates.value,
            "vector_i": self.vector_i.tolist(),
            "vector_j": self.vector_j.tolist(),
            "vector_k": self.vector_k.tolist(),
            # "transformation_matrix": self.transformation_matrix.tolist(),
            "origin": self.origin.tolist(),
            "axis": self.axis.tolist(),
            "vec": self.vec.tolist(),
        }
        outfile = os.path.join(path_to_results, "meshinfo_parameters.json")
        with open(outfile, "w") as infile:
            json.dump(json_data, infile)
        if self.data_mass is not None:
            self.data_mass.save_dataframe_to_hdf5(path_to_results)


class MeshInfoFile:
    def __init__(self, info: dict[int, MeshInfo | MeshInfoCyl]) -> None:
        """Parser of the meshinfo D1S-UNED file

        Parameters
        ----------
        info : dict[int, MeshInfo  |  MeshInfoCyl]
            contains the info related to the various meshes

        Attributes
        ----------
        info : dict[int, MeshInfo  |  MeshInfoCyl]
            contains the info related to the various meshes

        Examples
        --------
        parse a meshinfo file

        >>> from f4enix.output.meshinfo import MeshInfoFile
        ... mesh_info = MeshInfoFile.from_file('meshinfo')

        """
        self.info = info

    @classmethod
    def from_file(cls, meshinfofile: os.PathLike) -> MeshInfoFile:
        """Reads the meshinfo file and returns a dictionary where the keys are
        the mesh_id and the value are instances of DataMass.

        Parameters
        ----------
        meshinfofile : os.PathLike
            path to the file to be parsed

        Returns
        -------
        MeshInfoFile
            parsed meshinfo file
        """
        info = {}
        with open(meshinfofile, "r", encoding="utf-8") as infile:
            line = infile.readline()
            while line != "":
                if "Mesh tally number:" in line:
                    mesh_id = int(line.split()[-1])
                    mesh_info = _read_individual_mesh(infile)
                    info[mesh_id] = mesh_info
                line = infile.readline()
        return MeshInfoFile(info=info)

    # def get_first_meshinfo(self) -> MeshInfo | MeshInfoCyl:
    #     """Returns the meshinfo of the mesh with the lowest card id."""
    #     keys = list(self.keys())
    #     keys.sort()
    #     return self[keys[0]]


def _read_individual_mesh(infile):
    # skip 3 lines without relevant info
    for _ in range(3):
        infile.readline()
    words = infile.readline().split()
    if words[0] == "origin":
        # the mesh is cylindrical
        mesh_info = _read_header_cyl(words, infile)
    else:
        mesh_info = _read_header_cart(words, infile)
    mesh_info.data_mass = _read_voxels(infile, mesh_info)
    return mesh_info


def _read_header_cyl(words, infile) -> MeshInfoCyl:
    coordinates = CoordinateType.CYL
    origin = [float(x) for x in words[2:5]]
    axis = [float(x) for x in words[7:10]]
    vec = [float(x) for x in words[13:16]]
    words = infile.readline().split()
    vector_i = [float(x) for x in words[2:]]
    words = infile.readline().split()
    vector_j = [float(x) for x in words[2:]]
    # this line is buggy, it says revolutions but the units are in
    # radians, the two numbers may be merged ex: 0.000000006.28318531
    line = infile.readline()
    vector_k = re.findall(r"\d\.\d{8}", line)
    vector_k = [float(x) / (2 * np.pi) for x in vector_k]
    mesh_info = MeshInfoCyl(
        coordinates=coordinates,
        data_mass=None,
        vector_i=vector_i,
        vector_j=vector_j,
        vector_k=vector_k,
        origin=origin,
        axis=axis,
        vec=vec,
    )
    return mesh_info


def _read_header_cart(words, infile) -> MeshInfo:
    coordinates = CoordinateType.CART
    vector_i = [float(x) for x in words[2:]]
    words = infile.readline().split()
    vector_j = [float(x) for x in words[2:]]
    words = infile.readline().split()
    vector_k = [float(x) for x in words[2:]]
    mesh_info = MeshInfo(
        coordinates=coordinates,
        data_mass=None,
        vector_i=vector_i,
        vector_j=vector_j,
        vector_k=vector_k,
    )
    return mesh_info


def _read_voxels(infile, mesh_info) -> DataMass:
    data = {
        INDEX_KEY_VOXEL: [],
        INDEX_KEY_CELL: [],
        INDEX_KEY_MATERIAL: [],
        COLUMN_KEY_MASS_GRAMS: [],
    }
    number_of_voxels = (
        (len(mesh_info.vector_i) - 1)
        * (len(mesh_info.vector_j) - 1)
        * (len(mesh_info.vector_k) - 1)
    )
    for _i in range(number_of_voxels):
        voxel_id, volume, number_of_cells_inside = infile.readline().split()
        voxel_id = int(voxel_id)
        volume = float(volume)
        number_of_cells_inside = int(number_of_cells_inside)
        for _cell in range(number_of_cells_inside):
            cell_id, material_id, cell_mass = _read_cell_part(infile, volume)
            data[INDEX_KEY_VOXEL].append(voxel_id)
            data[INDEX_KEY_CELL].append(cell_id)
            data[INDEX_KEY_MATERIAL].append(material_id)
            data[COLUMN_KEY_MASS_GRAMS].append(cell_mass)
    if mesh_info.coordinates == CoordinateType.CYL:
        # WARNING: we divide the volume by 2 pi due to a bug in the D1S
        # generation of the meshinfo file for cyl coordinates.
        # It provides a volume calculated as if the theta coordinate units
        # were revolutions, but it is actually in radians.
        data[COLUMN_KEY_MASS_GRAMS] = np.array(data[COLUMN_KEY_MASS_GRAMS]) / (
            2 * np.pi
        )
    # the dataframe format has to be that of DataMass
    mass_dataframe = pd.DataFrame(data).set_index(
        [
            INDEX_KEY_VOXEL,
            INDEX_KEY_MATERIAL,
            INDEX_KEY_CELL,
        ]
    )
    data_mass = DataMass(mass_dataframe)
    return data_mass


def _read_cell_part(infile, volume):
    cell_id, density, material, volume_proportion = infile.readline().split()
    cell_id = int(cell_id)
    density = float(density)
    material = int(material)
    volume_proportion = float(volume_proportion)
    cell_volume_inside_voxel = volume * volume_proportion
    cell_mass = cell_volume_inside_voxel * density
    return cell_id, material, cell_mass
