from __future__ import annotations
import pyvista as pv
import numpy as np
import csv
import os
import logging

from copy import deepcopy
from tqdm import tqdm


ALLOWED_OUTPUT_FORMATS = ['point_cloud', 'ip_fluent', 'csv', 'vtk']


class PyVistaWrapper:
    def __init__(self, fn: os.PathLike, mesh: pv.DataSet) -> None:
        """This is a general wrapper for a pyvista object

        Parameters
        ----------
        fn : os.PathLike
            path to a .vtk file that will be read with pyvista and
            wrapped. It can also be interpreted as the name of the mesh.

        mesh : pv.DataSet
            pyvista representation of a grid.

        """
        self.filename = fn
        self.mesh = mesh
        # call this function to rewrite the mesh attributes after a change
        self._read_mesh_info()

    @classmethod
    def from_file(cls, fn: os.PathLike) -> PyVistaWrapper:
        """Generate the object from a .vtk file

        Parameters
        ----------
        fn : os.PathLike
            path to the .vtk file to be read

        Returns
        -------
        PyVistaWrapper
            initialized object
        """
        filename = os.path.basename(fn).split('.')[0]
        mesh = pv.read(fn)
        logging.info('{} read correctly'.format(fn))

        return cls(filename, mesh)

    def _read_mesh_info(self) -> None:
        self.centers = self.mesh.cell_centers().points
        self.points = self.mesh.points
        # List of arrays associated to cells
        self.cells_info = list(self.mesh.cell_data)
        # List of arrays associated to points
        self.points_info = list(self.mesh.point_data)
        self.n_coordinates = self.points.shape[1]  # Number of dimensions
        self.mesh_type = str(type(self.mesh)).split(".")[-1][:-2]
        logging.debug('info has been re-read on {}'.format(self.filename))

    def get_array_type(self, array_name: str) -> str:
        """get if an array is associated with either a point or a cells

        Parameters
        ----------
        array_name : str
            name of the array to be searched

        Returns
        -------
        str
            either 'cells' ors 'points'

        Raises
        ------
        KeyError
            raise a KeyError if the name of the array is not found
        """
        if array_name in self.cells_info:
            return "cells"
        elif array_name in self.points_info:
            return "points"
        else:
            raise KeyError('No array named {}'.format(array_name))

    def __str__(self) -> str:
        dimensions = (
            self.mesh.bounds[1] - self.mesh.bounds[0],
            self.mesh.bounds[3] - self.mesh.bounds[2],
            self.mesh.bounds[5] - self.mesh.bounds[4],
        )
        file = f"""
    Name: {self.filename}
        """
        formatted_mesh_bounds = ""
        for x in self.mesh.bounds:
            formatted_mesh_bounds += f" {x:.2f}"
        formatted_dimensions = ""
        for x in dimensions:
            formatted_dimensions += f" {x:.2f}"

        mesh = f"""
    Number of cells: {self.mesh.n_cells}
    Number of points: {self.mesh.n_points}
    Cells arrays: {str(self.cells_info)}
    Points arrays: {str(self.points_info)}
    Number of coordinates: {self.n_coordinates}
    Mesh bounds: {formatted_mesh_bounds}
    Mesh dimensions: {formatted_dimensions}
    Mesh type: {self.mesh_type}
        """
        return (file + "\n" + mesh)

    def __repr__(self) -> str:
        return self.__str__()

    def print_array_info(self, array_name: str) -> None:
        """print formatted information on the requested array.
        min and max val, integrals, averages. The output is differentiated if
        the requested array is a cell or point field.

        Parameters
        ----------
        array_name : str
            name of the array to be searched
        """

        if self.get_array_type(array_name) == "cells":  # For cells
            min_value = min(self.mesh[array_name])
            max_value = max(self.mesh[array_name])
            (
                integral_no_volume,
                average_no_volume,
                integral_volume,
                average_volume,
            ) = self._integral_and_average(array_name)
            print(
                f"""
                Minimum value: {min_value:.2e}
                Maximum value: {max_value:.2e}
                Integral without volume: {integral_no_volume:.2e}
                Integral with volume: {integral_volume:.2e}
                Average without volume: {average_no_volume:.2e}
                Average with volume: {average_volume:.2e}
                """
            )
        elif self.get_array_type(array_name) == "points":  # For points
            min_value = min(self.mesh[array_name])
            max_value = max(self.mesh[array_name])
            integral, average = self._integral_and_average(array_name)
            print(
                f"""
                Minimum value: {min_value:.2e}
                Maximum value: {max_value:.2e}
                Integral: {integral:.2e}
                Average: {average:.2e}
                """
            )
        else:
            print("This array doesn't belong to neither cells nor points")

    def _integral_and_average(self, array_name: str) -> tuple[float]:
        # for cells, the volume value is used to calculate the weight of each
        # value of the array
        if self.get_array_type(array_name) == "cells":
            integral_no_volume = sum(self.mesh[array_name])
            average_no_volume = integral_no_volume / self.mesh.n_cells
            cells_volume = abs(self.mesh.compute_cell_sizes()["Volume"])
            values = np.multiply(self.mesh[array_name], cells_volume)
            integral_volume = float(sum(values))
            average_volume = integral_volume / sum(cells_volume)

            return (integral_no_volume, average_no_volume, integral_volume,
                    average_volume)

        else:  # In meshtal.points_info
            integral = sum(self.mesh[array_name])
            average = integral / self.mesh.n_points

            return integral, average

    def translate(self, x: float = 0, y: float = 0, z: float = 0) -> None:
        """Translate the mesh in space. This action is performed in place.

        Parameters
        ----------
        x : float, optional
            Dx to be applied, by default 0
        y : float, optional
            Dy to be applied, by default 0
        z : float, optional
            Dz to be applied, by default 0
        """
        # convert the mesh to a suitable type if needed
        self._check_grid_type()

        self.mesh.translate((x, y, z), inplace=True)
        self._read_mesh_info()
        logging.info('Translation applied successfully.')

    def rotate(self, theta_x: float = 0, theta_y: float = 0,
               theta_z: float = 0) -> None:
        """Rotate a mesh in space. This action is performed in place.
        If more than one rotation is going to be applied, assume the the order
        of rotation applications will be x -> y -> z.

        Parameters
        ----------
        theta_x : float, optional
            _description_, by default 0
        theta_y : float, optional
            _description_, by default 0
        theta_z : float, optional
            _description_, by default 0
        """
        # convert the mesh to a suitable type if needed
        self._check_grid_type()

        self.mesh.rotate_x(theta_x, inplace=True)
        self.mesh.rotate_y(theta_y, inplace=True)
        self.mesh.rotate_z(theta_z, inplace=True)
        self._read_mesh_info()

        logging.info(f"Rotation applied successfully.")

    def scale(self, xyz:  float | int | list | tuple) -> None:
        """scale the grid geometry by a factor.
        Different factors can also be provided on the three dimensions.

        Parameters
        ----------
        xyz : float | int | list | tuple
            _description_
        """
        self.mesh.scale(xyz, inplace=True)
        self._read_mesh_info()

    def scalar_multiply(self, factor: float, array: str | list[str]) -> None:
        """scale a field or a list of field for a scalar factor. This
        modification is done inplace.

        Parameters
        ----------
        factor : float
            scalar factor to be applied
        array : str | list[str]
            field(s) to be modified.
        """
        if type(array) is str:
            array = [array]

        for name in array:
            self.mesh[name] = self.mesh[name] * factor

    def linear_combination(self, field_name: str, arrays: list[str],
                           coeff: list[int]) -> None:
        """perform a linear combination of different arrays in the mesh.
        The result is added to the mesh.

        Parameters
        ----------
        field_name : str
            name under which the new array should be stored
        arrays : list[str]
            list of arrays to be combined
        coeff : list[int]
            coefficients to be used for the linear combination
        """

        data = deepcopy(self.mesh[arrays[0]])*coeff[0]
        for name, c in zip(arrays[1:], coeff[1:]):
            data = data + c*deepcopy(self.mesh[name])

        self.mesh.add_field_data(data, field_name)

    def _check_grid_type(self) -> None:
        if (self.mesh_type == "StructuredGrid" or
            self.mesh_type == "UnstructuredGrid"):
            # These are admitted types
            mesh = self.mesh
            # The translate function does not work with RectilinearGrid Meshes.
            # So, a RectilinearGrid Mesh must be first converted to
            # StructuredGrid
        elif self.mesh_type == "RectilinearGrid":
            logging.warning(
                'Rectilinear grid was converted to structured grid')
            mesh = self.mesh.cast_to_structured_grid()
        else:
            TypeError('Cannot translate type {}'.format(self.mesh_type))

        self.mesh = mesh
        self._read_mesh_info()

    # WRITE A FILE IN A CHOSEN FORMAT
    def write_mesh(self, outpath: os.PathLike,
                   list_array_names: list[str] = None,
                   out_format: str = 'vtk') -> None:
        """Export the mesh to a file. vtk, csv, fluent (txt) and point cloud
        (txt) formats can be selected.

        Parameters
        ----------
        outpath : os.PathLike
            path to the output folder.
        list_array_names : list[str], optional
            arrays to be exported. The default is None, meaning that all the 
            available arrays will be used.
        out_format : str, optional
            output format. The allowed ones are ['point_cloud', 'ip_fluent',
            'csv', 'vtk']. Default is .vtk

        Raises
        ------
        KeyError
            raises KeyError if the output format is not allowed.
        """
        if list_array_names is None:
            list_array_names = list(self.mesh.array_names)

        file_name = f"{self.filename}_{out_format}"
        filepath = os.path.join(outpath, file_name)

        if out_format == 'vtk':
            if self.mesh_type == "StructuredGrid":
                ext = '.vts'
            elif self.mesh_type == "UnstructuredGrid":
                ext = '.vtu'
            elif self.mesh_type == "RectilinearGrid":
                ext = '.vtr'
            else:
                ext = '.vtk'

            self.mesh.save(filepath+ext)
            return

        # --- CSV writer ---
        elif out_format == "csv":
            new_name = filepath + '.csv'

            with open(new_name, "w", newline="") as outfile:
                writer = csv.writer(outfile)

                # TODO This may create some issues...
                values_type = self.get_array_type(list_array_names[0])
                if values_type == "cells":  # Take points or centers
                    f_points = self.centers
                else:  # Points
                    f_points = self.points

                for i in tqdm(range(len(f_points)),
                              unit=" Points", desc="Writing"):
                    csv_points = [
                        f"{f_points[i][0]:.3f}",
                        f" {f_points[i][1]:.3f}",
                        f" {f_points[i][2]:.3f}",
                    ]
                    for array_name in list_array_names:
                        csv_points.append(
                            f" {self.mesh[array_name][i]:.3f}")
                    writer.writerow(csv_points)

                logging.info(f"{new_name} created successfully!")
            return

        for array_name in list_array_names:
            values_type = self.get_array_type(list_array_names[0])
            if values_type == "cells":  # Take points or centers
                f_points = self.centers
            else:  # Points
                f_points = self.points
            values = self.mesh[array_name]

            # write depending on format
            # --- point cloud writer ---
            if out_format == "point_cloud":
                new_name = filepath + '.txt'
                with open(new_name, "w") as outfile:
                    outfile.write("x, y, z, value\n")
                    # TODO this can probably be optmized using
                    # outfile.writeline()
                    for i in tqdm(range(len(f_points)),
                                  unit=" Points", desc="Writing"):
                        outfile.write(f"{f_points[i][0]:.3f},")
                        outfile.write(f"{f_points[i][1]:.3f},")
                        outfile.write(f"{f_points[i][2]:.3f},")
                        outfile.write(f"{values[i]:.3f}\n")
                logging.info(f"{new_name} created successfully!")
                return

            # --- fluent writer ---
            elif out_format == "ip_fluent":
                new_name = filepath + '.txt'
                with open(new_name, "w") as outfile:
                    guion1 = "3"
                    n_coord = self.n_coordinates
                    n_values = str(len(f_points))
                    guion2 = "1"
                    uds = "uds-0"
                    beginning = f"{guion1}\n{n_coord}\n{n_values}\n{guion2}\n{uds}\n"
                    outfile.write(beginning)
                    outfile.write("(")
                    for i in tqdm(range(len(f_points)),
                                  unit=" x points", desc="Writing x"):
                        outfile.write(f"{f_points[i][0]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(range(len(f_points)),
                                  unit=" y points", desc="Writing y"):
                        outfile.write(f"{f_points[i][1]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(range(len(f_points)),
                                  unit=" z points", desc="Writing z"):
                        outfile.write(f"{f_points[i][2]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(range(len(f_points)),
                                  unit=" values", desc="Writing values"):
                        outfile.write(f"{values[i]:.3f}\n")

                    outfile.write(")\n")

                logging.info(f"{new_name} created successfully!")
                return

        raise KeyError(
            "Invalid format, these are the ones allowed: {}".format(ALLOWED_OUTPUT_FORMATS))

    def merge(self, grid) -> None:
        """Merge the current wrapper with a second one.

        Parameters
        ----------
        grid : PyVistaWrapper
            wrapper to be merged with the current one
        """
        self.mesh = self.mesh.merge(grid.mesh)
        # Regardless of the initial meshes,
        # the resulted mesh is an UnstructuredGrid
        logging.warning('Merge will cause an UnstructuredGrid to be created')
        self._read_mesh_info()
        new_name = self.filename[:-4] + "+" + grid.filename[:-4] + ".vtu"
        self.filename = new_name
        logging.info(f"Merge was successful'")
