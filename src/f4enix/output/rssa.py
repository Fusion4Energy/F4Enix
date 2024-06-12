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

import numpy as np
import os
from typing import BinaryIO, Dict
import matplotlib.pyplot as plt
import matplotlib.colors as colors

BYTE = np.byte
CHAR = np.char
INT = np.int32
FLOAT = np.float64
LONG = np.int64


class RSSA:
    def __init__(self, file: os.PathLike) -> None:
        """
        Representation of a RSSA file.

        Parameters
        ----------
        file : os.PathLike
            path to the rssa file to be parsed

        Attributes
        ----------
        self.parameters: dict
            {np1,  # Number of histories of the simulation, given as a negative number
            nrss,  # Number of tracks recorded
            nrcd,  # Number of values recorded for each particle, it should be 11
            njsw,  # Number of surfaces in JASW
            niss,  # Number of different histories that reached the SSW surfaces
            niwr,  # Number of cells in RSSA file
            mipts,  # Source particle type
            kjaq,  # Flag for macrobodies surfaces
            surfaces}
        self.tracks: np.ndarray
            each index of the array has 11 values
            0 a,  # History number of the particle, negative if uncollided
            1 b,  # Packed variable, the sign is the sign of the third direction cosine, starts with 8 = neutron, 16 = photon
            2 wgt,
            3 erg,
            4 tme,
            5 x,
            6 y,
            7 z,
            8 u,  # Particle direction cosine with X-axis
            9 v,  # Particle direction cosine with Y-axis, to calculate w (Z-axis) use the sign from b
            10 c  # Surface id
        filename : str
            name of the original file

        Examples
        --------
        Load a RSSA file and get some info

        >>> from f4enix.output.rssa import RSSA
        ... my_rssa = RSSA('small_cyl.w')
        ... print(my_rssa)
        RSSA file small_cyl was recorded using the following surfaces:
        Surface id: 1
        The surface type is a cylinder with a radius of 30.00
        The total amount of tracks recorded is 72083, of which 72083 were neutrons and 0 were photons.
        The simulation that produced this RSSA run 100000 histories
        The amount of independent histories that reached the RSSA surfaces was 70797.

        """
        self.filename = os.path.basename(file).split(".")[0]
        self.parameters, self.tracks = self._read_rssa(file)

    def __repr__(self) -> str:
        return self._get_info()

    def __str__(self) -> str:
        return self._get_info()

    @staticmethod
    def _read_rssa(filename: str) -> tuple[dict, np.array]:
        with open(filename, "rb") as infile:
            # This parameters hold information like the amount of histories or
            # the amount of tracks recorded
            parameters = _read_header(infile)
            tracks = _read_tracks(infile)

        return parameters, tracks

    def _calculate_grid_axes_cyl(
        self, z_int: int, theta_int: int, mask: np.ndarray = None
    ):
        """The axes are calculated without taking into account the type of
        particle"""
        if mask is None:
            mask = np.where(self.x)[0]

        z_min = self.z[mask].min()
        z_max = self.z[mask].max()

        thetas = np.angle(self.x[mask] + 1j * self.y[mask])  # in radians
        theta_min = np.min(thetas)
        theta_max = np.max(thetas)

        z_axis = np.linspace(z_min, z_max, z_int + 1)
        theta_axis = np.linspace(theta_min, theta_max, theta_int + 1)
        return z_axis, theta_axis

    def _generate_figures_current_cyl(
        self,
        particle: str = "n",
        z_int: int = 10,
        theta_int: int = 10,
        source_intensity: float = 1.7757e20,
        value_range: tuple[float, float] = None,
        mask: np.ndarray = None,
    ) -> tuple[plt.Figure, plt.Figure]:

        particle_mask = self._get_particle_mask(particle)
        if mask is not None:  # Apply to the mask a geom filter done earlier
            particle_mask = np.intersect1d(particle_mask, mask)

        # Create the empty grid
        z_axis, theta_axis = self._calculate_grid_axes_cyl(
            z_int, theta_int, particle_mask
        )
        grid_values = np.zeros((z_int, theta_int))
        error_grid = np.zeros((z_int, theta_int))

        # Calculate the indices at both axis of the grids for each track
        z_idx = _closest(z_axis, self.z[particle_mask])
        # in radians
        thetas = np.angle(self.x[particle_mask] + 1j * self.y[particle_mask])
        theta_idx = _closest(theta_axis, thetas)
        del thetas  # To relax a bit memory constraints

        # Populate both grids
        np.add.at(grid_values, (z_idx, theta_idx), self.wgt[particle_mask])
        np.add.at(error_grid, (z_idx, theta_idx), 1)

        # Normalize values
        # radius of the cylinder
        radius = np.linalg.norm([self.x[0], self.y[0]])
        extent = [
            radius * theta_axis[0],
            radius * theta_axis[-1],
            z_axis[0],
            z_axis[-1],
        ]  # used in the plotting
        area = abs(radius * (theta_axis[1] - theta_axis[0]) * (z_axis[1] - z_axis[0]))
        grid_values /= area
        grid_values /= abs(self.parameters["np1"])  # grid /= nps
        grid_values *= source_intensity
        # Give a relative error of 1 to empty voxels
        error_grid[np.where(error_grid == 0)] = 1
        error_grid = 1 / np.sqrt(error_grid)

        # Generate the figure values
        figure_values: plt.Figure = plt.figure()
        ax_values: plt.Axes = figure_values.add_subplot()
        ax_values.set_xlabel("Perimeter of the cylinder (cm)")
        ax_values.set_ylabel("Z (cm)")
        if particle == "n":
            ax_values.set_title("Neutron current through the surface [#/cm2/s]")
        else:
            ax_values.set_title("Photon current through the surface [#/cm2/s]")
        ax_values.imshow(grid_values, origin="lower", extent=extent)
        # Set the colors to log range
        if value_range is not None:
            log_max = int(np.log10(value_range[1]))
            log_min = int(np.log10(value_range[0]))
        else:
            # The +1 is needed so max=1234 => 10.000
            log_max = int(np.log10(grid_values.max())) + 1
            log_min = log_max - 10
        decades = log_max - log_min
        image_values = ax_values.images[0]
        # set the colormap and number of bins
        image_values.cmap = plt.get_cmap("jet", decades)
        # set the scale as log
        image_values.norm = colors.LogNorm(
            np.power(10.0, log_min), np.power(10.0, log_max)
        )
        _color_bar = figure_values.colorbar(image_values, orientation="horizontal")

        # Generate the errors figure
        figure_errors: plt.Figure = plt.figure()
        ax_errors: plt.Axes = figure_errors.add_subplot()
        ax_errors.set_xlabel("Perimeter of the cylinder (cm)")
        ax_errors.set_ylabel("Z (cm)")
        ax_errors.set_title("Relative error as 1/sqrt(N)")
        norm = colors.Normalize(0, 1)
        color_map = plt.get_cmap("jet", 10)
        ax_errors.imshow(
            error_grid, cmap=color_map, norm=norm, origin="lower", extent=extent
        )
        image_errors = ax_errors.images[0]
        figure_errors.colorbar(image_errors, orientation="horizontal")

        # Print information about the grid
        print(f"The area of a cell is {area:.2f}cm2")
        print(
            f"The resolution is {radius * (theta_axis[1] - theta_axis[0]):.2f}cm x {z_axis[1] - z_axis[0]:.2f}cm"
        )

        return figure_values, figure_errors

    def plot_cyl(
        self,
        particle: str = "n",
        z_int: int = 10,
        theta_int: int = 10,
        norm: float = 1,
        value_range: tuple[float, float] = None,
        x_range: tuple[float, float] = None,
        z_range: tuple[float, float] = None,
        outfolder: os.PathLike = None,
        filename: str = None,
    ) -> tuple[plt.Figure, plt.Figure]:
        """Plot the cylinder surface of RSSA.

        Parameters
        ----------
        particle : str, optional
            either neutron 'n' or photon 'p', by default 'n'
        z_int : int, optional
            z grid discretizations, by default 10
        theta_int : int, optional
            theta grid discretizations, by default 10
        norm : float, optional
            normalization to be applied to the results, by default 1
        value_range : tuple[float, float], optional
            min and max values to plot, by default None
        x_range : tuple[float, float], optional
            min x and min z to plot, by default None
        z_range : tuple[float, float], optional
            min z and max z to plot, by default None
        outfolder : os.PathLike, optional
            path to the output folder, default is None causing the plot not
            to be saved
        filename : str, optional
            ovverrides default name for the plot, by default None

        Returns
        -------
        tuple[plt.Figure, plt.Figure]
            Values and Error figures
        """
        # Filter tracks location if necessary
        mask = None
        if z_range is not None:
            mask1 = np.where(self.z > z_range[0])
            mask2 = np.where(self.z < z_range[1])
            mask = np.intersect1d(mask1, mask2)
            del mask1, mask2
        if x_range is not None:
            thetas = np.angle(self.x + 1j * self.y)  # in radians
            # radius of the cylinder
            radius = np.linalg.norm([self.x[0], self.y[0]])
            # Perimeter of the cylinder values, x values in the plot x-axis
            x_values = radius * thetas
            mask1 = np.where(x_values > x_range[0])
            mask2 = np.where(x_values < x_range[1])
            mask = np.intersect1d(mask, mask1)
            mask = np.intersect1d(mask, mask2)
            del mask1, mask2, thetas, radius, x_values

        figure_values, figure_errors = self._generate_figures_current_cyl(
            particle=particle,
            z_int=z_int,
            theta_int=theta_int,
            source_intensity=norm,
            value_range=value_range,
            mask=mask,
        )

        if outfolder is not None:
            if filename is None:
                filename = "{}_{}_z{}_theta{}".format(
                    self.filename, particle, z_int, theta_int
                )

            outval = os.path.join(outfolder, filename + ".jpeg")
            outerr = os.path.join(outfolder, filename + "_errors.jpeg")

            figure_values.savefig(outval, format="jpeg", dpi=1200)
            figure_errors.savefig(outerr + "_errors.jpeg", format="jpeg", dpi=1200)
        return figure_values, figure_errors

    def _get_info(self) -> str:
        info = f"RSSA file {self.filename} was recorded using the following surfaces:\n"
        for surface in self.parameters["surfaces"]:
            info += f'  Surface id: {surface["id"]}\n'

        sur_type = self.type
        if sur_type == "cyl":
            info += f"The surface type is a cylinder with a radius of {np.linalg.norm([self.x[0], self.y[0]]):.2f}\n"
        elif sur_type == "plane":
            info += f"The surface type is a plane...\n"

        n_tracks = self.tracks[self.mask_neutron_tracks].shape[0]
        p_tracks = self.tracks[self.mask_photon_tracks].shape[0]
        info += (
            f'The total amount of tracks recorded is {self.parameters["nrss"]}, of which {n_tracks} were neutrons'
            f" and {p_tracks} were photons.\n"
        )

        info += (
            f'The simulation that produced this RSSA run {np.abs(self.parameters["np1"])} histories\n'
            f'The amount of independent histories that reached the RSSA surfaces was {self.parameters["niss"]}.\n'
        )
        return info

    @property
    def mask_neutron_tracks(self) -> np.ndarray:
        """
        The neutron tracks
        """
        # Get all the bitarrays and don't pay attention to the sign
        bitarrays = np.abs(self.tracks[:, 1])
        # Neutrons start with 8 and photons with 16 followed by 1e8
        n_tracks = np.where(bitarrays < 9e8)[0]
        return n_tracks

    @property
    def mask_photon_tracks(self) -> np.array:
        """
        The photon tracks
        """
        # Get all the bitarrays and don't pay attention to the sign
        bitarrays = np.abs(self.tracks[:, 1])
        # Neutrons start with 8 and photons with 16 followed by 1e8
        p_tracks = np.where(bitarrays >= 9e8)[0]
        return p_tracks

    def _get_particle_mask(self, particle: str):
        if particle == "n":
            return self.mask_neutron_tracks
        elif particle == "p":
            return self.mask_photon_tracks
        else:
            raise ValueError(f"Particle was {particle}, not n or p...")

    @property
    def x(self) -> np.ndarray:
        """
        x positions
        """
        return self.tracks[:, 5]

    @property
    def y(self) -> np.ndarray:
        """
        y positions
        """
        return self.tracks[:, 6]

    @property
    def z(self) -> np.ndarray:
        """
        z positions
        """
        return self.tracks[:, 7]

    @property
    def wgt(self) -> np.ndarray:
        """
        weights
        """
        return self.tracks[:, 2]

    @property
    def energies(self) -> np.ndarray:
        """
        energies
        """
        return self.tracks[:, 3]

    @property
    def histories(self) -> np.ndarray:
        """
        histories
        """
        return np.abs(self.tracks[:, 0])

    @property
    def type(self):
        """
        RSSA surface type
        """
        # If there are more than 1 surf we cannot say if it is a cyl or a plane
        if len(self.parameters["surfaces"]) > 1:
            return "multiple"
        # Assume it is a cyl and calculate the radius of its tracks
        # if they are different it is a plane
        radius = np.sqrt(self.x**2 + self.y**2)
        if np.std(radius) < 1e-4:
            return "cyl"
        else:
            return "plane"


# --- Helper Functions ---
def _read_fortran_record(file: BinaryIO):
    count_1 = np.fromfile(file, INT, 1)[0]
    data = np.fromfile(file, np.byte, count_1)
    count_2 = np.fromfile(file, INT, 1)[0]
    if count_1 != count_2:
        raise ValueError(
            "The integers that go before and after the Fortran record are not equal..."
        )
    return data


def _read_header(file: BinaryIO) -> Dict:
    # First record
    data = _read_fortran_record(file)
    # The first line of the file with information like the code version, date and title
    format_record_id = data.tobytes().decode("UTF-8")
    if "d1suned" in format_record_id:
        # TODO: we could parse and store information like datetime and title
        _last_dump = np.frombuffer(data[-4:], INT)
        pass
    else:
        raise NotImplementedError(
            f"The code that generated this RSSA file has not been implemented"
            f" in this parser, see the code here: {format_record_id}..."
        )

    # Second record
    data = _read_fortran_record(file)
    np1 = np.frombuffer(data, LONG, 1, 0)[0]
    nrss = np.frombuffer(data, LONG, 1, 8)[0]
    nrcd = np.frombuffer(data, INT, 1, 16)[0]
    njsw = np.frombuffer(data, INT, 1, 20)[0]
    niss = np.frombuffer(data, LONG, 1, 24)[0]
    if nrcd != 11:
        raise NotImplementedError(
            f"The amount of values recorded for each particle should be 11 instead of {nrcd}..."
        )

    # Third record
    if np1 < 0:
        data = _read_fortran_record(file)
        niwr, mipts, kjaq = np.frombuffer(data, INT, 3)
    else:
        raise NotImplementedError(
            "The np1 value is not negative, as far as we understand it should be negative..."
        )

    # Fourth record
    surfaces = []
    for i in range(njsw):
        data = _read_fortran_record(file)
        surface = {"id": np.frombuffer(data, INT, 1, 0)[0]}
        if kjaq == 1:
            surface["info"] = np.frombuffer(data, INT, 1, 4)[0]
        else:
            surface["info"] = -1
        surface["type"] = np.frombuffer(data, INT, 1, 8)[0]
        surface["num_params"] = np.frombuffer(data, INT, 1, 12)[0]
        surface["params"] = np.frombuffer(data, INT, offset=16)
        surfaces.append(surface)

    # we read any extra records as determined by njsw+niwr...
    # no known case of their actual utility is known currently
    for j in range(njsw, njsw + niwr):
        _data = _read_fortran_record(file)
        raise NotImplementedError(
            "njsw + niwr values are bigger than njsw, behavior not explained"
        )

    # Summary record
    _data = _read_fortran_record(file)
    # Summary record not processed, its information does not interest us for now

    parameters = {
        "np1": np1,  # Number of histories of the simulation, given as a negative number
        "nrss": nrss,  # Number of tracks recorded
        "nrcd": nrcd,  # Number of values recorded for each particle, it should be 11
        "njsw": njsw,  # Number of surfaces in JASW
        "niss": niss,  # Number of different histories that reached the SSW surfaces
        "niwr": niwr,  # Number of cells in RSSA file
        "mipts": mipts,  # Source particle type
        "kjaq": kjaq,  # Flag for macrobodies surfaces
        "surfaces": surfaces,
    }
    return parameters


def _read_tracks(file: BinaryIO) -> np.ndarray:
    # Particle records
    # Read the whole remaining of the file at once, store all the bytes as a 1D numpy array
    data = np.fromfile(file, BYTE)
    # Reshape the array so each index holds the information of a single particle, we can do this because we
    #  know that the particle records have always the same length, 96 bytes
    data = data.reshape(-1, 96)
    # Remove the first and last 4 bytes, these are two integers that tell the record is 88 bytes long
    data = data[:, 4:-4]
    # Convert the array into a 1D array of float numbers instead of simply bytes
    data = np.frombuffer(data.flatten(), FLOAT)
    # Reshape the array so each index holds the information of a single particle, all the data is already converted
    #  from bytes to floats
    data = data.reshape(-1, 11)
    return data


def _closest(axis, values) -> np.array:
    """
    Given an axis (list of ordered equidistant values) and a value list,
    returns the indexes of the axis that are
     closest to the values.
    E.g. axis = [1, 2, 3, 4] value=[2.23, 4.98] result=[1, 2]
    """
    l0 = axis[0]
    lm = axis[1] - axis[0]
    idx = np.array((values - l0) / float(lm))
    idx = idx.astype(int)
    # Limit case where the value lies at the end of the axis
    idx[idx == len(axis) - 1] = len(axis) - 2

    return idx
