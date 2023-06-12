import os
import pyvista as pv
import numpy as np
from math import radians


class Atlas:
    def __init__(self, mesh: pv.PolyData, stl: pv.PolyData = None
                 ) -> None:
        super().__init__()
        self.mesh = mesh

        self.stl = stl

        if stl is not None:
            self._has_stl = True
        else:
            self._has_stl = False

        # Default legend settings
        self.legend_args = dict(
            title_font_size=18,
            label_font_size=18,
            shadow=True,
            # n_labels=7,
            italic=True,
            fmt="%.e",
            font_family="arial",
            vertical=True)

    def _get_plotter(self) -> pv.Plotter:
        # Initiate the plotter with all default actions if needed
        pl = pv.Plotter()
        pl.set_background('white')

        return pl

    def slice_on_axis(self, axis: str, n: int
                      ) -> tuple[list[pv.PolyData], list[pv.PolyData] | None]:
        """Creates a series of slice of the mesh distributed equally along the
        specified axis. Stl file will be sliced accordingly if present

        Parameters
        ----------
        axis : str
            either 'x', 'y' or 'z'
        n : int
            number of divisions

        Returns
        -------
        tuple[list[pv.PolyData], list[pv.PolyData]]
            list of mesh slices and stl slices. If no stl is provided the
            second list will be equal to None.
        """
        # Use the automatic slicing from pyvista
        mesh_slices = self.mesh.slice_along_axis(n, axis=axis)
        if self._has_stl:
            stl_slices = self._get_stl_slices(mesh_slices, axis=axis)
        else:
            stl_slices = None

        return mesh_slices, stl_slices

    def slice_toroidal(self, theta_increment: float,
                       center: list = [0, 0, 0]
                       ) -> tuple[list[pv.PolyData], list[pv.PolyData] | None]:
        """Create a series of slices of the mesh distributed toroidally equally
        around the vertical axis (z). Stl file will be sliced accordingly if
        present

        Parameters
        ----------
        theta_increment : float
            theta increments to be performed for each slide
        center : list, optional
            x, y, z point to be used as center for the
            slicing, by default [0, 0, 0].

        Returns
        -------
        tuple[list[pv.PolyData], list[pv.PolyData] | None]
            _description_
        """

        center = np.array(center)
        increment = radians(theta_increment)
        mesh_slices = []

        for theta in np.arange(0, np.pi, increment):
            normal = np.array(
                [np.cos(theta), np.sin(theta), 0.0]).dot(np.pi / 2.0)
            mesh_slices.append(self.mesh.slice(origin=center, normal=normal))

        if self._has_stl:
            stl_slices = self._get_stl_slices(mesh_slices, axis='z')
        else:
            stl_slices = None

        return mesh_slices, stl_slices

    def plot_slices(self, mesh_slices: list[pv.PolyData], array_name: str,
                    norm_axis: str,
                    stl_slices: list[pv.PolyData] = None,
                    min_max: tuple[float] = None,
                    log_scale: bool = True,
                    stl_rgb: list = [1, 1, 1],
                    n_colors: int = 256) -> None:

        for i, mesh_slice in enumerate(mesh_slices):
            pl = self._get_plotter()
            pl.add_mesh(mesh_slice, scalars=array_name,
                        scalar_bar_args=self.legend_args,
                        log_scale=log_scale,
                        # below_color='grey',
                        # above_color='purple',
                        clim=min_max, cmap='jet',
                        n_colors=n_colors)
            if stl_slices is not None:
                pl.add_mesh(stl_slices[i], color=stl_rgb)

    @staticmethod
    def _get_perpendicular_camera(mesh_slice: pv.PolyData, norm_axis: str,
                                  pl: pv.Plotter):
        # Initiate the plotter perpendicular to the plane
        # get the center of the slide
        mesh_slice

    def _get_stl_slices(self, mesh_slices: list[pv.PolyData], axis
                        ) -> list[pv.PolyData] | None:
        stl_slices = []
        # get the correspondent stl_slices
        for mesh_slice in mesh_slices:
            stl_slices.append(
                self.stl.slice(normal=axis, origin=mesh_slice.center))

        return stl_slices

    def _plot_single(self, name: str, xyz: list, cosines: list, minval: float,
                     maxval: float, legend: str, log_scale: str,
                     n_colours: int) -> None:
        pass
