import os
import pyvista as pv
import numpy as np
from math import radians, degrees

pv.set_plot_theme('document')


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
        pl = pv.Plotter(off_screen=True)
        pl.set_background('white')

        return pl

    def slice_on_axis(self, axis: str, n: int
                      ) -> list[tuple[str, pv.PolyData, pv.PolyData | None]]:
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
        list[tuple[str, pv.PolyData, pv.PolyData | None]]
            contains a list of (slice name, mesh slice, stl slice).
            the names are assigned according to slicing logic.
            In case of no stl assigned to the atlas, the stl slice will be
            equal to None.

        """
        # Use the automatic slicing from pyvista
        idxs = {'x': 0, 'y': 1, 'z': 2}

        # perform the slices
        slices = self.mesh.slice_along_axis(n, axis=axis)
        # get the corresponding stl slices if needed
        if self._has_stl:
            stl_slices = self._get_stl_slices(slices)

        outp = []
        # build the output
        for i, mslice in enumerate(slices):
            # TODO maybe a smarter function here
            name = 'P{} = {}'.format(axis, round(mslice.center[idxs[axis]], 1))
            if self._has_stl:
                stl_slice = stl_slices[i]
            else:
                stl_slice = None

            outp.append((name, mslice, stl_slice))

        return outp

    def slice_toroidal(self, theta_increment: float,
                       center: list = [0, 0, 0],
                       min_max_theta: tuple[float, float] = None,
                       ) -> list[tuple[str, pv.PolyData, pv.PolyData | None]]:
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
        min_max_theta : tuple(float, float), optional
            specify a minimum and max theta for the plots (degrees).
            This may help
            avoiding slices in empty regions in partial models if the mesh
            is not properly rotated. By default is None, meaning that the 
            slicing will start at theta = 0 and finish at theta = pi.

        Returns
        -------
        list[tuple[str, pv.PolyData, pv.PolyData | None]]
            contains a list of (slice name, mesh slice, stl slice).
            the names are assigned according to slicing logic.
            In case of no stl assigned to the atlas, the stl slice will be
            equal to None.
        """

        center = np.array(center)
        increment = radians(theta_increment)
        mesh_slices = []

        if min_max_theta is None:
            angles = np.arange(0, np.pi, increment)
        else:
            angles = np.arange(radians(min_max_theta[0]),
                               radians(min_max_theta[1]),
                               increment)

        for theta in angles:
            normal = np.array(
                [np.cos(theta), np.sin(theta), 0.0]).dot(np.pi / 2.0)
            mesh_slices.append(self.mesh.slice(origin=center, normal=normal))

        if self._has_stl:
            stl_slices = self._get_stl_slices(mesh_slices)

        outp = []
        # build the output
        for i, mslice in enumerate(mesh_slices):
            name = 'theta = {} deg'.format(round(degrees(angles[i]), 1))
            if self._has_stl:
                stl_slice = stl_slices[i]
            else:
                stl_slice = None

            outp.append((name, mslice, stl_slice))

        return outp

    def plot_slices(self,
                    slices: list[tuple[str, pv.PolyData, pv.PolyData | None]],
                    array_name: str,
                    outpath: os.PathLike,
                    min_max: tuple[float] = None,
                    log_scale: bool = True,
                    stl_rgb: str = 'white',
                    n_colors: int = 256) -> None:

        # Check that outpath exists
        if not os.path.exists(outpath):
            raise ValueError('{} does not exists'.format(outpath))

        for i, (name, mesh_slice, stl_slice) in enumerate(slices):
            pl = self._get_plotter()
            pl.add_mesh(mesh_slice, scalars=array_name,
                        scalar_bar_args=self.legend_args,
                        log_scale=log_scale,
                        below_color='grey',
                        above_color='purple',
                        clim=min_max, cmap='jet',
                        n_colors=n_colors)
            if stl_slice is not None:
                pl.add_mesh(stl_slice, color=stl_rgb)

            if i == 0:
                # ensure that all pictures will have the same bounds
                bounds = mesh_slice.bounds
            self._set_perpendicular_camera(mesh_slice, pl, bounds=bounds)
            filename = os.path.join(outpath, '{}.jpg'.format(name))
            pl.screenshot(filename)

    @staticmethod
    def _set_perpendicular_camera(mesh_slice: pv.PolyData,
                                  pl: pv.Plotter,
                                  bounds: list = None) -> None:
        # align camera: focus on center, position at center + normal
        center = mesh_slice.center
        pl.camera.focal_point = center
        pl.camera.position = center + mesh_slice.cell_normals[0]
        # reset camera to put entire mesh in view
        if bounds is None:
            pl.reset_camera()
        else:
            pl.reset_camera(bounds=bounds)

    def _get_stl_slices(self, mesh_slices: list[pv.PolyData]
                        ) -> list[pv.PolyData] | None:
        stl_slices = []
        # get the correspondent stl_slices
        for mesh_slice in mesh_slices:
            norm = mesh_slice.cell_normals[0]
            stl_slice = self.stl.slice(normal=norm, origin=mesh_slice.center)
            # give a very small translation in order to not be exactly
            # coincident, using the normal should be sufficient
            stl_slice.translate(norm, inplace=True)
            stl_slices.append(stl_slice)

        return stl_slices

    @staticmethod
    def _get_perpendicular_vector(v):
        if v[1] == 0 and v[2] == 0:
            if v[0] == 0:
                raise ValueError('zero vector')
            else:
                return np.cross(v, [0, 1, 0])

        return np.cross(v, [1, 0, 0])
