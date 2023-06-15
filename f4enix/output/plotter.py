import os
import pyvista as pv
import numpy as np
import logging
import docx
import win32com.client

from docx.shared import Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.oxml.ns import nsdecls
from docx.oxml import parse_xml

from pathlib import Path

from math import radians, degrees

pv.set_plot_theme('document')


class MeshPlotter:
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

            mesh_slice = self.mesh.slice(origin=center, normal=normal)

            # try to access its normal, if it has it, it is not empty
            try:
                mesh_slice.cell_normals[0]
            except KeyError:
                # it means that nothing can be sliced here because the
                # the slice is empty
                logging.warning(
                    'No slice can be done at theta={} deg'.format(degrees(theta)))
                continue

            mesh_slices.append(mesh_slice)

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
                    stl_color: str = 'white',
                    n_colors: int = 256) -> None:
        """Plot a series of slices to an outpath folder. The slices names
        are used as file names.

        Parameters
        ----------
        slices : list[tuple[str, pv.PolyData, pv.PolyData  |  None]]
            list of slices to be plotted. Usually produced with MeshPlotter
            methods.
        array_name : str
            name of the scalar array to be plotted
        outpath : os.PathLike
            path to the directory that will contain all the plots. This folder
            must be already created. Files with the same name will be
            overridden.
        min_max : tuple[float], optional
            min and max values to be set for the scalars in all plots
        log_scale : bool, optional
            if true, activate logarithmic scale, by default True
        stl_color : str, optional
            color for the stl, by default 'white'
        n_colors : int, optional
            number of discrete colors to be used in the legend, by default 256

        Raises
        ------
        ValueError
            if the path do not exists
        """

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
                pl.add_mesh(stl_slice, color=stl_color)

            if i == 0:
                # ensure that all pictures will have the same bounds
                bounds = mesh_slice.bounds

            self._set_perpendicular_camera(mesh_slice, pl, bounds=bounds)
            filename = os.path.join(outpath, '{}.png'.format(name))
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


class Atlas:

    def __init__(self, name: str = 'atlas',
                 landscape: str = True) -> None:
        """Class to handle the generation of the Atlas.

        Parameters
        ----------
        name : str, optional
            atlas name, by default 'atlas'
        landscape : str, optional
            if true the atlas will be produced in landscape orientation,
            by default True
        """

        self.name = name
        doc = docx.Document()
        doc.add_heading(name, level=0)
        self.doc = doc  # Word Document
        if landscape:
            self._change_orientation()

    def _insert_img(self, img: os.PathLike, width=Inches(7.5)) -> None:
        self.doc.add_picture(img, width=width)
        last_paragraph = self.doc.paragraphs[-1]
        last_paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER

    def build_from_root(self, root_path: os.PathLike) -> None:
        """Given a root folder, build an atlas with all the pictures found
        in the subfolders. Only one level of subfolders is supported at the
        moment.

        Parameters
        ----------
        root_path : os.PathLike
            path to the root folder.
        """

        # TODO only one level of folders admitted for the moment
        for folder in os.listdir(root_path):
            folderpath = os.path.join(root_path, folder)

            if not os.path.isdir(folderpath):
                continue

            self.doc.add_heading(folder, level=1)

            # This ensures correct order of the pictures that will be the same
            # of creation
            files = os.listdir(folderpath)
            paths = []
            for file in files:
                paths.append(os.path.join(folderpath, file))

            for filepath in sorted(paths, key=os.path.getmtime):
                heading = os.path.basename(filepath)[:-4]
                ext = os.path.basename(filepath)[-3:]
                if ext not in ['jpg', 'png', 'tif']:
                    continue

                self.doc.add_heading(heading, level=2)
                self._insert_img(filepath)

    def _change_orientation(self) -> docx.section.Section:
        current_section = self.doc.sections[-1]
        new_width, new_height = (current_section.page_height,
                                 current_section.page_width)
        new_section = self.doc.add_section()
        new_section.orientation = WD_ORIENT.LANDSCAPE
        new_section.page_width = new_width
        new_section.page_height = new_height

        return new_section

    def save(self, outpath: os.PathLike, pdfprint: bool = True) -> None:
        """
        Save word atlas and possibly export PDF

        Parameters
        ----------
        outpath : os.PathLike
            path to the folder where to save the atlas(es)
        pdfprint : Boolean, optional
            If True export also in PDF format

        Returns
        -------
        None.

        """
        outpath_word = os.path.join(outpath, self.name+'.docx')
        outpath_pdf = os.path.join(outpath, self.name+'.pdf')

        try:
            self.doc.save(outpath_word)
        except FileNotFoundError as e:
            print(' The following is the original exception:')
            print(e)
            print('\n it may be due to invalid characters in the file name')

        if pdfprint:
            in_file = outpath_word
            out_file = outpath_pdf

            try:
                word = win32com.client.Dispatch('Word.Application')
            except:
                raise NotImplementedError('Word not installed')

            try:
                doc = word.Documents.Open(in_file)
                doc.ExportAsFixedFormat(
                    OutputFileName=out_file,
                    # 17 = PDF output, 18=XPS output
                    ExportFormat=17,
                    OpenAfterExport=False,
                    # 0=Print (higher res), 1=Screen (lower res)
                    OptimizeFor=0,
                    # 0=No bookmarks,
                    # 1=Heading bookmarks only,
                    # 2=bookmarks match word bookmarks
                    CreateBookmarks=1,
                    DocStructureTags=True)
            finally:
                doc.Close()
                word.Quit()

    # @staticmethod
    # def _wrapper(paragraph, ptype):
    #     """
    #     Wrap a paragraph in order to add cross reference

    #     Parameters
    #     ----------
    #     paragraph : docx.Paragraph
    #         image to wrap.
    #     ptype : str
    #         type of paragraph to wrap

    #     Returns
    #     -------
    #     None.

    #     """
    #     if ptype == 'table':
    #         instruction = ' SEQ Table \\* ARABIC'
    #     elif ptype == 'figure':
    #         instruction = ' SEQ Figure \\* ARABIC'
    #     else:
    #         raise ValueError(ptype+' is not a supported paragraph type')

    #     run = run = paragraph.add_run()
    #     r = run._r
    #     fldChar = OxmlElement('w:fldChar')
    #     fldChar.set(qn('w:fldCharType'), 'begin')
    #     r.append(fldChar)
    #     instrText = OxmlElement('w:instrText')
    #     instrText.text = instruction
    #     r.append(instrText)
    #     fldChar = OxmlElement('w:fldChar')
    #     fldChar.set(qn('w:fldCharType'), 'end')
    #     r.append(fldChar)

    # @staticmethod
    # def _highlightCell(cell, color='FBD4B4'):
    #     shading_elm_1 = parse_xml(r'<w:shd {} w:fill="'.format(nsdecls('w')) +
    #                               color + r'"/>')
    #     cell._tc.get_or_add_tcPr().append(shading_elm_1)