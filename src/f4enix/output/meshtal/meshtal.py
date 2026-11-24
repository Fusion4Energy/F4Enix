from __future__ import annotations

import re
from copy import deepcopy
from typing import TYPE_CHECKING

import aci_tallies
import migjorn

from f4enix.core.constants import PathLike
from f4enix.input.MCNPinput import Input
from f4enix.output.meshtal.fmesh import COLUMN_LABELS, Fmesh

if TYPE_CHECKING:
    import pyvista as pv


class Meshtal:
    def __init__(
        self, filename: PathLike, filetype: str = "MCNP", inp: Input | None = None
    ):
        """Meshtal class to read and handle meshtal files.

        Parameters
        ----------
        filename : PathLike
            path to the meshtal file to be read.
        filetype : str, optional
            Output type of the file, by default "MCNP". Other options are "CDGS" and "CUV".
        inp : Input | None, optional
            MCNP input of a "MCNP" meshtal, see `from_d1s_meshtal`. Default is None.

        Attributes
        ----------
        filetype : str
            Type of the meshtal file, can be "MCNP", "CDGS", or "CUV".
        mesh : dict[int, Fmesh]
            Dictionary containing Fmesh objects indexed by their mesh ID.
        """
        self.filename = filename
        self.filetype = filetype
        self.mesh: dict[int, Fmesh] = {}
        if filetype == "MCNP":
            if inp is not None:
                inp = aci_tallies.Input.from_text(inp._model.to_source())
            self._file = aci_tallies.meshtal.open(filename, input=inp)
        elif filetype == "CUV":
            self._file = aci_tallies.cuv.open(filename)
        elif filetype == "CDGS":
            self._file = aci_tallies.cdgs.open(filename)
        else:
            raise KeyError(f"Unknown filetype {filetype!r}: use MCNP, CUV or CDGS")

    def __deepcopy__(self, memo: dict) -> Meshtal:
        # the file only reads: it cannot be copied, nor needs to
        memo[id(self._file)] = self._file
        memo[id(self)] = new = self.__class__.__new__(self.__class__)
        new.__dict__.update(deepcopy(self.__dict__, memo))
        return new

    @classmethod
    def from_d1s_meshtal(cls, filename: PathLike, inp: Input) -> Meshtal:
        """Read a D1SUNED meshtal binned by nuclide.

        D1SUNED writes the nuclide lists of the meshtal header (daughter or parent
        bins, ``etag`` or ``ttag`` = ``nuc``/``parent``) with a fixed width of 10
        characters per bin, so 8-digit ids such as metastable nuclides (e.g. Nb92m =
        41092900) overflow to '**********' and the file cannot be parsed. Each
        nuclide list is rebuilt from the ``emesh`` (or ``tmesh``)
        keyword of the corresponding FMESH card of the MCNP input (by
        aci-tallies, which also takes the exact mesh geometry from the input).
        ``eints``/``tints`` are not supported for the nuclide bins.

        Parameters
        ----------
        filename : PathLike
            path to the D1SUNED meshtal file.
        inp : Input
            MCNP input where the FMESH tallies are defined.

        Returns
        -------
        Meshtal
            meshtal object, ``readMesh()`` still has to be called.

        Examples
        --------
        >>> inp = Input.from_input("d1s_input.i")
        >>> meshtal = Meshtal.from_d1s_meshtal("meshtal", inp)
        >>> meshtal.readMesh()
        """
        return cls(filename, inp=inp)

    def create_filtered_mesh(
        self,
        mesh_id: int,
        binlabels: tuple[str],
        ebin: int | None = None,
        tbin: int | None = None,
    ) -> Fmesh:
        """Create a filtered mesh based on the specified mesh ID and bin labels. That
        is, it returns a new Fmesh object with the data filtered only according to the
        time and energy bins requested.

        Parameters
        ----------
        mesh_id : int
            ID of the mesh to filter.
        binlabels : tuple[str]
            labels for the bins to be used in the filtered mesh.
        ebin : int | None, optional
            Specific energy bin to extract, by default None
        tbin : int | None, optional
            specific time bin to extract, by default None

        Returns
        -------
        Fmesh
            Filtered Fmesh object containing only the specified bins.
        """
        return self.mesh[mesh_id].select_bins(ebin=ebin, tbin=tbin, binlabels=binlabels)

    def readMesh(
        self,
        mesh: int | list[int] | None = None,
        norm: str | None = None,
        cell_filters: list[int] | None = None,
    ) -> None:
        """Read the meshtal file and build the mesh dictionary.

        Parameters
        ----------
        mesh : int | list[int] | None, optional
            If None, all meshes are read. If an integer, only the mesh with that ID is
            read. If a list of integers, only the meshes with those IDs are read.
        norm : str | None, optional
            Normalization option for CUV meshes. Default is None. Valid options are
            'ctot' or 'celf'. See D1SUNED manual for additional details
        cell_filters : list[int] | None, optional
            List of cells with which voxels in CUV meshes are filtered. Default is None.

        """
        if not mesh:
            mesh = self.mesh_ids
        elif isinstance(mesh, int):
            mesh = [mesh]
        self.mesh = {id: self._read(id, norm, cell_filters) for id in mesh}

    def _read(self, id: int, norm: str | None, cells: list[int] | None) -> Fmesh:
        """Read a mesh, see `readMesh`"""
        if id not in self.mesh_ids:
            raise KeyError("bad tally entry")
        if self.filetype == "MCNP":
            if id in (failed := dict(self._file.failed)):
                raise ValueError(
                    f"Mesh tally {id} of {self.filename} cannot be read: {failed[id]}"
                )
            return Fmesh.from_aci(self._file[id], type="meshtal")
        if self.filetype == "CUV":
            # unknown options are ignored
            norm = norm if norm in ("ctot", "celf") else None
            tally = self._file.to_mesh(id, cells=cells, norm=norm)
            return Fmesh.from_aci(tally, type="CUV")
        info = self._file.info(id)
        return Fmesh.from_aci(
            self._file[id],
            type="CDGS",
            particle=None,
            comments=info.comment,
            cooling_time=info.cooling_time,
            strength=info.strength,
        )

    def write_all(self, outpath: PathLike, out_format: str = "vtk") -> None:
        """write all fmeshes to the outfolder in the specified format.

        Parameters
        ----------
        outpath : os.PathLike
            path to the output folder.
        out_format : str, optional
            output format. The allowed ones are ['point_cloud', 'ip_fluent',
            'csv', 'vtk']. Default is .vtk

        """
        for _, mesh in self.mesh.items():
            mesh.write(outpath, out_format=out_format)

    def collapse_grids(self, name_dict: dict[int, list[str]]) -> pv.DataSet:
        """If the all the fmeshes indicated in the dictionary are defined on
        the same
        structured grid, returns a grid onto which all the fmeshes are
        collapsed. That is, the returned grid will have all the field data that
        are stored in the different fmeshes.

        Parameters
        ----------
        name_dict : dict[int, list[str, str]]
            this dictionary is used to assign names to the fields that are
            added to the grid. The key of the dictionary is the number of the
            fmesh, while the list of strings are
            ['name of the values', 'name of the values statistical error'].
            For instance {104: ['neutron heating [W/cc]',
            'neutron heating [Rel err]']}

        Returns
        -------
        pv.DataSet
            returns the pyvista grid where all fields from the different
            Fmeshes have been added

        Raises
        ------
        RuntimeError
            if the fmeshes have different geometry if they do not contain
            only the default fields ['Value - Total', 'Error - Total'].
        """
        fmeshes = [self.mesh[key] for key in name_dict]
        if not all(fmeshes[0].sameMesh(fmesh) for fmesh in fmeshes):
            raise RuntimeError("the different fmeshes geom are not compatible")
        if any(tuple(fmesh.grid.array_names) != COLUMN_LABELS for fmesh in fmeshes):
            raise RuntimeError("no binning allowed for the collapse")
        grid = fmeshes[0].grid.copy()
        grid.clear_cell_data()
        for fmesh, names in zip(fmeshes, name_dict.values()):
            for new_name, name in zip(names, COLUMN_LABELS):
                grid[new_name] = fmesh.grid[name].copy()
        return grid

    def transform_fmesh(self, inp: Input) -> None:
        """Applies the transformation to the fmeshes in the meshtal object.
        Given an input object, rototranslate the fmeshes in the meshtal object
        according to their tr=... flag in the FMESH card, by using the corresponding
        transformation.

        Parameters
        ----------
        inp : Input
            MCNP input object including the fmesh cards and the transformations
        """
        transf_dict = {}
        pattern = r"tr\s*=\s*\d+"
        for card in inp._model.data_cards("fmesh"):
            if card.number not in self.mesh.keys():
                continue
            match = re.search(pattern, card.text, re.IGNORECASE)
            if match:
                tr_num = int(match.group().split("=")[-1].strip())
                transf_dict[card.number] = inp.transformations[f"TR{tr_num}"]
        self.transform_multiple_fmesh(transf_dict)

    def transform_multiple_fmesh(
        self, transf_dict: dict[int, migjorn.Transform]
    ) -> None:
        """Transform multiple fmeshes using a dict of migjorn Transform objects.

        Parameters
        ----------
        transf_dict : dict[int, migjorn.Transform]
            dictionary of fmesh numbers and transformation objects.
        """
        for key, transf in transf_dict.items():
            self.mesh[key].apply_transformation(transf)

    @property
    def mesh_ids(self) -> list[int]:
        """Returns the list of mesh IDs."""
        failed = [] if self.filetype == "CUV" else self._file.failed
        return sorted({*self._file.ids, *(id for id, _ in failed)})
