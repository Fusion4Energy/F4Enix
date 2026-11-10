from __future__ import annotations

import io
import re
from copy import deepcopy
from typing import TYPE_CHECKING

import migjorn

from f4enix.core.constants import (
    PAT_FMESH_BINS,
    PAT_FMESH_TAG,
    PAT_NUCLIDE_LIST,
    PathLike,
)
from f4enix.input.MCNPinput import Input
from f4enix.output.meshtal.aux_meshtal_functions import COLUMN_LABELS
from f4enix.output.meshtal.fmesh import Fmesh
from f4enix.output.meshtal.mesh_parser import (
    CDGSMeshParser,
    CUVMeshParser,
    MeshtalParser,
    Parser,
)

PARSER_SELECTOR: dict[str, type[Parser]] = {
    "MCNP": MeshtalParser,
    "CDGS": CDGSMeshParser,
    "CUV": CUVMeshParser,
}

if TYPE_CHECKING:
    import pyvista as pv


class Meshtal:
    def __init__(self, filename: PathLike, filetype: str = "MCNP"):
        """Meshtal class to read and handle meshtal files.

        Parameters
        ----------
        filename : PathLike
            path to the meshtal file to be read.
        filetype : str, optional
            Output type of the file, by default "MCNP". Other options are "CDGS" and "CUV".

        Attributes
        ----------
        filetype : str
            Type of the meshtal file, can be "MCNP", "CDGS", or "CUV".
        mesh : dict[int, Fmesh]
            Dictionary containing Fmesh objects indexed by their mesh ID.
        """
        self.filetype = filetype

        self._meshtal_parser: Parser = PARSER_SELECTOR[self.filetype](filename)
        self.mesh: dict[int, Fmesh] = {}

    @classmethod
    def from_d1s_meshtal(cls, filename: PathLike, inp: Input) -> "Meshtal":
        """Read a D1SUNED meshtal binned by nuclide.

        D1SUNED writes the nuclide lists of the meshtal header (daughter or parent
        bins, ``etag`` or ``ttag`` = ``nuc``/``parent``) with a fixed width of 10
        characters per bin, so 8-digit ids such as metastable nuclides (e.g. Nb92m =
        41092900) overflow to '**********' and the file cannot be parsed. Each
        nuclide list is rebuilt from the ``emesh`` (or ``tmesh``)
        keyword of the corresponding FMESH card of the MCNP input, then the
        corrected meshtal is parsed from memory (the original file is not
        modified). ``eints``/``tints`` are not supported for the nuclide bins.

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
        nuclide_lists = {
            str(card.number): _d1s_nuclide_lists(card.text)
            for card in inp._model.data_cards("fmesh")
        }

        lines = []
        tally = None
        n_list = 0
        with open(filename) as infile:
            for line in infile:
                if line.startswith(" Mesh Tally Number"):
                    tally = line[18:].strip()
                    n_list = 0
                match = PAT_NUCLIDE_LIST.match(line)
                if match:
                    bins = " ".join(nuclide_lists[tally][n_list])
                    line = f"{match.group()} {bins}\n"
                    n_list += 1
                lines.append(line)
        return cls(io.StringIO("".join(lines)))

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
        etbin_data = self.mesh[mesh_id].get_etbin_data(ebin=ebin, tbin=tbin)
        etbin_mesh = Fmesh(etbin_data, str(mesh_id), trsf=None, binlabels=binlabels)
        return etbin_mesh

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
            mesh_dict = {}
            for mesh_id in self._meshtal_parser.get_meshlist():
                mesh_dict[int(mesh_id)] = self._meshtal_parser.get_FMesh(
                    mesh_id,
                    norm,
                    cell_filters,
                )
            self.mesh = mesh_dict
        elif isinstance(mesh, int):
            self.mesh = {mesh: self._meshtal_parser.get_FMesh(mesh, norm, cell_filters)}
        elif isinstance(mesh, list):
            self.mesh = {
                m: self._meshtal_parser.get_FMesh(m, norm, cell_filters) for m in mesh
            }

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
        try:
            # check that the collapse is doable
            for i, key in enumerate(list(name_dict.keys())):
                fmesh = self.mesh[key]
                # Check they are all same size
                if i == 0:
                    try:
                        comparison = fmesh
                    except AttributeError:
                        # it means that the mesh has not been read
                        raise AttributeError(
                            "Please read the meshtal file first using readMesh()"
                        )
                else:
                    assert comparison.sameMesh(fmesh)
        except AssertionError:
            raise RuntimeError("the different fmeshes geom are not compatible")
        try:
            # check that the collapse is doable
            for i, key in enumerate(list(name_dict.keys())):
                fmesh = self.mesh[key]
                # check that there are only the two usual values, no binning
                assert tuple(fmesh.grid.array_names) == COLUMN_LABELS
        except AssertionError:
            raise RuntimeError("no binning allowed for the collapse")
        for i, key in enumerate(list(name_dict.keys())):
            fmesh = self.mesh[key]
            if i == 0:
                grid = deepcopy(fmesh.grid)
                for old_name, new_name in zip(grid.array_names, name_dict[key]):
                    grid.rename_array(old_name, new_name)
            else:
                # results just have to be added for the other fmeshes
                for array_name, id in zip(name_dict[key], COLUMN_LABELS):
                    grid[array_name] = fmesh.grid[id]
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
        return sorted(int(x) for x in self._meshtal_parser.get_meshlist())


def _d1s_nuclide_lists(card_text: str) -> list[list[str]]:
    """Nuclide lists of a D1SUNED FMESH card, in the order of the meshtal header."""
    text = "\n".join(line.split("$")[0] for line in card_text.splitlines())
    text = text.replace("&", " ")
    tags = {k.lower(): v.lower() for k, v in PAT_FMESH_TAG.findall(text)}
    values = {k.lower(): v.split() for k, v in PAT_FMESH_BINS.findall(text)}
    return [
        ["0.001"] + [f"{nuc}.{i:03d}" for i, nuc in enumerate(values[f"{axis}mesh"], 2)]
        for axis in ("t", "e")
        if tags.get(f"{axis}tag") in ("nuc", "parent")
    ]
