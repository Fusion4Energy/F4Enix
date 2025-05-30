import re
from copy import deepcopy

import pyvista as pv
from numjuggler import parser

from f4enix.constants import PathLike
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


class Meshtal:
    def __init__(self, filename: PathLike, filetype: str = "MCNP"):
        self.filetype = filetype

        self._meshtal_parser: Parser = PARSER_SELECTOR[self.filetype](filename)
        self.mesh: dict[int, Fmesh] = {}

    def create_filtered_mesh(
        self,
        mesh_id: int,
        binlabels: tuple[str],
        ebin: int | None = None,
        tbin: int | None = None,
    ) -> Fmesh:
        etbin_data = self.mesh[mesh_id].get_etbin_data(ebin=ebin, tbin=tbin)
        etbin_mesh = Fmesh(etbin_data, str(mesh_id), trsf=None, binlabels=binlabels)
        return etbin_mesh

    def readMesh(
        self,
        mesh: int | list[int] | None = None,
        norm: str | None = None,
        cell_filters: list[int] | None = None,
    ) -> None:
        """Read the meshtal file and build the mesh dictionary."""
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
        for key, data_card in inp.other_data.items():
            if "FMESH" in key:
                fmesh_num = data_card.name
                if fmesh_num not in self.mesh.keys():
                    continue
                match = False
                for line in data_card.lines:
                    # Search for the pattern in the line
                    match = re.search(pattern, line)
                    if match:
                        break
                if match:
                    transf_dict[fmesh_num] = inp.transformations[
                        "TR" + str(match.group().split("=")[-1])
                    ]
        self.transform_multiple_fmesh(transf_dict)

    def transform_multiple_fmesh(self, transf_dict: dict[int, parser.Card]) -> None:
        """Transforms multiple fmeshes in the meshtal object.
        Given a dictionary of fmeshes numbers and transformation cards,
        rototranslate the fmeshes in the meshtal object according to the associated
        transformation card in the dict


        Parameters
        ----------
        transf_dict : dict[int, parser.Card]
            dictionary of fmeshes numbers and transformation cards.
        """
        for key, transf in transf_dict.items():
            self.mesh[key].apply_transformation(transf)
