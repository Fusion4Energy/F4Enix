from f4enix.constants import PathLike
from f4enix.meshtal_2.mesh2vtk_2.Modules.FMesh import FMesh
from f4enix.meshtal_2.mesh2vtk_2.Modules.mesh_parser import (
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


class NewMesh:
    def __init__(self, fmesh: FMesh):
        self._fmesh = fmesh

    def print_info(self) -> str:
        return str(self._fmesh.get_info())


class NewMeshtal:
    def __init__(self, filename: PathLike, filetype: str = "MCNP"):
        self.filetype = filetype

        self._meshtal_parser: Parser = PARSER_SELECTOR[self.filetype](filename)
        self.mesh: dict[int, NewMesh] = self._build_mesh_dict()

    def _build_mesh_dict(self) -> dict[int, NewMesh]:
        """Build a dictionary of meshes from the meshtal file."""
        mesh_dict = {}
        for mesh_id in self._meshtal_parser.get_meshlist():
            fmesh = self._meshtal_parser.get_FMesh(mesh_id)
            mesh_dict[mesh_id] = NewMesh(fmesh)
        return mesh_dict
