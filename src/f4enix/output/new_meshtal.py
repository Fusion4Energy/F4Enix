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


class NewMeshtal:
    def __init__(self, filename: PathLike, filetype: str = "MCNP"):
        self.filetype = filetype

        self._meshtal_parser: Parser = PARSER_SELECTOR[self.filetype](filename)
        self.mesh: dict[int, FMesh] = {}

    def readMesh(
        self,
        mesh: int | list[int] | None = None,
        cell_filters: list[int] | None = None,
        norm: str | None = None,
    ) -> None:
        """Read the meshtal file and build the mesh dictionary."""
        if not mesh:
            self.mesh = self._build_mesh_dict()
        elif isinstance(mesh, int):
            self.mesh = {mesh: self._meshtal_parser.get_FMesh(mesh, norm, cell_filters)}
        elif isinstance(mesh, list):
            self.mesh = {
                m: self._meshtal_parser.get_FMesh(m, norm, cell_filters) for m in mesh
            }

    def _build_mesh_dict(self) -> dict[int, FMesh]:
        """Build a dictionary of meshes from the meshtal file."""
        mesh_dict = {}
        for mesh_id in self._meshtal_parser.get_meshlist():
            mesh_dict[int(mesh_id)] = self._meshtal_parser.get_FMesh(mesh_id)
        return mesh_dict
