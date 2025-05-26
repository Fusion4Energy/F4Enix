from pathlib import Path

from f4enix.constants import PathLike
from f4enix.meshtal_2.mesh2vtk_2.Modules.FMesh import FMesh
from f4enix.meshtal_2.mesh2vtk_2.Modules.mesh_parser import MeshtalParser


class NewMeshtal:
    def __init__(self, filename: PathLike, filetype: str = "MCNP"):
        self.filename: Path = Path(filename)
        self.filetype: str = filetype
