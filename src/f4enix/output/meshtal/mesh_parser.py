from abc import ABC, abstractmethod
from pathlib import Path

from f4enix.output.meshtal.aux_meshtal_functions import (
    _get_cdgsheader,
    _get_cdgsmesh_boundaries,
    _get_CDGSmesh_data,
    _get_CUVmesh_data,
    _get_header,
    _get_mesh_boundaries,
    _get_mesh_data,
    _scan_cdgsfile,
    _scan_cuvfile,
    _scan_meshfile,
    myOpen,
)
from f4enix.output.meshtal.fmesh import Fmesh, MeshData


class Parser(ABC):
    @abstractmethod
    def __init__(self, filename: str | Path): ...

    @abstractmethod
    def get_meshlist(self) -> tuple: ...

    @abstractmethod
    def get_FMesh(
        self, tally: int, norm: str | None = None, filter: list[int] | None = None
    ) -> Fmesh: ...


class MeshtalParser(Parser):
    """Read meshtal formatted file (support block, column or CuV format)."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.tallyPos = _scan_meshfile(self.fic)

    def get_meshlist(self) -> tuple:
        """Return list of mesh stored in the CUV file.

        Returns
        -------
        tuple
            A tuple containing the tally numbers of the meshes.
        """
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(
        self,
        tally: int,
        norm: str | None = None,
        filter: list[int] | None = None,
    ) -> Fmesh:
        """Returns Fmesh object representing a eshtally.

        Parameters
        ----------
        tally : int
            The mesh tally number to retrieve.
        norm : str | None, optional
            Normalization option for CUV meshes, not used.
        filter : list[int] | None, optional
            List of cells with which voxels in CUV meshes are filtered, not used.

        Returns
        -------
        Fmesh
            An Fmesh object containing the mesh data.
        """
        tally_str = str(tally)
        if tally_str not in self.tallyPos.keys():
            raise KeyError("bad tally entry")
        tallyPos, boundPos, dataPos = self.tallyPos[tally_str]
        geom, trsf, meshbins = _get_mesh_boundaries(self.fic, boundPos)

        nx1 = len(meshbins[0]) - 1
        nx2 = len(meshbins[1]) - 1
        nx3 = len(meshbins[2]) - 1
        ne = len(meshbins[3]) - 1 if meshbins[3].binbound else len(meshbins[3])
        nt = len(meshbins[4]) - 1 if meshbins[4].binbound else len(meshbins[4])

        if meshbins[3].totalbin:
            ne += 1
        if meshbins[4].totalbin:
            nt += 1

        shape = (nt, ne, nx3, nx2, nx1, 2)
        data = _get_mesh_data(self.fic, dataPos, geom, meshbins[4].explicit, shape)

        mesh = MeshData(geom, *meshbins, data)
        fm = Fmesh(mesh, tally_str, trsf)
        fm.type = "meshtal"
        fm.particle, fm.comments = _get_header(self.fic, tallyPos)

        return fm


class CUVMeshParser(Parser):
    """Read CUV formatted file."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.tallyPos = _scan_cuvfile(self.fic)

    def get_meshlist(self) -> tuple:
        """Return list of mesh stored in the CUV file.

        Returns
        -------
        tuple
            A tuple containing the tally numbers of the meshes.
        """
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(
        self, tally: int, norm: str | None = None, filter: list[int] | None = None
    ) -> Fmesh:
        """Returns Fmesh object representing a CUV meshtally.

        Parameters
        ----------
        tally : int
            The mesh tally number to retrieve.
        norm : str | None, optional
            Normalization option for CUV meshes.
        filter : list[int] | None, optional
            List of cells with which voxels in CUV meshes are filtered.

        Returns
        -------
        Fmesh
            An Fmesh object containing the mesh data.
        """
        tally_str = str(tally)
        if tally_str not in self.tallyPos.keys():
            raise KeyError("bad tally entry")
        tallyPos, boundPos, dataPos = self.tallyPos[tally_str]
        geom, trsf, meshbins = _get_mesh_boundaries(self.fic, boundPos, cuv=True)

        nx1 = len(meshbins[0]) - 1
        nx2 = len(meshbins[1]) - 1
        nx3 = len(meshbins[2]) - 1
        ne = len(meshbins[3]) - 1 if meshbins[3].binbound else len(meshbins[3])
        nt = len(meshbins[4]) - 1 if meshbins[4].binbound else len(meshbins[4])

        if meshbins[3].totalbin:
            ne += 1
        if meshbins[4].totalbin:
            nt += 1

        shape = (nt, ne, nx3, nx2, nx1, 2)
        data = _get_CUVmesh_data(self.fic, dataPos, shape, norm, filter)

        mesh = MeshData(geom, *meshbins, data)
        fm = Fmesh(mesh, tally_str, trsf)
        fm.type = "CUV"
        fm.particle, fm.comments = _get_header(self.fic, tallyPos)

        return fm


class CDGSMeshParser(Parser):
    """Read CDGS formatted file."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.srcmeshPos = _scan_cdgsfile(self.fic)

    def get_meshlist(self) -> tuple:
        """Return list of mesh stored in the CDGS file.

        Returns
        -------
        tuple
            A tuple containing the CDGS mesh ids.
        """
        return tuple(sorted(self.srcmeshPos.keys()))

    def get_FMesh(
        self, tally: int, norm: str | None = None, filter: list[int] | None = None
    ) -> Fmesh:
        """Return Fmesh object representing a CDGS source.
        The tally number is the source mesh number.

        Parameters
        ----------
        tally : int
            The source mesh number to retrieve.
        norm : str | None, optional
            Normalization option, not used in CDGS meshes.
        filter : list[int] | None, optional
            Filter options, not used in CDGS meshes.

        Returns
        -------
        Fmesh
            An Fmesh object containing the mesh data."""

        tally_str = str(tally)
        if tally_str not in self.srcmeshPos.keys():
            raise KeyError("bad tally entry")
        meshPos, dataPos = self.srcmeshPos[tally_str]
        geom, trsf, meshbins = _get_cdgsmesh_boundaries(self.fic, meshPos)

        nx1 = len(meshbins[0]) - 1
        nx2 = len(meshbins[1]) - 1
        nx3 = len(meshbins[2]) - 1
        ne = len(meshbins[3]) - 1 if meshbins[3].binbound else len(meshbins[3])
        nt = len(meshbins[4]) - 1 if meshbins[4].binbound else len(meshbins[4])

        if meshbins[3].totalbin:
            ne += 1
        if meshbins[4].totalbin:
            nt += 1

        shape = (nt, ne, nx3, nx2, nx1, 2)
        data = _get_CDGSmesh_data(self.fic, dataPos, shape)

        mesh = MeshData(geom, *meshbins, data)
        fm = Fmesh(mesh, tally_str, trsf)
        fm.type = "CDGS"
        fm.particle = None
        fm.cooling_time, fm.strength, fm.comments = _get_cdgsheader(self.fic, meshPos)

        return fm
