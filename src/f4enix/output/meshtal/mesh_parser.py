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
    def get_FMesh(self, tally: int, norm=None, filter=None) -> Fmesh: ...

    @abstractmethod
    def get_boundaries(self, tally): ...

    @abstractmethod
    def get_header(self, tally): ...


class MeshtalParser(Parser):
    """read meshtal formatted file (support block, column or CuV format)."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.tallyPos = _scan_meshfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of mesh stored in the meshtal"""
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(self, tally: int) -> Fmesh:
        """return Fmesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos, boundPos, dataPos = self.tallyPos[tally]
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
        fm = Fmesh(mesh, tally, trsf)
        fm.type = "meshtal"
        fm.particle, fm.comments = _get_header(self.fic, tallyPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        boundPos = self.tallyPos[tally][1]
        geom, trsf, meshbins = _get_mesh_boundaries(self.fic, boundPos)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.tallyPos[tally][0]
        particle, comments = _get_header(self.fic, tallyPos)
        return particle, comments


class CUVMeshParser(Parser):
    """read CUV formatted file."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.tallyPos = _scan_cuvfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of mesh stored in the CUV file"""
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(self, tally: int, norm=None, filter=None) -> Fmesh:
        """return Fmesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos, boundPos, dataPos = self.tallyPos[tally]
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
        fm = Fmesh(mesh, tally, trsf)
        fm.type = "CUV"
        fm.particle, fm.comments = _get_header(self.fic, tallyPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        boundPos = self.tallyPos[tally][1]
        geom, trsf, meshbins = _get_mesh_boundaries(self.fic, boundPos, cuv=True)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.tallyPos[tally][0]
        particle, comments = _get_header(self.fic, tallyPos)
        return particle, comments


class CDGSMeshParser(Parser):
    """read CDGS formatted file."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.srcmeshPos = _scan_cdgsfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of source meshes stored in the CCDGS file"""
        return tuple(sorted(self.srcmeshPos.keys()))

    def get_FMesh(self, tally: int, norm=None, filter=None) -> Fmesh:
        """return Fmesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.srcmeshPos.keys():
            print("bad tally entry")
            return
        meshPos, dataPos = self.srcmeshPos[tally]
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
        fm = Fmesh(mesh, tally, trsf)
        fm.type = "CDGS"
        fm.particle = None
        fm.cooling_time, fm.strength, fm.comments = _get_cdgsheader(self.fic, meshPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.srcmeshPos.keys():
            print("bad tally entry")
            return
        boundPos = self.srcmeshPos[tally][1]
        geom, trsf, meshbins = _get_cdgsmesh_boundaries(self.fic, boundPos, cuv=True)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.srcmeshPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.srcmeshPos[tally][0]
        particle, comments = _get_cdgsheader(self.fic, tallyPos)
        return particle, comments
