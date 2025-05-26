from pathlib import Path

from .FMesh import FMesh, MeshData
from .functions.read_functions import (
    get_cdgsheader,
    get_cdgsmesh_boundaries,
    get_CDGSmesh_data,
    get_CUVmesh_data,
    get_header,
    get_mesh_boundaries,
    get_mesh_data,
    scan_cdgsfile,
    scan_cuvfile,
    scan_meshfile,
)
from .functions.utils import myOpen


class MeshtalParser:
    """read meshtal formatted file (support block, column or CuV format)."""

    def __init__(self, filename: str | Path):
        self._filename = filename
        self._fic = myOpen(filename, "r")
        self._tallyPos = scan_meshfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of mesh stored in the meshtal"""
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(self, tally) -> FMesh:
        """return FMesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos, boundPos, dataPos = self.tallyPos[tally]
        geom, trsf, meshbins = get_mesh_boundaries(self.fic, boundPos)

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
        data = get_mesh_data(self.fic, dataPos, geom, meshbins[4].explicit, shape)

        mesh = MeshData(geom, *meshbins, data)
        fm = FMesh(mesh, tally, trsf)
        fm.type = "meshtal"
        fm.particle, fm.comments = get_header(self.fic, tallyPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        boundPos = self.tallyPos[tally][1]
        geom, trsf, meshbins = get_mesh_boundaries(self.fic, boundPos)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.tallyPos[tally][0]
        particle, comments = get_header(self.fic, tallyPos)
        return particle, comments


class CUVMeshParser:
    """read CUV formatted file."""

    def __init__(self, filename: str):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.tallyPos = scan_cuvfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of mesh stored in the CUV file"""
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(self, tally: int, norm=None, filter=None) -> FMesh:
        """return FMesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos, boundPos, dataPos = self.tallyPos[tally]
        geom, trsf, meshbins = get_mesh_boundaries(self.fic, boundPos, cuv=True)

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
        data = get_CUVmesh_data(self.fic, dataPos, shape, norm, filter)

        mesh = MeshData(geom, *meshbins, data)
        fm = FMesh(mesh, tally, trsf)
        fm.type = "CUV"
        fm.particle, fm.comments = get_header(self.fic, tallyPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        boundPos = self.tallyPos[tally][1]
        geom, trsf, meshbins = get_mesh_boundaries(self.fic, boundPos, cuv=True)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.tallyPos[tally][0]
        particle, comments = get_header(self.fic, tallyPos)
        return particle, comments


class CDGSMeshParser:
    """read CDGS formatted file."""

    def __init__(self, filename: str):
        self.filename = filename
        self.fic = myOpen(filename, "r")
        self.srcmeshPos = scan_cdgsfile(self.fic)

    def get_meshlist(self) -> tuple:
        """return list of source meshes stored in the CCDGS file"""
        return tuple(sorted(self.tallyPos.keys()))

    def get_FMesh(self, tally: int, norm=None, filter=None) -> FMesh:
        """return FMesh object with all data of the mesh tally number"""
        tally = str(tally)
        if tally not in self.srcmeshPos.keys():
            print("bad tally entry")
            return
        meshPos, dataPos = self.srcmeshPos[tally]
        geom, trsf, meshbins = get_cdgsmesh_boundaries(self.fic, meshPos)

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
        data = get_CDGSmesh_data(self.fic, dataPos, shape)

        mesh = MeshData(geom, *meshbins, data)
        fm = FMesh(mesh, tally, trsf)
        fm.type = "CDGS"
        fm.particle = None
        fm.cooling_time, fm.strength, fm.comments = get_cdgsheader(self.fic, meshPos)

        return fm

    def get_boundaries(self, tally):
        """return mesh information type, boundaries, ..."""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        boundPos = self.tallyPos[tally][1]
        geom, trsf, meshbins = get_cdgsmesh_boundaries(self.fic, boundPos, cuv=True)
        return geom, trsf, meshbins

    def get_header(self, tally):
        """return mesh header information"""
        tally = str(tally)
        if tally not in self.tallyPos.keys():
            print("bad tally entry")
            return
        tallyPos = self.tallyPos[tally][0]
        particle, comments = get_cdgsheader(self.fic, tallyPos)
        return particle, comments
