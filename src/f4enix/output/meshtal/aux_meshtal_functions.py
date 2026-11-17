from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from f4enix.output.meshtal.fmesh import Fmesh

import math
from typing import Any

import numpy as np
from vtkmodules.util import numpy_support
from vtkmodules.vtkCommonCore import VTK_FLOAT, vtkFloatArray, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid, vtkStructuredGrid

VALUE_LABEL = "Value - Total"
ERROR_LABEL = "Error - Total"
COLUMN_LABELS = (VALUE_LABEL, ERROR_LABEL)


class ExtraBin(np.ndarray):
    """class storing energy, time or specfic userbin"""

    def __new__(cls, data: np.ndarray, bintag: str, explicitbin: bool = True):
        obj = super().__new__(cls, data.shape)
        obj[:] = data[:]
        return obj

    def __init__(self, data: np.ndarray, bintag: str, explicitbin: bool = True):
        self._type = bintag
        self._explicit = explicitbin

        if bintag in ("erg", "tme"):
            self._binbound = True
            self._totalbin = len(data) != 2
        elif bintag in ("dtme", "nuc", "cel", "line"):
            self._binbound = False
            if bintag == "dtme":
                self._totalbin = False
            else:
                self._totalbin = len(data) != 1
        else:
            raise ValueError(
                f"Unknown bin type: {bintag}. Allowed types are 'erg', 'tme', 'dtme', 'nuc', 'cel', 'line'."
            )

        self._nvalue = len(data) - 1 if self._binbound else len(data)

    @property
    def nvalue(self) -> int:
        """
        number of intervals or discrte values
        """
        return self._nvalue

    @property
    def type(self) -> str:
        """
        type of data stored in the bin (erg, tme, nuc, cel, ...)
        """
        return self._type

    @property
    def explicit(self) -> bool:
        """
        True the bin is explcitly defined in the mesh definition.
        """
        return self._explicit

    @property
    def binbound(self) -> bool:
        """
        True if data are bin boundaries.
        """
        return self._binbound

    @property
    def totalbin(self) -> bool:
        """
        True if a total bin is included in data
        """
        return self._totalbin


def _makeVTKarray(array: np.ndarray) -> vtkFloatArray:
    return numpy_support.numpy_to_vtk(np.array(array).ravel(), array_type=VTK_FLOAT)


def _get_labels(
    nt: int, ne: int, labels: tuple[str, str] = COLUMN_LABELS
) -> list[list[tuple[str, str]]]:
    valstr, errstr = labels

    it_label = []
    if nt > 1 and ne > 1:
        for it in range(nt - 1):
            ie_label = []
            for ie in range(ne - 1):
                val_label = f"{valstr}_t{it + 1:03d}_e{ie + 1:03d}"
                err_label = f"{errstr}_t{it + 1:03d}_e{ie + 1:03d}"
                ie_label.append((val_label, err_label))
            val_label = f"{valstr}_t{it + 1:03d}_eTot"
            err_label = f"{errstr}_t{it + 1:03d}_eTot"
            ie_label.append((val_label, err_label))
            it_label.append(ie_label)

        ie_label = []
        for ie in range(ne - 1):
            val_label = f"{valstr}_tTot_e{ie + 1:03d}"
            err_label = f"{errstr}_tTot_e{ie + 1:03d}"
            ie_label.append((val_label, err_label))
        val_label = f"{valstr}_tTot_eTot"
        err_label = f"{errstr}_tTot_eTot"
        ie_label.append((val_label, err_label))
        it_label.append(ie_label)

    elif nt > 1:
        for it in range(nt - 1):
            val_label = f"{valstr}_t{it + 1:03d}"
            err_label = f"{errstr}_t{it + 1:03d}"
            ie_label = [(val_label, err_label)]
            it_label.append(ie_label)

        val_label = f"{valstr}_tTot"
        err_label = f"{errstr}_tTot"
        ie_label = [(val_label, err_label)]
        it_label.append(ie_label)

    elif ne > 1:
        ie_label = []
        for ie in range(ne - 1):
            val_label = f"{valstr}_e{ie + 1:03d}"
            err_label = f"{errstr}_e{ie + 1:03d}"
            ie_label.append((val_label, err_label))

        val_label = f"{valstr}_eTot"
        err_label = f"{errstr}_eTot"
        ie_label.append((val_label, err_label))
        it_label.append(ie_label)

    else:
        val_label = f"{valstr}"
        err_label = f"{errstr}"
        ie_label = [(val_label, err_label)]
        it_label.append(ie_label)

    return it_label


def _rectilinear_grid(
    mesh: "Fmesh", labels: tuple[str, str] = COLUMN_LABELS
) -> vtkRectilinearGrid:
    rgrid = vtkRectilinearGrid()
    rgrid.SetDimensions(mesh.nx1 + 1, mesh.nx2 + 1, mesh.nx3 + 1)

    bin1 = _makeVTKarray(mesh.x1bin)
    bin2 = _makeVTKarray(mesh.x2bin)
    bin3 = _makeVTKarray(mesh.x3bin)

    rgrid.SetXCoordinates(bin1)
    rgrid.SetYCoordinates(bin2)
    rgrid.SetZCoordinates(bin3)

    pdata = rgrid.GetCellData()

    bin_labels = _get_labels(mesh.nt, mesh.ne, labels)

    for it in range(mesh.nt):
        for ie in range(mesh.ne):
            vtkVal = _makeVTKarray(mesh.data[it, ie, :, :, :, 0])
            vtkErr = _makeVTKarray(mesh.data[it, ie, :, :, :, 1])
            vtkVal.SetName(bin_labels[it][ie][0])
            vtkErr.SetName(bin_labels[it][ie][1])
            pdata.AddArray(vtkVal)
            pdata.AddArray(vtkErr)
    return rgrid


def _structured_grid(
    mesh: "Fmesh", trsf=None, labels: tuple[str, str] = COLUMN_LABELS
) -> vtkStructuredGrid:
    if trsf is None:
        origin = (0.0, 0.0, 0.0)
    else:
        origin = trsf[0]

    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(mesh.nx1 + 1, mesh.nx2 + 1, mesh.nx3 + 1)

    pts = vtkPoints()
    pts.SetDataTypeToFloat()

    if mesh.geom == "cyl":
        if mesh.x3bin[1] > 0.5:
            mesh._x3bin = np.insert(mesh.x3bin, 1, 0.5, axis=0)
            mesh._nx3 += 1
            mesh._data = np.insert(mesh.data, 0, mesh.data[:, :, 0, :, :, :], axis=2)

    npts = (mesh.nx1 + 1) * (mesh.nx2 + 1) * (mesh.nx3 + 1)
    pts.SetNumberOfPoints(npts)

    n = 0
    if mesh.geom == "cyl":
        for t in mesh.x3bin:
            st = math.sin(t * 2 * math.pi)
            ct = math.cos(t * 2 * math.pi)
            for z in mesh.x2bin:
                z0 = z + origin[2]
                for r in mesh.x1bin:
                    x0 = r * ct + origin[0]
                    y0 = r * st + origin[1]
                    pts.InsertPoint(n, x0, y0, z0)
                    n += 1
    else:
        for z in mesh.x3bin:
            z0 = z + origin[2]
            for y in mesh.x2bin:
                y0 = y + origin[1]
                for x in mesh.x1bin:
                    pts.InsertPoint(n, x + origin[0], y0, z0)
                    n += 1

    sgrid.SetPoints(pts)

    pdata = sgrid.GetCellData()
    bin_labels = _get_labels(mesh.nt, mesh.ne, labels)

    for it in range(mesh.nt):
        for ie in range(mesh.ne):
            vtkVal = _makeVTKarray(mesh.data[it, ie, :, :, :, 0])
            vtkErr = _makeVTKarray(mesh.data[it, ie, :, :, :, 1])
            vtkVal.SetName(bin_labels[it][ie][0])
            vtkErr.SetName(bin_labels[it][ie][1])
            pdata.AddArray(vtkVal)
            pdata.AddArray(vtkErr)
    return sgrid


def _get_rotation_matrix(axs, vec):
    if vec is None:
        one = np.argmin(abs(axs))
        dummy = np.zeros(3)
        dummy[one] = 1.0
        vecX = np.cross(axs, dummy)
        vecX /= np.linalg.norm(vecX)
    else:
        vecX = vec

    vecY = np.cross(axs, vecX)
    rotmat = np.array((vecX, vecY, axs))
    return rotmat


class CustomTuple(tuple):
    origin: Any
    vec: Any
    axis: Any

    def __new__(cls, *args):
        return super().__new__(cls, args)
