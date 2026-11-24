from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from f4enix.output.meshtal.fmesh import Fmesh
    from f4enix.output.meshtal import Meshtal
# from os import PathLike
from f4enix.core.constants import PathLike

import math
from typing import Any

import numpy as np
import re
from vtkmodules.util import numpy_support
from vtkmodules.vtkCommonCore import VTK_FLOAT, vtkFloatArray, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid, vtkStructuredGrid

VALUE_LABEL = "Value - Total" # this is not universal, energy bins will have more keys
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


def opt_merge_arrays(A:np.ndarray, B:np.ndarray,
                     sigmaA:np.ndarray, sigmaB:np.ndarray):
    """
    Optimal merging is a weighted average which consideres
    uncertainties when asigning the weights.

    Arrays `A`, `B`, `sigmaA` and `sigmaB` must have the
    same shape.
    
    Returns tuple `(M, sigmaM)` where M is optimally merged
    array and sigmaM are standard deviations.
    """

    u = np.square(sigmaA) + np.square(sigmaB)
    u_nonzero = np.not_equal(u, 0)

    _div = np.zeros(shape=u.shape)
    np.divide(
        np.multiply(np.square(sigmaA), np.square(sigmaB)),
        u,
        where=u_nonzero,
        out=_div)
    opt_sigma = np.sqrt(_div)

    opt_value = np.zeros(shape=u.shape)
    np.divide(
            np.multiply(A, np.square(sigmaB)) +
            np.multiply(B, np.square(sigmaA)),
        u,
        where=u_nonzero,
        out=opt_value)

    # We need to determine where these arrays have values zero.
    # If one of the two arrays has value zero for some voxel,
    # we cannot optimaly merge them.
    # We should instead take the non zero value of the two.
    # If both arrays have a zero value at some voxel,
    # we shoud just return (0,0) as the optimaly merged value.

    # if element i is:
    # A_i == 0 and B_i == 0 --> M_i = 0
    # A_i != 0 and B_i == 0 --> M_i = A_i
    # A_i == 0 and B_i != 0 --> M_i = B_i
    # A_i != 0 and B_i != 0 --> M_i = merge(A_i, B_i)

    zeroA = np.equal(A, 0)
    zeroB = np.equal(B, 0)
    AnonzeroBzero 	= np.logical_and( np.logical_not(zeroA), zeroB)
    BnonzeroAzero 	= np.logical_and( np.logical_not(zeroB), zeroA)
    ABnonzero		= np.logical_and( np.logical_not(zeroA), np.logical_not(zeroB))

    shape = A.shape
    zeros = np.zeros(shape)	
    
    M = ( np.where(ABnonzero, opt_value, zeros)
        + np.where(AnonzeroBzero, A, zeros)
        + np.where(BnonzeroAzero, B, zeros)
        + zeros)
    sigmaM = ( np.where(ABnonzero, opt_sigma, zeros)
        + np.where(AnonzeroBzero, sigmaA, zeros)
        + np.where(BnonzeroAzero, sigmaB, zeros)
        + zeros)
    
    return (M, sigmaM)

#A# THIS IS WAY TO COMPLICATED ...
#A# def meshtal_to_ndarray(meshtal : "Meshtal"):
#A#     """
#A#     Convert `Meshtal` with mesh attribute to numpy array.
#A#     Numpy array will have the shape (M,N), where `N` is the number
#A#     of spatial voxels and M is the number of non spatial voxels
#A#     eg.: energy bins.
#A#     
#A#     Reads `Meshtal.mesh.grid.cell_data` dictionary and returns new
#A#     dictionary based on the keys of `cell_data`. 
#A#     
#A#     
#A#     Turns:
#A#     >>> mesh[tally_id].grid.cell_data = {
#A#     >>>     'Value - Total' : np.ndarray,
#A#     >>>     'Error - Total' : np.ndarray,
#A#     >>>     'Value - e001'  : np.ndarray,
#A#     >>>     'Error - e001'  : np.ndarray}
#A#     
#A#     Into:
#A#     >>> {tally_id : {
#A#     >>> 	'Value' : array([[Total],[e001]])
#A#     >>>     'Error' : array([[Total],[e001]])
#A#     >>> }}
#A#     """
#A# 
#A#     array_dict = {}
#A# 
#A#     for tally_id in meshtal.mesh.keys():
#A#         array_dict[tally_id] = {}
#A#         # print(f'{tally_id = }')
#A#         quantities = [VALUE_LABEL, ERROR_LABEL]
#A#         bin_matches = [re.match(rf'{VALUE_LABEL}(.*)', k) for k in
#A#                    meshtal.mesh[tally_id].grid.cell_data.keys()]
#A#         bins = np.unique([m.group(1) for m in bin_matches if m]) # bins can also be array ['']
#A#         # quantities = np.unique([m.group(1) for m in matches])
#A#         # bins = np.unique([m.group(2) for m in matches])
#A#         print(f'{quantities = }')
#A#         print(f'{bins = }')
#A# 
#A#         for q in quantities:
#A#             array_dict[tally_id][q] = np.array(
#A#                 [meshtal.mesh[tally_id].grid.cell_data[f'{q}{b}']
#A#                  for b in bins])
#A# 
#A#     return array_dict

#B# THIS IS WAY TOO COMPLICATED ...
#B# def ndarray_to_meshtal(array_dict : dict[np.ndarray], bin_keys : list = []):
#B#     """Reverse of meshtal_to_ndarray().
#B#     
#B#     The order of bin_keys matters! """
#B#     from f4enix.output.meshtal import Meshtal
#B#     meshtal = Meshtal.__new__(Meshtal)
#B#     meshtal.filetype = "MCNP"
#B#     meshtal.mesh = {}
#B#     for tally_id, quantity_dict in array_dict.items():
#B#         meshtal.mesh[tally_id] = pyvista.Grid ....
#B#         quantities = [VALUE_LABEL, ERROR_LABEL]
#B#         for q in quantities:
#B#             for ib, b in enumerate(bin_keys):
#B#                 meshtal.mesh[tally_id].cell_data[f'{q}{b}'] = quantity_dict[q][ib,:]



def init_meshtal(meshtal : PathLike | Meshtal, filetype : str = "MCNP"):
    """Return a `Meshtal` object with a mesh. Initialize from file or
    return back if it is already an instance of `Meshtal`."""
    from f4enix.output.meshtal import Meshtal

    if isinstance(meshtal, PathLike):
        mt = Meshtal(meshtal, filetype=filetype)
        mt.readMesh()
        return mt
    elif isinstance(meshtal, Meshtal):
        if meshtal.mesh:
            return meshtal
        else:
            meshtal.readMesh()
            return meshtal
    else:
        raise ValueError(f"accept_meshtal() expects meshtal as PathLike or Meshtal object. Unknown type provided: {type(meshtal)}")