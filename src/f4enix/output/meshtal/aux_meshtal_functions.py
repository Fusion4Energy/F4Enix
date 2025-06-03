from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from f4enix.output.meshtal.fmesh import Fmesh

import math
from pathlib import Path
from typing import Any

import numpy as np
import vtk
from vtk.util import numpy_support

VALUE_LABEL = "Value - Total"
ERROR_LABEL = "Error - Total"
COLUMN_LABELS = (VALUE_LABEL, ERROR_LABEL)


class myOpen:
    def __init__(self, filename: str | Path, mode: str):
        self._file = open(filename, mode)

    def _readline(self):
        return self._file.readline()

    def _tell(self):
        return self._file.tell()

    def _seek(self, pos: int):
        self._file.seek(pos)

    def _skipline(self, n: int = 1, verbose: bool = False):
        for i in range(n):
            line = self._file.readline()
            if verbose:
                print(f"skipped: {line}")

    def _get_multilines(self, nvalues: int):
        values = []
        while len(values) < nvalues:
            line = self._file.readline()
            values.extend(line.split())
        return np.array(values, dtype=np.double)


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


def _read_block(fic: myOpen, nl: int, array: np.ndarray) -> None:
    fic._skipline(2)
    for j in range(2):
        fic._skipline(2)
        for i in range(nl):
            line = fic._readline()
            array[i, :, j] = [float(x) for x in line.split()[1:]]
        fic._skipline()


def _get_block_data(
    fic: myOpen,
    position: int,
    type: str,
    explicit_time: bool,
    shape: tuple[int, int, int, int, int, int],
    dtme: bool = False,
) -> np.ndarray:
    fic._seek(position)
    if type == "ij":
        nc = shape[4]  # X
        nl = shape[3]  # Y
        nb = shape[2]  # Z
    elif type == "ik":
        nc = shape[4]  # X
        nl = shape[2]  # Z
        nb = shape[3]  # Y
    elif type == "jk":
        nc = shape[3]  # Y
        nl = shape[2]  # Z
        nb = shape[4]  # X

    ne = shape[1]
    nt = shape[0]
    blockdata = np.ndarray((nl, nc, 2))
    wrkshape = (nt, ne, nb, nl, nc, 2)
    data = np.ndarray(wrkshape)

    # for ttag != dtme outer loop over energies
    if not explicit_time:
        for ie in range(ne):
            fic._skipline(2)
            for ib in range(nb):
                _read_block(fic, nl, blockdata)
                data[0, ie, ib, :, :, :] = blockdata[:, :, :]
                fic._skipline()
            fic._skipline()

    else:
        for ie in range(ne):
            fic._skipline(3)
            for it in range(nt):
                fic._skipline(2)
                for ib in range(nb):
                    _read_block(fic, nl, blockdata)
                    data[it, ie, ib, :, :, :] = blockdata[:, :, :]
                    fic._skipline()
                fic._skipline()

    if type == "jk":
        data = np.moveaxis(data, 2, 4)
    elif type == "ik":
        data = np.swapaxes(data, 2, 3)

    return data


def _get_column_data(
    fic: myOpen, position: int, type: str, shape: tuple[int, int, int, int, int, int]
) -> np.ndarray:
    fic._seek(position)

    nt, ne, nx3, nx2, nx1, _ = shape
    if type == "col":
        ival = -2
        ierr = -1
    else:
        ival = -4
        ierr = -3

    fic._skipline(1)
    col_shape = (nt, ne, nx1, nx2, nx3, 2)
    data = np.ndarray(col_shape)
    for ie in range(ne):
        for it in range(nt):
            for i1 in range(nx1):
                for i2 in range(nx2):
                    for i3 in range(nx3):
                        strvalues = fic._readline().split()
                        val = float(strvalues[ival])
                        err = float(strvalues[ierr])
                        data[it, ie, i1, i2, i3, :] = (val, err)
    data = np.swapaxes(data, 2, 4)
    return data


def _read_cuv_cell_info(
    fic: myOpen, ncell: int
) -> tuple[list[tuple[int, float]], list[np.ndarray], list[np.ndarray]]:
    cell_info = []
    tot_value = []
    tot_error = []
    for ic in range(ncell):
        line = fic._readline().split()
        cell = int(float(line[0]))
        volf = float(line[1])
        cell_info.append((cell, volf))
        if len(line) == 4:
            tot_value.append(float(line[2]))
            tot_error.append(float(line[3]))

    addtot = len(tot_value) > 0

    val_tab = []
    err_tab = []
    for ic in range(ncell):
        line = fic._readline().split()
        values = np.array(line, dtype=np.double)
        line = fic._readline().split()
        errors = np.array(line, dtype=np.double)

        if addtot:
            np.append(values, tot_value[ic])
            np.append(errors, tot_error[ic])
        val_tab.append(values)
        err_tab.append(errors)

    return cell_info, val_tab, err_tab


def _read_cdgs_cell_info(
    fic: myOpen, ncell: int, ne: int
) -> tuple[list[tuple[int, float]], list[np.ndarray], list[np.ndarray]]:
    cell_info = []
    val_tab = []
    err_tab = []

    for ic in range(ncell):
        headline = fic._readline().split()
        cell = int(float(headline[0]))
        volf = float(headline[1])
        cell_info.append((cell, volf))

        values = fic._get_multilines(ne)
        errors = fic._get_multilines(ne)
        if ne != 1:
            nh = len(headline)
            if nh == 2:
                totval = np.sum(values)
                errval = 0.0
            elif nh == 3:
                totval = float(headline[2])
                errval = 0.0
            else:
                totval = float(headline[2])
                errval = float(headline[3])

            values = np.append(values, totval)
            errors = np.append(errors, errval)

        val_tab.append(values)
        err_tab.append(errors)

    return cell_info, val_tab, err_tab


def _select_cells(
    cell_info: list[tuple[int, float]],
    values: list[np.ndarray],
    errors: list[np.ndarray],
    filter: list[int],
) -> tuple[float, list[tuple[int, float]], list[np.ndarray], list[np.ndarray]]:
    new_cell_info = []
    new_values_tab = []
    new_errors_tab = []
    sumf = 0.0

    for cinfo, val, err in zip(cell_info, values, errors):
        if cinfo[0] in filter:
            new_cell_info.append(cinfo)
            new_values_tab.append(val)
            new_errors_tab.append(err)
            sumf += cinfo[1]
    return sumf, new_cell_info, new_values_tab, new_errors_tab


def _get_mean_value(
    cell_info: list[tuple[int, float]],
    val: list[np.ndarray],
    err: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    val_data = np.array(val, dtype=np.double)
    err_data = np.array(err, dtype=np.double)
    err_data = val_data * err_data
    tmp = np.array(cell_info)
    frac = tmp[:, 1]
    val_mean = np.matmul(frac, val_data)
    err_mean = np.where(val_mean != 0.0, np.matmul(frac, err_data) / val_mean, 0.0)
    return val_mean, err_mean


def _get_cuv_element(
    fic: myOpen, norm: str | None, filter: list[int] | None
) -> tuple[np.ndarray, np.ndarray]:
    line = fic._readline().split()
    volume = float(line[4])
    ncell = int(float(line[5]))
    cell_info, values, errors = _read_cuv_cell_info(fic, ncell)

    ne = len(values[0])

    if filter is not None:
        sumf, cell_info, values, errors = _select_cells(
            cell_info, values, errors, filter
        )
    else:
        sumf = 1.0

    if sumf == 0.0:
        values = np.zeros(ne)
        errors = np.zeros(ne)
        return values, errors

    if ncell > 1:
        values, errors = _get_mean_value(cell_info, values, errors)
    else:
        values = values[0]
        errors = errors[0]
    if norm == "ctot":
        values = values * volume
    elif norm == "celf":
        values = values / sumf
    else:
        raise ValueError(
            f"Unknown normalization option: {norm}. Allowed options are 'ctot' or 'celf'."
        )

    return values, errors


def _get_cdgs_element(
    fic: myOpen, ne: int, norm: str | None = None, filter: list[int] | None = None
) -> (
    tuple[int, np.ndarray, np.ndarray]
    | tuple[np.ndarray, np.ndarray]
    | tuple[None, None, None]
):
    line = fic._readline().split()
    if "end_source_data" == line[0] or "" == line[0]:
        return None, None, None

    index = int(float(line[0]))
    volume = float(line[-2])
    ncell = int(float(line[-1]))

    cell_info, values, errors = _read_cdgs_cell_info(fic, ncell, ne)

    if filter is not None:
        sumf, cell_info, values, errors = _select_cells(
            cell_info, values, errors, filter
        )
    else:
        sumf = 0.0
        for cf in cell_info:
            sumf += cf[1]

    if sumf == 0.0:
        values = np.zeros(ne)
        errors = np.zeros(ne)
        return values, errors

    if ncell > 1:
        values, errors = _get_mean_value(cell_info, values, errors)
    else:
        values = values[0]
        errors = errors[0]

    if norm == "total":
        values = values * volume
    elif norm == "fraction":
        values = values / sumf

    return index, values, errors


def _makeVTKarray(array: np.ndarray) -> vtk.vtkFloatArray:
    return numpy_support.numpy_to_vtk(np.array(array).ravel(), array_type=vtk.VTK_FLOAT)


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
) -> vtk.vtkRectilinearGrid:
    rgrid = vtk.vtkRectilinearGrid()
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
) -> vtk.vtkStructuredGrid:
    if trsf is None:
        origin = (0.0, 0.0, 0.0)
    else:
        origin = trsf[0]

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(mesh.nx1 + 1, mesh.nx2 + 1, mesh.nx3 + 1)

    pts = vtk.vtkPoints()
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
        dummy = np.array((0, 0, 0))
        dummy[one] = 1.0
        vecX = axs.cross(dummy)
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


def _get_transformation(line: str, cyl=True) -> tuple[np.ndarray, np.ndarray] | None:
    if cyl:
        if "VEC direction" in line:
            fmcnp6 = True
            org0 = line.index("at") + 2
            org1 = line.index("axis")

            axs0 = line.index("axis in") + 7
            axs1 = line.index("direction")

            vec0 = line.index("direction", axs1 + 1) + 9
        else:
            fmcnp6 = False
            org0 = line.index("at") + 2
            org1 = line.index(",")

            axs0 = line.index("axis in") + 7
            axs1 = line.index("direction")

        origin = np.array(line[org0:org1].split(), dtype=np.double)
        axis = np.array(line[axs0:axs1].split(), dtype=np.double)
        if fmcnp6:
            vec = np.array(line[vec0:].split(), dtype=np.double)
        else:
            vec = None

        axisZ = abs(axis.dot(np.array((0, 0, 1.0))) - 1.0) < 1e-12
        if vec is None:
            axisX = None
        else:
            axisX = abs(axis.dot(np.array((1.0, 0, 0))) - 1.0) < 1e-12

        if axisZ and (axisX is None or axisX):
            rotmat = None
        else:
            rotmat = _get_rotation_matrix(axis, vec)
    else:
        org0 = 0
        org1 = 44
        vecx0 = 45
        vecx1 = 87
        vecy0 = 88
        origin = np.array(line[org0:org1].split(), dtype=np.double)
        vecX = np.array(line[vecx0:vecx1].split(), dtype=np.double)
        vecY = np.array(line[vecy0:].split(), dtype=np.double)
        axisX = abs(vecX.dot(np.array((1.0, 0, 0.0))) - 1.0) < 1e-12
        axisY = abs(vecY.dot(np.array((0, 1.0, 0))) - 1.0) < 1e-12
        if axisX and axisY:
            rotmat = None
        else:
            rotmat = _get_rotation_matrix(vecX, vecY)

    if np.linalg.norm(origin) < 1e-12 and rotmat is None:
        return None
    if rotmat is None:
        rotmat = np.identity(3)
    result = CustomTuple(origin, rotmat)
    result.origin = origin
    result.vec = vec
    result.axis = axis
    return result


def _get_etbin_tag(line: str) -> str | None:
    if "Energy" in line:
        return "erg"
    elif "Decay time" in line:
        return "dtme"
    elif "Time" in line:
        return "tme"
    elif "Nuclide" in line:
        return "nuc"
    elif "Cell" in line:
        return "cel"


def _scan_meshfile(fic: myOpen) -> dict[str, tuple[int, int, int]]:
    tally = dict()
    while True:
        pos = fic._tell()
        line = fic._readline()[0:30]
        if " Mesh Tally Number" == line[0:18]:
            tnum = line[18:].strip()
            boundaries, data = _get_mesh_block(fic)
            block_pos = (pos, boundaries, data)
            tally[tnum] = block_pos
        elif line == "":
            break
    return tally


def _scan_cuvfile(fic: myOpen) -> dict[str, tuple[int, int, int]]:
    cuvmesh = dict()
    while True:
        pos = fic._tell()
        line = fic._readline()[0:30]
        if " Mesh Tally Number" == line[0:18]:
            tnum = line[18:].strip()
            boundaries, data = _skip_cuvdata(fic)
            cuvmesh[tnum] = (pos, boundaries, data)
        elif line == "":
            break
    return cuvmesh


def _scan_cdgsfile(fic: myOpen) -> dict[str, list]:
    srcmesh = dict()
    line = fic._readline().split()
    nmesh = int(float(line[1]))
    line = fic._readline().split()
    total_strength = float(line[1])

    for i in range(nmesh):
        pos = fic._tell()
        line = fic._readline()
        if "mesh_id" == line[0:7]:
            mshnum = line.split()[1]
            srcmesh[mshnum] = [pos]
            while True:
                pos = fic._tell()
                line = fic._readline()[0:11]
                if line == "source_data":
                    srcmesh[mshnum].append(pos)
                    break
            while True:
                line = fic._readline()[0:15]
                if line == "end_source_data":
                    break
    return srcmesh


def _skip_cuvdata(fic: myOpen) -> tuple[int, int]:
    line = fic._readline()[0:30]
    while line[0:22] != " Tally bin boundaries:":
        bound_pos = fic._tell()
        line = fic._readline()[0:30]
        if line == "":
            break

    geom, trsf, meshbins = _get_mesh_boundaries(fic)
    x1bin, x2bin, x3bin, ebin, tbin = meshbins
    n1 = len(x1bin) - 1
    n2 = len(x2bin) - 1
    n3 = len(x3bin) - 1
    nelemts = n1 * n2 * n3

    data_pos = fic._tell()
    fic._skipline(1)
    for i in range(nelemts):
        headline = fic._readline().split()
        index = int(float(headline[0]))
        ncel = int(float(headline[-1]))
        fic._skipline(3 * ncel)
    return (bound_pos, data_pos)


def _get_mesh_block(fic: myOpen) -> tuple[int, int]:
    blk_line = 0
    while blk_line < 2:
        line = fic._readline()
        pos = fic._tell()
        if line.strip(" ") == "\n":
            if blk_line == 0:
                bound_pos = pos
            elif blk_line == 1:
                data_pos = pos
            blk_line += 1
        elif line == "":
            break
    return (bound_pos, data_pos)


def _get_header(fic: myOpen, position: int) -> tuple[str, str | None]:
    particle_list = ("neutron", "photon", "electron", "proton", "deuteron")
    fic._seek(position)
    fic._skipline()
    comments = ""
    while True:
        line = fic._readline()
        if "mesh tally." in line:
            break
        comments += line
    if comments == "":
        comments = None

    part = ""
    for p in particle_list:
        if p in line:
            part = p
            break
    return part, comments


def _get_cdgsheader(fic: myOpen, position: int) -> tuple[float, float, str]:
    fic._seek(position)
    fic._skipline()
    comments = fic._readline()
    line = fic._readline().split()
    cooling = float(line[1])
    line = fic._readline().split()
    strength = float(line[1])

    return cooling, strength, comments


def _get_mesh_type(
    fic: myOpen, position: int, geom: str = "rec", explicit_time: bool = False
) -> str:
    fic._seek(position)
    line = fic._readline()

    if "Result" in line:
        return "col"
    elif "GEOMETRY" in line:
        return "oldCUV"

    if explicit_time:
        fic._skipline(4)
    else:
        fic._skipline()
    line = fic._readline()
    if geom == "rec":
        if "Z bin" in line:
            return "ij"
        elif "Y bin" in line:
            return "ik"
        elif "X bin" in line:
            return "jk"
        else:
            return "bad"
    else:
        if "Theta bin" in line:
            return "ij"
        elif "Z bin" in line:
            return "ik"
        elif "R bin" in line:
            return "jk"
        else:
            return "bad"


def _get_mesh_boundaries(
    fic: myOpen, position: int | None = None, cuv: bool | None = None
):
    if position is not None:
        fic._seek(position)
        fic._skipline()
    line = fic._readline()

    if "origin at" in line:
        trsf = _get_transformation(line)
        geom = "cyl"
        line = fic._readline()
    else:
        trsf = None
        geom = "rec"

    init = line.index(":") + 1
    x1bin = np.array(line[init:].split(), dtype=np.double)

    line = fic._readline()
    init = line.index(":") + 1
    x2bin = np.array(line[init:].split(), dtype=np.double)

    line = fic._readline()
    init = line.index(":") + 1
    x3bin = np.array(line[init:].split(), dtype=np.double)

    line = fic._readline()
    if "number:" in line:
        line = fic._readline()

    init = line.index(":") + 1

    bindata = np.array(line[init:].split(), dtype=np.double)
    bintag = _get_etbin_tag(line[0:init])
    tmpbin = ExtraBin(bindata, bintag)

    line = fic._readline()
    if line.strip(" ") == "\n":
        ebin = tmpbin
        bindata = np.array((0, 1e20), dtype=np.double)
        bintag = "tme"
        tbin = ExtraBin(bindata, bintag, explicitbin=False)
    else:
        tbin = tmpbin
        init = line.index(":") + 1
        bindata = np.array(line[init:].split(), dtype=np.double)
        bintag = _get_etbin_tag(line[0:init])
        ebin = ExtraBin(bindata, bintag)

    return geom, trsf, (x1bin, x2bin, x3bin, ebin, tbin)


def _get_cdgsmesh_boundaries(
    fic: myOpen, position: int
) -> tuple[
    str,
    tuple[np.ndarray, np.ndarray] | None,
    tuple[np.ndarray, np.ndarray, np.ndarray, ExtraBin, ExtraBin],
]:
    fic._seek(position)
    fic._skipline(4)
    line = fic._readline().split()
    ergbin = line[1] == "bins"
    line = fic._readline().split()
    nerg = int(float(line[1]))

    etag = "erg" if ergbin else "line"
    ebin = ExtraBin(fic._get_multilines(nerg), etag)
    tbin = ExtraBin(np.array((0, 1e20)), "tme", explicitbin=False)

    line = fic._readline().split()
    geom = line[1]

    line = fic._readline().split()
    nx1 = float(int(line[1]))
    nx2 = float(int(line[2]))
    nx3 = float(int(line[3]))

    line = fic._readline()
    trsf = _get_transformation(line, cyl=geom == "cyl")

    x1bin = fic._get_multilines(nx1)
    x2bin = fic._get_multilines(nx2)
    x3bin = fic._get_multilines(nx3)

    return geom, trsf, (x1bin, x2bin, x3bin, ebin, tbin)


def _get_mesh_data(
    fic: myOpen,
    position: int,
    geom: str,
    explicit_time: bool,
    shape: tuple[int, int, int, int, int, int],
) -> np.ndarray:
    type = _get_mesh_type(fic, position, geom, explicit_time)
    if type == "col" or type == "cf":
        data = _get_column_data(fic, position, type, shape)
    elif type == "oldCUV":
        # in mcnp5 version of d1suned CUV was written in meshtal files with CDGS format
        data = _get_CDGSmesh_data(fic, position, shape, oldCUV=True)
    else:
        data = _get_block_data(fic, position, type, explicit_time, shape)
    return data


def _get_CUVmesh_data(
    fic: myOpen,
    position,
    shape: tuple[int, int, int, int, int, int],
    norm: str | None = None,
    filter: list[int] | None = None,
) -> np.ndarray:
    fic._seek(position)
    fic._skipline()

    nt, ne, nx3, nx2, nx1, _ = shape
    nelement = nx1 * nx2 * nx3
    data = np.ndarray(shape)
    val_array = np.ndarray((nelement, ne))
    err_array = np.ndarray((nelement, ne))

    for it in range(nt):
        for i in range(nelement):
            val, err = _get_cuv_element(fic, norm, filter)
            val_array[i, :] = val[:]
            err_array[i, :] = err[:]

        val_array = np.transpose(val_array.reshape((nx1, nx2, nx3, ne)))
        err_array = np.transpose(err_array.reshape((nx1, nx2, nx3, ne)))
        data[it, :, :, :, :, 0] = val_array[:, :, :, :]
        data[it, :, :, :, :, 1] = err_array[:, :, :, :]
    return data


def _get_CDGSmesh_data(
    fic: myOpen,
    position: int,
    shape: tuple[int, int, int, int, int, int],
    oldCUV: bool = False,
) -> np.ndarray:
    fic._seek(position)
    if oldCUV:
        fic._skipline(4)
    else:
        fic._skipline()

    nt, ne, nx3, nx2, nx1, _ = shape
    nelement = nx1 * nx2 * nx3
    data = np.ndarray(shape)
    val_array = np.ndarray((nelement, ne))
    err_array = np.ndarray((nelement, ne))
    zero = np.zeros(ne)
    nextcell = True
    end_data = False
    ne_read = ne if ne == 1 else ne - 1
    for it in range(nt):
        for i in range(nelement):
            if nextcell:
                index, val, err = _get_cdgs_element(fic, ne_read)
                if index is None:
                    end_data = True
                else:
                    index -= 1

                nextcell = False
            if i == index:
                nextcell = not end_data
                val_array[i, :] = val[:]
                err_array[i, :] = err[:]
            else:
                val_array[i, :] = zero[:]
                err_array[i, :] = zero[:]

        val_array = np.transpose(val_array.reshape((nx1, nx2, nx3, ne)))
        err_array = np.transpose(err_array.reshape((nx1, nx2, nx3, ne)))
        data[it, :, :, :, :, 0] = val_array[:, :, :, :]
        data[it, :, :, :, :, 1] = err_array[:, :, :, :]
    return data
