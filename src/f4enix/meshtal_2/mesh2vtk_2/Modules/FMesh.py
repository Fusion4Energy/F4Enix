import numpy
import vtk

from .functions.utils import ExtraBin
from .functions.vtk_functions import rectilinear_grid, structured_grid


class MeshData:
    """class storing mesh and data information (value + error) from MC simulation"""

    def __init__(
        self,
        geom: str,
        x1bin: numpy.ndarray,
        x2bin: numpy.ndarray,
        x3bin: numpy.ndarray,
        ebin: ExtraBin = None,
        tbin: ExtraBin = None,
        data: numpy.ndarray = None,
    ):
        if geom not in ("rec", "cyl"):
            raise ("bad geometry label")
        self._geom = geom
        self._x1bin = x1bin
        self._x2bin = x2bin
        self._x3bin = x3bin

        self._nx1 = len(x1bin) - 1
        self._nx2 = len(x2bin) - 1
        self._nx3 = len(x3bin) - 1

        if ebin is None:
            self._ebin = ExtraBin(numpy.array((0, 100)), "erg")
        else:
            self._ebin = ebin
        self._ne = len(self._ebin) - 1 if self._ebin.binbound else len(self._ebin)
        if self._ebin.totalbin:
            self._ne += 1

        if tbin is None:
            self._tbin = ExtraBin(numpy.array((0, 1e20)), "tme")
        else:
            self._tbin = tbin
        self._nt = len(self._tbin) - 1 if self._tbin.binbound else len(self._tbin)
        if self._tbin.totalbin:
            self._nt += 1

        shape = (self._nt, self._ne, self._nx3, self._nx2, self._nx1, 2)
        self._data = numpy.ndarray(shape)

        if data is not None:
            if self._data.shape == data.shape:
                self._data = data[:, :, :, :, :, :]
            else:
                raise TypeError("Data dimension doesn't match expected dimension")

    def __add__(self, fmesh):
        return add_mesh(self, fmesh)

    def __mul__(self, factor):
        if not isinstance(factor, (int, float)):
            raise TypeError(f"MeshData can be multipled only by float or integer")
        return scale_mesh(self, factor)

    def __rmul__(self, factor):
        if not isinstance(factor, (int, float)):
            raise TypeError(f"MeshData can be multipled only by float or integer")
        return scale_mesh(self, factor)

    @property
    def geom(self) -> str:
        """
        Mesh geometry ('rec', 'cyl')
        """
        return self._geom

    @property
    def x1bin(self) -> numpy.ndarray:
        """
        X/R bins data (rec or cyl geometry)
        """
        return self._x1bin

    @property
    def x2bin(self) -> numpy.ndarray:
        """
        Y/Z bins data (rec or cyl geometry)
        """
        return self._x2bin

    @property
    def x3bin(self) -> numpy.ndarray:
        """
        Z/T bins data (rec or cyl geometry)
        """
        return self._x3bin

    @property
    def ebin(self) -> ExtraBin:
        """
        mcnp energy bin (or user bin assingned to mcnp energy bin)
        """
        return self._ebin

    @property
    def tbin(self) -> ExtraBin:
        """
        mcnp timebin (or user bin assingned to mcnp time bin)
        """
        return self._tbin

    @property
    def nx1(self) -> int:
        """
        number of intervals on X/R axis
        """
        return self._nx1

    @property
    def nx2(self) -> int:
        """
        number of intervals on Y/Z axis
        """
        return self._nx2

    @property
    def nx3(self) -> int:
        """
        number of intervals on Z/T axis
        """
        return self._nx3

    @property
    def ne(self) -> int:
        """
        number of energy intervals (or energy values if discrete) on energy grid (or defined user bin)
        """
        return self._ne

    @property
    def nt(self) -> int:
        """
        number of time intervals (or time values if discrete) on time grid (or defined user bin)
        """
        return self._nt

    @property
    def data(self) -> numpy.ndarray:
        """
        array with mesh values and statistical error
        """
        return self._data

    def add(self, fmesh):
        """add mesh value to the current mesh and return a new MeshData object"""
        return add_mesh(self, fmesh)

    def scale(self, factor: float):
        """multiply mesh values by the factor and return a new MeshData object"""
        return scale_mesh(self, factor)

    def get_etbin_data(self, ebin: int = None, tbin: int = None):
        """Return a sub mesh with only data of selected energy and time bins"""

        if ebin is None:
            ebin = 0
            new_earray = self._ebin
        else:
            if ebin < 0:
                if -ebin < self._ne + 1:
                    ebin += self._ne
                else:
                    raise TypeError("get_etbin_data: energy bin out of range")

            if ebin < self._ne:
                if self._ebin.totalbin and ebin == self._ne - 1:
                    if self._ebin.binbound:
                        new_earray = ExtraBin(
                            numpy.array((self._ebin[0], self._ebin[-1])),
                            self._ebin.type,
                        )
                    else:
                        new_earray = ExtraBin(numpy.array((-1,)), self._ebin.type)
                else:
                    if self._ebin.binbound:
                        new_earray = ExtraBin(
                            self._ebin[ebin : ebin + 2], self._ebin.type
                        )
                    else:
                        new_earray = ExtraBin(
                            self._ebin[ebin : ebin + 1], self._ebin.type
                        )
            else:
                raise TypeError("get_etbin_data: energy bin out of range")

        if tbin is None:
            tbin = 0
            new_tarray = self._tbin
        else:
            if tbin < 0:
                if -tbin < self._nt + 1:
                    tbin += self._nt
                else:
                    raise TypeError("get_etbin_data: time bin out of range")

            if tbin < self._nt:
                if self._tbin.totalbin and tbin == self._nt - 1:
                    if self._tbin.binbound:
                        new_tarray = ExtraBin(
                            numpy.array((self._tbin[0], self._tbin[-1])),
                            self._tbin.type,
                        )
                    else:
                        new_tarray = ExtraBin(numpy.array((-1,)), self._tbin.type)
                else:
                    if self._tbin.binbound:
                        new_tarray = ExtraBin(
                            self._tbin[tbin : tbin + 2], self._tbin.type
                        )
                    else:
                        new_tarray = ExtraBin(
                            self._tbin[tbin : tbin + 1], self._tbin.type
                        )
            else:
                raise TypeError("get_etbin_data: time bin out of range")

        bin_mesh = MeshData(
            self._geom, self._x1bin, self._x2bin, self._x3bin, new_earray, new_tarray
        )

        if bin_mesh.ne == 1 and bin_mesh.nt == 1:
            bin_mesh.data[0, 0, :, :, :, :] = self._data[tbin, ebin, :, :, :, :]
        elif bin_mesh.ne == 1:
            bin_mesh.data[:, 0, :, :, :, :] = self._data[:, ebin, :, :, :, :]
        elif bin_mesh.nt == 1:
            bin_mesh.data[0, :, :, :, :, :] = self._data[tbin, :, :, :, :, :]
        else:
            bin_mesh.data[:, :, :, :, :, :] = self._data[:, :, :, :, :, :]

        return bin_mesh


class FMesh(MeshData):
    """class storing all kind of mesh Fmesh, CuV, CDGS"""

    def __init__(self, mesh: MeshData, meshLabel: str, trsf=None):
        super().__init__(
            mesh.geom,
            mesh.x1bin,
            mesh.x2bin,
            mesh.x3bin,
            mesh.ebin,
            mesh.tbin,
            mesh.data,
        )

        self._trsf = trsf
        self._tally = meshLabel
        self._type = None
        self._comments = None
        self._particle = None

    def __add__(self, mesh2):
        return add_mesh(self, mesh2)

    def __mul__(self, factor):
        if not isinstance(factor, (int, float)):
            raise TypeError(f"Fmesh can be multipled only by float or integer")
        return scale_mesh(self, factor)

    def __rmul__(self, factor):
        if not isinstance(factor, (int, float)):
            raise TypeError(f"Fmesh can be multipled only by float or integer")
        return scale_mesh(self, factor)

    @property
    def trsf(self):
        """
        Mesh transformation
        """
        return self._trsf

    @property
    def tally(self) -> str:
        """
        Mesh label
        """
        return self._tally

    @property
    def type(self) -> str:
        """
        Mesh type (meshtal, CDGS, CUV, ...)
        """
        return self._type

    @type.setter
    def type(self, value: str):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"type should be string type not {type(value)}")
        self._type = value

    @property
    def comments(self) -> str:
        """
        Mesh comments
        """
        return self._comments

    @comments.setter
    def comments(self, value: str):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"comments should be string type not {type(value)}")
        self._comments = value

    @property
    def particle(self) -> str:
        """
        Transported particle
        """
        return self._particle

    @particle.setter
    def particle(self, value: str):
        if not isinstance(value, (str, type(None))):
            raise TypeError(f"particle should be string type not {type(value)}")
        self._particle = value

    def print_info(self) -> dict:
        info = {
            "tally": self.tally,
            "type": self.type,
            "comments": self.comments,
            "particle": self.particle,
            "x1bin": self.x1bin,
            "x2bin": self.x1bin,
            "x3bin": self.x1bin,
            "ebin": self.ebin,
            "tbin": self.tbin,
            "geometry": self.geom,
        }
        return info

    def convert2tally(self):
        raise NotImplementedError()

    def write_vtk(self, filename: str, binary: bool = False, binlabels=None) -> None:
        """Write Fmesh in vtk file"""

        if self._geom == "rec" and self._trsf is None:
            gW = vtk.vtkXMLRectilinearGridWriter()
            if binlabels:
                gd = rectilinear_grid(self, labels=binlabels)
            else:
                gd = rectilinear_grid(self)
            ext = ".vtr"
        else:
            gW = vtk.vtkXMLStructuredGridWriter()
            if binlabels:
                gd = structured_grid(self, trsf=self._trsf, labels=binlabels)
            else:
                gd = structured_grid(self, trsf=self._trsf)
            ext = ".vts"

        filename += ext
        gW.SetFileName(filename)
        if not binary:
            gW.SetDataModeToAscii()

        gW.SetInputData(gd)
        gW.Write()

    def sameMesh(self, other_mesh: MeshData) -> bool:
        if self.geom != other_mesh.geom:
            return False
        if self.data.shape != other_mesh.data.shape:
            return False

        mesh1_xbins = (self.x1bin, self.x2bin, self.x3bin)
        mesh2_xbins = (other_mesh.x1bin, other_mesh.x2bin, other_mesh.x3bin)

        for bin1, bin2 in zip(mesh1_xbins, mesh2_xbins):
            if numpy.any(abs(bin1 - bin2) > 1e-12):
                return False

        return True


def add_mesh(mesh1: MeshData, mesh2: MeshData) -> MeshData:
    if not same_mesh(mesh1, mesh2):
        print("not same mesh structure, cannot add values")
        return None

    else:
        fsum = MeshData(
            mesh1.geom, mesh1.x1bin, mesh1.x2bin, mesh1.x3bin, mesh1.ebin, mesh1.tbin
        )
        fsum.data[:, :, :, :, :, 0] = (
            mesh1.data[:, :, :, :, :, 0] + mesh2.data[:, :, :, :, :, 0]
        )
        sigma1 = mesh1.data[:, :, :, :, :, 0] * mesh1.data[:, :, :, :, :, 1]
        sigma2 = mesh2.data[:, :, :, :, :, 0] * mesh2.data[:, :, :, :, :, 1]
        val = fsum.data[:, :, :, :, :, 0]
        err = numpy.sqrt(sigma1 * sigma1 + sigma2 * sigma2)
        fsum.data[:, :, :, :, :, 1] = numpy.divide(err, val, where=val != 0)
    return fsum


def scale_mesh(mesh: MeshData, factor: float) -> MeshData:
    fscale = MeshData(
        mesh.geom, mesh.x1bin, mesh.x2bin, mesh.x3bin, mesh.ebin, mesh.tbin
    )
    fscale.data[:, :, :, :, :, :] = mesh.data[:, :, :, :, :, :]
    fscale.data[:, :, :, :, :, 0] = fscale.data[:, :, :, :, :, 0] * factor
    return fscale
