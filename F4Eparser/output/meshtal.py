import numpy as np
import vtk
import os
import time
from scipy.spatial.transform import Rotation as R
from io import open


# convert character to float
# able to convert no standard fortran exponent type 1.234-123
def dfloat(c: str) -> float:
    try:
        x = float(c)
    except:
        c = c.strip()
        m = float(c[0:-4])
        e = int(c[-4:])
        x = m * pow(10, e)
    return x


# return a list of each nth charaters of a string
def splitn(string: str, n: int) -> list[str]:
    return [string[i: i + n].strip() for i in range(0, len(string), n)]


# Salta lineas (definida en cien sitios)
def skipLines(f, n):
    for iskip in range(n):
        f.readline()
    return


# # read cell list from file for CUV filtering option
# def get_clist(fn: os.PathLike) -> list(int):
#     f = open(fn, "rt")
#     lines = f.readlines()
#     f.close()
#     clist = []
#     filter_mode = "accept"
#     for lin in lines:
#         line = lin.split("#")[0]

#         sc_num = line.count(":")
#         line = line.replace(":", " : ", sc_num)

#         val = line.split()
#         if "reject" in line:
#             filter_mode = "reject"
#             val.remove("reject")

#         scidx = 0
#         for i in range(sc_num):
#             sc = val.index(":", scidx)
#             low = int(val[sc - 1])
#             high = int(val[sc + 1])
#             del val[sc - 1: sc + 2]
#             scidx = sc
#             clist.extend(range(low, high + 1))

#         clist.extend(map(int, val))

#     clist.insert(0, filter_mode)
#     return clist


class Fmesh:
    dtype = np.float64
    cvarsCart = ("Z", "Y", "X")
    cvarsCyl = ("Theta", "Z", "R")
    IPT = ("neutron", "photon", "electron")

    # Reads from opened file
    def __init__(self, mshtl) -> None:
        """Fmesh object. It needs a Meshtal object to be initialized.

        Parameters
        ----------
        mshtl : Meshtal
            meshtal object where the Fmesh is located
        """
        self.filled = False
        self.ntally = -1
        self.part = "unknown"
        self.etag = "energy"
        self.type = "flux"
        self.dims = list()
        self.scaleFac = 1.0
        self.pos = 0
        self.startm = 0
        self.origin = np.zeros([4], self.dtype)
        self.meshtal = mshtl
        self.__format__ = "undef"
        self.comment = None
        self.dosecom = False
        self.usrbin = False
        self.readHead = False
        self.vec = None

        # Rotation matrix
        self.rotation = np.identity(3, self.dtype)
        return

    def __readMeshCom__(self, f) -> None:
        line = f.readline()
        # read comment part
        if line[0:5] == "     ":
            self.comment = line
            line = f.readline()
            while line[0:5] == "     ":
                self.comment += line
                line = f.readline()
            self.comment = self.comment.strip("\n")

        # read particle type
        for p in ("neutron", "photon", "electron"):
            if p in line:
                self.part = p
                break

        if "importance mesh" in line:
            self.type = "srcimp"
        elif "cell-under-voxel" in line:
            self.type = "cuv"

        # Read MCNP mesh comments
        line = f.readline()
        while "Tally bin boundaries:" not in line:
            if "response function" in line:
                self.dosecom = True
            if "Energy bin" in line:
                if "energy" not in line[10:-1]:
                    self.usrbin = True
            if "source tally" in line:
                self.type = "source"
            if "Tally target:" in line:
                self.target = int(line.split()[-1])
                self.__format__ = "impmlt"
            line = f.readline()

        if self.type == "srcimp" and self.__format__ == "undef":
            self.__format__ == "impij"

        if self.type == "cuv":
            self.__format__ = "cuv"

    def __readMeshDim__(self, f) -> None:
        mcnp5_Cyl = "  Cylinder origin at"
        mcnp6_Cyl = "               origin at"

        # Test if cylinder
        self.cart = True
        line = f.readline()
        if line.startswith(mcnp5_Cyl) or line.startswith(mcnp6_Cyl):
            ishft = 0
            if line.startswith(mcnp6_Cyl):
                ishft = 1
            vals = line.split()
            self.cart = False
            # zval may end with comma
            zval = vals[5 - ishft]
            if not zval[-1].isdigit():
                zval = zval[:-1]
            self.origin = np.array(
                (0.0, zval, vals[4 - ishft], vals[3 - ishft]), self.dtype
            )
            self.axis = np.array(
                (vals[8 - ishft], vals[9 - ishft], vals[10 - ishft]),
                self.dtype
            )
            if line.startswith(mcnp6_Cyl):
                self.vec = np.array((vals[13], vals[14], vals[15]), self.dtype)
            line = f.readline()
        # Read dimensions  including energy
        # (E,Z,Y,X) order in storage in cartesian
        # (E,Th,Z,R) order in cylindrical
        for n in range(3):
            i = line.index(":")
            self.dims.insert(0, np.array(line[i + 1:].split(), self.dtype))
            line = f.readline()
        self.ldims = [len(x) - 1 for x in self.dims]

        # process energy bins
        i = line.index(":")
        if self.usrbin:
            # decay time bins
            if "times" in line:
                self.etag = "times"
                self.dims.insert(0, np.array(line[i + 1:].split(), self.dtype))
                self.ener = self.dims[0]  # shallow copy
                self.ldims.insert(0, len(self.dims[0]))
            # source importance bins
            elif "Target tally" in line:
                self.etag = "tally"
                line = f.readline()
                i = line.index(":")
                self.dims.insert(0, np.array(line[i + 1:].split(), self.dtype))
                self.ener = self.dims[0]  # shallow copy
                self.ldims.insert(0, len(self.dims[0]))
            # cell or isotopes bin
            else:
                self.etag = line[0:i].split()[0].lower()
                if self.etag == "cells":
                    binlen = 12
                else:
                    binlen = 10
                binlist = splitn(line[i + 1:], binlen)[:-1]
                lastbin = binlist[-1]
                ib = lastbin.index(".")
                nbin = int(lastbin[ib + 1:])
                self.dims.insert(0, np.array(range(nbin), self.dtype))
                #        self.ener = np.array(line[i+1:].split(),self.dtype)
                self.ener = np.array(binlist, self.dtype)
                self.ldims.insert(0, nbin + 1)
        # energy bin
        else:
            if "number:" in line:
                i = line.index(":")
                ne = int(line[i + 1:])
                line = f.readline()
                i = line.index(":")
                if self.type not in ["srcimp", "cuv"]:
                    # check if onef or mltf
                    nelemts = (
                        self.ldims[0] * self.ldims[1] * self.ldims[2]
                    )  # dim of energy bin not yet inserted in ldims
                    # 6 columns 12 char each
                    # 2 block2 (flux, error)
                    nline = int(ne / 6)
                    nr = ne % 6
                    if nr != 0:
                        nline += 1
                    nblck = 2 * ((12 * ne) + nline)
                    pos = f.tell()
                    skipLines(f, 5)
                    onef = checkonef(f, nblck)
                    f.seek(pos)
                    if onef:
                        self.__format__ = "onef"
                    else:
                        self.__format__ = "mltf"

            self.dims.insert(0, np.array(line[i + 1:].split(), self.dtype))
            self.ener = self.dims[0]  # shallow copy
            nbin = 1
            if len(self.dims[0]) > 2:
                nbin = len(self.dims[0])
            self.ldims.insert(0, nbin)

        if self.__format__ == "undef":
            self.__format__ = "mcnp"
        self.startm = f.tell()

    def readMCNP(self, f) -> None:
        """Read the MCNP file.

        Parameters
        ----------
        f : file stream
            file stream to be read
        """
        # Proceed to read MCNP
        # Check if file is readable
        # There is one blank between fmesh

        if self.readHead:
            f.seek(self.startm)
        else:
            line = f.readline()
            if line == "":
                print("read error")
            self.ntally = int(line[-6:-1])

            # read mesh comments header
            self.__readMeshCom__(f)

            # read mesh dimensions and bins
            self.__readMeshDim__(f)

        line = f.readline()
        # Coordinates are origin-independent
        if self.cart:
            for i in range(1, 4):
                self.origin[i] = self.dims[i][0]
                self.dims[i] -= self.origin[i]
        # Memory storage for cell data
        # Reversed indeces to keep ordering in memory
        rshape = self.ldims[:]
        self.dat = np.zeros(rshape, self.dtype)
        self.err = np.zeros(rshape, self.dtype)
        # Data ordering in file
        line = f.readline()
        matFormat = "Rel Error" not in line

        # select correct keyword depending on the user bin
        if self.etag == "energy":
            eword = "Energy"
        elif self.etag == "cells":
            eword = "Cell"
        elif self.etag == "times":
            eword = "Time"
        else:
            eword = "Nuclide"
        energyCol = eword in line

        # All meshtal orders have energy as the least moving index, so we will
        # read variable of one energy bin, and then transpose to proper
        # ordering transposing takes memory

        # Matrix format, what order?
        if matFormat:
            skipLines(f, 1)
            if self.etag == "tally":
                self.tallyscore = [float(f.readline().split()[-1])]
            c = f.readline().split()[0]
            if self.cart:
                i = self.cvarsCart.index(c)
            else:
                i = self.cvarsCyl.index(c)
            if i == 0:  # ij
                iord = (1, 2, 3)  # Read order of storage
                itrn = (0, 1, 2)  # transpose of local storage
            elif i == 1:  # ik
                iord = (2, 1, 3)
                itrn = (1, 0, 2)
            elif i == 2:  # jk
                iord = (3, 1, 2)
                itrn = (1, 2, 0)
            else:
                print("Matrix transpose error")
        else:  # colformat x,y,z
            iord = (3, 2, 1)  # Z,Y,X para reshape
            itrn = (2, 1, 0)

        # Temporal variables for transposing
        # reshape of temp
        rshape = np.array([self.ldims[iord[i]] for i in range(3)])

        # ================================
        if matFormat:  # MATRIX FORMAT
            hlines = 3  # header lines to skip in first
            hlinesNext = 4  # Following hlines

            # Temporal variables for transposing
            xdat = np.zeros(rshape, self.dtype)
            xerr = np.zeros(rshape, self.dtype)

            for ie in range(self.ldims[0]):
                if ie > 0:
                    skipLines(f, 3)
                    if self.etag == "tally":
                        self.tallyscore.append(float(f.readline().split()[-1]))

                for ix1 in range(rshape[0]):
                    skipLines(f, hlines)
                    hlines = hlinesNext  # subsequent reads
                    for ix0 in range(rshape[1]):
                        xdat[ix1, ix0, :] = list(
                            map(dfloat, f.readline().split()[1:])
                        )
                    # use dfloat conversion in case of no standar fortan
                    # exponent
                    skipLines(f, 3)
                    for ix0 in range(rshape[1]):
                        xerr[ix1, ix0, :] = f.readline().split()[
                            1:
                        ]  # automatic conversion
                    skipLines(f, 2)
                # Copy of temp data
                self.dat[ie, :, :, :] = np.transpose(xdat, itrn)
                self.err[ie, :, :, :] = np.transpose(xerr, itrn)

        else:  # colformat
            n = np.prod(rshape)
            xdat = np.zeros(n, self.dtype)
            xerr = np.zeros(n, self.dtype)
            colDat = 31
            colErr = 43
            if energyCol:  # first column is energy
                colDat += 10
                colErr += 10
            for ie in range(self.ldims[0]):  # energy
                for ix in range(n):
                    line = f.readline()
                    try:
                        xdat[ix] = dfloat(line[colDat: colDat + 12])
                        xerr[ix] = line[colErr: colErr + 12]
                    except:
                        xdat[ix] = dfloat(line[colDat + 1: colDat + 13])
                        xerr[ix] = line[colErr + 1: colErr + 13]
                self.dat[ie, :, :, :] = np.transpose(xdat.reshape(rshape),
                                                     itrn)
                self.err[ie, :, :, :] = np.transpose(xerr.reshape(rshape),
                                                     itrn)

        # Removal of thick angular bin if any
        # only if:
        # Angular domain is a whole revolution (there may be a thick bin in
        # between)
        if not self.cart and abs(self.dims[1][-1] - self.dims[1][0] - 1.0) < 0.01:
            maxDth = 5.0 / self.ldims[1]
            for i in range(self.ldims[1]):
                if self.dims[1][i + 1] - self.dims[1][i] > maxDth:
                    self.ldims[1] -= 1
                    np.delete(self.dat, i, 1)
                    np.delete(self.err, i, 1)
                    # Reorder remaining bins
                    if i > 0 and i < self.ldims[1] - 1:
                        self.dims[1] = np.concatenate(
                            (self.dims[1][i + 1: -1] - 1.0, self.dims[1][0:i]),
                             axis=0)
                        self.dat = np.concatenate(
                            (self.dat[:, i + 1:, :, :], self.dat[:, : i - 1, :, :]),
                            axis=1,
                        )
                        self.err = np.concatenate(
                            (self.err[:, i + 1:, :, :], self.err[:, : i - 1, :, :]),
                            axis=1,
                        )
                    else:
                        np.delete(self.dims[1], i, 0)
                    # Leave loop
                    break

        self.filled = True

    # Read photonfile format of SRCIMP mesh (D1SUNED)
    def readSRCTYPE(self, f, cfilter=None):
        mcnp5_Cyl = "  Cylinder origin at"
        mcnp6_Cyl = "               origin at"
        self.etag = "energy"

        if self.readHead:
            f.seek(self.startm)
        else:
            line = f.readline()
            if line == "":
                print("Read error")
            self.ntally = int(line[-6:-1])

            # read mesh comments header
            self.__readMeshCom__(f)

            # read mesh dimensions and bins
            # erint("readDims")  # This is not defined!
            self.__readMeshDim__(f)

        # Coordinates are origin-independent
        if self.cart:
            for i in range(1, 4):
                self.origin[i] = self.dims[i][0]
                self.dims[i] -= self.origin[i]
        # Memory storage for cell data
        # Reversed indeces to keep ordering in memory
        if self.type == "cuv":
            skipLines(f, 2)
        else:
            if self.type == "srcimp":
                self.tallyscore = np.float(f.readline().split()[-1])
            skipLines(f, 5)

        nelemts = self.ldims[1] * self.ldims[2] * self.ldims[3]
        nbin = int(self.ldims[0])

        if nbin > 1:
            ne = nbin - 1
        else:
            ne = nbin
        rshape = (nbin, nelemts)

        # 6 columns 12 char each
        # 2 block2 (flux, error)
        nline = int(ne / 6)
        nr = ne % 6
        if nr != 0:
            nline += 1

        # check if onef or mltf
        if self.__format__ == "undef":
            nblck = 2 * ((12 * ne) + nline)
            pos = f.tell()
            onef = checkonef(f, nblck)
            f.seek(pos)
            if onef:
                self.__format__ = "onef"
            else:
                self.__format__ = "mltf"

        xdat = np.zeros(rshape, self.dtype)
        xerr = np.zeros(rshape, self.dtype)

        # source importance data format
        if self.type == "srcimp":
            for k in range(nelemts):
                vals = []
                errs = []
                tot = f.readline().split()[-1]
                for i in range(nline):
                    vals.extend(f.readline().split())
                for i in range(nline):
                    errs.extend(f.readline().split())

                if nbin > 1:
                    xdat[nbin - 1, k] = dfloat(tot)
                    xerr[nbin - 1, k] = 0.0
                for i in range(ne):
                    xdat[i, k] = dfloat(vals[i])
                    xerr[i, k] = float(errs[i])

        elif self.type == "cuv":
            # one/multiflux data format
            if cfilter is not None:
                if cfilter[0] == "reject":
                    accept = False
                else:
                    accept = True
            for k in range(nelemts):
                lvals = f.readline().split()
                nc = int(lvals[-1])
                voxvol = float(lvals[-2])

                celdata = [voxvol, [], []]
                cell = []

                # read cell data info
                for ic in range(nc):
                    vline = f.readline().split()
                    cline = int(vline[0])
                    cell.append(cline)
                    if (
                        cfilter is None
                        or (accept and cline in cfilter)
                        or (not accept and cline not in cfilter)
                    ):
                        celdata[1].append(float(vline[1]))
                        if nbin > 1:
                            celdata[2].append([float(vline[2]),
                                               float(vline[3])])

                voxvals = []
                voxerrs = []
                # read cell values and errors
                for ind, ic in enumerate(cell):
                    valc = []
                    errc = []
                    if (
                        cfilter is None
                        or (accept and ic in cfilter)
                        or (not accept and ic not in cfilter)
                    ):
                        for i in range(nline):
                            valc.extend(f.readline().split())
                        for i in range(nline):
                            errc.extend(f.readline().split())

                        vals = np.array([dfloat(x) for x in valc], self.dtype)
                        errs = np.array([dfloat(x) for x in errc], self.dtype)

                        if nbin > 1:
                            vals = np.append(vals, celdata[2][ind][0])
                            errs = np.append(errs, celdata[2][ind][1])

                        voxvals.append(vals)
                        voxerrs.append(errs)
                    else:
                        skipLines(f, 2 * nline)

                if len(voxvals) == 0:
                    vals = [0] * nbin
                    errs = [0] * nbin
                else:
                    # sum all bin if not total bin
                    vals, errs = sumCellInVox(
                        voxvals, voxerrs, celdata, Vmult=self.__mltopt__
                    )

                xdat[:, k] = vals[:]
                xerr[:, k] = errs[:]

        else:
            if cfilter is not None:
                if cfilter[0] == "reject":
                    accept = False
                else:
                    accept = True
            # one/multiflux data format
            for k in range(nelemts):
                lvals = f.readline().split()
                nc = int(lvals[-1])
                voxvol = float(lvals[-2])

                celfrac = [voxvol, []]
                cell = []
                for ic in range(nc):
                    vline = f.readline().split()
                    cline = int(vline[0])
                    cell.append(cline)
                    if (
                        cfilter is None
                        or (accept and cline in cfilter)
                        or (not accept and cline not in cfilter)
                    ):
                        celfrac[1].append(float(vline[1]))

                voxvals = []
                voxerrs = []

                if self.__format__ == "onef":
                    valc = []
                    errc = []
                    for i in range(nline):
                        valc.extend(f.readline().split())
                    for i in range(nline):
                        errc.extend(f.readline().split())

                    vals = np.array([dfloat(x) for x in valc], self.dtype)
                    errs = np.array(errc, self.dtype)
                    voxvals.append(vals)
                    voxerrs.append(errs)
                else:
                    for ic in cell:
                        valc = []
                        errc = []
                        if (
                            cfilter is None
                            or (accept and cline in cfilter)
                            or (not accept and cline not in cfilter)
                        ):
                            for i in range(nline):
                                valc.extend(f.readline().split())
                            for i in range(nline):
                                errc.extend(f.readline().split())

                            vals = np.array([dfloat(x) for x in valc],
                                            self.dtype)
                            errs = np.array(errc, self.dtype)
                            voxvals.append(vals)
                            voxerrs.append(errs)
                        else:
                            skipLines(f, 2 * nline)

                if len(voxvals) == 0:
                    vals = [0] * nbin
                    errs = [0] * nbin
                    tot = 0
                else:
                    # sum all bin if not total bin
                    vals, errs = sumElements(
                        voxvals,
                        voxerrs,
                        celfrac,
                        Vsum=self.__format__,
                        Vmult=self.__mltopt__,
                        corr=True,
                    )

                xdat[:, k] = vals[:]
                xerr[:, k] = errs[:]

        rshape = np.array([nbin, self.ldims[3], self.ldims[2], self.ldims[1]])
        self.dat = xdat.reshape(rshape)
        self.dat = self.dat.transpose(0, 3, 2, 1)
        self.err = xerr.reshape(rshape)
        self.err = self.err.transpose(0, 3, 2, 1)
        self.filled = True
        return

    #  end modifs

    def print_info(self) -> None:
        """print information related to the fmesh
        """
        if self.__format__ == "mltf":
            meshtype = "multiflux"
        else:
            meshtype = self.type

        if self.cart:
            geom = "rectangular"
        else:
            geom = "cylindrical"

        print(" Tally          : {}".format(self.ntally))
        if self.comment is not None:
            print(" Comments       : ")
            print(" {}".format(self.comment))
        print(" Particle       : {}".format(self.part))
        if self.__format__ != "mltf":
            print(" Mesh type      : {}".format(meshtype))
        else:
            print(
                " Mesh type      : {}, normalization {}".format(
                    meshtype, self.__mltopt__
                )
            )
        print(" Dose modif     : {}".format(self.dosecom))
        print(" Mesh geometry  : {}".format(geom))
        print(" Mesh origin    : {org[3]} {org[2]} {org[1]}".format(org=self.origin))
        if self.cart:
            print(
                " X dimensions   :{}".format(
                    self.__format_XYZ_Dim__(self.dims[3] + self.origin[3])
                )
            )
            print(
                " Y dimensions   :{}".format(
                    self.__format_XYZ_Dim__(self.dims[2] + self.origin[2])
                )
            )
            print(
                " Z dimensions   :{}".format(
                    self.__format_XYZ_Dim__(self.dims[1] + self.origin[1])
                )
            )
        else:
            print(" R dimensions   :{}".format(self.__format_XYZ_Dim__(self.dims[3])))
            print(" Z dimensions   :{}".format(self.__format_XYZ_Dim__(self.dims[2])))
            print(" T dimensions   :{}".format(self.__format_XYZ_Dim__(self.dims[1])))
        print(" Energy bins    :")
        self.print_EbinRange()

    def __format_XYZ_Dim__(self, vec: list, nval: int = 6) -> str:
        def format_XYZ_Dim_inter(vec, nval=6):
            interfound = False
            xtab = []
            dx0 = 0
            x1 = vec[0]
            i = 0
            for x2 in vec[1:]:
                dx = round(x2 - x1, 2)
                if dx == dx0:
                    i += 1
                    interfound = True
                else:
                    xtab.append([i, x1])
                    dx0 = dx
                    i = 0
                x1 = x2
            xtab.append([i, x2])

            if not interfound:
                return format_XYZ_Dim_long(vec)
            line = ""
            newline = "  {:10.3e}".format(xtab[0][1])

            for i, v in enumerate(xtab[1:]):
                newline += " {val[0]:4d}I {val[1]:10.3e}".format(val=v)
                if (i + 1) % nval == 0:
                    line += newline + "\n"
                    newline = "                 "

            if (i + 1) % nval != 0:
                line += newline
            else:
                line = line[:-1]
            return line

        def format_XYZ_Dim_long(vec: list, nval: int = 8) -> str:
            line = ""
            newline = ""
            for i, v in enumerate(vec):
                newline += " {:10.3e}".format(v)
                if (i + 1) % nval == 0:
                    line += newline + "\n"
                    newline = "                 "

            if (i + 1) % nval != 0:
                line += newline
            else:
                line = line[:-1]
            return line

        if len(vec) <= nval:
            return format_XYZ_Dim_long(vec)
        else:
            return format_XYZ_Dim_inter(vec)

    def print_EbinRange(self) -> None:
        """print the energy range bins
        """
        print("         flag : {}".format(self.etag))
        print("    bin index : bin range  ")
        if self.etag == "energy":
            nb = 0
            if len(self.ener) > 2:
                nb = len(self.ener) - 1
                for i in range(nb):
                    print(
                        "       {:4d}   :      {} - {}  MeV".format(
                            i, self.ener[i], self.ener[i + 1]
                        )
                    )

        elif self.etag in ["cells", "daughter", "parent"]:
            binlist = []
            ibin = 0
            for b in self.ener:
                decpart, intpart = np.modf(b)
                binnum = int(np.rint(abs(decpart) * 1000)) - 1
                if binnum == ibin:
                    binlist.append([int(intpart)])
                    ibin += 1
                else:
                    binlist[binnum].append(int(intpart))

            ib = 0
            print("       {:4d}   :  Other ".format(ib))
            for b in binlist[1:]:
                ib += 1
                line = "       {:4d}   :".format(ib)
                rang = 0
                for i, v in enumerate(b):
                    if i < len(b) - 1:
                        if b[i + 1] < 0:
                            rang = 1
                    if rang == 0:
                        line = line + " {},".format(v)
                    elif rang == 1:
                        line = line + " {}-".format(v)
                        rang = 2
                    else:
                        line = line + "{},".format(abs(v))
                        rang = 0
                print(line[:-1])

            ib += 1
            print("       {:4d}   :  Total ".format(ib))

        elif self.etag in ["times", "tally"]:
            for i, b in enumerate(self.ener):
                print("       {:4d}   :      {}".format(i, b))
        return

    # Checks whether it is the same mesh
    def sameMesh(self, xelf, checkErg: bool = False) -> bool:
        i1 = 1
        if checkErg:
            i1 = 0
            if self.etag != xelf.etag:
                return False

        if self.cart != xelf.cart:
            return False
        if self.ldims[i1:4] != xelf.ldims[i1:4]:
            return False
        if np.linalg.norm(self.origin - xelf.origin) > 0.001:
            return False
        if any(
            (np.linalg.norm(self.dims[i] - xelf.dims[i]) > 0.01 for i in range(i1, 4))
        ):
            return False
        if not self.cart:
            if np.linalg.norm(self.axis - xelf.axis) > 0.001:
                return False
            if self.vec is not None:
                if np.linalg.norm(self.vec - xelf.vec) > 0.001:
                    return False

        return True

    # Translate the mesh
    def translate(self, vec: list[float]) -> None:
        """translate the fmesh using a vector x, y, z.

        Parameters
        ----------
        vec : list[float]
            translation vector
        """
        self.origin += np.array(vec)

    # Salida en formato VTK
    # Escribe structured grid
    def getVTKsg(self) -> vtk.vtkStructuredGrid:
        """get the VTK structured grid format for the mesh

        Returns
        -------
        vtk.vtkStructuredGrid
            generated vtk object
        """
        import math

        kdims = self.dims
        pts = vtk.vtkPoints()
        pts.SetDataTypeToFloat()

        # Paso de celdas degeneradas en el centro
        # Nunca se considera el centro del cilindro
        if not self.cart and kdims[3][0] == 0.0:
            kdims[3][0] = kdims[3][1] * 0.01

        # Add an angular bin if first bin angle >0.5
        if not self.cart and kdims[1][1] > 0.5:
            self.ldims[1] += 1
            kdims[1] = np.insert(kdims[1], 1, 0.5, axis=0)
            self.dat = np.insert(self.dat, 0, self.dat[:, 0, :, :], axis=1)
            self.err = np.insert(self.err, 0, self.err[:, 0, :, :], axis=1)

        ptDims = np.array([len(kdims[3]), len(kdims[2]), len(kdims[1])])
        pts.SetNumberOfPoints(np.prod(ptDims))

        ps = []
        k = 0
        dpt = self.origin[:0:-1]  # origin en reversa
        for z in kdims[1]:
            if not self.cart:
                ang = 2.0 * math.pi * z
                sz = math.sin(ang)
                cz = math.cos(ang)
            for y in kdims[2]:
                for x in kdims[3]:
                    if self.cart:
                        pt = np.array((x, y, z))
                    else:
                        pt = np.array((x * cz, x * sz, y))  # por ahora axs = 0 0 1
                    ps.append(pt + dpt)
                    k += 1

        # cylindrical rotation
        ps -= dpt
        rotAxis = np.cross(self.axis, [0.0, 0.0, 1.0])
        if np.linalg.norm(rotAxis) == 0:
            rotM = R.from_rotvec([0, 0, 0])
        else:
            rotAxis = rotAxis / np.linalg.norm(rotAxis)
            ma = np.linalg.norm(self.axis)
            ang = -np.arccos(np.dot(self.axis, [0.0, 0.0, 1.0]) / ma)
            rotM = R.from_rotvec(rotAxis * ang)
        ps = rotM.apply(ps)
        ps += dpt
        for p in range(len(ps)):
            pts.InsertPoint(p, ps[p])

        # Puntos definidos
        sg = vtk.vtkStructuredGrid()
        sg.SetDimensions(*ptDims)

        sg.SetPoints(pts)
        sgcd = sg.GetCellData()
        # These slicing operations do not copy data
        # mySlice = [slice(-1,-2,-1),slice(None),slice(None),slice(None)]

        it = 0
        if self.etag not in ["times", "tally"]:
            # Dataset energia total
            it = 1
            sgcd.AddArray(
                makeVTKarray(self.dat[-1, :, :, :], "Value - Total", self.scaleFac)
            )
            sgcd.AddArray(makeVTKarray(self.err[-1, :, :, :], "Error - Total"))
        # Dataset other bins
        for ie in range(self.ldims[0] - it):
            sgcd.AddArray(
                makeVTKarray(
                    self.dat[ie, :, :, :], "ValueBin-{0:03d}".format(ie), self.scaleFac
                )
            )
            sgcd.AddArray(
                makeVTKarray(self.err[ie, :, :, :], "ErrorBin-{0:03d}".format(ie))
            )
        # Include weight windows in VTK file (now only one group)
        return sg

    def _writeVTKsg(self, ofn: os.PathLike) -> None:
        """Write the mesh to a vtk file

        Parameters
        ----------
        ofn : os.PathLike
            outpath for the file
        """

        # Escritura en disco
        off = vtk.vtkXMLStructuredGridWriter()
        off.SetFileName(ofn)
        # ASCII o binario (con o sin compresion)
        off.SetDataModeToAscii()
        # off.SetDataModeToBinary()
        # off.SetCompressorTypeToZLib()
        # Esto cambia con la version de VTK
        t = self.getVTKsg()
        self.meshtal.setVTKparams(t)
        if vtk.vtkVersion().GetVTKMajorVersion() >= 6:
            off.SetInputData(t)
        else:
            off.SetInput(t)
        off.Write()

    # Escribe rectilinear grid
    # No vale para cilindricas o rotadas
    def getVTKrg(self) -> vtk.vtkRectilinearGrid:
        """Get a rectilinear grid. This cannot be used for cylindrical or
        rotated grids

        Returns
        -------
        vtk.vtkRectilinearGrid
            generated rectilinear grid object
        """
        if not self.cart:
            print("Cylindrical meshtal cannot be plotted to RectangularGrid")

        if np.any(self.rotation != np.identity(3)):
            print("Rotated meshtal cannot be plotted to RectangularGrid")
            print("... but plotting anyway (no rotation)")

        xa = makeVTKarray(self.dims[3] + self.origin[3], "X (cm)")
        ya = makeVTKarray(self.dims[2] + self.origin[2], "Y (cm)")
        za = makeVTKarray(self.dims[1] + self.origin[1], "Z (cm)")

        rg = vtk.vtkRectilinearGrid()
        rg.SetDimensions(self.ldims[3] + 1, self.ldims[2] + 1,
                         self.ldims[1] + 1)
        rg.SetXCoordinates(xa)
        rg.SetYCoordinates(ya)
        rg.SetZCoordinates(za)
        rgcd = rg.GetCellData()
        it = 0
        if self.etag not in ["times", "tally"]:
            # Dataset energia total
            it = 1
            rgcd.AddArray(
                makeVTKarray(self.dat[-1, :, :, :], "Value - Total",
                             self.scaleFac)
            )
            rgcd.AddArray(makeVTKarray(self.err[-1, :, :, :], "Error - Total"))
        # Dataset other bins
        for ie in range(self.ldims[0] - it):
            rgcd.AddArray(
                makeVTKarray(
                    self.dat[ie, :, :, :], "ValueBin-{0:03d}".format(ie),
                    self.scaleFac
                )
            )
            rgcd.AddArray(
                makeVTKarray(self.err[ie, :, :, :],
                             "ErrorBin-{0:03d}".format(ie))
            )
        return rg

    def _writeVTKrg(self, ofn: os.PathLike) -> None:
        """Write the rectilinear grid

        Parameters
        ----------
        ofn : os.PathLike
            output vtk file path
        """
        # Escritura en disco
        off = vtk.vtkXMLRectilinearGridWriter()
        off.SetFileName(ofn)
        off.SetDataModeToAscii()
        # Esto cambia con la version de VTK
        t = self.getVTKrg()
        self.meshtal.setVTKparams(t)
        if vtk.vtkVersion().GetVTKMajorVersion() >= 6:
            off.SetInputData(t)
        else:
            off.SetInput(t)
        off.Write()

    def writeVTK(self, ofn: os.PathLike) -> None:
        """Write the fmesh to a vtk file

        Parameters
        ----------
        ofn : os.PathLike
            path to the vtk outfile
        """
        if self.cart:
            self._writeVTKrg(ofn)
        else:
            self._writeVTKsg(ofn)


class Meshtal:
    def __init__(self, fn: os.PathLike, filetype: str = "MCNP") -> None:
        """Class representing a parsed Meshtal file

        Parameters
        ----------
        fn : os.PathLike
            path to the meshtal file to be read
        filetype : str, optional
            type of file, by default "MCNP" which is the only file type
            implemented at the moment.
        """
        self.filename = fn
        # Parametros para VTK
        self.filetype = filetype
        self.f = open(fn, "rt")
        self.__readHeadMCNP__()
        self.mesh = self.__scanMCNP__()

        self.params = dict()
        self.params["creationTime"] = time.asctime()
        if "USER" in os.environ:
            self.params["author"] = os.environ["USER"]
        elif "USERNAME" in os.environ:
            self.params["author"] = os.environ["USERNAME"]
        if "HOSTNAME" in os.environ:
            self.params["host"] = os.environ["HOSTNAME"]
        elif "COMPUTERNAME" in os.environ:
            self.params["host"] = os.environ["COMPUTERNAME"][2:]
        self.params["path"] = os.getcwd()
        return

    def __scanMCNP__(self) -> dict[int, Fmesh]:
        mcnp5_Cyl = "  Cylinder origin at"
        mcnp6_Cyl = "               origin at"
        mesh = {}
        while True:
            line = self.f.readline()
            if line == "":
                break
            if "Mesh Tally Number" in line:
                t = Fmesh(self)
                ntally = int(line.split()[-1])
                t.ntally = ntally
                t.pos = self.f.tell() - len(line)

                t.__readMeshCom__(self.f)
                t.__readMeshDim__(self.f)
                t.readHead = True

                mesh[ntally] = t
        return mesh

    def readMesh(self, mesh: int | list[int] = None,
                 cell_filters: list[int] = None) -> None:
        """Parse a list of FMESHes

        Parameters
        ----------
        mesh : int | list[int], optional
            list of FMESH to be read, by default None means that all meshes
            will be parsed
        cell_filters : list[int], optional
            list of cells to be used as filters for the CuV approach,
            by default None

        Raises
        ------
        KeyError
            if the filetype is is not implemented
        """
        if self.filetype == "MCNP":
            # cycle on all meshes to be read
            if mesh is None:
                mesh = self.mesh.keys()
            elif type(mesh) is int:
                mesh = [mesh]

            for meshid in mesh:
                m = self.mesh[meshid]
                # move to the mesh location in the file
                self.f.seek(m.pos)
                if m.__format__ == "mcnp":
                    m.readMCNP(self.f)
                else:
                    # if m.__mltflt__ is not None:
                    #     flist = get_clist(m.__mltflt__)
                    # else:
                    #     flist = None
                    m.readSRCTYPE(self.f, cfilter=cell_filters)
        else:
            raise KeyError("This file type is not implemented: " +
                           str(self.filetype))

    def __readHeadMCNP__(self) -> None:
        vals = self.f.readline().split()
        self.code = vals[0]
        self.version = vals[2]
        self.probid = vals[-2] + " " + vals[-1]
        self.title = self.f.readline().strip()
        self.nps = int(
            float((self.f.readline().split()[-1]))
        )  # nps: int doesnt like decimals

    def writeVTK(self, ofn: os.PathLike) -> None:
        """write a .vtk file for all read fmeshes

        Parameters
        ----------
        ofn : os.PathLike
            path to the outfile name
        """
        mb = vtk.vtkMultiBlockDataSet()
        i = 0
        for k, m in self.mesh.items():
            if not m.filled:
                continue
            if m.cart:
                mb.SetBlock(i, m.getVTKrg())
            else:
                mb.SetBlock(i, m.getVTKsg())
            i += 1
        # A escribir
        off = vtk.vtkXMLMultiBlockDataWriter()
        self.setVTKparams(mb)
        off.SetInputData(mb)
        off.SetFileName(ofn)
        off.Write()

    def addVTKparams(self, xdict: dict) -> None:
        """Add more parameters to include in VTK

        Parameters
        ----------
        xdict : dict
            parameters to be added
        """
        self.params.update(xdict)

    def setVTKparams(self, xdat: vtk.vtkDataSet) -> None:
        """Modify VTK structure for including parameters

        Parameters
        ----------
        xdat : vtk.vtkDataSet
            vtk data set to be modified
        """
        # for x,v in self.params.iteritems():
        for x, v in self.params.items():
            t = vtk.vtkStringArray()
            t.SetName(x)
            t.InsertNextValue(v)
            xdat.GetFieldData().AddArray(t)

    def __repr__(self) -> str:
        print(self)

    def __str__(self) -> str:
        out = ''
        out = out + "\n Meshtally file : {}".format(self.filename)
        # tlist = self.mesh.keys()
        tlist = list(self.mesh.keys())
        tlist.sort()
        for t in tlist:
            m = self.mesh[t]
            if m.__format__ == "mltf":
                meshtype = "multiflux"
            else:
                meshtype = m.type
            if m.comment is not None:
                firstlinecom = m.comment.split("\n")[0].strip()
                out = out + (
                    "   Tally {} : {}  {} mesh   '{}' ".format(
                        t, m.part, meshtype, firstlinecom
                    )
                )
            else:
                out = out + (
                    "   Tally {} : {}  {} mesh ".format(t, m.part, meshtype))
        return out

    def print_info(self) -> None:
        """print information on the Meshtal
        """
        print(self)


# ================= END OF CLASS DEFINITIONS ======================


# This should be the same as above
def makeVTKarray(nArr, aName, sc=1.0):
    from vtk.util import numpy_support

    if sc == 1.0:
        vArr = numpy_support.numpy_to_vtk(nArr.ravel(), deep=1)
    else:
        # deep copy needed in both cases. Seems the multiplication does not
        # create space
        vArr = numpy_support.numpy_to_vtk(nArr.ravel() * sc, deep=1)
    vArr.SetName(aName)
    return vArr


# check if flux mesh from mcnpacab is onef or mltf format
def checkonef(f, lblck):
    ncel = 0
    while True:
        vals = f.readline().split()
        pos = f.tell()
        nvals = len(vals)
        if ncel > 1 and nvals == 9:
            return True

        if nvals == 9:
            ncel = int(vals[-1])
            pos = pos + 19 * ncel + lblck
            f.seek(pos)
        else:
            if nvals == 0:
                return True
            else:
                return False


def sumCellInVox(voxvals, voxerrs, celfrac, Vmult="none", corr=True,
                 nulcount=True):
    # sum the value in the voxel multiplied by the volume fraction
    valf = list(map(lambda x, y: x * y, voxvals, celfrac[1]))

    if nulcount:
        volOK = sum(celfrac[1])
    else:
        volOK = 0
        for i, x in enumerate(valf):
            if x[-1] != 0:
                volOK += celfrac[1][i]

    vals = sum(valf)
    if not corr:
        err2 = sum(map(lambda x, y: x * x * y * y, valf, voxerrs))
        errs = np.divide(
            np.sqrt(err2), vals, out=np.zeros_like(err2), where=(vals != 0)
        )
    else:
        err = sum(map(lambda x, y: abs(x * y), valf, voxerrs))
        errs = np.divide(err, vals, out=np.zeros_like(err), where=(vals != 0))

    # Default units are X per volume units
    if Vmult == "vtot":
        # multiply by the total volume of the voxel

        vals = celfrac[0] * vals  # here units are integral X
    elif Vmult == "celf":
        # multiply by the voxel volume fraction
        vals = (
            vals / volOK
        )  # here units are X per volume units (volume of no zero voxel)

    return vals, errs


def sumElements(voxvals, voxerrs, celfrac, Vsum="onef", Vmult="none", corr=False):
    if Vsum == "onef":
        # sum all value in the voxel
        vals = voxvals[0]
        errs = voxerrs[0]

        # insert value and error total over all energy bins
        xs = sum(vals)
        if xs != 0:
            ve = vals * errs
            if not corr:
                err2 = sum(ve * ve)
                errt = np.sqrt(err2) / xs
            else:
                errt = sum(abs(ve)) / xs
            errs = np.append(voxerrs[0], errt)
        else:
            errs = np.append(voxerrs[0], 0)
        vals = np.append(vals, xs)
    else:
        # sum the value in the voxel multiplied by the volume fraction
        if len(voxvals) > 1:
            volOK = 0
            valf = map(lambda x, y: x * y, voxvals, celfrac[1])

            # insert value and error total over all energy bins
            for i, x in enumerate(valf):
                xs = sum(x)
                if xs != 0:
                    volOK += celfrac[1][i]
                    ve = x * voxerrs[i]
                    if not corr:
                        err2 = sum(ve * ve)
                        errt = np.sqrt(err2) / xs
                    else:
                        errt = sum(abs(ve)) / xs
                    voxerrs[i] = np.append(voxerrs[i], errt)
                else:
                    voxerrs[i] = np.append(voxerrs[i], 0)
                valf[i] = np.append(x, xs)

            vals = sum(valf)

            if not corr:
                err2 = sum(map(lambda x, y: x * x * y * y, valf, voxerrs))
                errs = np.divide(
                    np.sqrt(err2), vals, out=np.zeros_like(err2),
                    where=(vals != 0)
                )
            else:
                err = sum(map(lambda x, y: abs(x * y), valf, voxerrs))
                errs = np.divide(err, vals, out=np.zeros_like(err),
                                 where=(vals != 0))
        else:
            volOK = celfrac[1][0]
            vals = voxvals[0]
            errs = voxerrs[0]

            # insert value and error total over all energy bins
            xs = sum(vals)
            if xs != 0:
                ve = vals * errs
                if not corr:
                    err2 = sum(ve * ve)
                    errt = np.sqrt(err2) / xs
                else:
                    errt = sum(abs(ve)) / xs
                errs = np.append(voxerrs[0], errt)
            else:
                errs = np.append(voxerrs[0], 0)
            vals = np.append(vals, xs)

    if Vmult == "vtot":
        # multiply by the total volume of the voxel
        vals = celfrac[0] * vals
    elif Vmult == "celf":
        # multiply by the voxel volume fraction
        vals = vals / volOK

    return vals, errs


def identical_mesh(m1, m2):
    part, mesh, mtype = True, True, True
    if m1.part != m2.part:
        part = False

    if not m1.sameMesh(m2):
        mesh = False

    if m1.type != m2.type:
        mtype = False

    return part, mesh, mtype


def scalemesh(m1, f1):
    mscal = Fmesh(m1.meshtal)
    mscal.filled = True
    mscal.dims = m1.dims
    mscal.ldims = m1.ldims
    mscal.cart = m1.cart
    mscal.part = m1.part
    mscal.etag = m1.etag
    mscal.type = m1.type
    mscal.origin = m1.origin
    mscal.rotation = m1.rotation
    if not m1.cart:
        mscal.axis = m1.axis
        mscal.vec = m1.vec

    mscal.dat = m1.dat * f1
    mscal.err = m1.err
    return mscal


# sum two meshes with the same dimemsions
def addmesh(m1, m2, f1=1.0, f2=1.0, corr=False):
    if m1.part != m2.part:
        print("Warning: particle type are different")

    if not m1.sameMesh(m2):
        print("mesh dimensions are not equal")
        return None

    if m1.type != m2.type:
        print("mesh type are different")
        return None

    msum = Fmesh(m1.meshtal)

    msum.filled = True
    msum.dims = m1.dims[:]
    msum.ldims = m1.ldims[:]
    msum.cart = m1.cart
    msum.part = m1.part
    msum.etag = m1.etag
    msum.type = m1.type
    msum.origin = m1.origin
    msum.rotation = m1.rotation
    if not m1.cart:
        msum.axis = m1.axis
        msum.vec = m1.vec

    msum.dat = m1.dat * f1 + m2.dat * f2

    if not corr:
        sf1 = f1 * f1
        sf2 = f2 * f2
        numerator = np.sqrt(sf1 * (m1.dat * m1.err) ** 2 + sf2 * (m2.dat * m2.err) ** 2)
        msum.err = np.divide(
            numerator, msum.dat, out=np.zeros_like(numerator), where=(msum.dat != 0)
        )
    else:
        numerator = abs(f1 * m1.dat * m1.err) + abs(f2 * m2.dat * m2.err)
        msum.err = np.divide(
            numerator, msum.dat, out=np.zeros_like(numerator), where=(msum.dat != 0)
        )

    return msum


# provide difference between two meshes
def diffmesh(m1, m2, absvalue=False, relative=False):
    if m1.part != m2.part:
        print("Warning: particle type are different")

    if not m1.sameMesh(m2):
        print("mesh dimensions are not equal")
        return None

    if m1.type != m2.type:
        print("mesh type are different")
        return None

    msum = Fmesh(m1.meshtal)

    msum.filled = True
    msum.dims = m1.dims[:]
    msum.ldims = m1.ldims[:]
    msum.cart = m1.cart
    msum.part = m1.part
    msum.etag = m1.etag
    msum.type = m1.type
    msum.origin = m1.origin
    msum.rotation = m1.rotation
    if not m1.cart:
        msum.axis = m1.axis
        msum.vec = m1.vec

    if absvalue:
        msum.dat = abs(m1.dat) - abs(m2.dat)
    else:
        msum.dat = m1.dat - m2.dat

    denominator = np.sqrt((m1.dat * m1.err) ** 2 + (m2.dat * m2.err) ** 2)
    msum.err = np.divide(
        msum.dat, denominator, out=np.zeros_like(msum.dat),
        where=(denominator != 0)
    )

    if relative:
        denominator = m1.dat + m2.dat
        msum.err = np.divide(
            2 * msum.dat,
            denominator,
            out=np.zeros_like(msum.dat),
            where=(denominator != 0),
        )

    return msum


# sum two meshes with the same dimemsions
def addbin(m1, binlist, flist=[], corr=False):
    nbins = len(binlist)
    if flist != []:
        if len(flist) != nbins:
            print("bin length and factor lenght is different")
            return None

    if m1.etag == "times":
        print("Times bin cannot be added")
        return None

    if nbins >= m1.ldims[0]:
        print("too many bin numbers")
        return None

    for b in binlist:
        if b >= m1.ldims[0] - 1:
            print("at least one bin index exceed mesh bin length")
            return None

    msum = Fmesh(m1.meshtal)
    msum.dtype = m1.dtype

    msum.filled = True
    msum.dims = m1.dims[:]
    msum.ldims = m1.ldims[:]
    msum.part = m1.part
    msum.etag = m1.etag
    msum.type = m1.type
    msum.origin = m1.origin
    msum.rotation = m1.rotation
    if not m1.cart:
        msum.axis = m1.axis
        msum.vec = m1.vec
    msum.cart = m1.cart

    msum.ldims[0] = 1
    msum.dims[0] = np.array([0], m1.dtype)

    msum.dat = np.zeros(msum.ldims, m1.dtype)
    msum.err = np.zeros(msum.ldims, m1.dtype)

    if flist == []:
        for b in binlist:
            msum.dat[0, :, :, :] += m1.dat[b, :, :, :]
            if not corr:
                msum.err[0, :, :, :] += (m1.err[b, :, :, :] * m1.dat[b, :, :, :]) ** 2
            else:
                msum.err[0, :, :, :] += abs(m1.err[b, :, :, :] * m1.dat[b, :, :, :])
    else:
        for i, b in enumerate(binlist):
            msum.dat[0, :, :, :] += m1.dat[b, :, :, :] * flist[i]
            if not corr:
                msum.err[0, :, :, :] += (
                    m1.err[b, :, :, :] * m1.dat[b, :, :, :] * flist[i]
                ) ** 2
            else:
                msum.err[0, :, :, :] += abs(
                    m1.err[b, :, :, :] * m1.dat[b, :, :, :] * flist[i]
                )

    if not corr:
        msum.err = np.divide(
            np.sqrt(msum.err),
            msum.dat,
            out=np.zeros_like(msum.err),
            where=(msum.dat != 0),
        )
    else:
        msum.err = np.divide(
            msum.err, msum.dat, out=np.zeros_like(msum.err),
            where=(msum.dat != 0)
        )

    return msum


# ================== END OF FUNCTION DEFINITIONS =================
