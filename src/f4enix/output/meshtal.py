"""This module is related to the parsing of meshtal files output from MCNP.

It includes parsing capabilities not only for classical MCNP fmesh, but also
for its variations introduced by the D1SUNED code, such as the Cell Under
Voxel (CuV) method.
"""

from __future__ import annotations

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or - 
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at: 
    https://eupl.eu/1.2/en/  
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""

import numpy as np
import vtk
import csv
import pyvista as pv
import os
import time
from scipy.spatial.transform import Rotation as R
from io import open
import logging
from tqdm import tqdm
from copy import deepcopy


ALLOWED_NORMALIZATIONS = ["vtot", "celf", None]
ALLOWED_OUTPUT_FORMATS = ["point_cloud", "ip_fluent", "csv", "vtk"]
ALLOWED_PARTICLES = [
    "Aneutron",
    "Alambda0",
    "Asigma+",
    "Asigma-",
    "Axi0",
    "Anu_e",
    "Anu_m",
    "Aproton",
    "Aomega-",
    "neutron",
    "photon",
    "electron",
    "mu_minus",
    "nu_e",
    "nu_m",
    "VALID",
    "proton",
    "lambda0",
    "sigma+",
    "sigma-",
    "xi0",
    "xi_minus",
    "omega",
    "mu_plus",
    "pi_plus",
    "pi_zero",
    "k_plus",
    "k0_short",
    "k0_long",
    "xi_plus",
    "deuteron",
    "triton",
    "helion",
    "alpha",
    "pi_minus",
    "k_minus",
    "heavyion",
]


# convert character to float
# able to convert no standard fortran exponent type 1.234-123
def _dfloat(c: str) -> float:
    try:
        x = float(c)
    except:
        c = c.strip()
        m = float(c[0:-4])
        e = int(c[-4:])
        x = m * pow(10, e)
    return x


# return a list of each nth charaters of a string
def _splitn(string: str, n: int) -> list[str]:
    return [string[i : i + n].strip() for i in range(0, len(string), n)]


# Salta lineas (definida en cien sitios)
def _skipLines(f, n):
    for _ in range(n):
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
    # IPT = ("neutron", "photon", "electron")

    # Reads from opened file
    def __init__(self, mshtl: Meshtal) -> None:
        """Object representing an Fmesh result.

        It needs a Meshtal object to be initialized. Once read, a
        PyVistaWrapper object will be assigned to it that can be used to
        manipulate and export the parsed data.

        Parameters
        ----------
        mshtl : Meshtal
            meshtal object where the Fmesh is located

        Attributes
        ----------
        grid: pv.DataSet
            pyvista grid of the parsed data.
        meshtal: Meshtal
            the same meshtal object provided as parameter.
        comment: str
            comment string for the tally provided in MCNP input
        part: Literal['neutron', 'photon', 'electron']
            particle tallied in the fmesh.
        type: str
            type of fmesh, e.g. 'cuv', 'source', etc.
        normalization: str
            type of normalization used for the CuV approach
        filled : bool
            if True, the Fmesh has been populated with data.
        cart: bool
            if True the mesh is cartesian, otherwise cylindrical.
        ntally: int
            tally number for the fmesh

        Examples
        --------
        For examples check out :py:class:`Meshtal` doc.

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

        # here is stored the vtk object once the Fmesh is read
        self.grid = None

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
        for p in ALLOWED_PARTICLES:
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
                (vals[8 - ishft], vals[9 - ishft], vals[10 - ishft]), self.dtype
            )
            if line.startswith(mcnp6_Cyl):
                self.vec = np.array((vals[13], vals[14], vals[15]), self.dtype)
            line = f.readline()
        # Read dimensions  including energy
        # (E,Z,Y,X) order in storage in cartesian
        # (E,Th,Z,R) order in cylindrical
        for n in range(3):
            i = line.index(":")
            self.dims.insert(0, np.array(line[i + 1 :].split(), self.dtype))
            line = f.readline()
        self.ldims = [len(x) - 1 for x in self.dims]

        # process energy bins
        i = line.index(":")
        if self.usrbin:
            # decay time bins
            if "times" in line:
                self.etag = "times"
                self.dims.insert(0, np.array(line[i + 1 :].split(), self.dtype))
                self.ener = self.dims[0]  # shallow copy
                self.ldims.insert(0, len(self.dims[0]))
            # source importance bins
            elif "Target tally" in line:
                self.etag = "tally"
                line = f.readline()
                i = line.index(":")
                self.dims.insert(0, np.array(line[i + 1 :].split(), self.dtype))
                self.ener = self.dims[0]  # shallow copy
                self.ldims.insert(0, len(self.dims[0]))
            # cell or isotopes bin
            else:
                self.etag = line[0:i].split()[0].lower()
                if self.etag == "cells":
                    binlen = 12
                else:
                    binlen = 10
                binlist = _splitn(line[i + 1 :], binlen)[:-1]
                lastbin = binlist[-1]
                ib = lastbin.index(".")
                nbin = int(lastbin[ib + 1 :])
                self.dims.insert(0, np.array(range(nbin), self.dtype))
                #        self.ener = np.array(line[i+1:].split(),self.dtype)
                self.ener = np.array(binlist, self.dtype)
                self.ldims.insert(0, nbin + 1)
        # energy bin
        else:
            if "number:" in line:
                i = line.index(":")
                ne = int(line[i + 1 :])
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
                    _skipLines(f, 5)
                    onef = _checkonef(f, nblck)
                    f.seek(pos)
                    if onef:
                        self.__format__ = "onef"
                    else:
                        self.__format__ = "mltf"

            self.dims.insert(0, np.array(line[i + 1 :].split(), self.dtype))
            self.ener = self.dims[0]  # shallow copy
            nbin = 1
            if len(self.dims[0]) > 2:
                nbin = len(self.dims[0])
            self.ldims.insert(0, nbin)

        if self.__format__ == "undef":
            self.__format__ = "mcnp"
        self.startm = f.tell()

    def _readMCNP(self, f) -> None:
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
            _skipLines(f, 1)
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
                    _skipLines(f, 3)
                    if self.etag == "tally":
                        self.tallyscore.append(float(f.readline().split()[-1]))

                for ix1 in range(rshape[0]):
                    _skipLines(f, hlines)
                    hlines = hlinesNext  # subsequent reads
                    for ix0 in range(rshape[1]):
                        xdat[ix1, ix0, :] = list(map(_dfloat, f.readline().split()[1:]))
                    # use dfloat conversion in case of no standar fortan
                    # exponent
                    _skipLines(f, 3)
                    for ix0 in range(rshape[1]):
                        xerr[ix1, ix0, :] = f.readline().split()[
                            1:
                        ]  # automatic conversion
                    _skipLines(f, 2)
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
                        xdat[ix] = _dfloat(line[colDat : colDat + 12])
                        xerr[ix] = line[colErr : colErr + 12]
                    except:
                        xdat[ix] = _dfloat(line[colDat + 1 : colDat + 13])
                        xerr[ix] = line[colErr + 1 : colErr + 13]
                self.dat[ie, :, :, :] = np.transpose(xdat.reshape(rshape), itrn)
                self.err[ie, :, :, :] = np.transpose(xerr.reshape(rshape), itrn)

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
                            (self.dims[1][i + 1 : -1] - 1.0, self.dims[1][0:i]), axis=0
                        )
                        self.dat = np.concatenate(
                            (self.dat[:, i + 1 :, :, :], self.dat[:, : i - 1, :, :]),
                            axis=1,
                        )
                        self.err = np.concatenate(
                            (self.err[:, i + 1 :, :, :], self.err[:, : i - 1, :, :]),
                            axis=1,
                        )
                    else:
                        np.delete(self.dims[1], i, 0)
                    # Leave loop
                    break

        self.filled = True
        # if it is filled it means that the vtk object can be created
        name = "{}_{}".format(self.meshtal.filename, self.ntally)
        if self.cart:
            self.grid = self._getVTKrg()
        else:
            self.grid = self._getVTKsg()

    # Read photonfile format of SRCIMP mesh (D1SUNED)
    def _readSRCTYPE(self, f, cfilter=None, norm=None):
        mcnp5_Cyl = "  Cylinder origin at"
        mcnp6_Cyl = "               origin at"
        self.etag = "energy"
        if norm is None:
            tag = "[/cc]"
        elif norm == "vtot":
            tag = "integral value"
        elif norm == "celf":
            tag = "[/cc-cell]"
        else:
            tag = ""

        self.tag = tag  # to be used to clarify output array

        self.normalization = norm

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
            _skipLines(f, 2)
        else:
            if self.type == "srcimp":
                self.tallyscore = np.float64(f.readline().split()[-1])
            _skipLines(f, 5)

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
            onef = _checkonef(f, nblck)
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
                for i in range(ne):
                    vals.extend(f.readline().split())
                for i in range(ne):
                    errs.extend(f.readline().split())

                if nbin > 1:
                    xdat[nbin - 1, k] = _dfloat(tot)
                    xerr[nbin - 1, k] = 0.0
                for i in range(ne):
                    xdat[i, k] = _dfloat(vals[i])
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
                            celdata[2].append([float(vline[2]), float(vline[3])])

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

                        vals = np.array([_dfloat(x) for x in valc], self.dtype)
                        errs = np.array([_dfloat(x) for x in errc], self.dtype)

                        if nbin > 1:
                            vals = np.append(vals, celdata[2][ind][0])
                            errs = np.append(errs, celdata[2][ind][1])

                        voxvals.append(vals)
                        voxerrs.append(errs)
                    else:
                        _skipLines(f, 2 * nline)

                if len(voxvals) == 0:
                    vals = [0] * nbin
                    errs = [0] * nbin
                else:
                    # sum all bin if not total bin
                    vals, errs = _sumCellInVox(voxvals, voxerrs, celdata, Vmult=norm)

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

                    vals = np.array([_dfloat(x) for x in valc], self.dtype)
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

                            vals = np.array([_dfloat(x) for x in valc], self.dtype)
                            errs = np.array(errc, self.dtype)
                            voxvals.append(vals)
                            voxerrs.append(errs)
                        else:
                            _skipLines(f, 2 * nline)

                if len(voxvals) == 0:
                    vals = [0] * nbin
                    errs = [0] * nbin
                    tot = 0
                else:
                    # sum all bin if not total bin
                    vals, errs = _sumElements(
                        voxvals,
                        voxerrs,
                        celfrac,
                        Vsum=self.__format__,
                        Vmult=norm,
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
        # if it is filled it means that the vtk object can be created
        # name = '{}_{}'.format(self.meshtal.filename, self.ntally)
        if self.cart:
            self.grid = self._getVTKrg()
        else:
            self.grid = self._getVTKrg()

    #  end modifs

    def print_info(self) -> None:
        """print information related to the fmesh"""
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
            print(" Mesh type      : {}".format(meshtype))
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
                return _format_XYZ_Dim_long(vec)
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

        def _format_XYZ_Dim_long(vec: list, nval: int = 8) -> str:
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
            return _format_XYZ_Dim_long(vec)
        else:
            return format_XYZ_Dim_inter(vec)

    def print_EbinRange(self) -> None:
        """print the energy range bins"""
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

    # # Translate the mesh
    # def translate(self, vec: list[float]) -> None:
    #     """translate the fmesh using a vector x, y, z.

    #     Parameters
    #     ----------
    #     vec : list[float]
    #         translation vector
    #     """
    #     self.origin += np.array(vec)

    # Salida en formato VTK
    # Escribe structured grid
    def _getVTKsg(self) -> pv.DataSet:
        """get the VTK structured grid format for the mesh

        Returns
        -------
        pv.DataSet
            generated vtk object
        """
        import math

        kdims = deepcopy(self.dims)
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
        if self.__format__ == "cuv":
            value_tag = "Value " + self.tag
        else:
            value_tag = "Value - Total"

        it = 0
        if self.etag not in ["times", "tally"]:
            # Dataset energia total
            it = 1
            sgcd.AddArray(
                _makeVTKarray(self.dat[-1, :, :, :], value_tag, self.scaleFac)
            )
            sgcd.AddArray(_makeVTKarray(self.err[-1, :, :, :], "Error - Total"))
        # Dataset other bins
        for ie in range(self.ldims[0] - it):
            sgcd.AddArray(
                _makeVTKarray(
                    self.dat[ie, :, :, :], "ValueBin-{0:03d}".format(ie), self.scaleFac
                )
            )
            sgcd.AddArray(
                _makeVTKarray(self.err[ie, :, :, :], "ErrorBin-{0:03d}".format(ie))
            )
        # Include weight windows in VTK file (now only one group)

        # TODO use pyvista to create the object
        # for the moment just wrap in a pyvista object
        return pv.wrap(sg)

    # def _writeVTKsg(self, ofn: os.PathLike) -> None:
    #     """Write the mesh to a vtk file

    #     Parameters
    #     ----------
    #     ofn : os.PathLike
    #         outpath for the file
    #     """
    #     try:
    #         name = 'tally_{}_{}.vts'.format(self.ntally, self.tag)
    #     except AttributeError:
    #         name = 'tally_{}.vts'.format(self.ntally)

    #     # Escritura en disco
    #     off = vtk.vtkXMLStructuredGridWriter()
    #     off.SetFileName(os.path.join(ofn, name))
    #     # ASCII o binario (con o sin compresion)
    #     off.SetDataModeToAscii()
    #     # off.SetDataModeToBinary()
    #     # off.SetCompressorTypeToZLib()
    #     # Esto cambia con la version de VTK
    #     t = self.getVTKsg()
    #     self.meshtal._setVTKparams(t)
    #     if vtk.vtkVersion().GetVTKMajorVersion() >= 6:
    #         off.SetInputData(t)
    #     else:
    #         off.SetInput(t)
    #     off.Write()

    # Escribe rectilinear grid
    # No vale para cilindricas o rotadas
    def _getVTKrg(self) -> pv.DataSet:
        """Get a rectilinear grid. This cannot be used for cylindrical or
        rotated grids

        Returns
        -------
        pv.DataSet
            generated rectilinear grid object
        """
        if not self.cart:
            logging.warning("Cylindrical meshtal cannot be plotted to RectangularGrid")

        if np.any(self.rotation != np.identity(3)):
            logging.warning("Rotated meshtal cannot be plotted to RectangularGrid")
            logging.warning("... but plotting anyway (no rotation)")

        xa = _makeVTKarray(self.dims[3] + self.origin[3], "X (cm)")
        ya = _makeVTKarray(self.dims[2] + self.origin[2], "Y (cm)")
        za = _makeVTKarray(self.dims[1] + self.origin[1], "Z (cm)")

        rg = vtk.vtkRectilinearGrid()
        rg.SetDimensions(self.ldims[3] + 1, self.ldims[2] + 1, self.ldims[1] + 1)
        rg.SetXCoordinates(xa)
        rg.SetYCoordinates(ya)
        rg.SetZCoordinates(za)
        rgcd = rg.GetCellData()
        it = 0
        if self.__format__ == "cuv":
            value_tag = "Value " + self.tag
        else:
            value_tag = "Value - Total"

        if self.etag not in ["times", "tally"]:
            # Dataset energia total
            it = 1
            rgcd.AddArray(
                _makeVTKarray(self.dat[-1, :, :, :], value_tag, self.scaleFac)
            )
            rgcd.AddArray(_makeVTKarray(self.err[-1, :, :, :], "Error - Total"))
        # Dataset other bins
        for ie in range(self.ldims[0] - it):
            rgcd.AddArray(
                _makeVTKarray(
                    self.dat[ie, :, :, :], "ValueBin-{0:03d}".format(ie), self.scaleFac
                )
            )
            rgcd.AddArray(
                _makeVTKarray(self.err[ie, :, :, :], "ErrorBin-{0:03d}".format(ie))
            )

        # TODO use pyvista to create the object
        # for the moment just wrap in a pyvista object
        return pv.wrap(rg)

    # def _writeVTKrg(self, ofn: os.PathLike) -> None:
    #     """Write the rectilinear grid

    #     Parameters
    #     ----------
    #     ofn : os.PathLike
    #         output vtk file path
    #     """
    #     try:
    #         name = 'tally_{}_{}.vtr'.format(self.ntally, self.tag)
    #     except AttributeError:
    #         name = 'tally_{}.vtr'.format(self.ntally)
    #     # Escritura en disco
    #     off = vtk.vtkXMLRectilinearGridWriter()
    #     off.SetFileName(os.path.join(ofn, name))
    #     off.SetDataModeToAscii()
    #     # Esto cambia con la version de VTK
    #     t = self._getVTKrg()
    #     self.meshtal._setVTKparams(t)
    #     if vtk.vtkVersion().GetVTKMajorVersion() >= 6:
    #         off.SetInputData(t)
    #     else:
    #         off.SetInput(t)
    #     off.Write()

    # def writeVTK(self, ofn: os.PathLike) -> None:
    #     """Write the fmesh to a vtk file

    #     Parameters
    #     ----------
    #     ofn : os.PathLike
    #         path to the vtk outfile folder
    #     """
    #     if not os.path.exists(ofn):
    #         raise IsADirectoryError('Directory does not exists {}'.format(ofn))

    #     try:
    #         name = 'tally_{}_{}'.format(self.ntally, self.tag)
    #     except AttributeError:
    #         name = 'tally_{}'.format(self.ntally)

    #     if self.cart:
    #         dataset = self._getVTKrg()
    #     else:
    #         dataset = self._getVTKrg()

    def write(
        self,
        outpath: os.PathLike,
        list_array_names: list[str] = None,
        out_format: str = "vtk",
        outfile: str = None,
    ) -> None:
        """Export the mesh to a file. vtk, csv, fluent (txt) and point cloud
        (txt) formats can be selected.

        Parameters
        ----------
        outpath : os.PathLike
            path to the output folder.
        list_array_names : list[str], optional
            arrays to be exported. The default is None, meaning that all the
            available arrays will be used.
        out_format : str, optional
            output format. The allowed ones are ['point_cloud', 'ip_fluent',
            'csv', 'vtk']. Default is .vtk
        outfile : str, optional
            name of the output file. If specified, overrides the default one.
            Do not include the extension of the file here. Default is None.

        Raises
        ------
        KeyError
            raises KeyError if the output format is not allowed.
        """
        if list_array_names is None:
            list_array_names = list(self.grid.array_names)

        if outfile is None:
            file_name = f"{self.meshtal.filename}_{self.ntally}_{out_format}"
        else:
            file_name = outfile

        # TODO either all cells or all point data are supported if not a vtk
        if out_format != "vtk":
            len_data = len(self.grid.cell_data)
            len_point = len(self.grid.point_data)
            if len_data > 0 and len_point == 0:
                f_points = self.grid.cell_centers().points
            elif len_point > 0:
                f_points = self.grid.points
            else:
                raise ValueError(
                    "mix between cell and point data is only supported for vtk"
                )

        filepath = os.path.join(outpath, file_name)
        mesh_type = str(type(self.grid)).split(".")[-1][:-2]

        if out_format == "vtk":
            if mesh_type == "StructuredGrid":
                ext = ".vts"
            elif mesh_type == "UnstructuredGrid":
                ext = ".vtu"
            elif mesh_type == "RectilinearGrid":
                ext = ".vtr"
            else:
                ext = ".vtk"

            self.grid.save(filepath + ext)
            return

        # --- CSV writer ---
        elif out_format == "csv":
            new_name = filepath + ".csv"

            with open(new_name, "w", newline="") as outfile:
                writer = csv.writer(outfile)

                # # TODO This may create some issues...
                # values_type = self.get_array_type(list_array_names[0])
                # if values_type == "cells":  # Take points or centers
                #     f_points = self.centers
                # else:  # Points
                #     f_points = self.points

                for i in tqdm(range(len(f_points)), unit=" Points", desc="Writing"):
                    csv_points = [
                        f"{f_points[i][0]:.3f}",
                        f" {f_points[i][1]:.3f}",
                        f" {f_points[i][2]:.3f}",
                    ]
                    for array_name in list_array_names:
                        csv_points.append(f" {self.grid[array_name][i]:.3f}")
                    writer.writerow(csv_points)

                logging.info(f"{new_name} created successfully!")
            return

        for array_name in list_array_names:
            # values_type = self.get_array_type(list_array_names[0])
            # if values_type == "cells":  # Take points or centers
            #     f_points = self.centers
            # else:  # Points
            #     f_points = self.points
            values = self.grid[array_name]

            # write depending on format
            # --- point cloud writer ---
            if out_format == "point_cloud":
                new_name = filepath + ".txt"
                with open(new_name, "w") as outfile:
                    outfile.write("x, y, z, value\n")
                    # TODO this can probably be optmized using
                    # outfile.writeline()
                    for i in tqdm(range(len(f_points)), unit=" Points", desc="Writing"):
                        outfile.write(f"{f_points[i][0]:.3f},")
                        outfile.write(f"{f_points[i][1]:.3f},")
                        outfile.write(f"{f_points[i][2]:.3f},")
                        outfile.write(f"{values[i]:.3f}\n")
                logging.info(f"{new_name} created successfully!")
                return

            # --- fluent writer ---
            elif out_format == "ip_fluent":
                new_name = filepath + ".txt"

                with open(new_name, "w") as outfile:
                    guion1 = "3"
                    n_coord = f_points.shape[1]  # self.n_coordinates
                    n_values = str(len(f_points))
                    guion2 = "1"
                    uds = "uds-0"
                    beginning = f"{guion1}\n{n_coord}\n{n_values}\n{guion2}\n{uds}\n"
                    outfile.write(beginning)
                    outfile.write("(")
                    for i in tqdm(
                        range(len(f_points)), unit=" x points", desc="Writing x"
                    ):
                        outfile.write(f"{f_points[i][0]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(
                        range(len(f_points)), unit=" y points", desc="Writing y"
                    ):
                        outfile.write(f"{f_points[i][1]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(
                        range(len(f_points)), unit=" z points", desc="Writing z"
                    ):
                        outfile.write(f"{f_points[i][2]:.3f}\n")

                    outfile.write(")\n")
                    outfile.write("(")

                    for i in tqdm(
                        range(len(f_points)), unit=" values", desc="Writing values"
                    ):
                        outfile.write(f"{values[i]:.3f}\n")

                    outfile.write(")\n")

                logging.info(f"{new_name} created successfully!")
                return

        raise KeyError(
            "Invalid format, these are the ones allowed: {}".format(
                ALLOWED_OUTPUT_FORMATS
            )
        )

    def _read_from_vtk(self, vtk_file: os.PathLike):
        # This is mostly used for quicker testing
        grid = pv.read(vtk_file)
        self.grid = grid


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

        Attributes
        ----------
        filename: str
            name of the file which is extracted by the fn path
        filetype: str
            same as the parameter
        mesh: dict[str, Fmesh]
            a dictionary containing all the Fmesh objects that are read from
            the meshtal file

        Examples
        --------
        Examples of usage of the Meshtal object

        >>> from f4enix.output.meshtal import Meshtal
        ... # Initialize the meshtal file
        ... file = 'cuvmsh'
        ... meshtal = Meshtal(file)
        ... print(meshtal)
        Meshtally file : cuvmsh   Tally 24 : neutron  cuv mesh   'Test cell under voxel'

        and then read all meshes available and access their associated
        PyVistaWrapper object

        >>> # if no mesh are specified all mesh are read
        ... meshtal.readMesh()
        ... # Once the mesh is read, the grid attribute gets filled
        ... print(type(meshtal.mesh[24].grid))
        ... meshtal.mesh[24].grid
        {24: <f4enix.output.meshtal.Fmesh at 0x1badb9d4310>}
        ...
        <class 'pyvista.core.grid.RectilinearGrid'>
        ...
        Name: cuvmsh_24
        ...
        Header	Data Arrays
        RectilinearGrid	Information
        N Cells	2293504
        N Points	2346125
        X Bounds	-1.700e+03, 1.700e+03
        Y Bounds	-1.700e+03, 1.700e+03
        Z Bounds	-1.360e+03, 1.740e+03
        Dimensions	137, 137, 125
        N Arrays	10
        Name	        Field	Type	N   	Min	       Max
        Value - Total	Cells	float64	1	0.000e+00	5.286e+08
        Error - Total	Cells	float64	1	0.000e+00	1.000e+00
        ValueBin-000	Cells	float64	1	0.000e+00	1.547e+05
        ErrorBin-000	Cells	float64	1	0.000e+00	1.000e+00
        ValueBin-001	Cells	float64	1	0.000e+00	3.497e+08
        ErrorBin-001	Cells	float64	1	0.000e+00	1.000e+00
        ValueBin-002	Cells	float64	1	0.000e+00	1.787e+08
        ErrorBin-002	Cells	float64	1	0.000e+00	1.000e+00
        ValueBin-003	Cells	float64	1	0.000e+00	0.000e+00
        ErrorBin-003	Cells	float64	1	0.000e+00	0.000e+00

        Cell Under Voxel (CuV) approach is supported where cell_filters and
        different normalization can be chosen. The fmesh can be exported to
        different formats

        >>> # Read a specific CuV fmesh filtering by a set of cells
        ... # and changing the normalization to get the integral result
        ... meshtal.readMesh(24, cell_filters=[1], norm='vtot')
        ... fmesh = meshtal.mesh[24]
        ... fmesh.print_info()
        ... outpath = 'folderpath'
        ... # it can be exported to different formats, extension will be
        ... # inferred automatically
        ... fmesh.write(outpath)
        ... fmesh.write(outpath, out_format='ip_fluent')
        Tally          : 24
        Comments       :
            Test cell under voxel
        Particle       : neutron
        Mesh type      : cuv
        Dose modif     : False
        Mesh geometry  : rectangular
        Mesh origin    : 0.0 0.0 0.0
        X dimensions   :   0.000e+00   50I  1.020e+02
        Y dimensions   :   0.000e+00   50I  1.020e+02
        Z dimensions   :   0.000e+00   50I  1.020e+02
        Energy bins    :
                flag : energy
            bin index : bin range
        Writing x: 100%|██████████| 132651/132651 [00:00<00:00, 204617.39 x points/s]
        Writing y: 100%|██████████| 132651/132651 [00:00<00:00, 217292.19 y points/s]
        Writing z: 100%|██████████| 132651/132651 [00:00<00:00, 227737.30 z points/s]
        Writing values: 100%|██████████| 132651/132651 [00:00<00:00, 770635.26 values/s]

        """
        logging.info("Loading Meshtal: {}".format(fn))
        self.filename = os.path.basename(fn).split(".")[0]
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

    def readMesh(
        self,
        mesh: int | list[int] = None,
        cell_filters: list[int] = None,
        norm: str = None,
    ) -> None:
        """Parse a list of FMESHes

        Parameters
        ----------
        mesh : int | list[int], optional
            list of FMESH to be read, by default None means that all meshes
            will be parsed
        cell_filters : list[int], optional
            list of cells to be used as filters for the CuV approach,
            by default None
        norm : str
            to be used only on CuV. Can be set either to "vtot" or "celf".
            Default is None which means that the normalzation will be done
            by dividing by the Voxel volume. 'vtot' means that the default
            (None) will be multiplied by the total voxel volume to obtain an
            integral value, while "celf" normalizes dividing by the volume of
            the cell fraction under the voxel.

        Raises
        ------
        KeyError
            if the filetype is is not implemented
        """
        logging.info("Reading fmesh: {}".format(mesh))
        if norm not in ALLOWED_NORMALIZATIONS:
            raise NameError("{} is not an allowed normalization keyword".format(norm))
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
                    m._readMCNP(self.f)
                else:
                    # if m.__mltflt__ is not None:
                    #     flist = get_clist(m.__mltflt__)
                    # else:
                    #     flist = None
                    m._readSRCTYPE(self.f, cfilter=cell_filters, norm=norm)
        else:
            raise KeyError("This file type is not implemented: " + str(self.filetype))

    def write_all(self, outpath: os.PathLike, out_format: str = "vtk") -> None:
        """write all fmeshes to the outfolder in the specified format.

        Parameters
        ----------
        outpath : os.PathLike
            path to the output folder.
        out_format : str, optional
            output format. The allowed ones are ['point_cloud', 'ip_fluent',
            'csv', 'vtk']. Default is .vtk

        """
        for _, mesh in self.mesh.items():
            mesh.write(outpath, out_format=out_format)

    def __readHeadMCNP__(self) -> None:
        # there is a different version that has a different header
        vals = self.f.readline()
        if "C ====" in vals:
            # this is the weird header that do not give you data
            self.code = None
            self.version = None
            self.probid = None
            self.title = None
        else:
            vals = vals.split()
            self.code = vals[0]
            self.version = vals[2]
            self.probid = vals[-2] + " " + vals[-1]
            self.title = self.f.readline().strip()

        self.nps = int(
            float((self.f.readline().split()[-1]))
        )  # nps: int doesnt like decimals

    def collapse_grids(self, name_dict: dict[int, list[str, str]]) -> pv.DataSet:
        """If the all the fmeshes indicated in the dictionary are defined on
        the same
        structured grid, returns a grid onto which all the fmeshes are
        collapsed. That is, the returned grid will have all the field data that
        are stored in the different fmeshes.

        Parameters
        ----------
        name_dict : dict[int, list[str, str]]
            this dictionary is used to assign names to the fields that are
            added to the grid. The key of the dictionary is the number of the
            fmesh, while the list of strings are
            ['name of the values', 'name of the values statistical error'].
            For instance {104: ['neutron heating [W/cc]',
            'neutron heating [Rel err]']}

        Returns
        -------
        pv.DataSet
            returns the pyvista grid where all fields from the different
            Fmeshes have been added

        Raises
        ------
        RuntimeError
            if the fmeshes have different geometry if they do not contain
            only the default fields ['Value - Total', 'Error - Total'].
        """

        ids = ["Value - Total", "Error - Total"]

        try:
            # check that the collapse is doable
            for i, key in enumerate(list(name_dict.keys())):
                fmesh = self.mesh[key]
                # Check they are all same size
                if i == 0:
                    try:
                        comparison = fmesh
                    except AttributeError:
                        # it means that the mesh has not been read
                        raise AttributeError(
                            "Please read the meshtal file first using readMesh()"
                        )
                else:
                    assert comparison.sameMesh(fmesh)
        except AssertionError:
            raise RuntimeError("the different fmeshes geom are not compatible")

        try:
            # check that the collapse is doable
            for i, key in enumerate(list(name_dict.keys())):
                fmesh = self.mesh[key]
                # check that there are only the two usual values, no binning
                assert fmesh.grid.array_names == ids
        except AssertionError:
            raise RuntimeError("no binning allowed for the collapse")

        for i, key in enumerate(list(name_dict.keys())):
            fmesh = self.mesh[key]
            if i == 0:
                grid = deepcopy(fmesh.grid)
                for old_name, new_name in zip(grid.array_names, name_dict[key]):
                    grid.rename_array(old_name, new_name)
            else:
                # results just have to be added for the other fmeshes
                for array_name, id in zip(name_dict[key], ids):
                    grid[array_name] = fmesh.grid[id]

        return grid

    # def writeVTK(self, ofn: os.PathLike) -> None:
    #     """write a .vtk file for all read fmeshes

    #     Parameters
    #     ----------
    #     ofn : os.PathLike
    #         path to the output folder
    #     """

    #     mb = vtk.vtkMultiBlockDataSet()
    #     i = 0
    #     for k, m in self.mesh.items():
    #         if not m.filled:
    #             continue
    #         if m.cart:
    #             mb.SetBlock(i, m._getVTKrg())
    #             extension = '.vtr'
    #         else:
    #             mb.SetBlock(i, m._getVTKsg())
    #             extension = '.vts'
    #         i += 1
    #     # A escribir
    #     outpath = os.path.join(ofn, self.filename+extension)
    #     off = vtk.vtkXMLMultiBlockDataWriter()
    #     self._setVTKparams(mb)
    #     off.SetInputData(mb)
    #     off.SetFileName(outpath)
    #     off.Write()

    # def _addVTKparams(self, xdict: dict) -> None:
    #     """Add more parameters to include in VTK

    #     Parameters
    #     ----------
    #     xdict : dict
    #         parameters to be added
    #     """
    #     self.params.update(xdict)

    # def _setVTKparams(self, xdat: vtk.vtkDataSet) -> None:
    #     """Modify VTK structure for including parameters

    #     Parameters
    #     ----------
    #     xdat : vtk.vtkDataSet
    #         vtk data set to be modified
    #     """
    #     # for x,v in self.params.iteritems():
    #     for x, v in self.params.items():
    #         t = vtk.vtkStringArray()
    #         t.SetName(x)
    #         t.InsertNextValue(v)
    #         xdat.GetFieldData().AddArray(t)

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        out = ""
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
                out = out + ("   Tally {} : {}  {} mesh ".format(t, m.part, meshtype))
        return out

    def print_info(self) -> None:
        """print information on the Meshtal"""
        print(self)


# ================= END OF CLASS DEFINITIONS ======================


# This should be the same as above
def _makeVTKarray(nArr, aName, sc=1.0):
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
def _checkonef(f, lblck):
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


def _sumCellInVox(voxvals, voxerrs, celfrac, Vmult="none", corr=True, nulcount=True):
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


def _sumElements(voxvals, voxerrs, celfrac, Vsum="onef", Vmult="none", corr=False):
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
                    np.sqrt(err2), vals, out=np.zeros_like(err2), where=(vals != 0)
                )
            else:
                err = sum(map(lambda x, y: abs(x * y), valf, voxerrs))
                errs = np.divide(err, vals, out=np.zeros_like(err), where=(vals != 0))
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


# def scalemesh(m1, f1):
#     mscal = Fmesh(m1.meshtal)
#     mscal.filled = True
#     mscal.dims = m1.dims
#     mscal.ldims = m1.ldims
#     mscal.cart = m1.cart
#     mscal.part = m1.part
#     mscal.etag = m1.etag
#     mscal.type = m1.type
#     mscal.origin = m1.origin
#     mscal.rotation = m1.rotation
#     if not m1.cart:
#         mscal.axis = m1.axis
#         mscal.vec = m1.vec

#     mscal.dat = m1.dat * f1
#     mscal.err = m1.err
#     return mscal


# # sum two meshes with the same dimemsions
# def addmesh(m1, m2, f1=1.0, f2=1.0, corr=False):
#     if m1.part != m2.part:
#         print("Warning: particle type are different")

#     if not m1.sameMesh(m2):
#         print("mesh dimensions are not equal")
#         return None

#     if m1.type != m2.type:
#         print("mesh type are different")
#         return None

#     msum = Fmesh(m1.meshtal)

#     msum.filled = True
#     msum.dims = m1.dims[:]
#     msum.ldims = m1.ldims[:]
#     msum.cart = m1.cart
#     msum.part = m1.part
#     msum.etag = m1.etag
#     msum.type = m1.type
#     msum.origin = m1.origin
#     msum.rotation = m1.rotation
#     if not m1.cart:
#         msum.axis = m1.axis
#         msum.vec = m1.vec

#     msum.dat = m1.dat * f1 + m2.dat * f2

#     if not corr:
#         sf1 = f1 * f1
#         sf2 = f2 * f2
#         numerator = np.sqrt(sf1 * (m1.dat * m1.err) ** 2 + sf2 * (m2.dat * m2.err) ** 2)
#         msum.err = np.divide(
#             numerator, msum.dat, out=np.zeros_like(numerator), where=(msum.dat != 0)
#         )
#     else:
#         numerator = abs(f1 * m1.dat * m1.err) + abs(f2 * m2.dat * m2.err)
#         msum.err = np.divide(
#             numerator, msum.dat, out=np.zeros_like(numerator), where=(msum.dat != 0)
#         )

#     return msum


# # provide difference between two meshes
# def diffmesh(m1, m2, absvalue=False, relative=False):
#     if m1.part != m2.part:
#         print("Warning: particle type are different")

#     if not m1.sameMesh(m2):
#         print("mesh dimensions are not equal")
#         return None

#     if m1.type != m2.type:
#         print("mesh type are different")
#         return None

#     msum = Fmesh(m1.meshtal)

#     msum.filled = True
#     msum.dims = m1.dims[:]
#     msum.ldims = m1.ldims[:]
#     msum.cart = m1.cart
#     msum.part = m1.part
#     msum.etag = m1.etag
#     msum.type = m1.type
#     msum.origin = m1.origin
#     msum.rotation = m1.rotation
#     if not m1.cart:
#         msum.axis = m1.axis
#         msum.vec = m1.vec

#     if absvalue:
#         msum.dat = abs(m1.dat) - abs(m2.dat)
#     else:
#         msum.dat = m1.dat - m2.dat

#     denominator = np.sqrt((m1.dat * m1.err) ** 2 + (m2.dat * m2.err) ** 2)
#     msum.err = np.divide(
#         msum.dat, denominator, out=np.zeros_like(msum.dat),
#         where=(denominator != 0)
#     )

#     if relative:
#         denominator = m1.dat + m2.dat
#         msum.err = np.divide(
#             2 * msum.dat,
#             denominator,
#             out=np.zeros_like(msum.dat),
#             where=(denominator != 0),
#         )

#     return msum


# # sum two meshes with the same dimemsions
# def addbin(m1, binlist, flist=[], corr=False):
#     nbins = len(binlist)
#     if flist != []:
#         if len(flist) != nbins:
#             print("bin length and factor lenght is different")
#             return None

#     if m1.etag == "times":
#         print("Times bin cannot be added")
#         return None

#     if nbins >= m1.ldims[0]:
#         print("too many bin numbers")
#         return None

#     for b in binlist:
#         if b >= m1.ldims[0] - 1:
#             print("at least one bin index exceed mesh bin length")
#             return None

#     msum = Fmesh(m1.meshtal)
#     msum.dtype = m1.dtype

#     msum.filled = True
#     msum.dims = m1.dims[:]
#     msum.ldims = m1.ldims[:]
#     msum.part = m1.part
#     msum.etag = m1.etag
#     msum.type = m1.type
#     msum.origin = m1.origin
#     msum.rotation = m1.rotation
#     if not m1.cart:
#         msum.axis = m1.axis
#         msum.vec = m1.vec
#     msum.cart = m1.cart

#     msum.ldims[0] = 1
#     msum.dims[0] = np.array([0], m1.dtype)

#     msum.dat = np.zeros(msum.ldims, m1.dtype)
#     msum.err = np.zeros(msum.ldims, m1.dtype)

#     if flist == []:
#         for b in binlist:
#             msum.dat[0, :, :, :] += m1.dat[b, :, :, :]
#             if not corr:
#                 msum.err[0, :, :, :] += (m1.err[b, :, :, :] * m1.dat[b, :, :, :]) ** 2
#             else:
#                 msum.err[0, :, :, :] += abs(m1.err[b, :, :, :] * m1.dat[b, :, :, :])
#     else:
#         for i, b in enumerate(binlist):
#             msum.dat[0, :, :, :] += m1.dat[b, :, :, :] * flist[i]
#             if not corr:
#                 msum.err[0, :, :, :] += (
#                     m1.err[b, :, :, :] * m1.dat[b, :, :, :] * flist[i]
#                 ) ** 2
#             else:
#                 msum.err[0, :, :, :] += abs(
#                     m1.err[b, :, :, :] * m1.dat[b, :, :, :] * flist[i]
#                 )

#     if not corr:
#         msum.err = np.divide(
#             np.sqrt(msum.err),
#             msum.dat,
#             out=np.zeros_like(msum.err),
#             where=(msum.dat != 0),
#         )
#     else:
#         msum.err = np.divide(
#             msum.err, msum.dat, out=np.zeros_like(msum.err),
#             where=(msum.dat != 0)
#         )

#     return msum


# ================== END OF FUNCTION DEFINITIONS =================
