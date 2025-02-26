"""
modified version of what could be found at:
https://github.com/kbat/mc-tools

This module deals with the parsing of MCNP MCTAL files

GNU Lesser General Public License v3.0
Permissions of this copyleft license are conditioned on making available
complete source code of licensed works and modifications under the same license
or the GNU GPLv3. Copyright and license notices must be preserved.
Contributors provide an express grant of patent rights. However, a larger work
using the licensed work through interfaces provided by the licensed work may be
distributed under different terms and without source code for the larger work.
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

import logging
import math
import os
import sys

import numpy as np
import pandas as pd


class Header:
    def __init__(self):
        """Header class. Contains a bunch of general information on the mctal
        file

        Attributes
        ----------
        kod : str
            name of the code that was used
        ver : str
            name of the version that was used
        probid : np.ndarray
            date and time when the problem was run
        knod : int
            dump number
        nps : int
            number of histories that were run
        rnr : int
            Number of pseudoradom numbers that were used
        title : str
            Problem identification line
        ntal : int
            number of tallies present in the file
        ntals : np.ndarray
            arrays of tally numbers
        npert : int
            number of perturbation
        """
        self.kod = ""  # Name of the code
        self.ver = ""  # Code version
        self.probid = np.array((), dtype=str)  # Date and time when the problem was run
        self.knod = 0  # The dump number
        self.nps = 0  # Number of histories that were run
        self.rnr = 0  # Number of pseudoradom numbers that were used
        self.title = ""  # Problem identification line
        self.ntal = 0  # Number of tallies
        self.ntals = np.array((), dtype=int)  # Array of tally numbers
        self.npert = 0  # Number of perturbations


class Tally:
    """This class is aimed to store all the information contained in a
    tally.
    """

    def __init__(self, tN: int) -> None:
        self.tallyNumber = tN  # Tally number
        self.typeNumber = 0  # Particle type number
        self.detectorType = None  # The type of detector tally where 0=none, 1=point, 2=ring, 3=pinhole radiograph,
        #     4=transmitted image radiograph (rectangular grid),
        #     5=transmitted image radiograph (cylindrical grid)
        # When negative, it provides the type of mesh tally
        self.radiograph = False  # Flag set to True is the tally is a radiograph tally.
        self.tallyParticles = np.array(
            (), dtype=int
        )  # List of 0/1 entries indicating which particle types are used by the tally
        self.tallyComment = np.array((), dtype=str)  # The FC card lines
        self.nCells = 0  # Number of cell, surface or detector bins
        self.mesh = False  # True if the tally is a mesh tally
        self.meshInfo = np.array(
            [0, 1, 1, 1], dtype=int
        )  # Mesh binning information in the case of a mesh tally
        self.nDir = 0  # Number of total vs. direct or flagged vs. unflagged bins
        self.nUsr = 0  # Number of user bins
        self.usrTC = None  # Total / cumulative bin in the user bins
        self.nSeg = 0  # Number of segment bins
        self.segTC = None  # Total / cumulative bin in the segment bins
        self.nMul = 0  # Number of multiplier bins
        self.mulTC = None  # Total / cumulative bin in the multiplier bins
        self.nCos = 0  # Number of cosine bins
        self.cosTC = None  # Total / cumulative bin in the cosine bins
        self.cosFlag = 0  # The integer flag of cosine bins
        self.nErg = 0  # Number of energy bins
        self.ergTC = None  # Total / cumulative bin in the energy bins
        self.ergFlag = 0  # The integer flag of energy bins
        self.nTim = 0  # Number of time bins
        self.timTC = None  # Total / cumulative bin in the time bins
        self.timFlag = 0  # The integer flag of time bins

        self.cells = np.array(())  # Array of cell     bin boundaries
        self.usr = np.array(())  # Array of user     bin boundaries
        self.seg = np.array(())  # Array of segments bin boundaries
        self.cos = np.array(())  # Array of cosine   bin boundaries
        self.erg = np.array(())  # Array of energy   bin boundaries
        self.tim = np.array(())  # Array of time     bin boundaries
        self.cora = np.array(
            ()
        )  # Array of cora     bin boundaries for mesh tallies (or lattices)
        self.corb = np.array(
            ()
        )  # Array of corb     bin boundaries for mesh tallies (or lattices)
        self.corc = np.array(
            ()
        )  # Array of corc     bin boundaries for mesh tallies (or lattices)

        self.tfc_jtf = np.array(())  # List of numbers in the tfc line
        self.tfc_dat = []  # Tally fluctuation chart data (NPS, tally, error, figure of merit)

        self.detectorTypeList = {
            -6: "smesh",
            -5: "cmesh",
            -4: "rmesh",
            # The line below duplicates the line above with short names for
            # tally naming during conversion.
            # See the function getDetectorType to see how this information
            # is used
            -3: "Spherical mesh tally",
            -2: "Cylindrical mesh tally",
            -1: "Rectangular mesh tally",
            0: "None",
            1: "Point",
            2: "Ring",
            3: "Pinhole radiograph",
            4: "Transmitted image rdiograph (rectangular grid)",
            5: "Transmitted image radiograph (cylindrical grid)",
            # The line below duplicates the line above, with short names
            # for tally naming during conversion.
            # See the function getDetectorType to see how this information
            # is used
            6: "pi",
            7: "tir",
            8: "tic",
        }

        self.particleListShort = {
            1: "Neutron",
            2: "Photon",
            3: "Neutron + Photon",
            4: "Electron",
            5: "Neutron + Electron",
            6: "Photon + Electron",
            7: "Neutron + Photon + Electron",
        }

        self.particleList = (
            #   1         2          3          4      5
            "Neutron",
            "Photon",
            "Electron",
            "Muon",
            "Tau",
            #         6                 7
            "Electron Neutrino",
            "Muon Neutrino",
            #     8            9          10         11         12
            "Tau Neutrino",
            "Proton",
            "Lambda 0",
            "Sigma +",
            "Sigma -",
            #   13           14          15          16
            "Cascade 0",
            "Cascade -",
            "Omega -",
            "Lambda c +",
            #   17              18              19         20
            "Cascade c +",
            "Cascade c 0",
            "Lambda b 0",
            "Pion +",
            #   21             22           23        24
            "Neutral Pion",
            "Kaon +",
            "K0 Short",
            "K0 Long",
            #  25   26     27       28     29      30        31
            "D +",
            "D 0",
            "D s +",
            "B +",
            "B 0",
            "B s 0",
            "Deuteron",
            #  32      33       34              35
            "Triton",
            "He3",
            "He4 (Alpha)",
            "Heavy ions",
        )

        self.binIndexList = ("f", "d", "u", "s", "m", "c", "e", "t", "i", "j", "k")

        self.isInitialized = False
        self.valsErrors = None  # Array of values and errors

    def _initializeValuesVectors(self) -> None:
        """This function initializes the 9-D matrix for the storage of
        values and errors."""

        nCells = self._getNbins("f")
        nCora = self._getNbins("i")
        nCorb = self._getNbins("j")
        nCorc = self._getNbins("k")
        nDir = self._getNbins("d")
        nUsr = self._getNbins("u")
        nSeg = self._getNbins("s")
        nMul = self._getNbins("m")
        nCos = self._getNbins("c")
        nErg = self._getNbins("e")
        nTim = self._getNbins("t")

        self.valsErrors = np.empty(
            (nCells, nDir, nUsr, nSeg, nMul, nCos, nErg, nTim, nCora, nCorb, nCorc, 2),
            dtype=float,
        )

        self.isInitialized = True

    #     def print(self, option=[]):
    #         """Tally printer
    #         """
    #         print("\033[1m[TALLY]\033[0m")
    #         print("Tally Number: %5d" % self.tallyNumber)

    #         print("Tally comment(s):")
    #         for comment in self.tallyComment:
    #             print("\t%s" % comment)

    #         mt = "Yes" if self.mesh else "No"
    #         print("Mesh tally: %s" % mt)
    #         rt = "Yes" if self.radiograph else "No"
    #         print("Radiograph Tally: %s" % rt)
    #         print("Detector type: %s" % self.getDetectorType())

    #         print("List of particles in tally:")
    #         for i, name in enumerate(self.getTallyParticles()):
    #             print(r"\t%2d - %s" % (i+1, name))
    #         if not self.mesh:
    #             print("Number of cells/point detectors/surfaces/macrobodies: %5d" % self.getNbins("f"))
    #         else:
    #                 mesh_tot = self.getNbins("i")*self.getNbins("j")*self.getNbins("k")

    #                 print ("Number of mesh tally bins: %5d" % mesh_tot)
    #                 print ("\tNumber of CORA bins: %5d" % self.getNbins("i"))
    #                 print ("\tNumber of CORB bins: %5d" % self.getNbins("j"))
    #                 print ("\tNumber of CORC bins: %5d" % self.getNbins("k"))
    #         print ("Number of tot vs. dir or flag vs. unflag bins: %5d" % self.getNbins("d"))
    #         print ("Number of user bins: %5d" % self.getNbins("u"))
    #         print ("Number of segments: %5d" % self.getNbins("s"))
    #         print ("Number of multipliers: %5d" % self.getNbins("m"))
    #         print ("Number of cosine bins: %5d" % self.getNbins("c"))
    #         print ("Number of energy bins: %5d" % self.getNbins("e"))
    #         print ("Number of time bins: %5d" % self.getNbins("t"))

    #         print ("Total values in the tally: %8d" % self.getTotNumber(False))

    # Not called
    # def getDetectorType(self, short=False) -> str:
    #     """Returns the type of the detector type used in the tally."""

    #     if not short:
    #         return self.detectorTypeList[self.detectorType]
    #     elif short and self.radiograph:
    #         return self.detectorTypeList[self.detectorType + 3]
    #     elif short and self.mesh:
    #         return self.detectorTypeList[self.detectorType - 3]

    # not called anywhere
    # def getTallyParticles(self) -> list[str]:
    #     """Returns the particles used in the tally.
    #        References can be found in Table 4-1 and page B-2 of the
    #        MCNPX manual."""

    #     particleNames = []

    #     if self.typeNumber > 0:
    #         particleNames.append(self.particleListShort[self.typeNumber])
    #     else:
    #         for i, _ in enumerate(self.particleList):
    #             try:
    #                 if self.tallyParticles[i] == 1:
    #                     particleNames.append(self.particleList[i])
    #             except IndexError:
    #                 # For some reasons there can be less than 35 particles
    #                 # listed. Skip in case.
    #                 pass
    #     return particleNames
    # This is not called anywhere
    # def getTotNumber(self, includeTotalBin=True) -> int:
    #     """Return the total number of bins."""

    #     nCells = self._getNbins("f", includeTotalBin)
    #     nCora = self._getNbins("i", includeTotalBin)
    #     nCorb = self._getNbins("j", includeTotalBin)
    #     nCorc = self._getNbins("k", includeTotalBin)
    #     nDir = self._getNbins("d", includeTotalBin)
    #     nUsr = self._getNbins("u", includeTotalBin)
    #     nSeg = self._getNbins("s", includeTotalBin)
    #     nMul = self._getNbins("m", includeTotalBin)
    #     nCos = self._getNbins("c", includeTotalBin)
    #     nErg = self._getNbins("e", includeTotalBin)
    #     nTim = self._getNbins("t", includeTotalBin)

    #     tot = (nCells * nDir * nUsr * nSeg * nMul * nCos * nErg * nTim *
    #            nCora * nCorb * nCorc)

    #     return tot

    def _insertCell(self, cN: int) -> bool:
        """Insert cell number."""

        if len(self.cells) <= self.nCells:
            self.cells = np.append(self.cells, cN)
            return True
        else:
            return False

    def _insertCorBin(self, axis: str, value) -> bool:
        """Insert cora/b/c values."""
        if axis == "a":
            if len(self.cora) <= self.meshInfo[1] + 1:
                self.cora = np.append(self.cora, value)
                return True
            else:
                return False
        if axis == "b":
            if len(self.corb) <= self.meshInfo[2] + 1:
                self.corb = np.append(self.corb, value)
                return True
            else:
                return False

        if axis == "c":
            if len(self.corc) <= self.meshInfo[3] + 1:
                self.corc = np.append(self.corc, value)
                return True
            else:
                return False

        raise ValueError("Invalid axis")

    def _insertUsr(self, uB: int) -> bool:
        """Insert usr bins."""

        if len(self.usr) <= self.nUsr:
            self.usr = np.append(self.usr, uB)
            return True
        else:
            return False

    def _insertSeg(self, sB: int) -> bool:
        """Insert seg bins."""

        if len(self.seg) <= self.nSeg:
            self.seg = np.append(self.seg, sB)
            return True
        else:
            return False

    def _insertCos(self, cB: int) -> bool:
        """Insert cosine bin."""

        if len(self.cos) <= self.nCos:
            self.cos = np.append(self.cos, cB)
            return True
        else:
            return False

    def _insertRadiograph(self, axis: str, rB: int) -> bool:
        """Insert radiograph coordinates on s and t-axis."""

        if axis == "s":
            if len(self.seg) <= self.nSeg + 1:
                self.seg = np.append(self.seg, rB)
                return True
            else:
                return False

        if axis == "t":
            if len(self.cos) <= self.nCos + 1:
                self.cos = np.append(self.cos, rB)
                return True
        else:
            return False

        raise ValueError("Invalid axis")

    def _insertErg(self, eB: int) -> bool:
        """Insert energy bin."""

        if len(self.erg) <= self.nErg:
            self.erg = np.append(self.erg, eB)
            return True
        else:
            return False

    def _insertTim(self, tB: int) -> bool:
        """Insert time bin."""

        if len(self.tim) <= self.nTim:
            self.tim = np.append(self.tim, tB)
            return True
        else:
            return False

    def _insertTfcJtf(self, jtf: int) -> bool:
        """Insert TFC jtf list."""

        if len(jtf) == 9:
            self.tfc_jtf = jtf
            return True
        else:
            return False

    def _insertTfcDat(self, dat) -> bool:
        """Insert TFC values."""

        if len(dat) <= 4:
            self.tfc_dat.append(dat)
            return True
        else:
            return False

    def _insertValue(self, c, d, u, s, m, a, e, t, f, i, j, k, val) -> None:
        """Insert tally value."""
        if self.isInitialized is False:
            self.initializeValuesVectors()

        self.valsErrors[c][d][u][s][m][a][e][t][f][i][j][k] = val

    def _getValue(self, f, d, u, s, m, c, e, t, i, j, k, v) -> float:
        """Return a value from tally."""
        return self.valsErrors[f][d][u][s][m][c][e][t][i][j][k][v]

    # _getAxis is not called anywhere, check with Davide
    # def _getAxis(self, axis: str):
    #     """Return an array containing the values of the axis bins.
    #     The desired axis is set by passing the
    #     corresponding letter as a function argument as
    #     defined in MCNPX manual (u,s,c,e,t) for the
    #     standard and (i,j,k) for mesh tallies axes (namely
    #     cora/b/c).
    #     """

    #     if axis == "u":
    #         if len(self.usr) != 0:
    #             return np.append([0], self.usr)

    #     if axis == "s":
    #         if len(self.seg) != 0:
    #             if self.radiograph:
    #                 return self.seg
    #             else:
    #                 first = self.seg[0] - 1.
    #                 return np.append([first], self.seg)

    #     if axis == "c":
    #         if len(self.cos) != 0:
    #             if self.radiograph:
    #                 return self.cos
    #             else:
    #                 first = -1.
    #                 return np.append([first], self.cos)

    #     if axis == "e":
    #         if len(self.erg) != 0:
    #             first = 0.0  # self.erg[0] - 1.
    #             return np.append([first], self.erg)

    #     if axis == "t":
    #         if len(self.tim) != 0:
    #             first = self.tim[0] - 1.
    #             return np.append([first], self.tim)

    #     if axis == "i":
    #         return self.cora

    #     if axis == "j":
    #         return self.corb

    #     if axis == "k":
    #         return self.corc

    #     return []

    def _getNbins(self, axis, inclTotBin=True) -> int:
        """Returns the number of bins relative to the desired axis.
        The correspondence is, as usual, (f,d,u,s,m,c,e,t) for standard 8D
        data, plus (i,j,k) for mesh tallies."""

        if axis == "f":
            nCells = 1 if self.nCells == 0 else self.nCells
            return nCells

        if axis == "i":
            return self.meshInfo[1]

        if axis == "j":
            return self.meshInfo[2]

        if axis == "k":
            return self.meshInfo[3]

        if axis == "d":
            nDir = 1 if self.nDir == 0 else self.nDir
            return nDir

        if axis == "u":
            nUsr = 1 if self.nUsr == 0 else self.nUsr
            nUsr = nUsr - 1 if self.usrTC == "t" and not inclTotBin else nUsr
            return nUsr

        if axis == "s":
            nSeg = 1 if self.nSeg == 0 else self.nSeg
            nSeg = nSeg - 1 if self.segTC == "t" and not inclTotBin else nSeg
            return nSeg

        if axis == "m":
            nMul = 1 if self.nMul == 0 else self.nMul
            nMul = nMul - 1 if self.mulTC == "t" and not inclTotBin else nMul
            return nMul

        if axis == "c":
            nCos = 1 if self.nCos == 0 else self.nCos
            nCos = nCos - 1 if self.cosTC == "t" and not inclTotBin else nCos
            return nCos

        if axis == "e":
            nErg = 1 if self.nErg == 0 else self.nErg
            nErg = nErg - 1 if self.ergTC == "t" and not inclTotBin else nErg
            return nErg

        if axis == "t":
            nTim = 1 if self.nTim == 0 else self.nTim
            nTim = nTim - 1 if self.timTC == "t" and not inclTotBin else nTim
            return nTim


class Mctal:
    def __init__(self, filepath: os.PathLike | str) -> None:
        """Object responsible for the parsing of MCNP mctal files.

        Parameters
        ----------
        filepath : os.PathLike | str
            path to the mctal file to be parsed

        Attributes
        ----------
        tallydata : dict[int, pd.DataFrame]
            dictionary that at each tally id associate a pandas dataframe
            containing the results. It supports multi-binning tallies.
        header : Header
            it is the parsed data of the mctal file. See the Header doc to
            understand how to access the data.

        Examples
        --------
        Parse the mctal file and access the data
        >>> # Import the mctal module
        ... from f4enix.output.mctal import Mctal
        ... # Parse the Mctal file
        ... file = 'mctal'
        ... mctal = Mctal(file)
        ... # get a summary of the min and max errors across tallies
        ... mctal.get_error_summary().sort_values(by='tally num')
            tally num	min rel error	max rel error
        0	    4	        NaN	            1.0000
        22	    6	        0.0007	        0.0272
        19	    14	        NaN	            1.0000
        1	    16	        0.0008	        0.0381
        """
        self.tallies = []
        self.thereAreNaNs = False
        self.header = Header()
        self.mctalFileName = os.path.basename(filepath)
        self.mctalFile = open(filepath, "r")

        # This variable will contain the read lines one by one, but it is
        # important to keep it global because the last value from getHeaders()
        # must be already available as first value for parseTally(). This will
        # also apply to successive calls to parseTally().
        self.line = None
        self.tallies = self._read()
        self.tallydata, self.totalbin = self._get_dfs()

    def get_error_summary(self, include_abs_err: bool = False) -> pd.DataFrame:
        """Return a dataframe containing a summary of the min and max errror
        registered in each tally. If both value and error are equal to zero,
        the errors will be set to NaN, since it means that nothing has been
        scored in the tally.

        Parameters
        ----------
        include_abs_err : bool, optional
            if True includes the absolute error in addition to the total one,
            by default is False

        Returns
        -------
        pd.DataFrame
            error summary
        """
        rows = []
        for tally, data in self.tallydata.items():
            min_error = data["Error"].min()
            min_idx = data["Error"].idxmin()
            max_error = data["Error"].max()
            max_idx = data["Error"].idxmax()

            min_val = data["Value"].iloc[min_idx]
            max_val = data["Value"].iloc[max_idx]

            # put a NaN, since no particle was tallied in the cell
            if min_error == 0 and min_val == 0:
                min_error = np.nan
            if max_error == 0 and max_val == 0:
                max_error = np.nan

            if include_abs_err:
                rows.append([tally, min_error, min_val, max_error, max_val])
            else:
                rows.append([tally, min_error, max_error])

        if include_abs_err:
            columns = [
                "tally num",
                "min rel error",
                "min abs err",
                "max rel error",
                "max abs err",
            ]
        else:
            columns = ["tally num", "min rel error", "max rel error"]

        df = pd.DataFrame(rows)
        df.columns = columns

        return df

    def _read(self) -> list[Tally]:
        """This function calls the functions getHeaders and parseTally in
        order to read the entier MCTAL file."""

        logging.info("\n\033[1;34m[Parsing file: %s...]\033[0m" % self.mctalFileName)

        self._getHeaders()
        self._getTallies()

        if self.thereAreNaNs:
            logging.info(
                "\n \033[1;30mThe MCTAL file contains one or more tallies with NaN values. Flagged.\033[0m\n",
                file=sys.stderr,
            )

        self.mctalFile.close()

        return self.tallies

    def _getHeaders(self) -> None:
        """This function reads the first lines from the MCTAL file. We call
        "header" what is written from the beginning to the first "tally"
        keyword."""

        self.line = self.mctalFile.readline().split()

        if len(self.line) == 7:
            self.header.kod = self.line[0]
            self.header.ver = self.line[1]
            pID_date = self.line[2]
            self.header.probid = np.append(self.header.probid, pID_date)
            pID_time = self.line[3]
            self.header.probid = np.append(self.header.probid, pID_time)
            self.header.knod = int(self.line[4])
            self.header.nps = int(self.line[5])
            self.header.rnr = int(self.line[6])

        elif len(self.line) == 3:
            self.header.knod = int(self.line[0])
            self.header.nps = int(self.line[1])
            self.header.rnr = int(self.line[2])

        self.header.title = self.mctalFile.readline().strip()

        self.line = self.mctalFile.readline().split()

        self.header.ntal = int(self.line[1])

        if self.header.ntal == 0:
            logging.error(
                "\n \033[1;31mNo tallies in this MCTAL file. Exiting.\033[0m\n",
                file=sys.stderr,
            )
            sys.exit(1)

        if len(self.line) == 4:
            self.header.npert = int(self.line[3])
            logging.error(
                "\n \033[1;31mMCTAL file with perturbation card. Not supported. Exiting.\033[0m\n",
                file=sys.stderr,
            )
            sys.exit(1)

        self.line = self.mctalFile.readline().split()

        while self.line[0].lower() != "tally":
            for char in self.line:
                self.header.ntals = np.append(self.header.ntals, int(char))
            self.line = self.mctalFile.readline().split()

    def _getTallies(self) -> None:
        """This function supervises the calls to parseTally() function.
        It will keep track of the position of cursor in the MCTAL file and stop execution when EOF is reached.
        """

        EOF = False
        while not EOF:
            EOF = self._parseTally()

    def _parseTally(self) -> bool:
        """This function parses an entire tally."""
        # TODO the linting is horrible here

        # The first line processed by this function is already in memory, either coming from the
        # last readline() in Header class or from the previous call to parseTally()

        tally = Tally(int(self.line[1]))

        logging.info(" \033[33mParsing tally: %5d\033[0m" % (tally.tallyNumber))

        tally.typeNumber = int(self.line[2])
        if self.line[3] != 0:
            tally.detectorType = int(self.line[3])

        if (
            tally.detectorType is not None
        ):  # check for None is needed for MCNP6 and F1 tally
            if tally.detectorType >= 3:
                tally.radiograph = True
            elif tally.detectorType <= -1:
                tally.mesh = True

        self.line = self.mctalFile.readline()
        # Some MCTAL files do not have the particle list
        # (e.g. produced by MCNP5) so we need to check if it's
        # present (line with particle list starts with space):
        if self.line[0] == " " and self.line[1] != " ":
            line = self.line.split()
            for p in line:
                if p.isdigit():
                    v = int(p)
                    tally.tallyParticles = np.append(tally.tallyParticles, v)
                else:
                    print(
                        "\n \033[1;31m Problem with particle list in tally %d in %s -> exiting.\033[0m\n"
                        % (tally.tallyNumber, self.mctalFileName),
                        file=sys.stderr,
                    )
                    print(self.line)
                    sys.exit(1)

        if len(tally.tallyParticles) != 0:
            self.line = self.mctalFile.readline()

        while self.line[0:5] == " " * 5 and self.line[0].lower() != "f":
            tally.tallyComment = np.append(tally.tallyComment, self.line[5:].rstrip())
            self.line = self.mctalFile.readline()

        line = self.line.split()

        if not tally.mesh:
            tally.nCells = int(line[1])
        else:
            tally.nCells = 1
            tally.meshInfo[0] = int(line[2])  # Unknown number
            tally.meshInfo[1] = int(line[3])  # number of cora bins
            tally.meshInfo[2] = int(line[4])  # number of corb bins
            tally.meshInfo[3] = int(line[5])  # number of corc bins

        # Fix for detector tallies
        if str(tally.tallyNumber)[-1] == "5":
            for k in range(tally.nCells):
                tally.cells = np.append(tally.cells, k + 1)

        self.line = self.mctalFile.readline()

        i = 0
        axisNumber = 0
        axisName = ("a", "b", "c")
        corsVals = (tally.meshInfo[1], tally.meshInfo[2], tally.meshInfo[3])

        while self.line[0].lower() != "d":  # CELLS
            if tally.mesh:
                for c in self.line.split():
                    if not tally._insertCorBin(axisName[axisNumber], float(c)):
                        raise IOError(
                            "Too many cells in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
                    i = i + 1
                    if i == (corsVals[axisNumber] + 1):
                        axisNumber = axisNumber + 1
                        i = 0
            elif "." in self.line and "E" not in self.line:
                for c in self.line.split():
                    if not tally._insertCell(float(c)):
                        raise IOError(
                            "Too many cells in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
            else:
                for c in self.line.split():
                    if not tally._insertCell(
                        int(c)
                    ):  # This means that for some reason you are trying to
                        # insert more cells than the number stated in f
                        raise IOError(
                            "Too many cells in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )

            self.line = self.mctalFile.readline()

        if tally.mesh is True:
            tally.cora = tally.cora[:-1]
            tally.corb = tally.corb[:-1]
            tally.corc = tally.corc[:-1]

        # DIR
        self.line = self.line.split()
        tally.nDir = int(self.line[1])

        while self.line[0].lower() != "u":
            self.line = self.mctalFile.readline()

        # USR
        self.line = self.line.split()
        if self.line[0].lower() == "ut":
            tally.usrTC = "t"
        if self.line[0].lower() == "uc":
            tally.usrTC = "c"
        tally.nUsr = int(self.line[1])

        # USR BINS
        self.line = self.mctalFile.readline()
        while self.line[0].lower() != "s":
            for u in self.line.split():
                if not tally._insertUsr(float(u)):
                    raise IOError(
                        "Too many user bins in the tally n. %d of %s"
                        % (tally.tallyNumber, self.mctalFileName)
                    )
            self.line = self.mctalFile.readline()

        # SEG
        self.line = self.line.split()
        if self.line[0].lower() == "st":
            tally.segTC = "t"
        if self.line[0].lower() == "sc":
            tally.segTC = "c"
        tally.nSeg = int(self.line[1])

        # SEGMENT BINS
        self.line = self.mctalFile.readline()
        while self.line[0].lower() != "m":
            for s in self.line.split():
                if tally.radiograph:
                    if not tally._insertRadiograph("s", float(s)):
                        raise IOError(
                            "Too many segment bins in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
                else:
                    if not tally._insertSeg(float(s)):
                        raise IOError(
                            "Too many segment bins in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
            self.line = self.mctalFile.readline()

        # MUL
        self.line = self.line.split()
        if self.line[0].lower() == "mt":
            tally.mulTC = "t"
        if self.line[0].lower() == "mc":
            tally.mulTC = "c"
        tally.nMul = int(self.line[1])

        while self.line[0].lower() != "c":
            self.line = self.mctalFile.readline()

        # COS
        self.line = self.line.split()
        if self.line[0].lower() == "ct":
            tally.cosTC = "t"
        if self.line[0].lower() == "cc":
            tally.cosTC = "c"
        tally.nCos = int(self.line[1])
        if len(self.line) == 3:
            tally.cosFlag = int(self.line[2])

        # COSINE BINS
        self.line = self.mctalFile.readline()
        while self.line[0].lower() != "e":
            for c in self.line.split():
                if tally.radiograph:
                    if not tally._insertRadiograph("t", float(c)):
                        raise IOError(
                            "Too many cosine bins in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
                else:
                    if not tally._insertCos(float(c)):
                        raise IOError(
                            "Too many cosine bins in the tally n. %d of %s"
                            % (tally.tallyNumber, self.mctalFileName)
                        )
            self.line = self.mctalFile.readline()

        if tally.radiograph is True:
            tally.seg = tally.seg[:-1]
            tally.seg = tally.seg[:-1]
            tally.seg = tally.seg[:-1]
            tally.cos = tally.cos[:-1]
            tally.cos = tally.cos[:-1]
            tally.cos = tally.cos[:-1]
        # ERG
        self.line = self.line.split()
        if self.line[0].lower() == "et":
            tally.ergTC = "t"
        if self.line[0].lower() == "ec":
            tally.ergTC = "c"
        tally.nErg = int(self.line[1])
        if len(self.line) == 3:
            tally.ergFlag = int(self.line[2])

        # ENERGY BINS
        self.line = self.mctalFile.readline()
        while self.line[0].lower() != "t":
            for e in self.line.split():
                if not tally._insertErg(float(e)):
                    raise IOError(
                        "Too many energy bins in the tally n. %d of %s"
                        % (tally.tallyNumber, self.mctalFileName)
                    )
            self.line = self.mctalFile.readline()

        # TIM
        self.line = self.line.split()
        if self.line[0].lower() == "tt":
            tally.timTC = "t"
        if self.line[0].lower() == "tc":
            tally.timTC = "c"
        tally.nTim = int(self.line[1])
        if len(self.line) == 3:
            tally.timFlag = int(self.line[2])

        # TIME BINS
        self.line = self.mctalFile.readline()
        while self.line.strip().lower() != "vals":
            for t in self.line.split():
                if not tally._insertTim(float(t)):
                    raise IOError(
                        "Too many time bins in the tally n. %d of %s"
                        % (tally.tallyNumber, self.mctalFileName)
                    )
            self.line = self.mctalFile.readline()

        # VALS
        f = 1
        Fld = []
        nFld = 0
        tally._initializeValuesVectors()

        nCells = tally._getNbins("f")
        nCora = tally._getNbins("i")
        nCorb = tally._getNbins("j")
        nCorc = tally._getNbins("k")
        nDir = tally._getNbins("d")
        nUsr = tally._getNbins("u")
        nSeg = tally._getNbins("s")
        nMul = tally._getNbins("m")
        nCos = tally._getNbins("c")
        nErg = tally._getNbins("e")
        nTim = tally._getNbins("t")

        for c in range(nCells):
            for d in range(nDir):
                for u in range(nUsr):
                    for s in range(nSeg):
                        for m in range(nMul):
                            for a in range(nCos):  # a is for Angle...forgive me
                                for e in range(nErg):
                                    for t in range(nTim):
                                        for k in range(nCorc):
                                            for j in range(nCorb):
                                                for i in range(nCora):
                                                    if (
                                                        f > nFld
                                                    ):  # f is for Field...again, forgive me
                                                        del Fld
                                                        del self.line
                                                        self.line = self.mctalFile.readline().strip()
                                                        Fld = self.line.split()
                                                        nFld = len(Fld) - 1
                                                        f = 0

                                                    if self.line[0:3] != "tfc":
                                                        # This needs to handle the bug like '8.23798-100'
                                                        try:
                                                            val = float(Fld[f])
                                                            err = float(Fld[f + 1])
                                                        except ValueError:
                                                            val = 0
                                                            err = 0

                                                        if math.isnan(
                                                            val
                                                        ) or math.isnan(err):
                                                            self.thereAreNaNs = True
                                                        tally._insertValue(
                                                            c,
                                                            d,
                                                            u,
                                                            s,
                                                            m,
                                                            a,
                                                            e,
                                                            t,
                                                            i,
                                                            j,
                                                            k,
                                                            0,
                                                            val,
                                                        )
                                                        tally._insertValue(
                                                            c,
                                                            d,
                                                            u,
                                                            s,
                                                            m,
                                                            a,
                                                            e,
                                                            t,
                                                            i,
                                                            j,
                                                            k,
                                                            1,
                                                            err,
                                                        )

                                                        f += 2

        del Fld

        if tally.mesh == False:
            # TFC JTF
            self.line = self.mctalFile.readline().strip().split()
            if self.line[0] != "tfc":
                raise IOError(
                    "There seem to be more values than expected in tally n. %d of %s"
                    % (tally.tallyNumber, self.mctalFileName)
                )

            del self.line[0]
            self.line = [int(i) for i in self.line]

            if not tally._insertTfcJtf(self.line):
                raise IOError("Wrong number of TFC jtf elements.")

            # TFC DAT
            self.line = self.mctalFile.readline().strip()
            while "tally" not in self.line and len(self.line) != 0:
                self.line = self.line.split()

                tfcDat = []

                tfcDat.append(int(self.line[0]))
                try:
                    val1 = float(self.line[1])
                except ValueError:
                    val1 = 0

                if math.isnan(val1) or math.isnan(float(self.line[2])):
                    self.thereAreNaNs = True

                tfcDat.append(val1)
                tfcDat.append(float(self.line[2]))
                if len(self.line) == 4:
                    if math.isnan(val1):
                        self.thereAreNaNs = True
                    tfcDat.append(float(self.line[3]))

                if not tally._insertTfcDat(tfcDat):
                    raise IOError(
                        "Wrong number of elements in TFC data line in the tally n. %d of %s"
                        % (tally.tallyNumber, self.mctalFileName)
                    )

                self.line = self.mctalFile.readline().strip()

        else:
            while "tally" not in self.line and len(self.line) != 0:
                self.line = self.mctalFile.readline().strip()

        self.tallies.append(tally)

        if self.line == "":
            self.line = self.line
            return True
        elif "tally" in self.line:
            self.line = self.line.split()
            return False

    def _get_dfs(
        self, collapse: bool = False
    ) -> tuple[dict[int, pd.DataFrame], dict[int, pd.DataFrame]]:
        """
        Retrieve and organize mctal data into a DataFrame.

        Parameters
        ----------
        collapse : bool
            collapse the Cell and segments binning in a single Cell-Segment
            one

        Returns
        -------
        tallydata : dict[int, pd.DataFrame]
            organized tally data.
        totalbin : dict[int, pd.DataFrame]
            organized tally data (only total bins).

        """
        tallydata = {}
        totalbin = {}

        for t in self.tallies:
            rows = []

            # --- Reorganize values ---
            # You cannot recover the following from the mctal
            # nDir = t._getNbins("d", False)
            nMul = t._getNbins("m", False)
            nSeg = t._getNbins("s", False)  # this can be used

            # Some checks for voids
            binnings = {
                "cells": t.cells,
                "user": t.usr,
                "dir": range(t.nDir),
                "segments": t.seg,
                "cosine": t.cos,
                "energy": t.erg,
                "time": t.tim,
                "cor A": t.cora,
                "cor B": t.corb,
                "cor C": t.corc,
            }

            # Cells may have a series of zeros, the last one may be for the
            # total
            cells = []
            # last_idx = len(binnings['cells'])-1
            for i, cell in enumerate(binnings["cells"]):
                if int(cell) == 0:
                    newval = "Input " + str(i + 1)
                    cells.append(newval)
                # Everything is fine, nothing to modify
                else:
                    # force it to be an int
                    cells.append(int(cell))
            binnings["cells"] = cells

            for name, binning in binnings.items():
                if len(binning) == 0:
                    binnings[name] = [np.nan]
            # Start iteration
            for f, fn in enumerate(binnings["cells"]):
                for d, dn in enumerate(binnings["dir"]):  # Unused
                    for u, un in enumerate(binnings["user"]):
                        for sn in range(1, nSeg + 1):
                            for m in range(nMul):  # (unused)
                                for c, cn in enumerate(binnings["cosine"]):
                                    for e, en in enumerate(binnings["energy"]):
                                        for nt, ntn in enumerate(binnings["time"]):
                                            for k, kn in enumerate(binnings["cor C"]):
                                                for j, jn in enumerate(
                                                    binnings["cor B"]
                                                ):
                                                    for i, ina in enumerate(
                                                        binnings["cor A"]
                                                    ):
                                                        val = t._getValue(
                                                            f,
                                                            d,
                                                            u,
                                                            sn - 1,
                                                            m,
                                                            c,
                                                            e,
                                                            nt,
                                                            i,
                                                            j,
                                                            k,
                                                            0,
                                                        )
                                                        err = t._getValue(
                                                            f,
                                                            d,
                                                            u,
                                                            sn - 1,
                                                            m,
                                                            c,
                                                            e,
                                                            nt,
                                                            i,
                                                            j,
                                                            k,
                                                            1,
                                                        )
                                                        rows.append(
                                                            [
                                                                fn,
                                                                dn,
                                                                un,
                                                                sn,
                                                                m,
                                                                cn,
                                                                en,
                                                                ntn,
                                                                ina,
                                                                jn,
                                                                kn,
                                                                val,
                                                                err,
                                                            ]
                                                        )

                # Only one total bin per cell is admitted
                val = t._getValue(f, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0)
                err = t._getValue(f, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1)
                if t.timTC is not None:
                    rows.append(
                        [fn, dn, un, sn, m, cn, en, "total", ina, jn, kn, val, err]
                    )
                    total = "Time"

                elif t.ergTC is not None:
                    rows.append(
                        [fn, dn, un, sn, m, cn, "total", ntn, ina, jn, kn, val, err]
                    )
                    total = "Energy"

                elif t.segTC is not None:
                    rows.append(
                        [fn, dn, un, "total", m, cn, en, ntn, ina, jn, kn, val, err]
                    )
                    total = "Segments"

                elif t.cosTC is not None:
                    rows.append(
                        [fn, dn, un, sn, m, "total", en, ntn, ina, jn, kn, val, err]
                    )
                    total = "Cosine"

                elif t.usrTC is not None:
                    rows.append(
                        [fn, dn, "total", sn, m, cn, en, ntn, ina, jn, kn, val, err]
                    )
                    total = "User"

            # --- Build the tally DataFrame ---
            columns = [
                "Cells",
                "Dir",
                "User",
                "Segments",
                "Multiplier",
                "Cosine",
                "Energy",
                "Time",
                "Cor C",
                "Cor B",
                "Cor A",
                "Value",
                "Error",
            ]
            df = pd.DataFrame(rows, columns=columns)

            # Default drop of multiplier and Dir
            # del df["Dir"]
            del df["Multiplier"]
            # --- Keep only meaningful binning ---
            # Drop NA
            df.dropna(axis=1, inplace=True)
            # Drop constant axes (if len is > 1)
            if len(df) > 1:
                for column in df.columns:
                    if column not in ["Value", "Error"]:
                        firstval = df[column].values[0]
                        # Should work as long as they are the exact same value
                        allequal = (df[column] == firstval).all()
                        if allequal:
                            del df[column]

            # Drop rows if they are exactly the same values
            # (untraced behaviour)
            df.drop_duplicates(inplace=True)

            # The double binning Surfaces/cells with segments can create
            # issues for JADE since if another binning is added
            # (such as energy) it is not supported. Nevertheless,
            # the additional segmentation can be quite useful and this can be
            # collapsed de facto in a single geometrical binning

            if (
                "Cells" in df.columns
                and "Segments" in df.columns
                and collapse
                and len(df)
            ) > 1:
                # Then we can collapse this in a single geometrical binning
                values = []
                for cell, segment in zip(df.Cells, df.Segments):
                    try:
                        val = str(int(cell)) + "-" + str(int(segment))
                    except ValueError:
                        # the cell or segment may not be a number
                        val = "{}-{}".format(cell, segment)
                    values.append(val)
                df["Cells-Segments"] = values
                # delete the collapsed columns
                del df["Cells"]
                del df["Segments"]

            # Sub DF containing only total bins
            try:
                dftotal = df[df[total] == "total"]
            except (KeyError, NameError):
                # KeyError : there is no total bin in df
                # NameError: total variable was not defined
                dftotal = None

            tallydata[t.tallyNumber] = df
            totalbin[t.tallyNumber] = dftotal

        return tallydata, totalbin
