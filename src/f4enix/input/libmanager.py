from __future__ import annotations

# -*- coding: utf-8 -*-
"""
This modules is related to handling nuclear data libraries.

It provides support for zaid and elements data and it used for different
operation such as the translation from a library to another of an MCNP input

"""

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


import json
import logging
import os
import re
import warnings

import numpy as np
import pandas as pd

import pandas as pd

from importlib.resources import files, as_file

from f4enix.input.xsdirpyne import Xsdir
import f4enix.resources as pkg_res
import f4enix.input.xsdirpyne as xs
from f4enix.input.xsdirpyne import Xsdir, OpenMCXsdir, SerpentXsdir
from pathlib import Path

XSDIR_CLASS = {
    "mcnp": Xsdir,
    "openmc": OpenMCXsdir,
    "serpent": SerpentXsdir,
    "d1s": Xsdir,
}
MCNP_TYPE_XSDIR = ["mcnp", "d1s", "serpent"]

# colors
CRED = "\033[91m"
CEND = "\033[0m"


MSG_DEFLIB = " The Default library {} was used for zaid {}"


class IsotopeDataParser:
    def __init__(self, isotopes_file: os.PathLike) -> None:
        # load the natural abundance file
        with isotopes_file as iso:
            abundances = pd.read_csv(iso, skiprows=2)
            abundances["idx"] = abundances["idx"].astype(str)
            abundances.set_index("idx", inplace=True)
            self.isotopes = abundances

    def get_formulazaid(self, formula):
        match = re.match(r"([a-z]+)([0-9]+)", formula, re.I)
        parts = match.groups()
        E, A = parts[0], int(parts[1])
        newiso = self.isotopes[self.isotopes["E"] == E]
        Z = newiso["Z"].values[0]
        zaid = "{0}{1:0>3}".format(Z, A)
        return zaid


class LibManager:

    # def __init__(self, xsdir_file, defaultlib='81c', activationfile=None,
    #             isotopes_file=None):
    def __init__(
        self,
        xsdir_path: os.PathLike | pd.DataFrame = None,
        defaultlib: str = None,
        activationfile: os.PathLike = None,
        isotopes_file: os.PathLike = None,
    ) -> None:
        """
        Object dealing with all complex operations that involves nuclear data.

        Parameters
        ----------
        xsdir_path : pd.DataFrame
            table related to libraries variables.
        defaultlib : str, optional
            lib suffix to be used as default in translation operations.
            If None, it reads from lib_df
        activationfile : str or path, optional
            path to the config file containing the reactions data for
            activation libraries. The default is None.
        isotopes_file : str or path, optional
            path to the isotopes files. If None (default) the file is searched
            in the current directory.

        Attributes
        ----------
        isotope_parser : IsotopeDataParser
            object dealing with the isotopes data.
        isotopes : pd.DataFrame
            contains the isotopes data.
        defaultlib : str
            lib suffix to be used as default in translation operations.
        data : dict[str, dict[str, Union[Xsdir, OpenMCXsdir, SerpentXsdir]]]
            contains the libraries data. first level keys are the codes, second
            level keys are the library suffixes. ultimate value is the xsdir
            object.
        codes : list
            list of codes available.
        libraries : dict[str, list[str]]
            contains the libraries available for each code.
        reactions : dict[str, pd.DataFrame]
            contains the reactions data for the different activation libraries.

        Returns
        -------
        None.

        """
        if xsdir_path is None:
            datapath = (
                os.path.join(os.getenv("DATAPATH", ""), "xsdir_mcnp6.2")
                if os.getenv("DATAPATH")
                else None
            )
            xsdir_path = Path(datapath)

        if xsdir_path is None:
            resources = files(pkg_res)
            xsdir_path = as_file(resources.joinpath("xsdir.txt"))
            xsdir = Xsdir(xsdir_path)

        try:
            xsdir_path = Path(xsdir_path)
        except TypeError:
            pass

        if isinstance(xsdir_path, os.PathLike):
            xsdir = Xsdir(xsdir_path)
            available_libs = list(set(np.array(xsdir.tablenames)[:, 1]))
            df_rows = []
            for lib in available_libs:
                df_rows.append([lib, "", "", xsdir_path])
            df_lib = pd.DataFrame(df_rows)
            df_lib.columns = ["Suffix", "Name", "Default", "MCNP"]
            xsdir_path = df_lib

        if isotopes_file is None:
            resources = files(pkg_res)
            iso_file = as_file(resources.joinpath("Isotopes.txt"))
        else:
            iso_file = as_file(Path(isotopes_file))

        self.isotope_parser = IsotopeDataParser(iso_file)
        self.isotopes = self.isotope_parser.isotopes

        # Convert all columns to lower case
        new_columns = []
        for column in xsdir_path:
            new_columns.append(column.lower())
        xsdir_path.columns = new_columns

        if defaultlib is None:
            try:
                self.defaultlib = xsdir_path[xsdir_path["default"] == "yes"][
                    "suffix"
                ].values[0]
            except IndexError:
                self.defaultlib = "81c"
        else:
            self.defaultlib = defaultlib

        self.data = {}
        self.codes = []
        xsdir_path.set_index("suffix", inplace=True)
        # Initilize the Xsdir object
        # self.XS = xs.Xsdir(xsdir_file)

        # this block of code needs to check the availability of the libraries.
        # Only libraries specified in the config file are checked, if paths
        # for the libraries are left empty, the library is not not checked and
        # it is not registered as available. If the path is not empty but
        # library is not found, a warning is raised, choice for interrupting the
        # session is left to the user.
        for code in xsdir_path.columns[2:]:
            code = code.lower()
            self.codes.append(code)
            self.data[code] = {}
            for library, row in xsdir_path.iterrows():
                path = row[code]
                # if the path is empty just ignore it
                if path is None or path == "":
                    logging.info("No path for %s library", library)
                    continue

                # if the path is not empty, check if the file exists
                # and if it does not, raise a warning since it may not be the
                # intended behaviour by the user
                if not os.path.exists(path):
                    logging.warning(
                        "Library %s for code %s not found at %s", library, code, path
                    )
                    # fatal_exception(path + " does not exist")

                if code in MCNP_TYPE_XSDIR:
                    xsdir = XSDIR_CLASS[code](path)
                    # verify that the library is actually in the xsdir
                    available_libs = set(np.array(xsdir.tablenames)[:, 1])
                    if library in available_libs:
                        self.data[code][library] = xsdir
                    else:
                        logging.warning(
                            "Library %s not present in MCNP XSDIR file: %s",
                            library,
                            path,
                        )

                elif code == "openmc":
                    self.data[code][library] = OpenMCXsdir(path, self, library)

                else:
                    raise ValueError(f"{code} code not implemented")

        # Identify different libraries installed. This is done checking H
        # libraries = self.check4zaid('1001')
        # libraries.extend(self.check4zaid('1000'))  # photons

        # """ Legacy library definition changed """
        # """
        # libraries = []
        # for table in self.XS:
        #     lib = table.name.split('.')[1]
        #     if lib not in libraries:
        #         libraries.append(lib)

        # self.libraries = libraries
        # """

        # libraries have now been checked at the source, they may be different
        # for each code
        libraries = {}
        for key, value in self.data.items():
            libraries[key] = []
            for lib, _ in value.items():
                libraries[key].append(lib)
        self.libraries = libraries

        # Load the activation reaction data if available
        reactions = {}
        if activationfile is not None:
            with as_file(Path(activationfile)) as infile:
                file = pd.ExcelFile(infile)
        else:
            resources = files(pkg_res)
            with as_file(resources.joinpath("activation_libs.xlsx")) as infile:
                file = pd.ExcelFile(infile)
        for sheet in file.sheet_names:
            # Load the df that also needs to be filled
            reactions[sheet] = file.parse(sheet).ffill()
            # translate the formula name to zaid

        # These are needed for faster operations
        newiso = self.isotopes.set_index(["E"])
        newiso = newiso.loc[~newiso.index.duplicated(keep="first")]
        self._newiso_byE = newiso.sort_index()

        newiso = self.isotopes.set_index(["Z"])
        newiso = newiso.loc[~newiso.index.duplicated(keep="first")]
        self._newiso_byZ = newiso.sort_index()

        self.reactions = reactions

    def check4zaid(self, zaid: str, code: str = "mcnp"):
        # Needs fixing
        """
        Check which libraries are available for the selected zaid and return it

        Parameters
        ----------
        zaid : str
            zaid string (e.g. 1001).
        code: str, optional
            Code to be looked up. Default is 'mcnp'

        Returns
        -------
        libraries : list
            list of libraries available for the zaid.

        """
        libraries = []
        if code in MCNP_TYPE_XSDIR:
            for lib in self.libraries[code]:
                xsdir = self.data[code][lib]
                if lib in xsdir.find_table(zaid, mode="default-fast"):
                    libraries.append(lib)
        else:
            raise NotImplementedError("{} not implemented yet".format(code))

        return libraries

    # def check_zaid(self, zaid, lib, code):
    #     XS = self._get_xs(code, lib)
    #     if isinstance(XS, str) or XS is None:
    #         return True
    #     elif len(XS.find_table(zaid, mode='default-fast')) > 0:
    #         return True
    #     else:
    #         return False

    def convertZaid(
        self, zaid: str, lib: str, code: str = "mcnp"
    ) -> dict[str, tuple[str, float, float]]:
        # Needs fixing
        """
        This methods will convert a zaid into the requested library

        modes:
            - 1to1: there is one to one correspondence for the zaid
            - natural: the zaids will be expanded using the natural abundance
            - absent: the zaid is not available in the library, a default one
              will be used or the natural one if available.

        Parameters
        ----------
        zaid : str
            zaid name (ex. 1001).
        lib : str
            library suffix (ex. 21c).
        code : str, optional
            code for which the translation is performed. default is MCNP

        Raises
        ------
        ValueError
            if the library is not available in the xsdir file or if there is
            no valid translation for the zaid.

        Returns
        -------
        translation : dict[str, tuple[str, float, float]]
            {zaidname:(lib,nat_abundance,Atomic mass)}.

        """
        # Check if library is available in Xsdir
        if lib not in self.libraries[code]:
            raise ValueError("Library " + lib + " is not available in xsdir file")

        zaidlibs = self.check4zaid(zaid, code)

        if code in MCNP_TYPE_XSDIR:
            XS = self.data[code][lib]
            # Natural zaid
            if zaid[-3:] == "000":
                # Check if zaid has natural info
                if XS.find_table(zaid + "." + lib, mode="exact"):
                    translation = {zaid: (lib, 1, 1)}  # mass not important

                else:  # Has to be expanded
                    translation = {}
                    reduced = self.isotopes[self.isotopes["Z"] == int(zaid[:-3])]
                    for idx, row in reduced.iterrows():
                        # zaid availability must be checked
                        if XS.find_table(idx + "." + lib, mode="exact"):
                            newlib = lib
                        elif self.data[code][self.defaultlib].find_table(
                            idx + "." + self.defaultlib, mode="exact"
                        ):
                            warnings.warn(MSG_DEFLIB.format(self.defaultlib, zaid))
                            newlib = self.defaultlib
                        else:
                            raise ValueError(
                                "No available translation for zaid :"
                                + zaid
                                + "It is needed for natural zaid expansion."
                            )

                        translation[idx] = (
                            newlib,
                            row["Mean value"],
                            row["Atomic Mass"],
                        )
            # 1to1
            elif lib in zaidlibs:
                translation = {zaid: (lib, 1, 1)}  # mass not important

            # No possible correspondence, natural or default lib has to be used
            else:
                # Check if the natural zaid is available
                natzaid = zaid[:-3] + "000"
                if XS.find_table(natzaid + "." + lib, mode="exact"):
                    translation = {natzaid: (lib, 1, 1)}  # mass not important
                # Check if default lib is available
                elif self.data[code][self.defaultlib].find_table(
                    zaid + "." + self.defaultlib, mode="exact"
                ):
                    warnings.warn(MSG_DEFLIB.format(self.defaultlib, zaid))
                    translation = {zaid: (self.defaultlib, 1, 1)}  # mass not imp
                else:
                    # Check if any zaid cross section is available
                    libraries = self.check4zaid(zaid, code)
                    # It has to be for the same type of particles
                    found = False
                    for library in libraries:
                        if library[-1] == lib[-1]:
                            found = True
                    # If found no lib is assigned
                    if found:
                        translation = {zaid: (None, 1, 1)}  # no masses
                    # If no possible translation is found raise error
                    else:
                        raise ValueError("No available translation for zaid :" + zaid)
        else:
            raise ValueError("Translation not required for code " + code)

        return translation

    def get_libzaids(self, lib: str, code: str = "mcnp"):
        # Needs fixing
        """
        Given a library, returns all zaids available

        Parameters
        ----------
        lib : str
            suffix of the library.
        code : str, optional
            code for which the zaids are recovered. default is MCNP

        Returns
        -------
        zaids : list[str]
            list of zaid names available in the library.

        """
        try:
            XS = self.data[code][lib]
        except KeyError:
            return []

        zaids = []

        if isinstance(XS, xs.Xsdir):
            for table in XS.find_zaids(lib):
                zaid = table.name.split(".")[0]
                if zaid not in zaids:
                    zaids.append(zaid)
        else:
            raise NotImplementedError("{} code is not yet implemented".format(code))

        return zaids

    def get_formulazaid(self, formula):
        return self.isotope_parser.get_formulazaid(formula)

    def get_zaidname(self, zaid: str) -> tuple[str, str]:
        """
        Given a zaid, its element name and formula are returned. E.g.,
        hydrogen, H1

        Parameters
        ----------
        zaid : str
            zaid number (e.g. 1001 for H1).

        Returns
        -------
        name : str
            element name (e.g. hydrogen).
        formula : str
            isotope name (e.g. H1).

        """
        if type(zaid) == str:
            splitted = zaid.split(".")
            elem = splitted[0][:-3]
            i = int(elem)
            isotope = splitted[0][-3:]

        else:
            i = int(zaid.element)
            isotope = zaid.isotope

        # newiso = self.isotopes.set_index("Z")
        # newiso = newiso.loc[~newiso.index.duplicated(keep="first")]

        name = self._newiso_byZ["Element"].loc[i]

        if int(isotope) > 0:
            formula = self._newiso_byZ["E"].loc[i] + "-" + str(int(isotope))
        else:
            formula = self._newiso_byZ["E"].loc[i]

        return name, formula

    def get_zaidnum(self, zaidformula: str) -> str:
        """
        Given a zaid formula return the correct number

        Parameters
        ----------
        zaidformula : str
            name of the zaid, e.g., H1.

        Returns
        -------
        zaidnum : str
            number of the zaid ZZZAA

        """
        # get the table and drop the duplicates
        # newiso = self.isotopes.set_index(["E"])
        # newiso = newiso.loc[~newiso.index.duplicated(keep="first")]
        # split the name
        patnum = re.compile(r"\d+")
        patname = re.compile(r"[a-zA-Z]+")
        try:
            num = patnum.search(zaidformula).group()
            name = patname.search(zaidformula).group()
        except AttributeError:
            raise ValueError("No correspondent zaid found for " + zaidformula)

        atomnumber = self._newiso_byE.loc[name, "Z"]

        zaidnum = "{}{:03d}".format(atomnumber, int(num))

        return zaidnum

    def select_lib(self, codes: list[str] = ["mcnp"]) -> str:
        """
        Prompt a library input selection with Xsdir availabilty check

        Parameters
        ----------
        code: list[str], optional
            code for which the library is selected. default is MCNP

        Returns
        -------
        lib : str
            Library to assess.

        """
        error = (
            CRED
            + """
 Error: {}
 The selected library is not available.
 """
            + CEND
        )
        # Add a counter to avoid falling in an endless loop
        i = 0
        while True:
            i += 1
            lib = input(" Select library (e.g. 31c or 99c-31c): ")
            # check that library is available in all requested codes

            if self._is_lib_available(lib, codes):
                break  # if the library is available for all codes, break loop

            elif lib[0] == "{":
                libs = json.loads(lib)
                # all libraries should be available
                tocheck = list(libs.values())
                tocheck.extend(list(libs.keys()))
                flag = True
                for val in tocheck:
                    if not self._is_lib_available(val, codes):
                        print(error.format(val))
                        flag = False
                if flag:
                    break

            elif "-" in lib:
                libs = lib.split("-")
                flag = True
                for val in libs:
                    if not self._is_lib_available(val, codes):
                        print(error.format(val))
                        flag = False
                if flag:
                    break

            elif lib == "back":
                break

            elif lib == "exit":
                break

            else:
                print(error.format(lib))

            if i > 20:
                raise ValueError("Too many wrong inputs")
        return lib

    def _is_lib_available(self, lib: str, codes: list[str] = ["mcnp"]) -> bool:
        flag_present = True
        for code in codes:
            if lib not in self.libraries[code]:
                flag_present = False
        return flag_present

    def get_zaid_mass(self, zaid: str) -> float:
        """
        Get the atomic mass of one zaid

        Parameters
        ----------
        zaid : matreader.Zaid
            Zaid to examinate.

        Returns
        -------
        m: float
            atomic mass.

        """
        try:
            m = self.isotopes["Atomic Mass"].loc[zaid.element + zaid.isotope]
        except KeyError:  # It means that it is a natural zaid
            # For a natural zaid the natural abundance mass is used
            df = self.isotopes.reset_index()
            df["Partial mass"] = df["Atomic Mass"] * df["Mean value"]
            masked = df.set_index("Z").loc[int(zaid.element)]
            m = masked["Partial mass"].sum()

        return float(m)

    def get_reactions(self, lib: str, parent: str) -> list[tuple[str, str]]:
        """
        get the reactions available for a specific zaid and parent nuclide

        Parameters
        ----------
        lib : str
            library suffix as in sheet name of the activation file.
        parent : str
            zaid number of the parent (e.g. 1001).

        Returns
        -------
        reactions : list[tuple[str, str]]
            contains tuple of (MT, daughter).

        """
        reactions = []
        try:
            df = self.reactions[lib].set_index("Parent")
            _, formula = self.get_zaidname(parent)
            formulazaid = formula.replace("-", "")  # eliminate the '-'
            # collect and provide as tuples
            subset = df.loc[formulazaid]
            try:
                for _, row in subset.iterrows():
                    MT = str(int(row["MT"]))
                    daughter = row["Daughter"]
                    daughter = self.get_zaidnum(daughter)
                    reactions.append((MT, daughter))

            except AttributeError:
                # then is not a DF but a Series
                MT = str(int(subset["MT"]))
                daughter = subset["Daughter"]
                daughter = self.get_zaidnum(daughter)
                reactions.append((MT, daughter))

        except KeyError:
            # library is not available or parent is not available
            pass

        return reactions