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
import re
import logging
import pandas as pd
import os
from importlib.resources import files, as_file

from f4enix.input.xsdirpyne import Xsdir
import f4enix.resources as pkg_res

# colors
CRED = '\033[91m'
CEND = '\033[0m'


MSG_DEFLIB = ' The Default library {} was used for zaid {}'


class LibManager:

    def __init__(self, xsdir_file: os.PathLike = None,
                 defaultlib: str = '81c', activationfile: os.PathLike = None,
                 isotopes_file: os.PathLike = None) -> None:
        """
        Object dealing with all complex operations that involves nuclear data

        Parameters
        ----------
        xsdir_file : str or path
            path to the MCNP xsdir reference file. If None (default) the
            default file shipped with the package is used.
        defaultlib : str, optional
            lib suffix to be used as default in translation operations.
            The default is '81c'.
        activationfile : str or path, optional
            path to the config file containing the reactions data for
            activation libraries. The default is None.
        isotopes_file : str or path, optional
            path to the isotopes files. If None (default) the default file
            shipped with the package is used.

        Returns
        -------
        None.

        Attributes
        ----------
        XS: Xsdir
            xsdir file representation
        isotopes: pd.DataFrame
            table containing different data for each isotopes such as formula,
            atomic number, natural abundance, atom mass, etc.
        defaultlib: str
            lib suffix to be used as default in translation operations.
        libraries: list[str]
            list of all the libraries present in the Xsdir file
        reactions: dict[str, pd.DataFrame]
            if an activation file was provided, this stores all the available
            decay reactions in d1s libraries.

        """
        # use both default files
        if xsdir_file is None and isotopes_file is None:
            # get the handlers for the default data files
            resources = files(pkg_res)
            XSDIR_FILE = as_file(resources.joinpath('xsdir.txt'))
            ISOTOPES_FILE = as_file(resources.joinpath('Isotopes.txt'))
            with XSDIR_FILE as xsdir, ISOTOPES_FILE as isotopes:
                abundances = pd.read_csv(isotopes, skiprows=2)
                self.XS = Xsdir(xsdir)

        # use only the default xsdir file
        elif xsdir_file is None:
            # get the handlers for the default data files
            resources = files(pkg_res)
            XSDIR_FILE = as_file(resources.joinpath('xsdir.txt'))
            with XSDIR_FILE as xsdir:
                self.XS = Xsdir(xsdir)
            abundances = pd.read_csv(isotopes_file, skiprows=2)

        # use only the default isotopes file
        elif isotopes_file is None:
            # get the handlers for the default data files
            resources = files(pkg_res)
            ISOTOPES_FILE = as_file(resources.joinpath('Isotopes.txt'))
            with ISOTOPES_FILE as isotopes:
                abundances = pd.read_csv(isotopes, skiprows=2)
            self.XS = Xsdir(xsdir_file)

        # no defaults
        else:
            self.XS = Xsdir(xsdir_file)
            abundances = pd.read_csv(isotopes_file, skiprows=2)

        abundances['idx'] = abundances['idx'].astype(str)
        abundances.set_index('idx', inplace=True)
        self.isotopes = abundances

        self.defaultlib = defaultlib

        # Identify different libraries installed. This is done checking H
        # libraries = self.check4zaid('1001')
        # libraries.extend(self.check4zaid('1000'))  # photons
        libraries = []
        for table in self.XS:
            lib = table.name.split('.')[1]
            if lib not in libraries:
                libraries.append(lib)

        self.libraries = libraries

        # Load the activation reaction data if available
        reactions = {}
        if activationfile is not None:
            file = pd.ExcelFile(activationfile)
        else:
            resources = files(pkg_res)
            with as_file(resources.joinpath('activation_libs.xlsx')) as infile:
                file = pd.ExcelFile(infile)
        for sheet in file.sheet_names:
            # Load the df that also needs to be filled
            reactions[sheet] = file.parse(sheet).ffill()

        # These are needed for faster operations
        newiso = self.isotopes.set_index(['E'])
        newiso = newiso.loc[~newiso.index.duplicated(keep='first')]
        self._newiso_byE = newiso.sort_index()

        newiso = self.isotopes.set_index(['Z'])
        newiso = newiso.loc[~newiso.index.duplicated(keep='first')]
        self._newiso_byZ = newiso.sort_index()

        self.reactions = reactions

    def check4zaid(self, zaid: str) -> list[str]:
        """
        Check which libraries are available for the selected zaid and return it

        Parameters
        ----------
        zaid : str
            zaid string (e.g. 1001).

        Returns
        -------
        libraries : list[str]
            list of libraries (tags) available for the zaid.

        """
        libraries = []
        for libname in self.XS.find_table(zaid, mode='default-fast'):
            libraries.append(libname)

        return libraries

    def convertZaid(self, zaid: str, lib: str
                    ) -> dict[str, tuple[str, float, float]]:
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
        if lib not in self.libraries:
            raise ValueError('Library '+lib+' is not available in xsdir file')

        zaidlibs = self.check4zaid(zaid)
        # Natural zaid
        if zaid[-3:] == '000':
            # Check if zaid has natural info
            if self.XS.find_table(zaid+'.'+lib, mode='exact'):
                translation = {zaid: (lib, 1, 1)}  # mass not important

            else:  # Has to be expanded
                translation = {}
                reduced = self.isotopes[self.isotopes['Z'] == int(zaid[:-3])]
                for idx, row in reduced.iterrows():
                    # zaid availability must be checked
                    if self.XS.find_table(idx+'.'+lib, mode='exact'):
                        newlib = lib
                    elif self.XS.find_table(idx+'.'+self.defaultlib,
                                            mode='exact'):
                        logging.warning(
                            MSG_DEFLIB.format(self.defaultlib, zaid))
                        newlib = self.defaultlib
                    else:
                        raise ValueError('No available translation for zaid :' +
                                         zaid+'It is needed for natural zaid expansion.')

                    translation[idx] = (newlib, row['Mean value'],
                                        row['Atomic Mass'])
        # 1to1
        elif lib in zaidlibs:
            translation = {zaid: (lib, 1, 1)}  # mass not important

        # No possible correspondence, natural or default lib has to be used
        else:
            # Check if the natural zaid is available
            natzaid = zaid[:-3]+'000'
            if self.XS.find_table(natzaid+'.'+lib, mode='exact'):
                translation = {natzaid: (lib, 1, 1)}  # mass not important
            # Check if default lib is available
            elif self.XS.find_table(zaid+'.'+self.defaultlib, mode='exact'):
                logging.warning(
                    MSG_DEFLIB.format(self.defaultlib, zaid))
                translation = {zaid: (self.defaultlib, 1, 1)}  # mass not imp
            else:
                # Check if any zaid cross section is available
                libraries = self.check4zaid(zaid)
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
                    raise ValueError('No available translation for zaid :' +
                                     zaid)

        return translation

    def get_libzaids(self, lib: str) -> list[str]:
        """
        Given a library, returns all zaids available

        Parameters
        ----------
        lib : str
            suffix of the library.

        Returns
        -------
        zaids : list[str]
            list of zaid names available in the library.

        """
        zaids = []

        for table in self.XS.find_zaids(lib):
            zaid = table.name.split('.')[0]
            if zaid not in zaids:
                zaids.append(zaid)

        return zaids

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
            splitted = zaid.split('.')
            elem = splitted[0][:-3]
            i = int(elem)
            isotope = splitted[0][-3:]

        else:
            i = int(zaid.element)
            isotope = zaid.isotope

        # newiso = self.isotopes.set_index('Z')
        # newiso = newiso.loc[~newiso.index.duplicated(keep='first')]

        name = self._newiso_byZ['Element'].loc[i]
        formula = self._newiso_byZ['E'].loc[i]+'-'+str(int(isotope))

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

        # newiso = self.isotopes.set_index(['E'])
        # newiso = newiso.loc[~newiso.index.duplicated(keep='first')]

        # split the name
        patnum = re.compile(r'\d+')
        patname = re.compile(r'[a-zA-Z]+')
        try:
            num = patnum.search(zaidformula).group()
            name = patname.search(zaidformula).group()
        except AttributeError:
            raise ValueError('No correspondent zaid found for '+zaidformula)

        atomnumber = self._newiso_byE.loc[name, 'Z']

        zaidnum = "{}{:03d}".format(atomnumber, int(num))

        return zaidnum

    def select_lib(self) -> str:
        """
        Prompt a library input selection with Xsdir availabilty check

        Returns
        -------
        lib : str
            Library to assess.

        """
        error = CRED+'''
 Error: {}
 The selected library is not available.
 '''+CEND
        # Add a counter to avoid falling in an endless loop
        i = 0
        while True:
            i += 1
            lib = input(' Select library (e.g. 31c or 99c-31c): ')
            if lib in self.libraries:
                break

            elif lib[0] == '{':
                libs = json.loads(lib)
                # all libraries should be available
                tocheck = list(libs.values())
                tocheck.extend(list(libs.keys()))
                flag = True
                for val in tocheck:
                    if val not in self.libraries:
                        print(error.format(val))
                        flag = False
                if flag:
                    break

            elif '-' in lib:
                libs = lib.split('-')
                flag = True
                for val in libs:
                    if val not in self.libraries:
                        print(error.format(val))
                        flag = False
                if flag:
                    break

            else:
                print(error.format(lib))

            if i > 20:
                raise ValueError('Too many wrong inputs')
        return lib

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
            m = self.isotopes['Atomic Mass'].loc[zaid.element+zaid.isotope]
        except KeyError:  # It means that it is a natural zaid
            # For a natural zaid the natural abundance mass is used
            df = self.isotopes.reset_index()
            df['Partial mass'] = df['Atomic Mass']*df['Mean value']
            masked = df.set_index('Z').loc[int(zaid.element)]
            m = masked['Partial mass'].sum()

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
            df = self.reactions[lib].set_index('Parent')
            _, formula = self.get_zaidname(parent)
            formulazaid = formula.replace('-', '')  # eliminate the '-'
            # collect and provide as tuples
            subset = df.loc[formulazaid]
            try:
                for _, row in subset.iterrows():
                    MT = str(int(row['MT']))
                    daughter = row['Daughter']
                    daughter = self.get_zaidnum(daughter)
                    reactions.append((MT, daughter))

            except AttributeError:
                # then is not a DF but a Series
                MT = str(int(subset['MT']))
                daughter = subset['Daughter']
                daughter = self.get_zaidnum(daughter)
                reactions.append((MT, daughter))

        except KeyError:
            # library is not available or parent is not available
            pass

        return reactions
