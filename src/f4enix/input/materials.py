"""
This module is related to the parsing of a material section of a MCNP input

The information is organized as follows:

MatCardsList
     Material
         Submaterial
               Element
                   Zaid

"""

from __future__ import annotations

import migjorn

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

import copy
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from contextlib import contextmanager
from decimal import Decimal
import pandas as pd

from f4enix.core.constants import (
    AVOGADRO_NUMBER,
    PAT_COMMENT,
    PAT_MAT,
    PAT_MX,
    PAT_SPACE,
)
from f4enix.input.libmanager import LibManager
from f4enix.core.irradiation import Nuclide

LM = LibManager()


# -------------------------------------
# == CLASSES FOR MATERIAL READING ==
# -------------------------------------
class Zaid:
    def __init__(
        self,
        fraction: str | float,
        nuclide: Nuclide,
    ) -> None:
        """
        Object representing a Zaid

        Parameters
        ----------
        fraction : str/float
            fraction of the zaid.
        nuclide : Nuclide
            nuclide object containing the element, isotope and library.

        Returns
        -------
        None.

        All the init parameters are saved as attributes. In addition:

        Attributes
        ----------
        name: str
            AAZZZ[lib]

        """
        self.fraction = float(fraction)
        self.nuclide = nuclide
        self.fullname = nuclide.write_to_formula().split(".")[0]
        self.name = nuclide.write_to_int_string()

    @property
    def element(self) -> int:
        return int(str(self.nuclide.zaid)[:-3])

    @property
    def isotope(self) -> int:
        return self.nuclide.isotope

    @property
    def library(self) -> str | None:
        return self.nuclide.lib

    @classmethod
    def from_string(cls, string: str) -> Zaid:
        """
        Generate a zaid object from an MCNP string

        Parameters
        ----------
        string : str
            original MCNP string.

        Returns
        -------
        Zaid
            created zaid.

        """
        # Divide fraction from zaid
        patSpacing = re.compile(r"[\s\t]+")
        items = patSpacing.split(string)

        # ZAID
        nuclide = Nuclide.from_int_string(items[0])

        # identify fraction
        fraction = items[1]

        return cls(fraction, nuclide)

    def to_text(
        self, abundance: float | None = None, elem_mass_fraction: float | None = None
    ) -> str:
        """
         Get the zaid string ready for MCNP material card

        Returns
        -------
        str
            zaid string.

        """
        fraction = "{:.6E}".format(Decimal(self.fraction))
        # Add INFO
        if abundance is None:
            abundance = ""
        else:
            abundance = "%s" % float("%.5g" % float(abundance))

        if elem_mass_fraction is None:
            mass_fraction = ""
        else:
            mass_fraction = "%s" % float("%.5g" % float(elem_mass_fraction * 100))

        abundance = "AB(%) " + abundance
        weight = "WEIGHT(%) " + mass_fraction
        inline_comm = "    $ " + self.fullname
        args = (self.name, fraction, inline_comm, weight, abundance)

        return "{0:>15} {1:>18} {2:<12} {3:<10} {4:<10}".format(*args)


class Element:
    def __init__(self, zaidList: list[Zaid]) -> None:
        """
        Generate an Element object starting from a list of zaids.
        It will collapse multiple instance of a zaid into a single one

        Parameters
        ----------
        zaidList : list
            list of zaids constituting the element.

        Returns
        -------
        None.

        Attributes
        ----------
        Z : str
            element str notation (AA) read from the zaids.
        zaids: list[Zaid]

        """
        zaids = {}
        for zaid in zaidList:
            # If already in dic sum the fractions
            if zaid.name in zaids.keys():
                zaids[zaid.name] = zaids[zaid.name] + zaid.fraction
            else:
                zaids[zaid.name] = zaid.fraction

        zaidList = []
        for name, fraction in zaids.items():
            zaidList.append(Zaid.from_string(name + " " + str(fraction)))

        self.Z = zaid.element
        self.name = zaid.nuclide.element
        self.zaids = zaidList

    def _get_abundances(self) -> dict[str, float]:
        """
        Update zaids abundance and mass fraction in the submaterial.

        Parameters
        ----------
        mass_fraction : float
            mass fraction of the element in the submaterial.

        Returns
        -------
        dict[str, float]
            dictionary of zaid abundances.

        """
        abundances = {}
        tot_fraction = 0
        for zaid in self.zaids:
            tot_fraction = tot_fraction + zaid.fraction

        for zaid in self.zaids:
            ab = zaid.fraction / tot_fraction * 100
            abundances[zaid.name] = ab
        return abundances

    def get_fraction(self) -> float:
        """
        Get the sum of the fraction of the zaids composing the element

        Returns
        -------
        fraction : float
            element fraction.

        """
        fraction = 0
        for zaid in self.zaids:
            fraction = fraction + zaid.fraction

        return fraction


class Material:
    # init method for zaid
    def __init__(
        self,
        name: str,
        zaids: list[Zaid],
        header: str = None,
        additional_keys: list[str] = None,
        mx_cards: list = None,
    ) -> None:
        """
        Generate a Material Object starting from a list of Zaid. Usually this kind of objects are generated
        directly reading a full material card, and rarely instanciated directly
        with the __init__ method.

        Parameters
        ----------
        name : str
            if the first submaterial, the name is the name of the material
            (e.g. m1).
        zaids : list[Zaid]
            list of zaids composing the material.
        header : str, optional
            Header of the submaterial. The default is None.
        additional_keys : list[str], optional
            list of additional keywords in the submaterial. The default is
            None.
        mx_cards : list, optional
            list of mx_cards in the material if present. The default is None.

        Returns
        -------
        None.

        Attributes
        ----------
        zaids: list[Zaid]
            list of zaids in the material
        elements: list[Element]
            list of elements in the material
        header: str
            comment in the MCNP input file that is the header of the material
        additional_keys: list[str]
            list of additional keys that may be part of the material
        """

        # List of zaids object of the submaterial
        self._zaids = zaids
        self.name = name.strip()
        self.elements, self._zaids = self._collapse_zaids()
        self.header = header

        # Additional keys as plib,hlib etc.
        if additional_keys is None:
            additional_keys = []
        self.additional_keys = additional_keys

        if mx_cards is None:
            self.mx_cards = []
        else:
            self.mx_cards = mx_cards

    @property
    def zaids(self) -> list[Zaid]:
        return self._zaids

    @zaids.setter
    def zaids(self, value: list[Zaid]) -> None:
        self._zaids = value
        self.elements, self._zaids = self._collapse_zaids()

    @classmethod
    def from_migjorn(cls, material: migjorn.Material) -> Material:
        """
        Generate a material from a migjorn material

        Parameters
        ----------
        material : migjorn.Material
            migjorn material object.

        Returns
        -------
        Material
            generated material.

        """
        zaids = []
        for zaid, fraction in material.entries:
            zaid = Zaid(fraction, Nuclide.from_int_string(zaid))
            zaids.append(zaid)

        # TODO: missing header property in mig material
        header = f"C temporary header: {material.id}"

        return cls(f"M{material.id}", zaids=zaids, header=header)

    @classmethod
    def from_zaids(
        cls,
        zaids: Sequence[tuple[str | int, float]],
        libman: LibManager,
        lib: str,
        name: str = "",
        mat_id: int = 1,
    ) -> Material:
        """Generate a material giving a list of zaids or elements.

        Parameters
        ----------
        zaids : list[tuple[str | int, float]]
            list of (zaid, fraction) couples, e.g., (1001, -0.1). If negative,
            fractions are intended as
            mass fraction, if positive as atom fractions. To define a material
            using elements, one can use natural zaids notation, e.g.,
            (1000, -1). The translation operation to the requested library
            will take care of expanding such zaids using natural abundances.
        libman : LibManager
            The usual library manager needed for XS operations.
        lib : str
            suffix of the library to apply to the material, e.g., 31c.
        name : str, optional
            this will be put in the header of the material card as a comment,
            by default ''.
        mat_id : int, optional
            material id, by default 1.

        Returns
        -------
        Material
            F4Enix material object
        """
        zaid_list = []
        for zaid, fraction in zaids:
            # first assume it was given as str
            try:
                zaid = libman.get_zaidnum(str(zaid))
            # else assume it was given as normal zaid
            except ValueError:
                zaid = str(zaid)
            zaid_list.append(Zaid(fraction, Nuclide.from_int_string(zaid)))

        material = cls(f"M{mat_id}", zaids=zaid_list, header=f"C {name}")
        material.translate(lib, libman)
        return material

    @classmethod
    def from_text(cls, text: list[str] | str) -> Material:
        """
        Generate a submaterial from MCNP input text

        Parameters
        ----------
        text : list[str] | str
            Original text of the MCNP input.

        Returns
        -------
        Material
            generated submaterial.

        """
        # Get a list of string splitting on newlines if a simple
        # string is provided
        if type(text) is str:
            text = text.splitlines()
        # As a first thing, let's be sure that no nasty "\r" special characters are
        # present in the text
        text = [line.replace("\r", "") for line in text]
        text = [line.replace("\n", "") for line in text]

        in_header = True
        header = ""
        zaidList = []
        additional_keys_list = []
        for line in text:
            material_starts = PAT_MAT.match(line)
            zaids = None
            additional_keys = None
            start = 0

            # Determine if still in header + parse of first line
            if material_starts:
                in_header = False
                name = material_starts.group()
                start = material_starts.end()

            if in_header:
                header = header + line
                continue

            # parse the material body
            zaids, additional_keys = _readLine(line[start:])
            if zaids is not None:
                zaidList.extend(zaids)
            if additional_keys is not None:
                additional_keys_list.extend(additional_keys)

        return cls(
            name,
            zaidList,
            header=header,
            additional_keys=additional_keys_list,
        )

    def _collapse_zaids(self) -> tuple[list[Element], list[Zaid]]:
        """
        Organize zaids into their elements and collapse mutiple istances

        Returns
        -------
        tuple[list[Element], list[Zaid]]

        """
        elements = {}
        for zaid in self.zaids:
            if zaid.element not in elements.keys():
                elements[zaid.element] = [zaid]
            else:
                elements[zaid.element].append(zaid)

        elemList = []
        for _, zaids in elements.items():
            elemList.append(Element(zaids))

        collapsed_zaids = []
        for elem in elemList:
            for zaid in elem.zaids:
                collapsed_zaids.append(zaid)

        return elemList, collapsed_zaids

    def to_text(self) -> str:
        """
        Write to text in MNCP format the submaterial

        Returns
        -------
        str
            formatted submaterial text.

        """
        # get additional data for the comments
        fake_mat = copy.deepcopy(self)
        fake_mat.switch_fraction("atom", LM)
        fake_mat.switch_fraction("mass", LM)
        df = fake_mat.get_info()

        if self.header is not None:
            text = self.header + "\n"
        else:
            text = ""

        text = text + self.name.upper() + "\n"

        for elem in self.elements:
            abundances = elem._get_abundances()
            elem_mass_fraction = (
                df.groupby("Element")["Element Mass Fraction"].mean().loc[elem.name]
            )
            for zaid in elem.zaids:
                text = (
                    text
                    + zaid.to_text(
                        abundance=abundances[zaid.name],
                        elem_mass_fraction=elem_mass_fraction,
                    )
                    + "\n"
                )

        # Add additional keys
        if len(self.additional_keys) > 0:
            text = text + "\t"
            for key in self.additional_keys:
                text = text + " " + key
        # Add mx cards
        for mx in self.mx_cards:
            for line in mx.lines:
                line = line.strip("\n")
                text = text + "\n" + line.upper()

        return text.strip("\n")

    def translate(self, newlib: dict | str, lib_manager: LibManager) -> None:
        """
        This method implements the translation logic of JADE. All zaids are
        translated accordingly to the newlib specified.

        Parameters
        ----------
        newlib : dict | str
            There are a few ways that newlib can be provided:

            1) str (e.g. 31c), the new library to translate to will be the
            one indicated;

            2) dic (e.g. {'98c' : '99c', '31c: 32c'}), the new library is
            determined based on the old library of the zaid

            3) dic (e.g. {'98c': [list of zaids], '31c': [list of zaids]}),
            the new library to be used is explicitly stated depending
            on the zaidnum.
        lib_manager : LibManager
            Object handling libraries operation.

        Returns
        -------
        None.

        """
        newzaids = []
        for zaid in self.zaids:
            # Implement the capability to translate to different libraries
            # depending on the starting one
            if type(newlib) == dict:
                # Check for which kind of dic it is
                if type(list(newlib.values())[0]) == str:
                    # The assignment is based on old lib
                    try:
                        newtag = newlib[zaid.library]
                    except KeyError:
                        # the zaid should have been assigned to a library
                        raise ValueError(
                            """
 Zaid {} was not assigned to any library""".format(zaid.name)
                        )

                else:
                    # The assignment is explicit, all libs need to be searched
                    newtag = None
                    for lib, zaids in newlib.items():
                        if zaid.nuclide.zaid in zaids:
                            newtag = lib
                            break
                    # Check that a library has been actually found
                    if newtag is None:
                        # the zaid should have been assigned to a library
                        raise ValueError(
                            """
 Zaid {} was not assigned to any library""".format(zaid.name)
                        )
            else:
                newtag = newlib

            # if it is a dosimetry library, the translation needs to be ignored
            if zaid.library in lib_manager.dosimetry_lib:
                # fake a 1to1 translation where the original suffix is retained
                translation = {zaid.nuclide.zaid: (zaid.library, 1, 1)}
            else:
                try:
                    translation = lib_manager.convertZaid(
                        str(zaid.nuclide.zaid), newtag
                    )
                except ValueError:
                    # No Available translation was found, ignore zaid
                    # Only video warning, to propagate to the log would be too much
                    print(
                        "  WARNING: no available translation was found for "
                        + zaid.name
                        + ".\n  The zaid has been ignored. "
                    )
                    continue

            # Check if it is  atomic or mass fraction
            if float(zaid.fraction) < 0:
                ref_mass = 0
                for key, item in translation.items():
                    ref_mass = ref_mass + item[1] * item[2]

                for key, item in translation.items():
                    fraction = str(item[1] * item[2] / ref_mass * zaid.fraction)
                    library = item[0]

                    newzaids.append(
                        Zaid(
                            fraction,
                            Nuclide.from_int_string(f"{key}.{library}"),
                        )
                    )

            else:
                for key, item in translation.items():
                    fraction = str(item[1] * zaid.fraction)
                    library = item[0]

                    newzaids.append(
                        Zaid(
                            fraction,
                            Nuclide.from_int_string(f"{key}.{library}"),
                        )
                    )

        self.zaids = newzaids

    def get_tot_fraction(self) -> float:
        """
        Returns the total material fraction
        """
        fraction = 0
        for zaid in self.zaids:
            fraction = fraction + zaid.fraction

        return fraction

    def switch_fraction(
        self, ftype: str, lib_manager: LibManager, inplace: bool = True
    ) -> Material | None:
        """
        Switch between atom or mass fraction for the material card.
        If the material is already switched the command is ignored.

        Parameters
        ----------
        ftype : str
            Either 'mass' or 'atom' to chose the type of switch.
        lib_manager : libmanager.LibManager
            Handles zaid data.
        inplace : bool
            if True the densities of the isotopes are changed inplace,
            otherwise a copy of the material is provided. DEFAULT is True

        Raises
        ------
        KeyError
            if ftype is not either 'atom' or 'mass'.

        Returns
        -------
        material : Material | None
            The material with switched fractions if inplace is False, otherwise None.

        """
        # Get total fraction
        totf = self.get_tot_fraction()

        if totf < 0 and ftype == "mass" or totf > 0 and ftype == "atom":
            if not inplace:
                return copy.deepcopy(self)
            else:
                # The switch is not needed, return None
                return None

        norm = 0
        for zaid in self.zaids:
            atom_mass = lib_manager.get_zaid_mass(zaid)
            if ftype == "atom":
                norm = norm + (-1 * zaid.fraction / atom_mass)
            elif ftype == "mass":
                norm = norm + (zaid.fraction * atom_mass)
            else:
                raise KeyError(ftype + " is not a valid key error [atom, mass]")

        if inplace:
            for zaid in self.zaids:
                atom_mass = lib_manager.get_zaid_mass(zaid)
                if ftype == "atom":
                    zaid.fraction = (-1 * zaid.fraction / atom_mass) / norm
                else:
                    zaid.fraction = (-1 * zaid.fraction * atom_mass) / norm
            self.elements, self._zaids = self._collapse_zaids()
            return None
        else:
            new_zaids = []
            for zaid in self.zaids:
                atom_mass = lib_manager.get_zaid_mass(zaid)
                newz = copy.deepcopy(zaid)
                if ftype == "atom":
                    newz.fraction = (-1 * zaid.fraction / atom_mass) / norm
                else:
                    newz.fraction = (-1 * zaid.fraction * atom_mass) / norm
                new_zaids.append(newz)
            mat = copy.deepcopy(self)
            mat.zaids = new_zaids
            return mat

    def _get_info_df(self) -> pd.DataFrame:
        """
        Returns DataFrame containing the different fractions of the elements
        and zaids.

        Returns
        -------
        table of information of fractions of elements and zaids in the material

        """
        dic_zaids = {
            "Element": [],
            "Isotope": [],
            "Zaid Fraction": [],
            "Elem Fraction": [],
        }
        for elem in self.elements:
            fraction = elem.get_fraction()
            # dic_element['Element'].append(elem.Z)
            for zaid in elem.zaids:
                dic_zaids["Elem Fraction"].append(fraction)
                elementname = zaid.nuclide.element
                dic_zaids["Element"].append(elementname)
                dic_zaids["Isotope"].append(
                    zaid.fullname + " [" + str(zaid.nuclide.zaid) + "]"
                )
                dic_zaids["Zaid Fraction"].append(zaid.fraction)

        df_zaids = pd.DataFrame(dic_zaids)

        return df_zaids

    def get_info(self) -> pd.DataFrame:
        """Get information on the fraction of the different elements and zaids contained
        in the materials.

        Returns
        -------
        df_complete: pd.DataFrame
            detailed dataframe
        df_elem: pd.DataFrame
            dataframe grouped at element level
        """
        material_atom = self.switch_fraction("atom", LM, inplace=False)
        material_mass = self.switch_fraction("mass", LM, inplace=False)

        df = self._get_info_df()
        df_a = material_atom._get_info_df()
        df_m = material_mass._get_info_df()

        df["Material"] = self.name
        df["Atom Fraction"] = df_a["Zaid Fraction"] / df_a["Zaid Fraction"].sum()
        df["Mass Fraction"] = df_m["Zaid Fraction"] / df_m["Zaid Fraction"].sum()
        df["Element Mass Fraction"] = (
            df_m["Elem Fraction"]
            / df_m.groupby("Element")["Elem Fraction"].mean().sum()
        )  # normalize the elemental fractions
        df["Element Atom Fraction"] = (
            df_a["Zaid Fraction"]
            / df_a.groupby("Element")["Zaid Fraction"].mean().sum()
        )  # normalize the elemental fractions

        df.set_index(["Material", "Element", "Isotope"], inplace=True)

        return df

    def scale_fractions(self, norm_factor: float) -> None:
        """
        Scale the zaids fractions using a normalizing factor

        Parameters
        ----------
        norm_factor : float
            scaling factor.

        Returns
        -------
        None.

        """
        for zaid in self.zaids:
            zaid.fraction = zaid.fraction * norm_factor

        self.elements, self._zaids = self._collapse_zaids()

    def get_tad(self, density: int | float, lib_manager: LibManager) -> float:
        """Return the total atom density of the material given the mass density.

        Parameters
        ----------
        density : int | float
            mass density in g/cm^3.
        lib_manager : LibManager
            Library manager for the conversion.

        Returns
        -------
        float
            atom density of the material in barn^-1 cm^-1.
        """
        self.switch_fraction("atom", lib_manager, inplace=True)
        mass_number = 0
        tot_fraction = self.get_tot_fraction()

        for zaid in self.zaids:
            mass_number = (
                mass_number
                + zaid.fraction
                / tot_fraction
                * lib_manager.isotopes["Atomic Mass"].loc[zaid.name.split(".")[0]]
            )

        return density / mass_number * AVOGADRO_NUMBER * 1e-24

    def get_density(self, tad: int | float, lib_manager: LibManager) -> float:
        """Return the density of the material given the total atom density.

        Parameters
        ----------
        tad : int | float
            total atom density in barn^-1 cm^-1.
        lib_manager : LibManager
            Library manager for the conversion.

        Returns
        -------
        float
            mass density of the material in g/cm^3.
        """
        self.switch_fraction("mass", lib_manager, inplace=True)
        mass_number_rec = 0
        tot_fraction = self.get_tot_fraction()

        for zaid in self.zaids:
            mass_number_rec = (
                mass_number_rec
                + zaid.fraction
                / tot_fraction
                / lib_manager.isotopes["Atomic Mass"].loc[zaid.name.split(".")[0]]
            )
        mass_number = 1 / mass_number_rec

        return tad * 1e24 / AVOGADRO_NUMBER * mass_number

    def fractions_to_atom_densities(
        self, lib_manager: LibManager, density: float
    ) -> None:
        """Given a specific mass density of the material, replace the atom fractions
        in the material with un-normalized atom fractions which correspond to the
        zaid atom densitiy in the material.

        This is obtained multiplying the normalized atom fraction of each zaid by the
        total atomic density of the material.

        Parameters
        ----------
        lib_manager : LibManager
            Library manager for the conversion.
        density : float
            mass density of the material in g/cm^3.
        """
        # force a double switch towards atom fractions
        self.switch_fraction("mass", lib_manager, inplace=True)
        self.switch_fraction("atom", lib_manager, inplace=True)
        tad = self.get_tad(density, lib_manager)
        for zaid in self.zaids:
            zaid.fraction = zaid.fraction * tad
        self.elements, self._zaids = self._collapse_zaids()


# Support function for Submaterial
def _readLine(string: str) -> tuple[list[Zaid], list[str] | None]:
    patSpacing = re.compile(r"[\s\t]+")
    patComment = re.compile(r"\$")
    patnumber = re.compile(r"\d+")

    pieces = patSpacing.split(string.strip())

    # kill comment section
    i = 0
    for i, piece in enumerate(pieces):
        if patComment.match(pieces[i]) is not None:
            del pieces[i:]
            break

    i = 0
    zaids = []
    additional_keys = None
    while True:
        try:
            # Check if it is zaid or keyword
            if patnumber.match(pieces[i]) is None or pieces[i] == "":
                additional_keys = pieces[i:]
                break
            else:
                zaidstring = pieces[i] + " " + pieces[i + 1]
                zaid = Zaid.from_string(zaidstring)
                zaids.append(zaid)

            i = i + 2

        except IndexError:
            break

    return zaids, additional_keys


class MatCardsList(Sequence):
    def __init__(self, materials: list[Material]) -> None:
        """
        Object representing the list of materials included in an MCNP input.
        This class is a child of the Sequence base class.

        Parameters
        ----------
        materials : list[Material]
            list of materials.

        Returns
        -------
        None.

        Attributes
        ----------
        materials: list[Material]
            list of the material objects that compose the material section
        matdic: dict[str, Material]
            dictionary of the material objects that compose the material
            section. This attribute is used mostly as internal object, the
            recommended way to access materials is to use directly the
            "dictionary" capabilities of the MatCardsList object.

        Examples
        --------

        >>> from f4enix.input.inputAPI import MatCardsList
        ... # initialize from file
        ... materials = MatCardsList.from_input('inputfile.i')
        ... # get a specific material
        ... mat1 = materials['m1']

        """
        self.materials = materials
        # Build also the dictionary
        self.matdic = self._compute_dic()

    def __len__(self) -> int:
        return len(self.materials)

    def __repr__(self) -> str:
        return str(self.matdic)

    def __str__(self) -> str:
        return str(self.matdic)

    def __getitem__(self, key: str | int) -> Material:
        if type(key) is int:
            return self.materials[key]
        else:
            return self.matdic[key.upper()]

    def append(self, material: Material) -> None:
        self.materials.append(material)
        self.matdic = self._compute_dic()

    def extend(self, materials: list[Material]) -> None:
        if type(materials) is not list:
            raise TypeError("'materials' should be a list of materials")
        # manually extend to be sure there are no duplicates
        for material in materials:
            if material.name.upper() not in self.matdic.keys():
                self.materials.append(material)
        self.matdic = self._compute_dic()

    def remove(self, item: Material) -> None:
        self.materials.remove(item)  # TODO this should get the key instead
        self.matdic = self._compute_dic()

    def _compute_dic(self) -> dict[str, Material]:
        matdic = {}
        for material in self.materials:
            matdic[material.name.upper()] = material

        return matdic

    @classmethod
    def from_input(cls, inputfile: os.PathLike) -> "MatCardsList":
        """Parse material cards from an MCNP input file using migjorn.

        Parameters
        ----------
        inputfile : os.PathLike
            MCNP input file containing the material section.

        Returns
        -------
        MatCardsList
        """
        model = migjorn.Model.from_file(str(inputfile))
        return cls.from_migjorn(model)

    @classmethod
    def from_migjorn(cls, model: migjorn.Model) -> "MatCardsList":

        # build the materials
        materials = []
        for mat in model.materials:
            materials.append(Material.from_migjorn(mat))

        groups: dict[str, list[tuple[str, str]]] = {}
        for mat in model.materials:
            groups.setdefault(str(mat.id), []).append(("M", mat.text))

        mat_card_list = cls(materials)

        for card in model.data_cards:
            if PAT_MX.match(card.name):
                mat_id = card.name.upper().replace("MX", "M")
                mat_card_list[mat_id].mx_cards.append(card)

        return mat_card_list

    def to_text(self) -> str:
        """
        return text of the material cards in order

        Returns
        -------
        str
            material card list MCNP formatted text.

        """
        text = ""
        for material in self.materials:
            text = text + "\n" + material.to_text()

        return text.strip("\n")

    def translate(self, newlib: str | dict, lib_manager: LibManager) -> None:
        """
        This method allows to translate the material cards to another library.
        The zaid are collapsed again to get the new elements

        Parameters
        ----------
        newlib : dict or str
            There are a few ways that newlib can be provided:

            1) str (e.g. 31c), the new library to translate to will be the
            one indicated;

            2) dic (e.g. {'98c' : '99c', '31c: 32c'}), the new library is
            determined based on the old library of the zaid

            3) dic (e.g. {'98c': [list of zaids], '31c': [list of zaids]}),
            the new library to be used is explicitly stated depending
            on the zaidnum.

        lib_manager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        None.

        """
        for material in self.materials:
            material.translate(newlib, lib_manager)

    def get_info(
        self, lib_manager: LibManager, zaids: bool = False
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Get the material informations in terms of fraction and composition
        of the material card

        Parameters
        ----------
        lib_manager : libmanager.LibManager
            To handle element name recovering.
        zaids : bool, optional
            Consider or not the zaid level. The default is False.

        Returns
        -------
        df : pd.DataFrame
            Raw infos on the fractions.
        df_elem : pd.DataFrame
            processed info for the element: normalized fraction added both for
            material and submaterial.

        """
        df_list = []
        df_elem_list = []
        for mat in self.materials:
            df, df_elem = mat.get_info(lib_manager, zaids=zaids)
            df_list.append(df)
            df_elem_list.append(df_elem)

        df = pd.concat(df_list)
        df_elem = pd.concat(df_elem_list)

        return df, df_elem

    def generate_material(
        self,
        materials: list[str],
        percentages: list[float],
        newlib: str,
        libmanager: LibManager,
        fractiontype="atom",
        mat_name="M1",
    ) -> Material:
        """
        Starting from an MCNP input, materials contained in its material list
        can be used to generate a new material combining them.

        Parameters
        ----------
        materials : list[str]
            list of materials to mix (e.g. ['m1', 'M2']).
        percentages : list[float]
            percentages associated to the source materials in the new materials
            (e.g. [0.1, 0.9)]. Their are intended as atom or mass fraction
            depending on the fractiontype that is specified.
        newlib : str
            library for the new material.
        fractiontype : str, optional
            type of fraction to use in the new material (either 'atom' or
            'mass'. The default is 'atom'.
        mat_name : str, optional
            Material card name of the new generated material. the default is
            'M1'
        Returns
        -------
        Material
            Newly created material

        """
        mat_name = mat_name.upper()

        if re.match(r"^M\d{1,7}$", mat_name) is None:
            print("\nMaterial name not valid, set to M1\n")
            mat_name = "M1"

        # Translate to requested lib
        self.translate(newlib, libmanager)

        # Collect all submaterials
        zaids = []
        main_header = ""
        for materialname, percentage in zip(materials, percentages):
            materialname = materialname.upper()
            percentage_str = str(round(float(percentage) * 100, 2)) + "%"
            main_header = f"{main_header}C Material: {materialname} Percentage: {percentage_str} ({fractiontype})\n"
            material = copy.deepcopy(self[materialname])
            # Ensure materials have the requested fraction type
            material.switch_fraction(fractiontype, libmanager)

            # Scale fractions
            totfraction = material.get_tot_fraction()

            # normalized & scaled
            norm_factor = float(percentage) / totfraction
            if fractiontype == "mass":
                norm_factor = -norm_factor
            material.scale_fractions(norm_factor)

            zaids.extend(material.zaids)

        # Generate new material and matlist
        newmat = Material(mat_name, zaids=zaids, header=main_header)

        return newmat

    def fractions_to_atom_densities(
        self, lib_manager: LibManager, density: float
    ) -> None:
        """Given a specific mass density of the material, replace the atom fractions
        in the material with un-normalized atom fractions which correspond to the
        zaid atom densitiy in the material.

        This is obtained multiplying the normalized atom fraction of each zaid by the
        total atomic density of the material.

        Parameters
        ----------
        lib_manager : LibManager
            Library manager for the conversion.
        density : float
            mass density of the material in g/cm^3.
        """
        for material in self.materials:
            material.fractions_to_atom_densities(lib_manager, density)


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
