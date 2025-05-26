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
from numjuggler import parser as par

from f4enix.constants import AVOGADRO_NUMBER, PAT_COMMENT, PAT_MAT, PAT_MX
from f4enix.input.libmanager import LibManager


def indent(elem, level: int = 0) -> None:
    """Indent an XML element and its children to make the XML structure more human-readable.

    Parameters
    ----------
    elem : xml.etree.ElementTree.Element
        XML element to be indented
    level : int, optional
        Current level of indentation, by default 0
    """

    # Create the indentation string based on the specified level
    i = "\n" + level * "  "
    if len(elem):
        # If the element has child elements
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


# -------------------------------------
# == CLASSES FOR MATERIAL READING ==
# -------------------------------------
class Zaid:
    def __init__(
        self,
        fraction: str | float,
        element: str,
        isotope: str,
        library: str,
        ab: str = "",
        fullname: str = "",
        elem_mass_fraction: float | None = None,
    ) -> None:
        """
        Object representing a Zaid

        Parameters
        ----------
        fraction : str/float
            fraction of the zaid.
        element : str
            element part of the zaid (AA).
        isotope : str
            isotope part of the zaid (ZZZ).
        library : str
            library suffix (e.g. 99c).
        ab : str, optional
            abundance of the zaid in the material. The default is ''.
        fullname : str, optional
            formula name (e.g. H1). The default is ''.
        elem_mass_fraction : float, optional
            mass fraction of the element in the submaterial. The default is None

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
        self.element = element
        self.isotope = isotope
        self.library = library
        self.ab = ab
        self.fullname = fullname
        self.elem_mass_fraction = elem_mass_fraction

        if self.library is None:
            self.name = self.element + self.isotope
        else:
            self.name = self.element + self.isotope + "." + self.library

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
        pieces = items[0].split(".")
        # Try to identify library if present
        try:
            library = pieces[1]
        except IndexError:
            library = None

        # Identify element and isotope
        element = pieces[0][:-3]
        isotope = pieces[0][-3:]

        # identify fraction
        fraction = items[1]

        return cls(fraction, element, isotope, library)

    def to_text(self) -> str:
        """
         Get the zaid string ready for MCNP material card

        Returns
        -------
        str
            zaid string.

        """
        fraction = "{:.6E}".format(Decimal(self.fraction))
        if self.library is None:
            zaidname = self.element + self.isotope
        else:
            zaidname = self.element + self.isotope + "." + self.library
        # Add INFO
        try:
            abundance = "%s" % float("%.5g" % float(self.ab))
        except ValueError:
            abundance = ""
        try:
            mass_fraction = "%s" % float("%.5g" % float(self.elem_mass_fraction * 100))
        except TypeError:
            mass_fraction = ""

        abundance = "AB(%) " + abundance
        weight = "WEIGHT(%) " + mass_fraction
        inline_comm = "    $ " + self.fullname
        args = (zaidname, fraction, inline_comm, weight, abundance)

        return "{0:>15} {1:>18} {2:<12} {3:<10} {4:<10}".format(*args)

    def to_xml(self, libmanager: LibManager, submaterial: SubMaterial) -> None:
        """Generate XML content for a nuclide within a material.

        Parameters
        ----------
        libmanager :
            libmanager
        submaterial :
            The XML element for the material where the nuclide content will be added.
        """
        nuclide = self.get_fullname(libmanager).replace("-", "")
        if self.fraction < 0.0:
            ET.SubElement(
                submaterial, "nuclide", name=nuclide, wo=str(abs(self.fraction))
            )
        else:
            ET.SubElement(
                submaterial, "nuclide", name=nuclide, ao=str(abs(self.fraction))
            )

    def get_fullname(self, libmanager: LibManager) -> str:
        """
        Get the formula name of the zaid (e.g. H1)

        Parameters
        ----------
        libmanager : libmanager.LibManager
            libmanager handling the libraries operations.

        Returns
        -------
        formula : str
            zaid formula name.

        """
        name, formula = libmanager.get_zaidname(self)

        return formula


#    def update_info(self,ab,fullname):
#        """
#        Update zaid info
#        """
#        self.additional_info['ab'] = ab
#        self.additional_info['fullname'] = fullname


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
        self.zaids = zaidList

    def update_zaidinfo(self, libmanager: LibManager, mass_fraction: float) -> None:
        """
        Update zaids infos through a libmanager. Info are the formula name and
        the abundance in the material.

        Parameters
        ----------
        libmanager : libmanager.LibManager
            libmanager handling the libraries operations.
        mass_fraction : float
            mass fraction of the element in the submaterial.

        Returns
        -------
        None.

        """
        tot_fraction = 0
        for zaid in self.zaids:
            tot_fraction = tot_fraction + zaid.fraction

        for zaid in self.zaids:
            fullname = zaid.get_fullname(libmanager)
            ab = zaid.fraction / tot_fraction * 100
            #            zaid.update_info(ab,fullname)
            zaid.ab = ab
            zaid.fullname = fullname
            zaid.elem_mass_fraction = mass_fraction

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


class SubMaterial:
    # init method for zaid
    def __init__(
        self,
        name: str,
        zaidList: list[Zaid],
        elemList: list[Element] = None,
        header: str = None,
        additional_keys: list[str] = None,
    ) -> None:
        """
        Generate a SubMaterial Object starting from a list of Zaid and
        eventually Elements list. Usually this kind of objects are generated
        directly reading a full material card, and rarely instanciated directly
        with the __init__ method.

        Parameters
        ----------
        name : str
            if the first submaterial, the name is the name of the material
            (e.g. m1).
        zaidList : list[Zaid]
            list of zaids composing the submaterial.
        elemList : list[Element], optional
            list of elements composing the submaterial. The default is None.
        header : str, optional
            Header of the submaterial. The default is None.
        additional_keys : list[str], optional
            list of additional keywords in the submaterial. The default is
            None.

        Returns
        -------
        None.

        Attributes
        ----------
        zaidList: list[Zaid]
            list of zaids in the sub-material
        elements: list[Element]
            list of elements in the sub-material
        header: str
            comment in the MCNP input file that is the header of the submat
        additional_keys: list[str]
            list of additional keys that may be part of the material

        """

        # List of zaids object of the submaterial
        self.zaidList = zaidList

        # Name of the material
        if name is not None:
            self.name = name.strip()  # Be sure to strip spaces
        else:
            self.name = None

        # List of elements in material
        if elemList is None:
            self._collapse_zaids()
        else:
            self.elements = elemList

        # Header of the submaterial
        self.header = header

        # Additional keys as plib,hlib etc.
        if additional_keys is None:
            additional_keys = []
        self.additional_keys = additional_keys

    @classmethod
    def from_text(cls, text: str) -> SubMaterial:
        """
        Generate a submaterial from MCNP input text

        Parameters
        ----------
        text : list[str]
            Original text of the MCNP input.

        Returns
        -------
        SubMaterial
            generated submaterial.

        """
        # Useful patterns
        patSpacing = re.compile(r"[\s\t]+")
        patComment = PAT_COMMENT
        patName = PAT_MAT
        searchHeader = True
        header = ""
        zaidList = []
        additional_keys_list = []
        for line in text:
            zaids = None
            additional_keys = None
            # Header MUST be at the top of the text block
            if searchHeader:
                # Get header
                if patComment.match(line) is None:
                    searchHeader = False
                    # Special treatment for first line
                    try:
                        name = patName.match(line).group()
                    except AttributeError:
                        # There is no material name
                        name = None

                    pieces = patSpacing.split(line)
                    if len(pieces) > 1:
                        # CASE1: only material name+additional spacing
                        if pieces[1] == "":
                            pass  # no more actions for this line
                        # CASE2: material name + zaids or only zaids
                        else:
                            if name is None:
                                start = 0
                            else:
                                start = patName.match(line).end()
                            zaids, additional_keys = _readLine(line[start:])
                    # CASE3: only material name and no spacing
                    else:
                        pass  # no more actions for this line
                else:
                    header = header + line
                    continue
            else:
                zaids, additional_keys = _readLine(line)

            if zaids is not None:
                zaidList.extend(zaids)

            if additional_keys is not None:
                additional_keys_list.extend(additional_keys)

        return cls(
            name,
            zaidList,
            elemList=None,
            header=header[:-1],
            additional_keys=additional_keys_list,
        )

    def _collapse_zaids(self) -> None:
        """
        Organize zaids into their elements and collapse mutiple istances

        Returns
        -------
        None.

        """
        elements = {}
        for zaid in self.zaidList:
            if zaid.element not in elements.keys():
                elements[zaid.element] = [zaid]
            else:
                elements[zaid.element].append(zaid)

        elemList = []
        for _, zaids in elements.items():
            elemList.append(Element(zaids))

        self.elements = elemList

    def to_text(self) -> str:
        """
        Write to text in MNCP format the submaterial

        Returns
        -------
        str
            formatted submaterial text.

        """
        if self.header is not None:
            text = self.header + "\n"
        else:
            text = ""
        # if self.name is not None:
        #     text = text+'\n'+self.name
        if self.elements is not None:
            for elem in self.elements:
                for zaid in elem.zaids:
                    text = text + zaid.to_text() + "\n"
        else:
            for zaid in self.zaidList:
                text = text + zaid.to_text() + "\n"

        # Add additional keys
        if len(self.additional_keys) > 0:
            text = text + "\t"
            for key in self.additional_keys:
                text = text + " " + key

        return text.strip("\n")

    def to_xml(self, libmanager: LibManager, material: Material) -> None:
        """Generate XML content for a material and add it to a material tree.

        Parameters
        ----------
        libmanager :
            libmanager handling the libraries operations.
        material :
            The XML tree where the material content will be added.
        """

        # matid = id
        # matname = str(self.name)
        # matdensity = str(abs(density))
        # if density < 0:
        #    density_units = "g/cc"
        # else:
        #    density_units = "atom/b-cm"
        # submaterial = ET.SubElement(material_tree, "material", id=matid, name=matname)
        # ET.SubElement(submaterial, "density", value=matdensity, units=density_units)
        if self.elements is not None:
            for elem in self.elements:
                for zaid in elem.zaids:
                    zaid.to_xml(libmanager, material)
        else:
            for zaid in self.zaidList:
                zaid.to_xml(libmanager, material)

    def translate(
        self, newlib: dict | str, lib_manager: LibManager, code: str = "mcnp"
    ) -> None:
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
        for zaid in self.zaidList:
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
                    zaidnum = zaid.element + zaid.isotope
                    for lib, zaids in newlib.items():
                        if zaidnum in zaids:
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
                translation = {zaid.element + zaid.isotope: (zaid.library, 1, 1)}
            else:
                try:
                    translation = lib_manager.convertZaid(
                        zaid.element + zaid.isotope, newtag, code
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
                    element = str(key)[:-3]
                    isotope = str(key)[-3:]
                    library = item[0]

                    newzaids.append(Zaid(fraction, element, isotope, library))

            else:
                for key, item in translation.items():
                    fraction = str(item[1] * zaid.fraction)
                    element = str(key)[:-3]
                    isotope = str(key)[-3:]
                    library = item[0]

                    newzaids.append(Zaid(fraction, element, isotope, library))

        self.zaidList = newzaids
        self._collapse_zaids()

    def _update_info(
        self, lib_manager: LibManager, element_mass_fractions: pd.Series
    ) -> None:
        """
        This methods allows to update the in-line comments for every zaids
        containing additional information

        Parameters
        ----------
        lib_manager : libmanager.LibManager
            Library manager for the conversion.
        element_mass_fractions : pd.Series
            Series containing the mass fractions of the elements in the submaterial.

        Returns
        -------
        None.

        """
        self._collapse_zaids()  # To be sure to have adjourned elements

        for elem in self.elements:
            fullname = Zaid.from_string(elem.Z + "000 -1").get_fullname(lib_manager)
            element_name = fullname.split("-")[0]
            mass_fraction = element_mass_fractions[element_name]
            elem.update_zaidinfo(lib_manager, mass_fraction)

        # TODO
        # Here the zaidlist of the submaterial should be adjourned or the next
        # collapse zaid will cancel the informations. If update info is used
        # as last operations there are no problems.

    def get_info(self, lib_manager: LibManager) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Returns DataFrame containing the different fractions of the elements
        and zaids

        Parameters
        ----------
        lib_manager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        df_el : pd.DataFrame
            table of information of the submaterial on an elemental level.
        df_zaids : pd.DataFrame
            table of information of the submaterial on a zaid level.

        """
        # dic_element = {'Element': [], 'Fraction': []}
        # dic_zaids = {'Element': [], 'Zaid': [], 'Fraction': []}
        dic_element = {"Element": [], "Fraction": []}
        dic_zaids = {"Element": [], "Isotope": [], "Fraction": []}
        for elem in self.elements:
            fraction = elem.get_fraction()
            # dic_element['Element'].append(elem.Z)
            dic_element["Fraction"].append(fraction)
            for zaid in elem.zaids:
                fullname = zaid.get_fullname(lib_manager)
                elementname = fullname.split("-")[0]
                # dic_zaids['Element'].append(elem.Z)
                # dic_zaids['Zaid'].append(zaid.isotope)
                dic_zaids["Element"].append(elementname)
                dic_zaids["Isotope"].append(
                    fullname + " [" + str(zaid.element) + str(zaid.isotope) + "]"
                )
                dic_zaids["Fraction"].append(zaid.fraction)

            dic_element["Element"].append(elementname)

        df_el = pd.DataFrame(dic_element)
        df_zaids = pd.DataFrame(dic_zaids)

        return df_el, df_zaids

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
        for zaid in self.zaidList:
            zaid.fraction = zaid.fraction * norm_factor

        self._collapse_zaids()


# Support function for Submaterial
def _readLine(string: str) -> tuple[list[Zaid], list[str] | None]:
    patSpacing = re.compile(r"[\s\t]+")
    patComment = re.compile(r"\$")
    patnumber = re.compile(r"\d+")

    pieces = patSpacing.split(string)
    # kill first piece if it is void
    if pieces[0] == "":
        del pieces[0]
    # kill last piece if it is void
    if pieces[-1] == "":
        del pieces[-1]

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


class Material:
    def __init__(
        self,
        zaids: list[Zaid],
        elem: list[Element],
        name: str,
        submaterials: list[SubMaterial] = None,
        mx_cards: list = None,
        header: str = None,
        density: float = None,
    ) -> None:
        """
        Object representing an MCNP material

        Parameters
        ----------
        zaids : list[zaids]
            zaids composing the material.
        elem : list[elem]
            elements composing the material.
        name : str
            name of the material (e.g. m1).
        submaterials : list[Submaterials], optional
            list of submaterials composing the material. The default is None.
        mx_cards : list, optional
            list of mx_cards in the material if present. The default is None.
        header : str, optional
            material header. The default is None.
        density : float, optional
            material density, used for OpenMC materials. Default is None.
        Returns
        -------
        None.

        Attributes
        ----------
        all __init__ parameters are stored as attributes.

        """

        self.zaids = zaids
        self.elem = elem
        self.submaterials = submaterials
        self.name = name.strip()
        if mx_cards is None:
            self.mx_cards = []
        else:
            self.mx_cards = mx_cards
        self.header = header
        self.density = density

        # Adjust the submaterial and headers reading
        try:
            # The first submaterial header is actually the material header.
            # If submat is void it has to be deleted (only header), otherwise
            # it means it has no header
            submat = submaterials[0]
            if len(submat.zaidList) == 0:  # Happens in reading from text
                # self.header = submat.header
                del self.submaterials[0]
            # else:
            #     self.submaterials[0].header = None
        except IndexError:
            self.header = None

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
            zaid_list.append(Zaid(fraction, zaid[:-3], zaid[-3:], None))

        submat = SubMaterial("", zaid_list)
        submat.translate(lib, libman)
        material = cls(
            None, None, f"M{mat_id}", submaterials=[submat], header=f"C {name}"
        )
        material._update_info(libman)
        return material

    @classmethod
    def from_text(cls, text: list[str]) -> Material:
        """
        Create a material from MCNP formatted text

        Parameters
        ----------
        text : list[str]
            Transport code formatted text representing the material.

        Returns
        -------
        matreader.Material
            material object created.

        """
        # split the different submaterials
        patC = PAT_COMMENT
        pat_matHeader = PAT_MAT
        inHeader = True
        subtext = []
        submaterials = []

        # As a first thing, let's be sure that no nasty "\r" special characters are
        # present in the text
        text = [line.replace("\r", "") for line in text]

        for line in text:
            checkComment = patC.match(line)
            checkHeaderMat = pat_matHeader.match(line)

            if checkHeaderMat is not None:
                header = "".join(subtext)
                subtext = []

            if inHeader:
                subtext.append(line)
                if checkComment is None:  # The end of the header is found
                    inHeader = False
            else:
                if checkComment is None:  # Still in the material
                    subtext.append(line)
                else:  # a new header starts
                    submaterials.append(SubMaterial.from_text(subtext))
                    inHeader = True
                    subtext = [line]

        submaterials.append(SubMaterial.from_text(subtext))

        return cls(
            None, None, submaterials[0].name, submaterials=submaterials, header=header
        )

    def to_text(self) -> str:
        """
        Write the material to MCNP formatted text

        Returns
        -------
        str
            MCNP formatte text representing the material.

        """
        if self.density is not None:
            if self.header is not None:
                text = (
                    self.header.strip("\n")
                    + "\n"
                    + self.name.lower().strip("\n")
                    + " "
                    + str(self.density)
                )
            else:
                text = self.name.lower() + " " + str(self.density)
        else:
            if self.header is not None:
                text = self.header.strip("\n") + "\n" + self.name.upper().strip("\n")
            else:
                text = self.name.upper()
        if self.submaterials is not None:
            for submaterial in self.submaterials:
                text = text + "\n" + submaterial.to_text()
            # Add mx cards
            for mx in self.mx_cards:
                for line in mx.lines:
                    line = line.strip("\n")
                    text = text + "\n" + line.upper()
        else:
            text = "  Not supported yet, generate submaterials first"
            pass  # TODO

        return text.strip("\n")

    def to_xml(self, libmanager: LibManager, material_tree: ET.Element) -> None:
        """Generate XML content for a material and its submaterials.

        Parameters
        ----------
        libmanager :
            libmanager
        material_tree :
            The XML element for the material where content will be added.
        """
        matid = re.sub("[^0-9]", "", str(self.name))
        matname = str(self.name)
        matdensity = str(abs(self.density))
        if self.density < 0:
            density_units = "g/cc"
        else:
            density_units = "atom/b-cm"
        material = ET.SubElement(material_tree, "material", id=matid, name=matname)
        ET.SubElement(material, "density", value=matdensity, units=density_units)
        if self.submaterials is not None:
            for submaterial in self.submaterials:
                submaterial.to_xml(libmanager, material)

    def translate(
        self,
        newlib: dict | str,
        lib_manager: LibManager,
        code: str = "mcnp",
        update: bool = None,
    ) -> None:
        """
        This method allows to translate all submaterials to another library

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
            object handling all libraries operations.
        code : str, optional
            Monte Carlo code for material format. The default is 'mcnp'.
        update : bool, optional
            if True, material infos are updated. The default is True.

        Returns
        -------
        None.
        """
        for submat in self.submaterials:
            submat.translate(newlib, lib_manager, code)

        self._update_info(lib_manager)

    def get_tot_fraction(self) -> float:
        """
        Returns the total material fraction
        """
        fraction = 0
        for submat in self.submaterials:
            for zaid in submat.zaidList:
                fraction = fraction + zaid.fraction

        return fraction

    def add_mx(self, mx_cards: list) -> None:
        """
        Add a list of mx_cards to the material
        """
        self.mx_cards.append(mx_cards)

    def _update_info(self, lib_manager: LibManager) -> None:
        """
        This methods allows to update the in-line comments for every zaids
        containing additional information

        lib_manager: (LibManager) Library manager for the conversion
        """
        fake_mat = copy.deepcopy(self)
        fake_mat.switch_fraction("atom", lib_manager)
        fake_mat.switch_fraction("mass", lib_manager)
        _, df_elem = fake_mat.get_info(lib_manager)

        for i, submaterial in enumerate(self.submaterials):
            element_fractions = df_elem.loc[self.name, i + 1]["Sub-Material Fraction"]
            submaterial._update_info(lib_manager, element_fractions)

    def switch_fraction(
        self, ftype: str, lib_manager: LibManager, inplace: bool = True
    ) -> list[SubMaterial]:
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
        submaterials : list[SubMaterial]
            list of the submaterials where fraction have been switched

        """
        # Get total fraction
        totf = self.get_tot_fraction()
        new_submats = []

        if ftype == "atom":  # mass2atom switch
            if totf < 0:  # Check if the switch must be effectuated
                # x_n = (x_m/m)/sum(x_m/m)
                # get sum(x_m/m)
                norm = 0
                for submat in self.submaterials:
                    for zaid in submat.zaidList:
                        atom_mass = lib_manager.get_zaid_mass(zaid)
                        norm = norm + (-1 * zaid.fraction / atom_mass)

                for submat in self.submaterials:
                    new_zaids = []
                    new_submat = copy.deepcopy(submat)
                    for zaid in submat.zaidList:
                        atom_mass = lib_manager.get_zaid_mass(zaid)
                        if inplace:
                            zaid.fraction = (-1 * zaid.fraction / atom_mass) / norm
                        else:
                            newz = copy.deepcopy(zaid)
                            newz.fraction = (-1 * zaid.fraction / atom_mass) / norm
                            new_zaids.append(newz)
                    new_submat.zaidList = new_zaids
                    # new_submat._update_info(lib_manager)
                    # adjourn the element list
                    new_submat._collapse_zaids()
                    submat._collapse_zaids()
                    new_submats.append(new_submat)
            else:
                new_submats = self.submaterials

        elif ftype == "mass":  # atom2mass switch
            if totf > 0:  # Check if the switch must be effectuated
                # x_n = (x_m*m)/sum(x_m*m)
                # get sum(x_m*m)
                norm = 0
                for submat in self.submaterials:
                    for zaid in submat.zaidList:
                        atom_mass = lib_manager.get_zaid_mass(zaid)
                        norm = norm + (zaid.fraction * atom_mass)

                for submat in self.submaterials:
                    new_zaids = []
                    new_submat = copy.deepcopy(submat)
                    for zaid in submat.zaidList:
                        atom_mass = lib_manager.get_zaid_mass(zaid)
                        if inplace:
                            zaid.fraction = (-1 * zaid.fraction * atom_mass) / norm
                        else:
                            newz = copy.deepcopy(zaid)
                            newz.fraction = (-1 * zaid.fraction * atom_mass) / norm
                            new_zaids.append(newz)
                    new_submat.zaidList = new_zaids
                    # adjourn the element list
                    new_submat._collapse_zaids()
                    submat._collapse_zaids()
                    # new_submat._update_info(lib_manager)
                    new_submats.append(new_submat)
            else:
                new_submats = self.submaterials

        else:
            raise KeyError(ftype + " is not a valid key error [atom, mass]")

        # self._update_info(lib_manager)

        return new_submats

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
        for submat in self.submaterials:
            for zaid in submat.zaidList:
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
        for submat in self.submaterials:
            for zaid in submat.zaidList:
                mass_number_rec = (
                    mass_number_rec
                    + zaid.fraction
                    / tot_fraction
                    / lib_manager.isotopes["Atomic Mass"].loc[zaid.name.split(".")[0]]
                )
        mass_number = 1 / mass_number_rec

        return tad * 1e24 / AVOGADRO_NUMBER * mass_number

    def get_info(self, lib_manager: LibManager, zaids: bool = False):
        """Get information on the fraction of the different elements and zaids contained
        in the materials.

        Parameters
        ----------
        lib_manager : LibManager
            library manager to handle lib operations
        zaids : bool, optional
            If true the info is output also at zaid level, by default False

        Returns
        -------
        df_complete: pd.DataFrame
            detailed dataframe
        df_elem: pd.DataFrame
            dataframe grouped at element level
        """
        infos = []
        complete_infos = []
        submats_atom = self.switch_fraction("atom", lib_manager, inplace=False)
        submats_mass = self.switch_fraction("mass", lib_manager, inplace=False)
        i = 0
        for submat, submat_a, submat_m in zip(
            self.submaterials, submats_atom, submats_mass
        ):
            dic_el, dic_zaids = submat.get_info(lib_manager)
            dic_el_a, dic_zaids_a = submat_a.get_info(lib_manager)
            dic_el_m, dic_zaids_m = submat_m.get_info(lib_manager)

            if zaids:
                dic = dic_zaids
                dic_a = dic_zaids_a
                dic_m = dic_zaids_m
            else:
                dic = dic_el
                dic_a = dic_el_a
                dic_m = dic_el_m

            dic["Material"] = self.name
            dic["Submaterial"] = i + 1
            infos.append(dic)

            c_dic = copy.deepcopy(dic)
            c_dic["Atom Fraction"] = dic_a["Fraction"]
            c_dic["Mass Fraction"] = dic_m["Fraction"]
            complete_infos.append(c_dic)

            i = i + 1

        df = pd.concat(infos)
        df_complete = pd.concat(complete_infos)
        del df_complete["Fraction"]

        if zaids:
            df.set_index(
                ["Material", "Submaterial", "Element", "Isotope"], inplace=True
            )
            df_complete.set_index(
                ["Material", "Submaterial", "Element", "Isotope"], inplace=True
            )

        else:
            df.set_index(["Material", "Submaterial", "Element"], inplace=True)
            df_complete.set_index(["Material", "Submaterial", "Element"], inplace=True)

        # Additional df containing normalized element fraction of submaterial
        # and material

        # Get total fractions
        df_elem = df.groupby(["Material", "Submaterial", "Element"]).sum()
        df_sub = df.groupby(["Material", "Submaterial"]).sum()
        df_mat = df.groupby(["Material"]).sum()

        # Compute percentages
        sub_percentage = []
        mat_percentage = []
        for idx, row in df_elem.iterrows():
            matID = idx[0]
            elemID = idx[1]
            sub_percentage.append(
                row["Fraction"] / df_sub["Fraction"].loc[(matID, elemID)]
            )
            mat_percentage.append(row["Fraction"] / df_mat["Fraction"].loc[matID])

        df_elem["Sub-Material Fraction"] = sub_percentage
        df_elem["Material Fraction"] = mat_percentage

        return df_complete, df_elem

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
        for submat in self.submaterials:
            for zaid in submat.zaidList:
                zaid.fraction = zaid.fraction * tad
            submat._collapse_zaids()
        self._update_info(lib_manager)


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
    def from_input(cls, inputfile: os.PathLike) -> MatCardsList:
        """
        This method use the numjuggler parser to help identify the mcards in
        the input. Then the mcards are parsed using the classes defined in this
        module

        Parameters
        ----------
        inputfile : os.PathLike
            MCNP input file containing the material section.

        Returns
        -------
        MatCardsList
            new material card list generated.

        """
        matPat = PAT_MAT
        mxPat = PAT_MX
        commentPat = PAT_COMMENT
        # Using parser the data cards are extracted from the input.
        # Comment section are interpreted as cards by the parser
        with suppress_stdout():
            # Suppress output from tab replacing
            cards = par.get_cards_from_input(inputfile)
            cardsDic = par.get_blocks(cards)
        datacards = cardsDic[5]

        materials = []
        previous_lines = [""]
        mx_cards = []
        mx_found = False

        for datacard in datacards:
            lines = datacard.lines

            # Check if it is a material card
            if matPat.match(lines[0]) is not None:
                # Check if previous card is the header
                if commentPat.match(previous_lines[0]):
                    previous_lines.extend(lines)
                    material = Material.from_text(previous_lines)
                else:
                    material = Material.from_text(lines)

                materials.append(material)

            # Check if the current is an mx cards
            if mxPat.match(lines[0]) is not None:
                mx_cards.append(lines)
                mx_found = True

            # If not Add mx cards if previous one was an mx
            elif mx_found:
                materials[-1].add_mx(mx_cards)
                mx_cards = []
                mx_found = False

            else:
                mx_found = False

            previous_lines = lines

        # If material is last datacard
        if mx_found:
            materials[-1].add_mx(mx_cards)

        return cls(materials)

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

    def to_xml(self, libmanager: LibManager) -> str:
        """Generate an XML representation of materials and return it as a string.

        Parameters
        ----------
        libmanager :
            libmanager

        Returns
        -------
        str
            The XML representation of materials as a string.
        """

        # Create XML element to represent the collection of materials.
        material_tree = ET.Element("materials")

        for material in self.materials:
            material.to_xml(libmanager, material_tree)

        # Apply indentation to the generated XML data.
        indent(material_tree)

        return ET.tostring(material_tree, encoding="unicode", method="xml")

    def translate(
        self, newlib: str | dict, lib_manager: LibManager, code: str = "mcnp"
    ) -> None:
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

            # Rebuild elements
            for submat in material.submaterials:
                submat._collapse_zaids()

        # this is a lazy fix, if performance issues are encountered the
        # update_zaid_info method of submaterials should be looked at
        self.update_info(lib_manager)

    def update_info(self, lib_manager: LibManager) -> None:
        """
        This methods allows to update the in-line comments for every zaids
        containing additional information

        Parameters
        ----------
        lib_manager : libmanager.Libmanager
            Library manager for the conversion.

        Returns
        -------
        None.

        """
        for mat in self.materials:
            mat._update_info(lib_manager)

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

        if re.match(r"^M\d{1,7}$", mat_name) is None:
            print("\nMaterial name not valid, set to M1\n")
            mat_name = "M1"

        # Translate to requested lib
        self.translate(newlib, libmanager)

        # Collect all submaterials
        submaterials = []
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
            current_submaterials = []
            for j, submat in enumerate(material.submaterials):
                # normalized & scaled
                norm_factor = float(percentage) / totfraction
                if fractiontype == "mass":
                    norm_factor = -norm_factor
                submat.scale_fractions(norm_factor)

                if submat.header is None:
                    comment = "C no submat header"
                else:
                    comment = submat.header
                # Add info to the header in order to back-trace the generation
                submat.header = f"C {materialname}, submaterial {j + 1}\n{comment}"

                # Drop additional keys if present
                submat.additional_keys = []
                current_submaterials.append(submat)

            # Change the header of the first submaterial to include the mat. 1
            new_sub_header = (
                str(material.header).strip("\n") + "\n" + current_submaterials[0].header
            ).strip("\n")

            current_submaterials[0].header = new_sub_header
            submaterials.extend(current_submaterials)

        # Generate new material and matlist
        newmat = Material(
            None, None, mat_name, submaterials=submaterials, header=main_header
        )
        newmat._update_info(libmanager)

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
