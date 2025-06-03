"""
This module is related to the parsing and manipulation of MCNP input files.

The parser is built on the numjuggler python module.

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

import json
import logging
import os
import re
from copy import deepcopy
from typing import Sequence

import numpy as np
import pandas as pd
from numjuggler import likefunc as lf
from numjuggler import parser

from f4enix.constants import (
    PAT_CARD_KEY,
    PAT_COMMENT,
    PAT_F_TR_CARD_KEY,
    PAT_FMESH_KEY,
    PAT_NP,
    UNION_INTERSECT_SYMBOLS,
)
from f4enix.input.auxiliary import debug_file_unicode
from f4enix.input.d1suned import IrradiationFile, Reaction, ReactionFile
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import MatCardsList, Material

PAT_MT = re.compile(r"m[tx]\d+", re.IGNORECASE)
PAT_BLANK_LINE = re.compile(r"\n[\s\t]*\n")
ADD_LINE_FORMAT = "         {}\n"


class CardsDict(dict):
    def __init__(self, card_type: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.card_type = card_type

    def update(self, *args, **kwargs):
        """Override the update method to customize how entries are added."""
        for key, value in dict(*args, **kwargs).items():
            self[key] = value

    def __setitem__(self, key, value):
        """Override the __setitem__ method to customize how entries are added."""

        # Ensure the key is a string
        if not isinstance(key, str):
            raise ValueError("Keys must be strings")

        # Ensure the value is an instance of parser.Card
        if not isinstance(value, parser.Card):
            raise ValueError("Values must be instances of parser.Card")

        # Ensure that the key matches the card number
        if "*" in key:
            num = key[1:]
        else:
            num = key
        if num != str(value.name):
            logging.debug(
                f"The card number {value.name} does not match the key {key}, it will be overridden"
            )
            value.name = int(num)
            value._set_value_by_type(self.card_type, int(num))

        # Call the superclass method to actually add the entry
        super().__setitem__(key, value)


class Input:
    def __init__(
        self,
        cells: dict[str, parser.Card],
        surfs: dict[str, parser.Card],
        materials: MatCardsList,
        transformations: dict[str, parser.Card],
        other_data: dict[str, parser.Card],
        tally_keys: list[int],
        fmesh_keys: list[int],
        header: str,
    ) -> None:
        """Class representing an MCNP input file.

        cells, surfaces, materials and transformations are handled explicitly.
        All other datacards are treated generically for the moment being.
        The parsing is built on the numjuggler python module.

        Parameters
        ----------
        cells: dict[str, parser.Card]
            cleaned numjuggler cards for each cells in the input. keys are the
            number of the cells.
        surfs: dict[str, parser.Card]
            cleaned numjuggler cards for each surface in the input. keys are
            the number of the surfaces.
        materials: MatCardsList
            material cards section of the input.
        transformations: dict[str, parser.Card]
            list of transformation datacards (i.e. TRn)
        other_data: dict[str, parser.Card]
            list of all remaining datacards that are treated in a generic way
        tally_keys: list[int]
            ids of the tallies available in the input
        fmesh_keys: list[int]
            ids of the fmeshes available in the input
        header: str
            header of the input file

        Attributes
        ----------
        cells: dict[str, parser.Card]
            cleaned numjuggler cards for each cells in the input. keys are the
            number of the cells.
        surfs: dict[str, parser.Card]
            cleaned numjuggler cards for each surface in the input. keys are
            the number of the surfaces.
        materials: MatCardsList
            material cards section of the input.
        transformations: dict[str, parser.Card]
            list of transformation datacards (i.e. TRn)
        other_data: dict[str, parser.Card]
            list of all remaining datacards that are treated in a generic way
        tally_keys: list[int]
            ids of the tallies available in the input
        fmesh_keys: list[int]
            ids of the fmeshes available in the input
        header: str
            header of the input file

        Examples
        --------
        The most common way to initiliaze an Input object is from a MCNP input
        file

        >>> from f4enix.input.MCNPinput import Input
        ... # Read the input file
        ... inp = Input.from_input(inpfile)

        >>> inp.cells
        ... # inp.surfs
        {'1': <numjuggler.parser.Card at 0x26535f02a70>,
         '2': <numjuggler.parser.Card at 0x26535f03130>,
         ...
         '128': <numjuggler.parser.Card at 0x26535f48f70>}

        The input can be translated to another library and rewritten to a file

        >>> from f4enix.input.libmanager import LibManager
        ... # Initialize a default nuclear data libraries manager
        ... libmanager = LibManager()
        ... # Translate the input to another library
        ... inp.translate('21c', libmanager)
        ... inp.write(outfile_path)

        Retrieve (and possibly modify) different cards in the input

        >>> print('Cell number 1:')
        ... print(inp.get_cells_by_id([1]))
        ... print('Surfaces number 10 and 20:')
        ... print(inp.get_surfs_by_id([10, 20]))
        ... print('All cells with material M1')
        ... print(inp.get_cells_by_matID(1))
        ... # use this method only if the previous ones
        ... # are not enough
        ... print('Generic way to obtain cards')
        ... print(inp._get_cards_by_id(['SDEF', 'IMP:N,P'], inp.other_data))
        Cell number 1:
        {'1': <numjuggler.parser.Card object at 0x0000026535F02A70>}
        Surfaces number 10 and 20:
        {'10': <numjuggler.parser.Card object at 0x0000026535F491B0>,
        '20': <numjuggler.parser.Card object at 0x0000026535F49390>}
        All cells with material M1
        {'22': <numjuggler.parser.Card object at 0x00000265360E7850>}
        Generic way to obtain cards
        {'SDEF': <numjuggler.parser.Card object at 0x0000026535F4B190>,
        'IMP:N,P': <numjuggler.parser.Card object at 0x0000026535F4A890>}

        Extract a subset of cells depending on some condition and from that
        dump out a minimal working MCNP input that includes all necessary
        surfaces, transformations and materials.

        >>> # --- Extract cells based on material ---
        ... inp = Input.from_input(inpfile)
        ... cells_ids = []  # store here the cells to be extracted
        ... selected_mat = 11
        ... for key, cell in inp.cells.items():
        ...     # get the material of the cell using numjuggler API
        ...     mat = cell._get_value_by_type('mat')
        ...     if mat == selected_mat:
        ...         cells_ids.append(key)
        ... print(cells_ids)
        ... # extract the cells subset to a file
        ... inp.extract_cells(cells_ids, 'outfile.i')

        Get useful summary of the material section of the input in a
        pandas.DataFrame object.

        >>> from f4enix.input.MCNPinput import Input
        ... from f4enix.input.libmanager import LibManager
        ... # Initialize a default nuclear data libraries manager
        ... libmanager = LibManager(xsdir_file='my_xsdir')
        ... inp = Input.from_input(inpfile)
        ... inp.materials.get_info(libmanager, complete=True)
        (                              Atom Fraction  Mass Fraction
        Material Submaterial Element
        m1       1           H             0.021630      -0.019046
                 2           C             0.018920      -0.198524
                 3           N             0.002060      -0.025207
                 4           O             0.027060      -0.378226
                 5           Mg            0.001190      -0.025268
        ...                                     ...            ...
        M29      1           Nb            0.000005      -0.000100
                             Co            0.000041      -0.000500
                             Fe            0.055451      -0.648443
        M30      1           H             0.063398      -0.112140
                             O             0.031622      -0.887860
        [333 rows x 2 columns],
                                    Fraction  Sub-Material Fraction  \
        Material Submaterial Element
        M24      1           Ag       0.000016               0.000016
                             Al       0.000631               0.000631
                             As       0.000023               0.000023
                             Au       0.000009               0.000009
                             Bi       0.000163               0.000163
        ...                                ...                    ...
        m7       2           O        0.033428               1.000000
        m8       1           Be       0.002970               0.034219
        ...
                             Cu                0.944769
                             Ni                0.021012
        m9       1           Be                1.000000
        [333 rows x 3 columns])


        """
        self.cells = cells
        self.surfs = surfs
        self.materials = materials
        self.transformations = transformations
        self.other_data = other_data
        self.tally_keys = tally_keys
        self.fmesh_keys = fmesh_keys
        self.header = header

        # # store also all the original cards for compatibility
        # # with some numjuggler modes
        # cells.extend(surfs)
        # cells.extend(data)

    @property
    def cells(self) -> dict[str, parser.Card]:
        """Get the cells attribute."""
        return self._cells

    @cells.setter
    def cells(self, value: dict[str, parser.Card]) -> None:
        """Set the cells attribute."""
        # enforce that the keys are all corresponding to the cell number
        # in case of difference, the dictionary key takes precedence
        for key, card in value.items():
            if key != str(card.name):
                logging.debug(
                    f"The cell number {card.name} is different from Key {key} and will be replaced"
                )
                card.name = int(key)
                card._set_value_by_type("cel", int(key))

        self._cells = CardsDict("cel", value)

    @property
    def surfs(self) -> dict[str, parser.Card]:
        """Get the surfs attribute."""
        return self._surfs

    @surfs.setter
    def surfs(self, value: dict[str, parser.Card]) -> None:
        """Set the surfs attribute."""
        # enforce that the keys are all corresponding to the surface number
        # in case of difference, the dictionary key takes precedence
        for key, card in value.items():
            if "*" in key:
                num = key[1:]
            else:
                num = key
            if num != str(card.name):
                logging.debug(
                    f"The surface number {card.name} is different from Key {key} and will be replaced"
                )
                card.name = int(num)
                card._set_value_by_type("sur", int(num))

        self._surfs = CardsDict("sur", value)

    @classmethod
    def from_input(cls, inputfile: os.PathLike | str) -> Input:
        """Generate an Input object from an MCNP input text file using
        numjuggler as a parser

        Parameters
        ----------
        inputfile : os.PathLike | str
            input file

        Returns
        -------
        Input
            Input object
        """
        cells, surfaces, data, header = _get_input_arguments(inputfile)
        cells = cls._to_dict(cells)
        surfs = cls._to_dict(surfaces)

        (
            materials,
            transformations,
            other_data,
        ) = cls._parse_data_section(cls, data)

        # get a list of the tally keys
        tally_keys = []
        fmesh_keys = []
        for key, card in other_data.items():
            if card.dtype == "Fn":
                tally_keys.append(card.name)
            elif PAT_FMESH_KEY.match(key):
                fmesh_keys.append(card.name)

        return cls(
            cells,
            surfs,
            materials,
            transformations,
            other_data,
            tally_keys,
            fmesh_keys,
            header,
        )

    def write(self, outfilepath: os.PathLike | str, wrap: bool = False) -> None:
        """write the input to a file

        Parameters
        ----------
        outfilepath : os.PathLike | str
            path to the output file
        wrap : bool
            if true the text is wrapped at 80 char. May cause slowdowns
        """
        logging.info("Writing to {}".format(outfilepath))

        self.write_blocks(
            outfilepath,
            wrap,
            self.cells,
            self.surfs,
            self.materials,
            self.header,
            self.transformations,
            self.other_data,
        )

        logging.info("File was written correctly")

    def merge(self, other_inp: Input, ensure_updated_dicts: bool = False) -> None:
        """Merge the input with another input object

        Parameters
        ----------
        other_inp : Input
            input to be merged
        ensure_updated_dicts : bool
            if True, the keys of both inputs are updated. This may not be needed
            in certain applications and can be set to false if it is sure that
            the keys are already up to date. Default is False.
        """
        if ensure_updated_dicts:
            self._update_card_keys()
            other_inp._update_card_keys()

        self._safe_dict_update(self.cells, other_inp.cells)
        self._safe_dict_update(self.surfs, other_inp.surfs)
        self.materials.extend(other_inp.materials.materials)
        self._safe_dict_update(self.transformations, other_inp.transformations)
        # I do not want to stop merging for 2 SDEF cards for instance
        self.other_data.update(other_inp.other_data)
        # ignore duplicated keys for these for the moment being
        self.tally_keys.extend(other_inp.tally_keys)
        self.fmesh_keys.extend(other_inp.fmesh_keys)

    @staticmethod
    def _safe_dict_update(original_dict: dict, ext_dict: dict):
        for key, val in ext_dict.items():
            if key in original_dict.keys():
                raise KeyError("Duplicated card entry: " + key)
            original_dict[key] = val

    def _update_card_keys(self) -> None:
        """This function is pretty costly but it allows to ensure that cards
        values and cards keys are consistent. It useful for instance after
        a renumbering operation
        """
        newcells = {}
        newsurfs = {}
        newdata = {}
        newtrans = {}
        for cards, newset in zip(
            [self.cells, self.surfs, self.other_data, self.transformations],
            [newcells, newsurfs, newdata, newtrans],
        ):
            for card in cards.values():
                key = _get_card_key(card)
                newset[key] = card
        # update all dicts
        self.cells = newcells
        self.surfs = newsurfs
        self.other_data = newdata
        self.transformations = newtrans

    def renumber(
        self,
        cells: int | None = None,
        surfs: int | None = None,
        universes: int | None = None,
        translations: int | None = None,
        renum_all: int | None = None,
        update_keys: bool = False,
    ) -> None:
        """Renumber cards of the input files. Either all cards can be renumbered
        using the renum_all parameter or onnly a subset between cells, surfs,
        universes and translations. Materials are not supported for the time
        being.

        Parameters
        ----------
        cells : int, optional
            offset to be used for cells, by default None
        surfs : int, optional
            offset to be used for surfaces, by default None
        universes : int, optional
            offset to be used for universes, by default None
        translations : int, optional
            offset to be used for translations, by default None
        renum_all : int, optional
            offset to be used for all supported cards, by default None. It will
            trump all other specified offsets.
        update_keys : bool, optional
            if True, the keys of the different dictionaries are updated according
            to the new numbering. This is useful if more operations are needed
            on the file. If the file will be directly written after renumbering
            though, this costly operation is unnecessary. Default is False.
        """
        if renum_all is not None:
            # assign the same value to all the other variables
            cells = surfs = universes = translations = renum_all
        # for the moment do not care about logging the changes
        # see here https://github.com/travleev/numjuggler/blob/e7659eb5abe54d84e3982e8fa0775ad5caf3a04a/numjuggler/main.py#L1423
        maps = {}
        for offset, card_type in zip(
            [cells, surfs, universes, translations],
            ["cel", "sur", "u", "tr"],
        ):
            if offset is not None:
                # renumber only requested cards
                maps[card_type] = lf.LikeFunction()
                maps[card_type].default = lf.add_func(int(offset))
                # do not modify zero numbers (important for material
                # numbers)
                maps[card_type].mappings[lf.Range(0)] = lf.const_func(0)
                maps[card_type].doc = f"Function for {offset} from command line"
        for cards in [self.cells, self.surfs, self.other_data, self.transformations]:
            for card in cards.values():
                card.apply_map(maps)

        if update_keys:
            self._update_card_keys()

    def translate(self, newlib: str | dict, libmanager: LibManager) -> None:
        """
        Translate the input to another library

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

        libmanager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        None.

        """

        try:
            if isinstance(newlib, str) and newlib[0] == "{":
                # covert the dic
                newlib = json.loads(newlib)
        except KeyError:
            # It is already a dict, pass
            pass

        self.update_zaidinfo(libmanager)
        self.materials.translate(newlib, libmanager)

    def update_zaidinfo(self, lib_manager: LibManager):
        """
        This methods allows to update the in-line comments for every zaids
        containing additional information

        Parameters
        ----------
        lib_manager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        None.

        """

        self.materials.update_info(lib_manager)

    @staticmethod
    def set_cell_void(cell: parser.Card) -> None:
        if cell.ctype == 3 and cell.get_m() != 0:
            cell.hidden["~"][0] = ""
            cell._set_value_by_type("mat", 0)
            cell._Card__m = 0  # necessary for the get val
        else:
            logging.warning(f"cell {cell.name} is either already void or not a cell")

    @staticmethod
    def _print_cards(cards: dict[str, parser.Card], wrap: bool = False) -> list[str]:
        text = []
        for _, card in cards.items():
            text_candidate = card.card(wrap=wrap).strip("\n") + "\n"
            # delete all '\r' special characters
            text_candidate = text_candidate.replace("\r", "")
            # avoid blank lines
            text_candidate = PAT_BLANK_LINE.sub("\n", text_candidate)
            text.append(text_candidate)
        return text

    @staticmethod
    def _to_dict(cards: list[parser.Card]) -> dict[str, parser.Card]:
        new_cards = {}
        flag_add = False
        for card in cards:
            card.get_values()
            try:
                key = card.name
                key = _get_card_key(card)
                if key in new_cards.keys():
                    raise KeyError("Duplicated card entry: " + key)

            except AttributeError:
                # This means that this is a fake card just made by comments
                # it should be merged with the following card

                # this is true for comments, but sometimes it happens also
                # with real cards due to bugs in numjuggler

                # let's first check if it is a comment
                if PAT_COMMENT.match(card.lines[0]) is not None:
                    comment = card.lines
                    flag_add = True
                    continue
                # and then if it is a proper card
                else:
                    key = card.card().split()[0].upper()
                    if key in new_cards.keys():
                        raise KeyError("Duplicated card entry: " + key)

            if flag_add:
                comment.extend(card.lines)
                card.lines = comment
                card.get_input()
                card.get_values()
                flag_add = False

            new_cards[key] = card

        return new_cards

    @staticmethod
    def _get_cards_by_id(ids: list[str], cards: dict) -> dict[str, parser.Card]:
        selected_cards = {}
        for id_card in ids:
            try:
                selected_cards[id_card] = cards[id_card]
            except KeyError:
                try:
                    # sometimes it may be with an asterisk
                    selected_cards["*" + id_card] = cards["*" + id_card]
                except KeyError:
                    raise KeyError("The card is not available in the input")

        return selected_cards

    def get_cells_by_id(
        self, ids: Sequence[int | str], make_copy: bool = False
    ) -> dict[str, parser.Card]:
        """given a list of cells id return a dictionary of such cells

        Parameters
        ----------
        ids : list[int | str]
            cells id to be extracted
        make_copy : bool
            if True, makes a deepcopy of the cells instead of working on the
            original ones. Default is False.

        Returns
        -------
        dict
            extracted cells
        """
        str_ids = []
        for id in ids:
            str_ids.append(str(id))
        if make_copy:
            return deepcopy(self._get_cards_by_id(str_ids, self.cells))
        else:
            return self._get_cards_by_id(str_ids, self.cells)

    def get_surfs_by_id(self, ids: list[int]) -> dict[str, parser.Card]:
        """given a list of surfaces id return a dictionary of such surfaces

        Parameters
        ----------
        ids : list[int]
            cells id to be extracted

        Returns
        -------
        dict
            extracted surfaces
        """
        str_ids = []
        for id in ids:
            str_ids.append(str(id))
        return self._get_cards_by_id(str_ids, self.surfs)

    def get_materials_subset(self, ids: list[str] | str) -> MatCardsList | Material:
        """given a list of material ids generate a new MatCardsList with
        the requested subset

        Parameters
        ----------
        ids : Union(list[str]), str)
            ids of the materials to put into the subset

        Returns
        -------
        MatCardsList | Material
            new materials subset. A single material is returned if only one was
            requested
        """
        if type(ids) is str:
            return self.materials[ids.upper()]
        else:
            materials = []
            for id_mat in ids:
                materials.append(self.materials[id_mat.upper()])
            return MatCardsList(materials)

    def get_data_cards(self, ids: list[str] | str) -> dict[str, parser.Card]:
        """Get a tranformation card or an other data card by its key.

        For the moment, transformation cards mixed with other data cards is
        not supported.

        Parameters
        ----------
        ids : list[str] | str
            keys of the cards to retrieve

        Returns
        -------
        dict[str, parser.Card]
            retrieved cards
        """
        if isinstance(ids, str):
            ids2use = [ids]
        else:
            ids2use = ids

        try:
            cards = self._get_cards_by_id(ids2use, self.other_data)
        except KeyError:
            cards = self._get_cards_by_id(ids2use, self.transformations)

        return cards

    def _parse_data_section(
        self, cardslist: list[parser.Card]
    ) -> tuple[MatCardsList, dict[str, parser.Card], dict[str, parser.Card]]:
        # first of all correct numjuggler parser
        cards = self._to_dict(cardslist)

        materials = []  # store here the materials
        transformations = {}  # store translations
        other_data = {}  # store here the other data cards

        for key, card in cards.items():
            key = self._clean_card_name(key)
            try:
                if card.values[0][1] == "mat":
                    if PAT_MT.match(card.lines[0]) or PAT_MT.match(card.lines[-1]):
                        # mt or mx cards should be added to the previous
                        # material
                        materials[-1].add_mx(card)
                    else:
                        materials.append(Material.from_text(card.lines))

                elif card.dtype == "TRn":
                    transformations[key] = card

                else:
                    other_data[key] = card

            except IndexError:
                # this means that there were no values
                other_data[key] = card

        return MatCardsList(materials), transformations, other_data

    def extract_cells(
        self,
        cells: list[int],
        outfile: os.PathLike | str,
        renumber_offsets: dict | None = None,
        keep_universe: bool = True,
        extract_fillers: bool = True,
    ):
        """given a list of cells, dumps a minimum MCNP working file.

        The file will includes all the requested cells, defined surfaces,
        materials and translations.

        Parameters
        ----------
        cells : list[int]
            desired list of cells
        outfile : os.PathLike | str
            path to the file where the MCNP input needs to be dumped
        renumber_offsets : dict, optional
            apply the self.renumber() function to the extracted input.
            the dict will be passed as keyargs to the function.
            Default is None.
        keep_universe: bool
            If True keeps the 'U=' key in the cell cards, otherwise that is
            removed. Default is True.
        extract_fillers: bool
            if True extract also the cells belonging to a universe that is
            used in a 'FILL=' keyword. This happens recursively. default is
            True.
        """
        logging.info("write MCNP reduced input")

        if renumber_offsets is not None or not keep_universe:
            make_copy = True
        else:
            make_copy = False

        cells_cards, surfs, materials = self._extraction_function(
            cells, keep_universe, extract_fillers, make_copy=make_copy
        )

        newinput = Input(
            cells_cards,
            surfs,
            materials,
            deepcopy(self.transformations),
            deepcopy(self.other_data),
            self.tally_keys,
            self.fmesh_keys,
            self.header,
        )

        if renumber_offsets is not None:
            newinput.renumber(**renumber_offsets)

        newinput.write(outfile)

    @staticmethod
    def write_blocks(
        file: os.PathLike,
        wrap: bool,
        cells_cards: dict[str, parser.Card],
        surfs: dict[str, parser.Card],
        materials: MatCardsList,
        header: list[str] | None = None,
        trans: dict[str, parser.Card] | None = None,
        other_data: dict[str, parser.Card] | None = None,
    ):
        """Writes F4Enix dicts of cells, surfaces and data cards.
        The method receives cells, surfaces, materials F4Enix dicts and
        optionally header, transformation and other data F4Enix dicts and
        prints the MCNP input

        Parameters
        ----------
        file : os.PathLike
            path of the MCNP input that will be printed
        wrap : bool
            flag to check if the input should be wrapped to 80 characters per
            line
        cells_cards : dict[str, parser.Card]
            F4Enix dict of cells
        surfs : dict[str, parser.Card]
            F4Enix dict of surfaces
        materials : MatCardsList
           MatCardsList object including the materials objects to be printed
        header : list[str], optional
            list of lines of header of MCNP input, by default None
        trans : dict[str, parser.Card], optional
            F4Enix dict of transformations, by default None
        other_data : dict[str, parser.Card], optional
            fEnix dict of MCNP data cards, by default None
        """

        with open(file, "w") as outfile:
            # Add the header lines
            if header is not None:
                for line in header:
                    # remove special '\r' characters
                    newline = line.replace("\r", "")
                    outfile.write(newline)
            else:
                outfile.write("C\n")
            # Add the cells
            outfile.writelines(Input._print_cards(cells_cards, wrap=wrap))
            # Add a break
            outfile.write("\n")
            # Add the surfaces
            outfile.writelines(Input._print_cards(surfs, wrap=wrap))
            # Add a break
            outfile.write("\n")
            # Add materials
            if trans is not None:
                outfile.writelines(Input._print_cards(trans))
            if materials is not None and len(materials.matdic) > 0:
                outfile.write(materials.to_text() + "\n")
            # other data is not mandatory to be written
            if other_data is not None:
                outfile.writelines(Input._print_cards(other_data, wrap=wrap))

    def _extraction_function(
        self,
        cells: list[int],
        keep_universe: bool = True,
        extract_fillers: bool = True,
        make_copy: bool = False,
    ) -> tuple[dict[str, parser.Card], dict[str, parser.Card], MatCardsList]:
        logging.info("Collecting the cells, surfaces, materials and transf.")
        cset = set(cells)

        # first, get all surfaces needed to represent the cn cell.
        sset = set()  # surfaces
        mset = set()  # material
        # tset = set()  # transformations
        self._collect_hash_uni(cset, extract_fillers)

        # sort the set
        cset = list(cset)
        try:
            cset.sort()
        except TypeError:
            # if the list is not sortable, it means that it is a list of strings
            # and we need to convert it to integers
            logging.warning(
                "The list of cell ids are not int (strings?). Trying to convert it."
            )
            cset = [int(c) for c in cset]
            cset.sort()

        # create a copy if modifications are needed on the cells
        if make_copy:
            cells_cards = self.get_cells_by_id(cset, make_copy=True)
        else:
            cells_cards = self.get_cells_by_id(cset)

        # Get all surfaces and materials
        for i, (_, cell) in enumerate(cells_cards.items()):
            for v, t in cell.values:
                if t == "sur":
                    sset.add(v)
                elif t == "mat":
                    if int(v) != 0:  # void material is not defined in a card
                        mset.add("M" + str(v))

            if not keep_universe and cell.values[0][0] in cells:
                Input.remove_u(cells_cards[_])

        # Do not bother for the moment in selecting also the transformations

        #                     elif t == 'tr':
        #                         tset.add(v)

        # # final run: for all cells find surfaces, materials, etc.
        # for key, surf in self.surfs.items():
        #     if key in sset:
        #         # surface card can refer to tr
        #         for v, t in surf.values:
        #             if t == 'tr':
        #                 tset.add(v)
        # order surfaces
        sset = list(sset)
        sset.sort()
        # get surfaces dict
        surfs = self.get_surfs_by_id(sset)
        # get materials dict
        materials = self.get_materials_subset(mset)
        # this needs to be always a mat list, even if only one material is available
        if isinstance(materials, Material):
            materials = MatCardsList([materials])

        return cells_cards, surfs, materials

    def _collect_hash_uni(self, cset: set, extract_fillers: bool):
        # duplicate the final set and work on a dynamic set that contains only
        # new cells at each loop
        cell_set = deepcopy(cset)

        # next runs: find all other cells:
        again = True
        while again:
            again = False
            new_set = set()
            uni_set = set()
            # loop over cells to extract
            for cell_num in cell_set:
                c = self.cells[str(cell_num)]
                # get hash cells in the cells that have to be extracted
                cref = c.get_refcells()
                # add the hash cells to extraction list
                new_set |= cref
                # collect universes in the definition of cells
                if extract_fillers:
                    fill = c.get_f()
                    if fill is not None:
                        uni_set.add(fill)
            # if one wants to extract also lower levels, loop over universes
            # and collect their cells
            if extract_fillers:
                for _, c in self.cells.items():
                    if c.get_u() in uni_set:
                        new_set.add(c.values[0][0])
            # get the new set with the cells to be checked
            cell_set = new_set - cell_set
            # check if loop is to be repeated
            if cell_set:
                again = True
                cset |= cell_set

    def extract_universe(
        self,
        universe: int,
        outfile: os.PathLike | str,
        renumber_offsets: dict | None = None,
        keep_universe: bool = False,
    ):
        """Dumps a minimum MCNP working file that
        includes all the cells, surfaces, materials and
        translations of the universe. The resulting file doesn't have the universe
        keyword in the cell definitions

        Parameters
        ----------
        universe : int
            universe id to be extracted
        outfile : os.PathLike | str
            path to the file where the MCNP input needs to be dumped
        renumber_offsets : dict, optional
            apply the self.renumber() function to the extracted input.
            the dict will be passed as keyargs to the function.
            Default is None.
        keep_universe : bool
            determines if the u=... card should be kept or not in cells'
            definitions. Defult is False.
        """
        cell_ids_to_extract = []
        for cell_id, cell in self.cells.items():
            cell_universe = cell.get_u()

            if cell_universe == universe:
                cell_ids_to_extract.append(cell.values[0][0])

        self.extract_cells(
            cells=cell_ids_to_extract,
            outfile=outfile,
            renumber_offsets=renumber_offsets,
            keep_universe=keep_universe,
        )

    @staticmethod
    def _clean_card_name(key: str) -> str:
        # this is to clean cases like:
        # *TR1 -> TR1
        # F6:N,P -> F6
        key = key.upper()
        try:
            newkey = PAT_F_TR_CARD_KEY.search(key).group()
            # handle the TF edge case
            if "TF" in key.upper():
                newkey = "T" + newkey
        except AttributeError:
            logging.debug("the following key was not cleaned: " + key)
            newkey = key

        return newkey

    def get_cells_by_matID(
        self, matID: int | str, deepcopy_flag: bool = True
    ) -> dict[str, parser.Card]:
        """Given a material ID return a dictionary {key, card} of all
        the cells to which that material is assigned to.

        The cells that are returned are deepcopies of the original ones.

        Parameters
        ----------
        matID : int | str
            material ID to filter the cells
        deepcopy_flag: bool
            if False, the cells are not copied. Default is True

        Returns
        -------
        dict[int, parser.Card]
            cells to which the material is assigned to
        """
        logging.debug("get cells for material {} requested".format(matID))
        filtered_cells = {}
        for key, cell in self.cells.items():
            if cell._get_value_by_type("mat") == int(matID):
                if deepcopy_flag:
                    filtered_cells[key] = deepcopy(cell)
                else:
                    filtered_cells[key] = cell
        return filtered_cells

    def scale_densities(self, factor: float) -> None:
        """Scale the density values of all cells by the same factor. Void
        cells are ignored. Resulting density will be equal to
        original_density*factor

        Parameters
        ----------
        factor : float
            scaling factors for the densities
        """
        for _, cell in self.cells.items():
            if not cell._get_value_by_type("mat") == 0:
                density = cell.get_d()
                newdensity = "{:.5e}".format(density * factor)
                cell.set_d(newdensity)

    def get_densities_range(self) -> pd.DataFrame:
        """Return a DataFrame listing for all material the minimum and maximum
        density values encountered in the input cells.

        Returns
        -------
        pd.DataFrame
            Range of densities used for each material.
        """
        lm = LibManager()
        materials = []
        for _, cell in self.cells.items():
            mat_id = cell.get_m()
            if mat_id == 0:
                continue  # ignore void cells
            dens = float(cell.get_d())

            # if the density is negative (mass) leave it as it is,
            # if atomic fraction is used, convert it to mass first
            if dens < 0:
                dens = abs(dens)
            else:
                dens = self.materials[f"M{mat_id}"].get_density(dens, lm)

            materials.append([mat_id, dens])

        df = pd.DataFrame(materials, columns=["Material ID", "Density"])
        min_vals = df.groupby("Material ID")["Density"].min()
        max_vals = df.groupby("Material ID")["Density"].max()
        summary = pd.DataFrame()
        summary["Min density [g/cc]"] = min_vals
        summary["Max density [g/cc]"] = max_vals
        summary["Material ID"] = max_vals.index
        return summary.set_index("Material ID").sort_index()

    def get_cells_summary(self) -> pd.DataFrame:
        """Get a summary of infos for each cell

        A DataFrame is returned where for each cell is listed the material,
        density, universe and filler is present.

        Returns
        -------
        pd.DataFrame
            Summary of cells info
        """
        rows = []
        for key, cell in self.cells.items():
            row = {"cell": int(key)}
            row["material"] = cell.get_m()
            row["density"] = cell.get_d()
            row["universe"] = cell.get_u()
            row["filler"] = cell.get_f()
            rows.append(row)

        df = pd.DataFrame(rows)
        return df.set_index("cell").sort_index()

    def _get_tally_cards(
        self,
        idx: int,
    ) -> list[str]:
        keys = []
        pat = re.compile(r"F[a-zA-Z]*{}$".format(idx))
        for key, _ in self.other_data.items():
            if pat.match(key) is not None:
                keys.append(key)
        return keys

    def _retrieve_input(self, tag: str) -> str:
        # get the FC comment excluding the FC tag
        comment_line = self.other_data[tag].input[0]
        inp = comment_line.replace(tag, "").strip()
        inp = inp.replace(tag.lower(), "").strip()
        return inp

    def get_tally_summary(self, fmesh: bool = False) -> pd.DataFrame:
        """Get a summary of the tallies defined in the input

        Both normal tallies and fmeshes can be requested. For each tally the
        number, particle, description and multiplier are listed (if available)

        Parameters
        ----------
        fmesh : bool, optional
            if True produced a summary for the fmehses instead of for the
            normal tallies, by default False

        Returns
        -------
        pd.DataFrame
            summary info on defined tallies

        """

        if fmesh:
            tag_tally = "FMESH"
            tallies = self.fmesh_keys
        else:
            tag_tally = "F"
            tallies = self.tally_keys

        rows = []
        for key in tallies:
            desc = pd.NA
            particle = pd.NA
            multiplier = None
            card_keys = self._get_tally_cards(key)
            for aux_key in card_keys:
                if aux_key[:2] == "FC":
                    desc = self._retrieve_input(aux_key)
                elif aux_key == tag_tally + str(key):
                    line = self.other_data[aux_key].input[0]
                    particle = PAT_NP.search(line).group().upper().strip(":")
                elif aux_key[:2] == "FM":
                    multiplier = self._retrieve_input(aux_key).split()

            row = {"Tally": key, "Particle": particle, "Description": desc}

            if multiplier is not None:
                row["Normalization"] = multiplier[0]
                if len(multiplier) > 1:
                    row["Other multipliers"] = multiplier[1:]
                else:
                    row["Other multipliers"] = pd.NA
            else:
                row["Normalization"] = pd.NA
                row["Other multipliers"] = pd.NA

            rows.append(row)

        return pd.DataFrame(rows).set_index("Tally").sort_index()

    def replace_material(
        self,
        new_mat_id: int,
        new_density: str,
        old_mat_id: int,
        u_list: list[int] | None = None,
    ) -> None:
        """Replace a material and density in the input with other values.

        Parameters
        ----------
        new_mat_id : int
            id of the new material (0 for void)
        new_density : str
            new value for the density (including sign)
        old_mat_id : int
            id of the material to be replaced
        u_list : list[int]
            change the material only if cells belong to one of the universes
            in the list. By default is None, all cells are affected.
        """

        if new_mat_id < 0 or new_mat_id < 0:
            raise ValueError("Wrong values for the material ids")
        for _, cell in self.cells.items():
            in_universe = False
            # check if universe is a parameter
            if u_list is not None:
                if cell.get_u() in u_list:
                    in_universe = True
            else:
                # always in universe, default is always
                in_universe = True

            # If the material needs change and in universe
            if cell.get_m() == old_mat_id and in_universe:
                self.replace_cell_material(cell, new_mat_id, new_density)

    @staticmethod
    def replace_cell_material(cell: parser.Card, new_mat_id: int, new_density: str):
        old_mat_id = cell.get_m()
        # Void needs to be handle in a specific way
        if old_mat_id == 0 and new_mat_id == 0:
            logging.warning("Replacing void with void")
        if old_mat_id == 0:
            Input.add_material_to_void_cell(cell, new_mat_id, new_density)
        elif new_mat_id == 0:
            Input.set_cell_void(cell)
        else:
            cell._set_value_by_type("mat", new_mat_id)
            cell._Card__m = new_mat_id  # necessary for the get val
            cell.set_d(new_density)

    @staticmethod
    def add_material_to_void_cell(
        cell: parser.Card,
        new_mat_id: int,
        new_density: str,
    ) -> None:
        """Sets a material (and density) to a void cell.

        Parameters
        ----------
        cell : parser.Card
            cell to be modified.
        new_mat_id : int
            id of the new material.
        new_density : str
            new value for the density (including sign).
        """

        if int(float(new_density)) == 0 or new_mat_id <= 0:
            raise ValueError("Wrong values for the new material and density")

        if cell.ctype == 3 and cell.get_m() == 0:
            cell.hidden["~"].insert(0, str(new_density))
            cell._set_value_by_type("mat", new_mat_id)
            cell._Card__m = new_mat_id  # necessary for the get val
            cell._Card__d = float(new_density)
            # Introduce parentheses before the third word in the first row
            first_row = cell.input[0].split()
            first_row.insert(2, r"~")
            cell.input[0] = " ".join(first_row)
        else:
            logging.warning(f"cell {cell.name} is not a void cell")

    @staticmethod
    def add_cell_fill_u(
        cell: parser.Card,
        param: str,
        param_value: int,
        inplace: bool = True,
    ) -> parser.Card:
        """Adds a u=/fill= keyword parameter to a cell that doesn't have it.

        Parameters
        ----------
        cell : parser.Card
            numjuggler cell card to which the surface will be added
        param : str
            can be 'u' or 'fill', it tells the parameter to be added to the cell.
        param_value : int
            the value of the parameter 'u' or 'fill' to be added to the cell.
        inplace: bool
            if False a deepcopy is created. By default is True.

        Returns
        -------
        parser.Card
            numjuggler card of the modified cell
        """

        valid_params = ["u", "fill"]
        if param.lower() not in [p.lower() for p in valid_params]:
            raise ValueError(f"Invalid parameter {param}")
        if cell.get_u() and param.lower() == "u":
            raise ValueError("Cell already has a universe defined")
        if cell.get_f() and param.lower() == "fill":
            raise ValueError("Cell already has a filler defined")

        cell._Card__f = -1
        cell._Card__u = -1
        len_param = len(str(param_value))

        new_cell = Input._add_symbol_to_cell_input(
            cell,
            f"{param.lower()}={{:<{len_param}}}",
            inplace=inplace,
            parentheses=False,
        )

        for k in range(len(new_cell.values) - 1, -1, -1):
            if new_cell.values[k][1] == "sur" or new_cell.values[k][1] == "cel":
                break

        new_cell.values.insert(k + 1, (param_value, param.lower()))

        return new_cell

    @staticmethod
    def add_surface(
        cell: parser.Card,
        add_surface: int,
        new_cell_num: int | None = None,
        mode: str = "intersect",
        inplace: bool = True,
    ) -> parser.Card:
        """Adds a surface to cell's definition as union or intersection.

        Parameters
        ----------
        cell : parser.Card
            numjuggler cell card to which the surface will be added
        add_surface : int
            the surface number to be added to cell's definition. It should
            include the sign.
        new_cell_num : int, optional
            new number of the cell after the addition of the surface to cell's
            definition, by default None. If a new number is specified, the
            modifications are done on a copy of the original cell, otherwise
            these are done inplace.
        mode : str, optional
            can be 'union' or 'intersect', it tells the operation with which the
            surface is added to cell's definition, by default 'intersect'
        inplace: bool
            if False a deepcopy is created. By default is True.

        Returns
        -------
        parser.Card
            numjuggler card of the modified cell
        """
        if add_surface >= 0:
            template_addition = (
                UNION_INTERSECT_SYMBOLS[mode]
                + "{:<"
                + str(len(str(add_surface)))
                + "} "
            )
        else:
            template_addition = (
                UNION_INTERSECT_SYMBOLS[mode]
                + "-{:<"
                + str(len(str(add_surface)) - 1)
                + "} "
            )

        new_cell = Input._add_symbol_to_cell_input(
            cell, template_addition, inplace=inplace, new_cell_num=new_cell_num
        )

        for k in range(len(new_cell.values) - 1, -1, -1):
            if new_cell.values[k][1] == "sur" or new_cell.values[k][1] == "cel":
                break

        new_cell.values.insert(k + 1, (abs(add_surface), "sur"))

        return new_cell

    @staticmethod
    def hash_cell(
        cell: parser.Card,
        hash_id: int,
        new_cell_num: int = None,
        inplace: bool = True,
    ) -> parser.Card:
        """Hash a cell to a new cell number and update the hash dictionary.

        Parameters
        ----------
        cell : parser.Card
            numjuggler cell card to be hashed
        hash_id : int
            id of the hash cell.
        new_cell_num : int, optional
            new cell number to which the cell will be hashed. By default is None
        inplace : bool, optional
            if False a deepcopy is created. By default is True.
        """
        template_addition = "#{:<" + str(len(str(hash_id)) - 1) + "} "
        new_cell = Input._add_symbol_to_cell_input(
            cell, template_addition, inplace=inplace, new_cell_num=new_cell_num
        )
        # add the hash cell to the cell values
        for k in range(len(new_cell.values) - 1, -1, -1):
            if new_cell.values[k][1] == "sur" or new_cell.values[k][1] == "cel":
                break

        new_cell.values.insert(k + 1, (hash_id, "cel"))

        return new_cell

    def hash_multiple_cells(self, hash_dict: dict[int, list[int]]) -> None:
        """all keys in the hash dict correspond to hash cells to be added
        to the cells inlcuded in the hash dict value.

        Parameters
        ----------
        hash_dict : dict[int, list[int]]
            info on the cells to hash and with what

        """
        for hash_id, cells in hash_dict.items():
            for cell_num in cells:
                cell = self.cells[str(cell_num)]
                self.hash_cell(cell, hash_id, inplace=True)

    def cells_union(
        self,
        cell_num_list: list[str],
        new_cell_num: int = None,
    ) -> None:
        """Given a list of cells, it creates a new cell that is the union of
        the cells in the list, starting from the first cell in the list.
        The resulting cell will be put in the dict of cells
        of the input (with the new number if provided), while the old cells will
        be deleted. If a renumbering was done, the input must be re-read.

        Parameters
        ----------
        cell_list : list[parser.Card]
            list of cells to be united
        new_cell_num : int, optional
            new number of the union cell, by default None (old number is kept).
        """

        cell_list = []
        for cell_num in cell_num_list:
            cell_list.append(self.cells[cell_num])

        if new_cell_num is None:
            cel = cell_list[0]
            new_cell_num = cel.values[0][0]

        new_cell = cell_list[0]

        for cell in cell_list[1:]:
            new_cell = Input._add_symbol_to_cell_input(
                new_cell,
                " : (" + cell.get_geom().replace("\n", " ") + ")",
                inplace=False,
                parentheses=True,
                new_cell_num=new_cell_num,
            )
            cell_lines = new_cell.card(wrap=True).splitlines(keepends=True)
            new_cell = parser.Card(cell_lines, 3, new_cell.pos)
            new_cell.get_values()

        for cell_num in cell_num_list:
            self.cells.pop(cell_num)

        self.cells[str(new_cell_num)] = new_cell

    @staticmethod
    def _add_symbol_to_cell_input(
        cell: parser.Card,
        add_symbol: str,
        inplace: bool = True,
        parentheses: bool = True,
        new_cell_num: int | None = None,
    ) -> parser.Card:
        if inplace:
            new_cell = cell
        else:
            new_cell = deepcopy(cell)

        # Introduce parentheses before the third word in the first row
        first_row = new_cell.input[0].split()

        if parentheses:
            if new_cell.get_m() == 0:
                first_row[1] = first_row[1] + " ("
            else:
                first_row[2] = first_row[2] + " ("

        new_cell.input[0] = " ".join(first_row)

        # Check all rows if there are letters in the row
        for i, line in enumerate(new_cell.input):
            row = line.split()
            keywords = False
            for m, words in enumerate(row):
                if any(c.isalpha() for c in words):
                    param_cards_idx = m
                    keywords = True
                    break
            if not keywords:
                param_cards_idx = len(row)
            if keywords or (not keywords and i == len(new_cell.input) - 1):
                if parentheses:
                    row.insert(param_cards_idx, ") " + add_symbol)
                else:
                    row.insert(param_cards_idx, add_symbol)
                new_cell.input[i] = " ".join(row)

                if new_cell.input[i][:5] != "     " and i != 0:
                    new_cell.input[i] = "     " + new_cell.input[i]
                break

        # renumber the cell if requested
        if new_cell_num is not None:
            new_cell.name = new_cell_num
            new_cell._set_value_by_type("cel", new_cell_num)

        return new_cell

    @staticmethod
    def remove_u(cell: parser.Card) -> None:
        """given a cell, it removes the universe option from its definition.

        Parameters
        ----------
        cell : parser.Card
            cell from which the universe has to be removed

        """
        # initialize new input list
        new_input = []
        # remove universe option from input template
        for input_part in cell.input:
            new_input.append(re.sub(r"[uU]=\{:<\d+\}", "", input_part))

        # assign new input to cell
        cell.input = new_input
        # remove value associated to the universe in 'values'
        for b, (t, v) in enumerate(cell.values):
            if v == "u":
                cell.values.pop(b)
                break
        # reset universe private value (i know this is not a good practice, tbd)
        cell._Card__u = None

    def add_F_tally(
        self,
        tally_ID: int,
        particles: list[str],
        cells: Sequence[int | str],
        add_total: bool = False,
        energies: Sequence[float] | np.ndarray | None = None,
        description: str | None = None,
        multiplier: str | None = None,
        add_SD: bool = True,
    ):
        """Add a F-tally to the input

        Parameters
        ----------
        tally_ID : int
            tally ID
        particles : list[str]
            particle to tally
        cells : list[int | str]
            list of cells to tally. One can also add a string with a complex cell
            definition, like unions of cells and cells in a universe,
            e.g. "((1 2 3) < 10)"
        add_total : bool, optional
            if True adds the total tally, by default False
        energies : list[float], optional
            list of energies for the tally, by default None
        description : str, optional
            description of the tally, by default None
        multiplier : str, optional
            multiplier of the tally, by default None
        add_SD : bool, optional
            if True adds the SD 1 card, by default True
        """
        # Add the tally card
        tally_ID = int(tally_ID)
        # add the IDs to the tally ids
        self.tally_keys.append(tally_ID)

        particles_str = particles[0]
        if len(particles) > 1:
            for particle in particles[1]:
                particles_str += "," + particle
        # add description if available
        if description is not None:
            line = [f"FC{tally_ID} {description}\n"]
            self.other_data[f"FC{tally_ID}"] = parser.Card(line, 5, -1)
        # tally main body
        lines = [f"F{tally_ID}:{particles_str}\n"]
        lines.append("     ")
        # --- add cells ---
        if add_total:
            # ensure it is a list
            cells = list(cells)
            cells.append("T")
        for cell in cells:
            lines[1] += str(cell) + " "
        lines[1] += "\n"
        self.other_data[f"F{tally_ID}"] = parser.Card(lines, 5, -1)
        card_lines = (
            self.other_data[f"F{tally_ID}"].card(wrap=True).splitlines(keepends=True)
        )
        self.other_data[f"F{tally_ID}"] = parser.Card(card_lines, 5, -1)
        # add energies if requested
        if energies is not None:
            lines = [f"E{tally_ID}\n"]
            lines.append("     ")
            for energy in energies:
                lines[1] += f"{energy:.4e} "
            lines[1] += "\n"
            self.other_data[f"E{tally_ID}"] = parser.Card(lines, 5, -1)
            card_lines = (
                self.other_data[f"E{tally_ID}"]
                .card(wrap=True)
                .splitlines(keepends=True)
            )
            self.other_data[f"E{tally_ID}"] = parser.Card(card_lines, 5, -1)
        # add SD if requested
        if add_SD:
            repetitions = len(cells) - 1
            if repetitions != 0:
                lines = [f"SD{tally_ID} 1 {repetitions}R\n"]
            else:
                lines = [f"SD{tally_ID} 1\n"]
            self.other_data[f"SD{tally_ID}"] = parser.Card(lines, 5, -1)
        # add multiplier if available
        if multiplier is not None:
            line = [f"FM{tally_ID} {multiplier}\n"]
            self.other_data[f"FM{tally_ID}"] = parser.Card(line, 5, -1)

    def add_stopCard(self, nps: int | float = 1e7):
        """
        This does not append a stop card anymore but simply changes the NPS card. If
        a float is provided, it is converted to an integer.

        Parameters
        ----------
        nps : int | float
            number of particles to simulate. Default 1e7

        Returns
        -------
        None.

        """

        line = "NPS " + str(int(nps)) + " \n"

        card = parser.Card([line], 5, -1)
        self.other_data["NPS"] = card

    def check_range(self, range: list[int], who: str = "cell") -> bool:
        """Check if the provided range is not within the used index, i.e., if the range
        is free. Both cells and surfaces can be checked.

        Parameters
        ----------
        range : list[int]
            list of indices to be checked
        who : str, optional
            either 'surf' or 'cell', by default 'cell'

        Returns
        -------
        bool
            True if the provided range is free, False otherwise

        Raises
        ------
        ValueError
            only surf or cell are accepted as who
        """
        if who == "cell":
            used_index = self.cells.keys()
        elif who == "surf":
            used_index = self.surfs.keys()
        else:
            raise ValueError("who can be only 'cell' or 'surf'")
        # check that the provided range is not within the used index
        # need to do the check for strings since there may be asterisks
        for i in range:
            for j in used_index:
                index = j.replace("*", "")
                index = index.replace("+", "")
                if str(i) == index:
                    return False
        return True

    def delete_fill_cards(self) -> None:
        """Delete all fill cards from the input cells."""
        for key, cell in self.cells.items():
            if cell.get_f() is not None:
                # remove the fill card if present
                cell.remove_fill()

                # Be sure to have correct fields generated
                cell_lines = cell.card(wrap=True).splitlines(keepends=True)
                new_cell = parser.Card(cell_lines, 3, cell.pos)
                new_cell.get_values()
                self.cells[key] = new_cell


class D1S_Input(Input):
    def __init__(
        self,
        *args,
        irrad_file: IrradiationFile | None = None,
        reac_file: ReactionFile | None = None,
    ) -> None:
        """Children of the :py:class:`Input`.

        it includes also the reaction and irradiation files necessary for a
        D1S-UNED run and defines additional methods related to them.

        Parameters
        ----------
        args : list
            list of arguments to be passed to the super constructor.
        irrad_file : IrradiationFile, optional
            irradiation file object, by default None
        reac_file : ReactionFile, optional
            readtion file object, by default None

        Attributes
        ----------
        irrad_file : IrradiationFile
            irradiation file object
        reac_file : ReactionFile
            readtion file object

        Examples
        --------
        translate the input defining an activation and transport library.
        the reaction file will be used to identify to which isotopes the
        activation library has to be assigned.

        >>> from f4enix.input.MCNPinput import D1S_Input
        ... from f4enix.input.libmanager import LibManager
        ... d1s_inp = D1S_Input.from_input('d1stest.i', irrad_file='irr_test',
        ...                                reac_file='reac_fe')
        ... d1s_inp.smart_translate('99c', '00c', LibManager())
        """

        super().__init__(*args)
        self.irrad_file = irrad_file
        self.reac_file = reac_file

    @classmethod
    def from_input(
        cls,
        inputfile: os.PathLike,
        irrad_file: os.PathLike | None = None,
        reac_file: os.PathLike | None = None,
    ) -> D1S_Input:
        """Generate a D1S-UNED input file.

        this includes also the reaction and irradiation files.

        Parameters
        ----------
        inputfile : os.PathLike
            path to the MCNP input (D1S)
        irrad_file : os.PathLike, optional
            path to the irradiation file, by default None (no file associated)
        reac_file : os.PathLike, optional
            path to the reaction file, by default None (no file associated)

        Returns
        -------
        D1S_Input
            generated D1S_Input object
        """
        cells, surfaces, data, header = _get_input_arguments(inputfile)
        cells = cls._to_dict(cells)
        surfs = cls._to_dict(surfaces)

        (
            materials,
            transformations,
            other_data,
        ) = cls._parse_data_section(cls, data)

        # get a list of the tally keys
        tally_keys = []
        fmesh_keys = []
        for key, card in other_data.items():
            if card.dtype == "Fn":
                tally_keys.append(card.name)
            elif PAT_FMESH_KEY.match(key):
                fmesh_keys.append(card.name)

        if irrad_file is not None:
            newirrad_file = IrradiationFile.from_text(irrad_file)
        else:
            newirrad_file = None

        if reac_file is not None:
            newreac_file = ReactionFile.from_text(reac_file)
        else:
            newreac_file = None

        return cls(
            cells,
            surfs,
            materials,
            transformations,
            other_data,
            tally_keys,
            fmesh_keys,
            header,
            irrad_file=newirrad_file,
            reac_file=newreac_file,
        )

    def get_potential_paths(self, libmanager: LibManager, lib: str) -> list[Reaction]:
        """Given an activation library, return a list of all possible reactions
        paths foreseen by the libmanager that can originate from the material
        section of the input.

        Parameters
        ----------
        libmanager : LibManager
            Handlers of cross sections operations
        lib : str
            activation library to be used (e.g. 99c)

        Returns
        -------
        list[Reaction]
            list of Reaction objects describing all possible paths included in
            LibManager
        """
        reactions = []
        for material in self.materials.materials:
            for submat in material.submaterials:
                for zaid in submat.zaidList:
                    parent = zaid.element + zaid.isotope
                    zaidreactions = libmanager.get_reactions(lib, parent)
                    # if len(zaidreactions) > 0:
                    #     # it is a parent only if reactions are available
                    #     parentlist.append(parent)
                    for MT, daughter in zaidreactions:
                        reactions.append((parent, MT, daughter))
                        # daughterlist.append(daughter)

        reactions = list(set(reactions))
        reactions.sort()
        # --- Build the reactions and reaction file ---
        reaction_list = []
        for parent, MT, daughter in reactions:
            parent = parent + "." + lib
            # Build a comment
            _, parent_formula = libmanager.get_zaidname(parent)
            if daughter[-3:] == "900":
                # Then it is metastable
                _, daughter_formula = libmanager.get_zaidname(daughter[:-3])
                daughter_formula += "m"
            else:  # ground state
                _, daughter_formula = libmanager.get_zaidname(daughter)
            comment = "{} -> {}".format(parent_formula, daughter_formula)

            rx = Reaction(parent, MT, daughter, comment=comment)
            reaction_list.append(rx)

        return reaction_list

    def get_reaction_file(
        self, libmanager: LibManager, lib: str, set_as_attribute: bool = True
    ) -> ReactionFile:
        """
        Get a reaction file suitable for the input.

        The reaction file is built selecting from all the possible reaction
        paths that can originate in the model due to its material cards only
        the reactions that lead to a daughter listed in the irradiation file.
        By default this is added as the reac_file for the input.

        Parameters
        ----------
        libmanager : LibManager
            Object handling all cross-sections related operations.
        lib : str
            library suffix to be used.
        set_as_attribute: bool
            if True (default) the reactionfile is saved as the self.reac_file

        Returns
        -------
        ReactionFile
            Object representing the react file for D1S.

        Raises
        ------
        ValueError
            if no irradiation files have been assigned yet to the input

        """
        # irrad file is necessary for this operation
        # recover all available daughters
        if self.irrad_file is None:
            raise ValueError("irrad_file attribute cannot be None for this operation")
        else:
            available_daughters = self.irrad_file.get_daughters()

        # Recover all possible reactions
        reactions = self.get_potential_paths(libmanager, lib)

        # perform the selection
        selected_reactions = []
        for reaction in reactions:
            if reaction.daughter in available_daughters:
                # add the reaction to the one to use
                selected_reactions.append(reaction)

        reac_file = ReactionFile(selected_reactions)

        if set_as_attribute:
            self.reac_file = reac_file

        return reac_file

    def smart_translate(
        self,
        activation_lib: str,
        transport_lib: str,
        libmanager: LibManager,
        fix_natural_zaid: bool = False,
    ) -> None:
        """
        Translate the input to another library without relying on old libs.

        Both the activation and transport libraries are changed. The reaction
        file and PIKMT cards are also translated.
        Differently from self.translate(),
        libraries are not modified based on old ones but based on the reaction
        file. This is, to all parent zaids listed in the reactions will be
        assigned the activation_lib, to all others the transport_lib.

        Parameters
        ----------
        activation_lib : str
            library to be used for activation, e.g., 99c
        activation_lib : dict[str, str]
            library to be used for activation, e.g., 31c
        libmanager : LibManager
            Library manager for the conversion.
        fix_natural_zaid: bool
            if True, and additional initial translation with the transport lib
            is done in order to expand the natural zaids. If the transport lib
            do not expand the natural zaid there may be issues.

        Returns
        -------
        None

        """
        if self.reac_file is None:
            raise ValueError("reac_file cannot be None for this operation")

        if fix_natural_zaid:
            # get a first translation to avoid issues with old natural zaids
            self.translate(transport_lib, libmanager)

        active_zaids = []
        transp_zaids = []

        for reaction in self.reac_file.reactions:
            # strip the lib from the parent
            parent = reaction.parent.split(".")[0]
            active_zaids.append(parent)
            reaction.change_lib(activation_lib)

        # Now check for the remaing materials in the input to be assigned
        # to transport
        for material in self.materials.materials:
            for submaterial in material.submaterials:
                for zaid in submaterial.zaidList:
                    zaidnum = zaid.element + zaid.isotope
                    if zaidnum not in active_zaids and zaidnum not in transp_zaids:
                        transp_zaids.append(zaidnum)

        newlib = {activation_lib: active_zaids, transport_lib: transp_zaids}

        # trigger translation of the PIKMT card
        self.add_PIKMT_card()

        # Translate the input with the new lib
        self.materials.translate(newlib, libmanager)

    def add_PIKMT_card(self) -> None:
        """
        Add a PIKMT card to the input file. If a PKMT card is already present
        this will be overridden.

        Returns
        -------
        None.

        """
        if self.reac_file is None:
            logging.warning("No reaction file has been assigned to the input")
            return
        key = "PIKMT"
        lines = [key + "\n"]
        for parent in self.reac_file.get_parents():
            lines.append("         {}    {}\n".format(parent, 0))

        card = parser.Card(lines, 5, -1)
        self.other_data[key] = card  # should override other PKMT cards

    def add_track_contribution(
        self, tallykey: str, zaids: list[str], who: str = "parent"
    ):
        """
        Given a list of zaid add the FU bin in the requested tallies in order
        to collect the contribution of them to the tally.

        Parameters
        ----------
        tallykey : str
            ID of the tally onto which to operate (e.g. F4).
        zaids : list[str]
            list of zaid numbers of the parents/daughters (e.g. 1001).
        who : str, optional
            either 'parent' or 'daughter' specifies the types of zaids to
            be tracked. The default is 'parent'.

        Raises
        ------
        ValueError
            check for admissible who parameter.

        Returns
        -------
        bool
            return True if lines were added correctly

        """
        card = self.other_data[tallykey]
        num = str(_get_num_tally(tallykey))

        card.lines.append("FU" + num + " 0\n")

        if who == "parent":
            for zaid in zaids:
                card.lines.append(ADD_LINE_FORMAT.format("-" + str(zaid)))
        elif who == "daughter":
            for zaid in zaids:
                card.lines.append(ADD_LINE_FORMAT.format(zaid))
        else:
            raise ValueError(who + ' is not an admissible "who" parameters')
        card.get_input()


def _get_input_arguments(inputfile: os.PathLike | str) -> tuple:
    name = os.path.basename(inputfile).split(".")[0]

    # Get the blocks using numjuggler parser
    logging.info("Reading file: {}".format(name))
    jug_cards = parser.get_cards_from_input(inputfile)
    try:
        jug_cardsDic = parser.get_blocks(jug_cards)
    except UnicodeDecodeError as e:
        logging.error("The file contains unicode errors, scan initiated")
        txt = debug_file_unicode(inputfile)
        logging.error("The following error where encountered: \n" + txt)
        raise e

    logging.debug("Reading has finished")

    # Parse the different sections
    header = jug_cardsDic[2][0].lines
    cells = jug_cardsDic[3]
    surfaces = jug_cardsDic[4]
    data = jug_cardsDic[5]

    return cells, surfaces, data, header


def _get_num_tally(key: str) -> int:
    patnum = re.compile(r"\d+")
    try:
        num = patnum.search(key).group()
    except AttributeError:
        # The pattern was not found
        raise ValueError(key + " is not a valid tally ID")

    return int(num)


def _get_card_key(card: parser.Card) -> str:
    # it may happen that comments have been added to the card. The first non
    # comment line is the one to be used
    d = "\n"
    lines = [e + d for e in card.card().split(d) if e]
    for line in lines:
        if PAT_COMMENT.match(line) is None:
            return line.split()[0].upper()

    # if this point is reached the name has not been found
    raise ValueError(f"No key was found for card {card.card()}")
