import os
import logging
import json
import re
from numjuggler import parser
import shutil

from f4eparser.input.materials import MatCardsList, Material
from f4eparser.input.libmanager import LibManager
from f4eparser.input.auxiliary import debug_file_unicode
from copy import deepcopy


PAT_MT = re.compile(r'm[tx]\d+', re.IGNORECASE)


class Input:
    def __init__(
            self,
            cells: list[parser.Card],
            surfs: list[parser.Card],
            data: list[parser.Card],
            header: list = None) -> None:

        self.cells = self._to_dict(cells)
        self.surfs = self._to_dict(surfs)

        (self.materials, self.transformations,
         self.other_data) = self._parse_data_section(data)

        self.header = header

        # store also all the original cards for compatibility
        # with some numjuggler modes
        cells.extend(surfs)
        cells.extend(data)

    @classmethod
    def from_input(cls, inputfile: os.PathLike):
        """ Generate an Input object from an MCNP input text file using
        numjuggler as a parser

        Parameters
        ----------
        inputfile : os.PathLike
            input file

        Returns
        -------
        Input
            Input object
        """
        name = os.path.basename(inputfile).split(".")[0]

        # Get the blocks using numjuggler parser
        logging.info('Reading file: {}'.format(name))
        jug_cards = parser.get_cards_from_input(inputfile)
        try:
            jug_cardsDic = parser.get_blocks(jug_cards)
        except UnicodeDecodeError as e:
            logging.error('The file contains unicode errors, scan initiated')
            txt = debug_file_unicode(inputfile)
            logging.error('The following error where encountered: \n'+txt)
            raise e

        logging.debug('Reading has finished')

        # Parse the different sections
        header = jug_cardsDic[2][0].lines
        cells = jug_cardsDic[3]
        surfaces = jug_cardsDic[4]
        data = jug_cardsDic[5]

        return cls(cells, surfaces, data, header=header)

    def write(self, outfilepath: os.PathLike) -> None:
        """write the input to a file

        Parameters
        ----------
        outfilepath : os.PathLike
            path to the output file
        """
        logging.info('Writing to {}'.format(outfilepath))

        with open(outfilepath, 'w') as outfile:

            # Add the header lines
            for line in self.header:
                outfile.write(line)
            # Add the cells
            self._write_cards(self.cells, outfile)
            # Add a break
            outfile.write('\n')
            # Add the surfaces
            self._write_cards(self.surfs, outfile)
            # Add a break
            outfile.write('\n')
            # Add the material section
            outfile.write(self.materials.to_text())
            # Add the rest of the data cards
            self._write_cards(self.transformations, outfile)
            # Add a break
            outfile.write('\n')
            self._write_cards(self.other_data, outfile)
            # Add a break
            outfile.write('\n')

        logging.info('File was written correctly')

    def translate(self, newlib: str, libmanager: LibManager) -> None:
        """
        Translate the input to another library

        Parameters
        ----------
        newlib : str
            suffix of the new lib to translate to.
        libmanager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        None.

        """

        try:
            if newlib[0] == '{':
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
    def _write_cards(cards: dict[str, parser.Card], outfile) -> None:
        for _, card in cards.items():
            for line in card.card(wrap=True):
                outfile.write(line)

    @staticmethod
    def _to_dict(cards: list[parser.Card]) -> dict[str, parser.Card]:
        new_cards = {}
        flag_add = False
        for card in cards:
            card.get_values()
            try:
                key = card.name
            except AttributeError:
                # This means that this is a fake card just made by comments
                # it should be merged with the following card
                comment = card.lines
                flag_add = True
                continue

            if flag_add:
                comment.extend(card.lines)
                card.lines = comment
                card.get_input()
                card.get_values()
                flag_add = False

            new_cards[key] = card

        return new_cards

    @staticmethod
    def _get_cards_by_id(ids: list, cards: dict) -> dict:
        selected_cards = {}
        for id_card in ids:
            selected_cards[id_card] = cards[id_card]

        return selected_cards

    def get_cells_by_id(self, ids: list[int]) -> dict:
        """given a list of cells id return a dictionary of such cells

        Parameters
        ----------
        ids : list[int]
            cells id to be extracted

        Returns
        -------
        dict
            extracted cells
        """
        return self._get_cards_by_id(ids, self.cells)

    def get_surfs_by_id(self, ids: list[int]) -> dict:
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
        return self._get_cards_by_id(ids, self.surfs)

    def get_materials_subset(self, ids: list[str] | str) -> MatCardsList:
        """given a list of material ids generate a new MatCardsList with
        the requested subset

        Parameters
        ----------
        ids : Union(list[str]), str)
            ids of the materials to put into the subset

        Returns
        -------
        MatCardsList
            new materials subset
        """
        if type(ids) is str:
            return self.materials[ids.upper()]
        else:
            materials = []
            for id_mat in ids:
                materials.append(self.materials[id_mat.upper()])
            return MatCardsList(materials)

    def _parse_data_section(self, cards: list[parser.Card]
                            ) -> tuple[MatCardsList,
                                       list[parser.Card],
                                       list[parser.Card]]:

        # first of all correct numjuggler parser
        cards = self._to_dict(cards)

        materials = []  # store here the materials
        transformations = {}  # store translations
        other_data = {}  # store here the other data cards

        for key, card in cards.items():
            try:
                if card.values[0][1] == 'mat':
                    if (PAT_MT.match(card.lines[0]) or
                            PAT_MT.match(card.lines[-1])):
                        # mt or mx cards should be added to the previous
                        # material
                        materials[-1].add_mx(card)
                    else:
                        materials.append(Material.from_text(card.lines))

                elif card.dtype == 'TRn':
                    transformations[key] = card

                else:
                    other_data[key] = card

            except IndexError:
                # this means that there were no values
                other_data[key] = card

        return MatCardsList(materials), transformations, other_data

    def extract_cells(self, cells: list[int], outfile: os.PathLike):
        """given a list of cells, dumps a minimum MCNP working file that
        includes all the requested cells, defined surfaces, materials and
        translations.

        Parameters
        ----------
        cells : list[int]
            desired list of cells
        outfile : os.PathLike
            path to the file where the MCNP input needs to be dumped
        """
        logging.info('Collecting the cells, surfaces, materials and transf.')
        cset = set(cells)

        # first, get all surfaces needed to represent the cn cell.
        sset = set()  # surfaces
        mset = set()  # material
        # tset = set()  # transformations

        # next runs: find all other cells:
        again = True
        while again:
            again = False
            for key, c in self.cells.items():
                if key in cset:
                    cref = c.get_refcells()
                    if cref.difference(cset):
                        again = True
                        cset = cset.union(cref)

        # Get all surfaces and materials
        cells = self.get_cells_by_id(cset)
        for _, cell in cells.items():
            for v, t in cell.values:
                if t == 'sur':
                    sset.add(v)
                elif t == 'mat':
                    mset.add('M'+str(v))

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
        logging.info('write MCNP reduced input')
        with open(outfile, 'w') as outfile:
            # Add the header lines
            for line in self.header:
                outfile.write(line)
            # Add the cells
            self._write_cards(cells, outfile)
            # Add a break
            outfile.write('\n')
            # Add the surfaces
            surfs = self.get_surfs_by_id(sset)
            self._write_cards(surfs, outfile)
            # Add a break
            outfile.write('\n')
            # Add materials
            materials = self.get_materials_subset(mset)
            outfile.write(materials.to_text())
            outfile.write('\n')
            self._write_cards(self.transformations, outfile)

        logging.info('input written correctly')
