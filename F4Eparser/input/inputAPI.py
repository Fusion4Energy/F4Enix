import os
import logging
import json
from numjuggler import parser

from f4eparser.input.materials import MatCardsList, Material
from f4eparser.input.libmanager import LibManager


class Input:
    def __init__(
            self,
            cells: list[parser.Card],
            surfs: list[parser.Card],
            data: list[parser.Card],
            header: list = None) -> None:

        self.cells = self._to_dict(cells)
        self.surfs = self._to_dict(surfs)
        self.materials, self.other_data = self._parse_data_section(data)
        self.header = header

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
        jug_cardsDic = parser.get_blocks(jug_cards)
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
            self._write_cards(self.other_data, outfile)

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
            for line in card.card():
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

    def _parse_data_section(self, cards: list[parser.Card]
                            ) -> tuple[MatCardsList, list[parser.Card]]:

        # first of all correct numjuggler parser
        cards = self._to_dict(cards)

        materials = []  # store here the materials
        other_data = {}  # store here the other data cards

        for key, card in cards.items():
            try:
                if card.values[0][1] == 'mat':
                    materials.append(Material.from_text(card.lines))
                else:
                    other_data[key] = card

            except IndexError:
                # this means that there were no values
                other_data[key] = card

        return MatCardsList(materials), other_data
