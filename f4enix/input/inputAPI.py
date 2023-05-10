"""
This module is related to the parsing and manipulation of MCNP input files.

The parser is built on the numjuggler python module.

"""
import os
import logging
import json
import re
from numjuggler import parser
import shutil

from f4enix.input.materials import MatCardsList, Material
from f4enix.input.libmanager import LibManager
from f4enix.input.auxiliary import debug_file_unicode
from f4enix.constants import PAT_COMMENT, PAT_CARD_KEY
from copy import deepcopy


PAT_MT = re.compile(r'm[tx]\d+', re.IGNORECASE)


class Input:
    def __init__(
            self,
            cells: list[parser.Card],
            surfs: list[parser.Card],
            data: list[parser.Card],
            header: list = None) -> None:
        """Class representing an MCNP input file.

        cells, surfaces, materials and transformations are handled explicitly.
        All other datacards are treated generically for the moment being.
        The parsing is built on the numjuggler python module.

        Parameters
        ----------
        cells : list[parser.Card]
            list of numjuggler.parser.Card objects containing MCNP cells.
        surfs : list[parser.Card]
            list of numjuggler.parser.Card objects containing MCNP surfaces.
        data : list[parser.Card]
            list of numjuggler.parser.Card objects containing MCNP data cards.
        header : list, optional
            list of strings that compose the MCNP header, by default None

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
        transformations: list[parser.Card]
            list of transformation datacards (i.e. TRn)
        other_data: list[parser.Card]
            list of all remaining datacards that are treated in a generic way

        Examples
        --------
        The most common way to initiliaze an Input object is from a MCNP input
        file

        >>> from f4enix.input.inputAPI import Input
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

        >>> from f4enix.input.inputAPI import Input
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

        self.cells = self._to_dict(cells)
        self.surfs = self._to_dict(surfs)

        (self.materials, self.transformations,
         self.other_data) = self._parse_data_section(data)

        self.header = header

        # # store also all the original cards for compatibility
        # # with some numjuggler modes
        # cells.extend(surfs)
        # cells.extend(data)

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

    def write(self, outfilepath: os.PathLike, wrap: bool = False) -> None:
        """write the input to a file

        Parameters
        ----------
        outfilepath : os.PathLike
            path to the output file
        wrap : bool
            if true the text is wrapped at 80 char. May cause slowdowns
        """
        logging.info('Writing to {}'.format(outfilepath))

        with open(outfilepath, 'w') as outfile:

            # Add the header lines
            for line in self.header:
                outfile.write(line)
            # Add the cells
            outfile.writelines(self._print_cards(self.cells, wrap=wrap))
            # Add a break
            outfile.write('\n')
            # Add the surfaces
            outfile.writelines(self._print_cards(self.surfs, wrap=wrap))
            # Add a break
            outfile.write('\n')
            # Add the material section (they exit without the \n)
            outfile.write(self.materials.to_text()+'\n')
            # Add the translations
            outfile.writelines(self._print_cards(self.transformations,
                                                 wrap=wrap))
            # Add the rest of the datacards
            outfile.writelines(self._print_cards(self.other_data, wrap=wrap))
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
    def _print_cards(cards: dict[str, parser.Card],
                     wrap: bool = False) -> list[str]:
        text = []
        for _, card in cards.items():
            text.append(card.card(wrap=wrap).strip('\n')+'\n')
        return text

    @staticmethod
    def _to_dict(cards: list[parser.Card]) -> dict[str, parser.Card]:
        new_cards = {}
        flag_add = False
        for card in cards:

            card.get_values()
            try:
                key = card.name
                key = card.card().split()[0].upper()
                if key in new_cards.keys():
                    raise KeyError('Duplicated card entry: '+key)

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
                        raise KeyError('Duplicated card entry: '+key)

            if flag_add:
                comment.extend(card.lines)
                card.lines = comment
                card.get_input()
                card.get_values()
                flag_add = False

            new_cards[key] = card

        return new_cards

    @staticmethod
    def _get_cards_by_id(ids: list[str], cards: dict) -> dict:
        selected_cards = {}
        for id_card in ids:
            try:
                selected_cards[id_card] = cards[id_card]
            except KeyError:
                try:
                    # sometimes it may be with an asterisk
                    selected_cards['*'+id_card] = cards['*'+id_card]
                except KeyError:
                    raise KeyError('The card is not available in the input')

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
        str_ids = []
        for id in ids:
            str_ids.append(str(id))
        return self._get_cards_by_id(str_ids, self.cells)

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
        str_ids = []
        for id in ids:
            str_ids.append(str(id))
        return self._get_cards_by_id(str_ids, self.surfs)

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
            key = self._clean_card_name(key)
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
        # make sure these are str
        cset = []
        for cell in cells:
            cset.append(str(cell))
        cset = set(cset)

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
                    if int(v) != 0:  # void material is not defined in a card
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
            outfile.writelines(self._print_cards(cells))
            # Add a break
            outfile.write('\n')
            # Add the surfaces
            surfs = self.get_surfs_by_id(sset)
            outfile.writelines(self._print_cards(surfs))
            # Add a break
            outfile.write('\n')
            # Add materials
            materials = self.get_materials_subset(mset)
            outfile.write(materials.to_text()+'\n')
            outfile.writelines(self._print_cards(self.transformations))

        logging.info('input written correctly')

    @staticmethod
    def _clean_card_name(key: str) -> str:
        # this is to clean cases like:
        # *TR1 -> TR1
        # F6:N,P -> F6
        try:
            newkey = PAT_CARD_KEY.search(key).group()
        except AttributeError:
            logging.debug('the following key was not cleaned: '+key)
            newkey = key

        return newkey

    def get_cells_by_matID(self, matID: int,
                           deepcopy_flag: bool = True
                           ) -> dict[str, parser.Card]:
        """Given a material ID return a dictionary {key, card} of all
        the cells to which that material is assigned to.

        The cells that are returned are deepcopies of the original ones.

        Parameters
        ----------
        matID : int
            material ID to filter the cells
        deepcopy_flag: bool
            if False, the cells are not copied. Default is True

        Returns
        -------
        dict[int, parser.Card]
            cells to which the material is assigned to
        """
        logging.debug('get cells for material {} requested'.format(matID))
        filtered_cells = {}
        for key, cell in self.cells.items():
            if cell._get_value_by_type('mat') == int(matID):
                if deepcopy_flag:
                    filtered_cells[key] = deepcopy(cell)
                else:
                    filtered_cells[key] = cell
        return filtered_cells
