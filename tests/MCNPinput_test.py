import os
import numpy as np
import pytest
from copy import deepcopy
from importlib.resources import files, as_file
from numjuggler import parser
import f4enix.resources as pkg_res
import tests.resources.input as input_res
import tests.resources.libmanager as lib_res

from f4enix.input.MCNPinput import Input
from f4enix.input.libmanager import LibManager


resources_inp = files(input_res)
resources_lib = files(lib_res)
resources_pkg = files(pkg_res)

# INP_EX_PATH = os.path.join(resources_inp, 'test_exceptions.i')
# DIS_INP_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest.i')
# DIS_NOPKMT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_noPKMT.i')
# DIS_GETREACT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_getreact.i')

# IRRAD_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_irrad')
# REACT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_react')


class TestInput:
    with as_file(resources_inp.joinpath('test.i')) as FILE1:
        testInput = Input.from_input(FILE1)

    with as_file(resources_inp.joinpath('various_bugs.i')) as file:
        bugInput = Input.from_input(file)
    # exceptInput = InputFile.from_text(INP_EX_PATH)
    with (as_file(resources_lib.joinpath('Activation libs.xlsx')) as ACTIVATION_FILE,
          as_file(resources_lib.joinpath('xsdir')) as XSDIR_FILE,
          as_file(resources_pkg.joinpath('Isotopes.txt')) as ISOTOPES_FILE):

        lm = LibManager(XSDIR_FILE, activationfile=ACTIVATION_FILE,
                        isotopes_file=ISOTOPES_FILE)

    def test_from_input(self):
        inp = deepcopy(self.testInput)
        self._check_macro_properties(inp)

    # def test_jt60_bug(self, tmpdir):
    #     with as_file(resources_inp.joinpath('jt60.i')) as file:
    #         # MT AND MX CARDS
    #         jt60_input = Input.from_input(file)
    #     # check that writing and re-reading does not change anything
    #     outfile = tmpdir.mkdir('sub').join('jt_tmp.i')
    #     jt60_input.write(outfile)

    #     inp1 = Input.from_input(outfile)
    #     inp1.write(outfile)
    #     print(inp1.materials.matdic)
    #     print(inp1.get_materials_subset('m14').to_text())
    #     inp2 = Input.from_input(outfile)
    #     print(inp2.get_materials_subset('m14').to_text())
    #     outfile2 = tmpdir.mkdir('sub2').join('jt_tmp2.i')
    #     inp2.write(outfile2)

    #     with open(outfile, 'r') as infile1, open(outfile2, 'r') as infile2:
    #         for line1, line2 in zip(infile1, infile2):
    #             assert line1 == line2

    def test_write(self, tmpdir):
        # read
        inp = deepcopy(self.testInput)
        # write
        outfile = tmpdir.mkdir('sub').join('tempfile.i')
        inp.write(outfile)
        # re-read
        inp2 = Input.from_input(outfile)
        self._check_macro_properties(inp2)

        # test if translations are rewritten correctly
        inp = deepcopy(self.bugInput)
        outfile = tmpdir.mkdir('sub2').join('tempfile2.i')
        inp.write(outfile)
        inp2 = Input.from_input(outfile)
        _ = inp2.transformations['TR1']

        assert True

    def _check_macro_properties(self, inp: Input):
        # check some macro properties
        assert inp.header[0].strip('\n') == 'This is the header'
        assert len(inp.cells) == 128
        assert len(inp.surfs) == 129
        assert len(inp.materials) == 25

    def test_update_zaidinfo(self):
        newinput = deepcopy(self.testInput)
        newinput.update_zaidinfo(self.lm)
        assert True

    def test_translate(self):
        # The test for a correct translation of material card is already done
        # in materials. here we only check that it goes trough without errors
        newinput = deepcopy(self.testInput)
        newinput.translate('00c', self.lm)
        newinput = deepcopy(self.testInput)
        newinput.translate('{"31c": "00c", "70c": "81c"}', self.lm)
        assert True

    def test_get_cells_by_id(self):
        cards = self.testInput.get_cells_by_id([1, 2])
        cards = self.testInput.get_cells_by_id(['1', '2'])
        assert True

    def test_get_surfs_by_id(self):
        cards = self.testInput.get_surfs_by_id([1, 2])
        cards = self.testInput.get_surfs_by_id(['1', '2'])
        assert True

    def test_get_materials_subset(self):
        materials = 'm23'
        _ = self.testInput.get_materials_subset(materials)
        materials = ['m22', 'M30']
        _ = self.testInput.get_materials_subset(materials)
        assert True

    # def test_print_cards(self):
    #     newinput = deepcopy(self.testInput)
    #     print(newinput._print_cards(newinput.cells))
    #     assert False

    def test_extract_cells(self, tmpdir):
        newinput = deepcopy(self.testInput)
        cells = [23, 24, 25, 31]
        outfile = tmpdir.mkdir('sub').join('extract.i')
        newinput.extract_cells(cells, outfile)
        # re-read
        inp2 = Input.from_input(outfile)
        assert len(inp2.cells) == 5
        assert len(inp2.surfs) == 10
        assert len(inp2.materials) == 3

    def test_duplicated_nums(self):
        # There was a bug reading material 101
        self.bugInput.get_materials_subset(['m101'])
        assert True

    def test_missing_data_cards(self):
        _ = self.bugInput.other_data['SP2']
        _ = self.bugInput.transformations['TR1']
        _ = self.bugInput.other_data['CUT:N']

        assert True

    def test_clean_card_name(self):
        cardnames = ['*TR1', 'f6:n,p']
        expected = ['TR1', 'f6']

        for name, exp in zip(cardnames, expected):
            assert Input._clean_card_name(name) == exp

    @pytest.mark.parametrize('flag', [True, False])
    def test_get_cells_by_matID(self, flag):
        newinput = deepcopy(self.testInput)
        cells = newinput.get_cells_by_matID(13, deepcopy_flag=flag)
        for filtered, expected in zip(cells.keys(), range(2, 22)):
            assert filtered == str(expected)

    def test_scale_densities(self):
        newinput = deepcopy(self.testInput)
        d1 = newinput.cells['49'].get_d()
        d2 = newinput.cells['52'].get_d()
        d2 = newinput.cells['53'].get_d()
        newinput.scale_densities(0.33333333333)
        assert newinput.cells['49'].get_d() == 0.0412067
        assert newinput.cells['52'].get_d() == 0
        assert newinput.cells['53'].get_d() == -2.60000
