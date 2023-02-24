import os
import numpy as np
from copy import deepcopy
from importlib.resources import files, as_file
import mcparser.resources as pkg_res
import tests.resources.input as input_res
import tests.resources.libmanager as lib_res

from mcparser.input.inputAPI import Input
from mcparser.input.libmanager import LibManager


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
    # exceptInput = InputFile.from_text(INP_EX_PATH)
    with (as_file(resources_lib.joinpath('Activation libs.xlsx')) as ACTIVATION_FILE,
          as_file(resources_lib.joinpath('xsdir')) as XSDIR_FILE,
          as_file(resources_pkg.joinpath('Isotopes.txt')) as ISOTOPES_FILE):

        lm = LibManager(XSDIR_FILE, activationfile=ACTIVATION_FILE,
                        isotopes_file=ISOTOPES_FILE)

    def test_from_input(self):
        inp = deepcopy(self.testInput)
        self._check_macro_properties(inp)

    def test_write(self, tmpdir):
        # read
        inp = deepcopy(self.testInput)
        # write
        outfile = tmpdir.mkdir('sub').join('tempfile.i')
        inp.write(outfile)
        # re-read
        inp2 = Input.from_input(outfile)
        self._check_macro_properties(inp2)

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
