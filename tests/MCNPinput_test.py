import os
import numpy as np
import pytest
from copy import deepcopy
from importlib.resources import files, as_file
from numjuggler import parser
import f4enix.resources as pkg_res
import tests.resources.input as input_res
import tests.resources.libmanager as lib_res
import pandas as pd

from f4enix.input.MCNPinput import Input, D1S_Input
from f4enix.input.libmanager import LibManager
from f4enix.input.d1suned import ReactionFile, IrradiationFile


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
    with as_file(resources_inp.joinpath("test.i")) as FILE1:
        testInput = Input.from_input(FILE1)

    with as_file(resources_inp.joinpath("various_bugs.i")) as file:
        bugInput = Input.from_input(file)
    # exceptInput = InputFile.from_text(INP_EX_PATH)
    with (
        as_file(resources_lib.joinpath("Activation libs.xlsx")) as ACTIVATION_FILE,
        as_file(resources_lib.joinpath("xsdir")) as XSDIR_FILE,
        as_file(resources_pkg.joinpath("Isotopes.txt")) as ISOTOPES_FILE,
    ):
        lm = LibManager(
            XSDIR_FILE, activationfile=ACTIVATION_FILE, isotopes_file=ISOTOPES_FILE
        )

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
        outfile = tmpdir.mkdir("sub").join("tempfile.i")
        inp.write(outfile)
        # re-read
        inp2 = Input.from_input(outfile)
        self._check_macro_properties(inp2)

        # test if translations are rewritten correctly
        inp = deepcopy(self.bugInput)
        outfile = tmpdir.mkdir("sub2").join("tempfile2.i")
        inp.write(outfile)
        inp2 = Input.from_input(outfile)
        _ = inp2.transformations["TR1"]

        assert True

    def _check_macro_properties(self, inp: Input):
        # check some macro properties
        assert inp.header[0].strip("\n") == "This is the header"
        assert len(inp.cells) == 128
        assert len(inp.surfs) == 129
        assert len(inp.materials) == 25
        assert len(inp.tally_keys) == 7
        assert len(inp.fmesh_keys) == 5

    def test_update_zaidinfo(self):
        newinput = deepcopy(self.testInput)
        newinput.update_zaidinfo(self.lm)
        assert True

    def test_translate(self):
        # The test for a correct translation of material card is already done
        # in materials. here we only check that it goes trough without errors
        newinput = deepcopy(self.testInput)
        newinput.translate("00c", self.lm)
        newinput = deepcopy(self.testInput)
        newinput.translate('{"31c": "00c", "70c": "81c"}', self.lm)
        assert True

        # let's check also that abundances info is correctly added
        assert "$ H-1    AB(%) 99.988" in newinput.materials.materials[0].to_text()

    def test_get_cells_by_id(self):
        cards = self.testInput.get_cells_by_id([1, 2])
        cards = self.testInput.get_cells_by_id(["1", "2"])
        assert True

    def test_get_surfs_by_id(self):
        cards = self.testInput.get_surfs_by_id([1, 2])
        cards = self.testInput.get_surfs_by_id(["1", "2"])
        assert True

    def test_get_materials_subset(self):
        materials = "m23"
        _ = self.testInput.get_materials_subset(materials)
        materials = ["m22", "M30"]
        _ = self.testInput.get_materials_subset(materials)
        assert True

    def test_get_data_cards(self):
        _ = self.testInput.get_data_cards("SDEF")

        try:
            self.testInput.get_data_cards("adas")
            assert False
        except KeyError:
            assert True

    def test_get_cells_summary(self):
        df = self.testInput.get_cells_summary()
        assert len(df) == 128
        assert len(df.columns) == 4

    # def test_print_cards(self):
    #     newinput = deepcopy(self.testInput)
    #     print(newinput._print_cards(newinput.cells))
    #     assert False

    def test_extract_cells(self, tmpdir):
        newinput = deepcopy(self.testInput)
        cells = [23, 24, 25, 31]
        outfile = tmpdir.mkdir("sub").join("extract.i")
        newinput.extract_cells(cells, outfile, renumber_from=1)
        # re-read
        inp2 = Input.from_input(outfile)
        assert len(inp2.cells) == 5
        assert len(inp2.surfs) == 10
        assert len(inp2.materials) == 3
        assert list(inp2.cells.keys()) == ["1", "2", "3", "4", "5"]
        assert inp2.cells["3"].values[-2][0] == 1

        with as_file(resources_inp.joinpath("test_1.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        outfile = os.path.join(os.path.dirname(outfile), "extract_fillers.i")

        mcnp_input.extract_cells([50], outfile, extract_fillers=True, renumber_from=500)

        # re-read
        result = Input.from_input(outfile)

        assert len(result.cells) == 7
        assert mcnp_input.cells["10"].values[0][0] == 10

    def test_extract_universe(self, tmpdir):
        with as_file(resources_inp.joinpath("test_universe.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        outfile = tmpdir.mkdir("sub").join("extract_universe.i")

        universe = 125
        mcnp_input.extract_universe(universe, outfile)

        # re-read
        result = Input.from_input(outfile)

        assert len(result.cells) == 3
        assert len(result.surfs) == 2
        assert len(result.materials) == 1
        for _, cell in result.cells.items():
            assert cell.get_u() is None
        assert mcnp_input.cells["21"].get_u() == universe

    def test_duplicated_nums(self):
        # There was a bug reading material 101
        self.bugInput.get_materials_subset(["m101"])
        assert True

    def test_missing_data_cards(self):
        _ = self.bugInput.other_data["SP2"]
        _ = self.bugInput.transformations["TR1"]
        _ = self.bugInput.other_data["CUT:N"]

        assert True

    def test_clean_card_name(self):
        cardnames = ["*TR1", "f6:n,p"]
        expected = ["TR1", "f6"]

        for name, exp in zip(cardnames, expected):
            assert Input._clean_card_name(name) == exp

    @pytest.mark.parametrize("flag", [True, False])
    def test_get_cells_by_matID(self, flag):
        newinput = deepcopy(self.testInput)
        cells = newinput.get_cells_by_matID(13, deepcopy_flag=flag)
        for filtered, expected in zip(cells.keys(), range(2, 22)):
            assert filtered == str(expected)

    def test_scale_densities(self):
        newinput = deepcopy(self.testInput)
        d1 = newinput.cells["49"].get_d()
        d2 = newinput.cells["52"].get_d()
        d2 = newinput.cells["53"].get_d()
        newinput.scale_densities(0.33333333333)
        assert newinput.cells["49"].get_d() == 0.0412067
        assert newinput.cells["52"].get_d() == 0
        assert newinput.cells["53"].get_d() == -2.60000

    @pytest.mark.parametrize(
        ["id", "expected"],
        [[94, ["FC94", "F94", "FM94"]], [214, ["FC214", "FMESH214", "FM214"]]],
    )
    def test_get_tally_cards(self, id, expected):
        keys = self.testInput._get_tally_cards(id)
        assert keys == expected

    def test_get_tally_summary(self):
        summary = self.testInput.get_tally_summary()
        assert len(summary) == 7
        assert summary.loc[194].values.tolist() == [
            "N",
            "T in Li pt2 appm/FPY",
            "3.8566e10",
            ["25", "205"],
        ]
        assert summary.loc[204].values.tolist() == ["N", pd.NA, pd.NA, pd.NA]

        summary = self.testInput.get_tally_summary(fmesh=True)
        assert len(summary) == 5
        assert summary.loc[224].values.tolist() == [
            "P",
            "FMESH Photon Heating [MeV/cc/n_s]",
            "-1",
            ["0", "-5", "-6"],
        ]

    def test_set_cell_void(self):
        newinput = deepcopy(self.testInput)
        Input.set_cell_void(newinput.cells["49"])
        assert (
            newinput.cells["49"].card()
            == "49   0     -128 129 48  -49               $imp:n,p=1\n"
        )

    def test_replace_material(self):
        with as_file(resources_inp.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.replace_material(10, "-2", 4)
        assert newinp.cells["21"].get_m() == 10
        assert newinp.cells["21"].get_d() == -2

        newinp.replace_material(0, "10", 10, u_list=[125])
        assert newinp.cells["21"].get_m() == 0

    def test_add_surface(self):
        newinput = deepcopy(self.testInput)
        sur = 180
        new_cell = Input.add_surface(
            newinput.cells["27"], -sur, None, "intersect", False
        )
        assert not any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in newinput.cells["27"].values
        )
        assert not newinput.cells["27"].input[0].split()[-1] == r"-{:<3}"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[-1] == r"-{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[3] == "(-128"
        assert new_cell.card(wrap=False, comment=False).split()[-1] == str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[-2][-1] == ")"

        new_cell = Input.add_surface(
            newinput.cells["27"], -sur, None, "intersect", True
        )
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in newinput.cells["27"].values
        )
        assert newinput.cells["27"].input[0].split()[-1] == r"-{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[3] == "(-128"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[-1] == r"-{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[-1] == str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[-2][-1] == ")"
        assert newinput.cells["27"].card(wrap=False, comment=False).split()[-1] == str(
            -sur
        )

        with as_file(resources_inp.joinpath("test_universe.i")) as FILE:
            mcnp_input = Input.from_input(FILE)
        new_cell = Input.add_surface(mcnp_input.cells["1"], sur, None, "union", False)
        assert new_cell.values[0][0] == 1
        assert not any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in mcnp_input.cells["1"].values
        )
        assert not mcnp_input.cells["1"].input[0].split()[-3] == r":{:<3}"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[-3] == r":{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[-3] == ":" + str(sur)
        assert new_cell.card(wrap=False, comment=False).split()[-4] == ")"
        assert (
            mcnp_input.cells["1"].card(wrap=False, comment=False).split()[-2]
            == "imp:n=1"
        )
        assert new_cell.card(wrap=False, comment=False).split()[2] == "(-1"

        new_cell = Input.add_surface(mcnp_input.cells["1"], -sur, 50, "union", False)
        assert new_cell.values[0][0] == 50
        assert mcnp_input.cells["1"].values[0][0] == 1
        assert not any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in mcnp_input.cells["1"].values
        )
        assert not mcnp_input.cells["1"].input[0].split()[-3] == r":-{:<3}"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[-3] == r":-{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[-3] == ":" + str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[-4] == ")"
        assert (
            mcnp_input.cells["1"].card(wrap=False, comment=False).split()[-2]
            == "imp:n=1"
        )
        assert new_cell.card(wrap=False, comment=False).split()[2] == "(-1"
        assert new_cell.card(wrap=False, comment=False).split()[0] == "50"
        assert mcnp_input.cells["1"].card(wrap=False, comment=False).split()[0] == "1"

        sur = 5555
        new_cell = Input.add_surface(mcnp_input.cells["22"], -sur, 50, "union", True)
        assert new_cell.values[0][0] == 50
        assert mcnp_input.cells["22"].values[0][0] == 50
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in mcnp_input.cells["22"].values
        )
        assert mcnp_input.cells["22"].input[0].split()[4] == r":-{:<4}"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[4] == r":-{:<4}"
        assert new_cell.card(wrap=False, comment=False).split()[4] == ":" + str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[3] == ")"
        assert (
            mcnp_input.cells["22"].card(wrap=False, comment=False).split()[-1]
            == "U=125"
        )
        assert new_cell.card(wrap=False, comment=False).split()[2] == "(-22"
        assert new_cell.card(wrap=False, comment=False).split()[0] == "50"
        assert mcnp_input.cells["22"].card(wrap=False, comment=False).split()[0] == "50"


class TestD1S_Input:
    with (
        as_file(resources_inp.joinpath("d1stest.i")) as inp_file,
        as_file(resources_inp.joinpath("d1stest_irrad")) as irrad_file,
        as_file(resources_inp.joinpath("d1stest_react")) as react_file,
    ):
        inp = D1S_Input.from_input(
            inp_file, reac_file=react_file, irrad_file=irrad_file
        )

    with (
        as_file(resources_lib.joinpath("Activation libs.xlsx")) as ACTIVATION_FILE,
        as_file(resources_lib.joinpath("xsdir")) as XSDIR_FILE,
        as_file(resources_pkg.joinpath("Isotopes.txt")) as ISOTOPES_FILE,
    ):
        lm = LibManager(
            XSDIR_FILE, activationfile=ACTIVATION_FILE, isotopes_file=ISOTOPES_FILE
        )

    def test_smart_translate(self):
        # This test needs to be improved
        with (
            as_file(resources_inp.joinpath("d1stest_irrad_st")) as irrad_file,
            as_file(resources_inp.joinpath("d1stest_react_st")) as react_file,
        ):
            react_file = ReactionFile.from_text(react_file)
            irrad_file = IrradiationFile.from_text(irrad_file)

        newinp = deepcopy(self.inp)
        newinp.irrad_file = irrad_file
        newinp.reac_file = react_file

        activation_lib = "98c"
        transport_lib = "00c"
        newinp.smart_translate(
            activation_lib, transport_lib, self.lm, fix_natural_zaid=True
        )

        translation = newinp.materials.to_text()

        assert translation.count("98c") == 4
        assert translation.count("00c") == 145
        assert newinp.reac_file.reactions[0].parent == "24050.98c"

    def test_add_PKMT_card(self):
        with as_file(resources_inp.joinpath("d1stest_noPKMT.i")) as inp_file:
            newinp = D1S_Input.from_input(inp_file)
        newinp.reac_file = self.inp.reac_file

        newinp.add_PIKMT_card()
        card = newinp.other_data["PIKMT"]
        assert len(card.lines) == 17

    def test_get_reaction_file(self):
        with (
            as_file(resources_inp.joinpath("d1stest_getreact.i")) as inp_file,
            as_file(resources_inp.joinpath("d1stest_irrad_getreact")) as irr_file,
        ):
            newinp = D1S_Input.from_input(inp_file, irrad_file=irr_file)

        lib = "99c"
        reacfile = newinp.get_reaction_file(self.lm, lib)
        assert ["24050"] == reacfile.get_parents()

    def test_get_potential_paths(self):
        reaction_list = self.inp.get_potential_paths(self.lm, "98c")
        assert len(reaction_list) == 32

    @pytest.mark.parametrize(["who", "sign"], [["parent", "-"], ["daughter", ""]])
    def test_add_track_contribution(self, tmpdir, who, sign):
        zaids = ["1001", "1002"]
        tallyID = "F124"

        # --- Test parents---
        inp = deepcopy(self.inp)
        inp.add_track_contribution(tallyID, zaids, who=who)
        # dump and reread the input
        tmpfile = os.path.join(tmpdir, "tmp.i")
        inp.write(tmpfile)
        newinp = D1S_Input.from_input(tmpfile)
        # get the new injected card
        card = newinp.other_data["FU124"]
        for line, exp in zip(
            card.lines[-3:], ["FU124 0", sign + "1001", sign + "1002"]
        ):
            assert line.strip() == exp
