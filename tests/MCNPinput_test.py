import os
from copy import deepcopy
from importlib.resources import as_file, files

import numpy as np
import pandas as pd
import pytest

import f4enix.resources as pkg_res
import tests.resources.input as input_res
import tests.resources.libmanager as lib_res
from f4enix.input.d1suned import IrradiationFile, ReactionFile
from f4enix.input.libmanager import LibManager
from f4enix.input.MCNPinput import D1S_Input, Input

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

    def test_check_range(self):
        inp = self.testInput
        assert not inp.check_range([1, 2])
        assert inp.check_range(range(1000, 10010))
        assert not inp.check_range([1, 1e4])

        assert not inp.check_range([89], who="surf")
        assert inp.check_range([360], who="surf")
        assert not inp.check_range([1360], who="surf")

    def test_cell_property(self):
        inp = deepcopy(self.testInput)
        # verify that the name and values have been changed
        inp.cells = {"1": inp.cells["1"], "2": inp.cells["1"]}
        assert inp.cells["2"].name == 2
        assert inp.cells["2"].values[0] == (2, "cel")

        # adding a new value should still force the keys and names to be
        # the same
        inp.cells["3"] = deepcopy(inp.cells["1"])
        assert inp.cells["3"].name == 3
        assert inp.cells["3"].values[0] == (3, "cel")

        # no random stuff can be added to the dictionary
        with pytest.raises(ValueError):
            inp.cells[1] = inp.cells["1"]
        with pytest.raises(ValueError):
            inp.cells["5"] = 1

        # verify that also the dictionary update works as expected
        inp.cells.update({"4": deepcopy(inp.cells["1"]), "2": deepcopy(inp.cells["1"])})
        assert inp.cells["4"].name == 4
        assert inp.cells["4"].values[0] == (4, "cel")
        assert inp.cells["2"].name == 2
        assert inp.cells["2"].values[0] == (2, "cel")

    def test_surf_property(self):
        inp = deepcopy(self.testInput)
        # verify that the name and values have been changed
        inp.surfs = {"1": inp.surfs["1"], "2": inp.surfs["1"]}
        assert inp.surfs["2"].name == 2
        assert inp.surfs["2"].values[0] == (2, "sur")

        # adding a new value should still force the keys and names to be
        # the same
        inp.surfs["3"] = deepcopy(inp.surfs["1"])
        assert inp.surfs["3"].name == 3
        assert inp.surfs["3"].values[0] == (3, "sur")

    def test_from_input(self):
        inp = deepcopy(self.testInput)
        self._check_macro_properties(inp)

        inp = deepcopy(self.bugInput)
        assert inp.other_data["WWE:N"]
        assert inp.other_data["WWN1:N"]

    def test_hash_cell(self):
        cell = self.testInput.cells["2"]
        hashed_cell = Input.hash_cell(cell, 12, inplace=False)
        assert hashed_cell.card() == "2 13 7.2058E-02 ( -128 129 1 -2 ) #12 \n"

    def test_hash_multiple_cells(self):
        inp = deepcopy(self.testInput)
        inp.hash_multiple_cells({12: [2, 3, 4]})
        assert inp.cells["2"].card() == "2 13 7.2058E-02 ( -128 129 1 -2 ) #12 \n"

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

    def test_renumber(self, tmpdir):
        with as_file(resources_inp.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)
        testInput.renumber(renum_all=100, update_keys=True)
        testInput.write(tmpdir.join("renum.i"))
        # check that the update keys have worked
        assert testInput.transformations["TR101"]
        assert testInput.cells["101"]
        # check some
        newinp = Input.from_input(tmpdir.join("renum.i"))
        assert newinp.cells["101"].get_f() == 225
        assert "122 0      -122 imp:n=1" in newinp.cells["122"].card()
        assert newinp.cells["122"].get_u() == 225
        assert newinp.transformations["TR101"]

    def test_add_material(self, tmpdir):
        with as_file(resources_inp.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)
        testInput.add_material_to_void_cell(testInput.cells["22"], 10, -1.1)
        testInput.add_material_to_void_cell(testInput.cells["99"], 94, 1.1)
        testInput.add_material_to_void_cell(testInput.cells["21"], 90, 1.1)
        assert testInput.cells["22"].get_m() == 10
        assert testInput.cells["22"].get_d() == -1.1
        assert testInput.cells["99"].get_m() == 94
        assert testInput.cells["99"].get_d() == 1.1
        assert testInput.cells["21"].get_m() == 4
        assert testInput.cells["21"].get_d() == -1.0

        testInput.write(tmpdir.join("new_mat.i"))
        testInput = Input.from_input(tmpdir.join("new_mat.i"))
        assert testInput.cells["22"].get_m() == 10
        assert testInput.cells["22"].get_d() == -1.1
        assert testInput.cells["99"].get_m() == 94
        assert testInput.cells["99"].get_d() == 1.1
        assert testInput.cells["21"].get_m() == 4
        assert testInput.cells["21"].get_d() == -1.0

    def test_add_cell_fill_u(self, tmpdir):
        with as_file(resources_inp.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)

        new = testInput.add_cell_fill_u(testInput.cells["99"], "U", 50, inplace=False)
        assert testInput.cells["99"].get_u() == None
        assert new.get_u() == 50
        new = testInput.add_cell_fill_u(testInput.cells["99"], "u", 50, inplace=False)
        assert testInput.cells["99"].get_u() == None
        assert new.get_u() == 50

        testInput.add_cell_fill_u(testInput.cells["99"], "U", 50, inplace=True)
        assert testInput.cells["99"].get_u() == 50

        new = testInput.add_cell_fill_u(
            testInput.cells["22"], "FILL", 250, inplace=False
        )
        assert testInput.cells["22"].get_f() == None
        assert new.get_f() == 250
        new = testInput.add_cell_fill_u(
            testInput.cells["22"], "fill", 250, inplace=False
        )
        assert testInput.cells["22"].get_f() == None
        assert new.get_f() == 250

        testInput.add_cell_fill_u(testInput.cells["22"], "FILL", 250, inplace=True)
        assert testInput.cells["22"].get_f() == 250

        testInput.write(tmpdir.join("new_fill.i"))
        testInput = Input.from_input(tmpdir.join("new_fill.i"))
        assert testInput.cells["22"].get_f() == 250
        assert testInput.cells["99"].get_u() == 50

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

    def test_merge(self):
        with as_file(resources_inp.joinpath("test_1.i")) as FILE1:
            inp1 = Input.from_input(FILE1)
        with as_file(resources_inp.joinpath("test_universe.i")) as FILE2:
            inp2 = Input.from_input(FILE2)

        # this should not be allowed due to duplicate surf
        dest = deepcopy(inp1)
        try:
            dest.merge(inp2)
            assert False
        except KeyError:
            assert True

        # renumber and try again
        dest = deepcopy(inp1)
        inp2.renumber(renum_all=1000, update_keys=True)
        dest.merge(inp2)
        assert len(dest.cells) == 14
        assert len(dest.surfs) == 9

    def _check_macro_properties(self, inp: Input):
        # check some macro properties
        assert inp.header[0].strip("\n").strip("\r") == "This is the header"
        assert len(inp.cells) == 128
        assert len(inp.surfs) == 130
        assert len(inp.materials) == 25
        assert len(inp.tally_keys) == 7
        assert len(inp.fmesh_keys) == 5

    def test_update_zaidinfo(self):
        newinput = deepcopy(self.testInput)
        newinput.update_zaidinfo(self.lm)
        assert True

    def test_update_card_keys(self):
        # test a bug
        inp = deepcopy(self.bugInput)
        inp._update_card_keys()

    def test_translate(self):
        # The test for a correct translation of material card is already done
        # in materials. here we only check that it goes trough without errors
        newinput = deepcopy(self.testInput)
        newinput.translate("00c", self.lm)
        newinput = deepcopy(self.testInput)
        newinput.translate('{"31c": "00c", "70c": "81c"}', self.lm)
        assert True

        # let's check also that abundances info is correctly added
        assert (
            "$ H-1    WEIGHT(%) 100.0 AB(%) 99.988"
            in newinput.materials.materials[0].to_text()
        )

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
        renumber_offsets = {"cells": 1}
        newinput.extract_cells(cells, outfile, renumber_offsets=renumber_offsets)
        # re-read
        inp2 = Input.from_input(outfile)
        assert len(inp2.cells) == 5
        assert len(inp2.surfs) == 10
        assert len(inp2.materials) == 3
        assert list(inp2.cells.keys()) == ["16", "24", "25", "26", "32"]

        with as_file(resources_inp.joinpath("test_1.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        outfile = os.path.join(os.path.dirname(outfile), "extract_fillers.i")

        renumber_offsets = {"cells": 500}
        mcnp_input.extract_cells(
            [50], outfile, extract_fillers=True, renumber_offsets=renumber_offsets
        )

        # re-read
        result = Input.from_input(outfile)

        assert len(result.cells) == 7
        assert mcnp_input.cells["10"].values[0][0] == 10

        # test extract without renumbering
        result.extract_cells([550], outfile)

        # test extract with strings instead of ints
        result.extract_cells(["550"], outfile)

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

    def test_missing_data_cards(self, tmpdir):
        # Check that all these data do not go missing after a rewrite
        with as_file(resources_inp.joinpath("various_bugs.i")) as file:
            inp = Input.from_input(file)
        inp.write(tmpdir.join("bug.i"))

        reinp = Input.from_input(tmpdir.join("bug.i"))
        _ = reinp.other_data["SP2"]
        _ = reinp.transformations["TR1"]
        _ = reinp.other_data["CUT:N"]
        _ = reinp.other_data["WWN1:P"]
        _ = reinp.other_data["WWN1:N"]
        _ = reinp.other_data["F96"]
        _ = reinp.other_data["F30004"]
        _ = reinp.other_data["TF30004"]

        assert True

    @pytest.mark.parametrize(
        ["cardname", "clean_name"],
        [
            ["*TR1", "TR1"],
            ["f6:n,p", "F6"],
            ["WWN1:n", "WWN1:N"],
            ["WWE:n", "WWE:N"],
            ["TF300004", "TF300004"],
        ],
    )
    def test_clean_card_name(self, cardname, clean_name):
        assert Input._clean_card_name(cardname) == clean_name

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

    def test_add_F_tally(self):
        newinput = deepcopy(self.testInput)
        cells = range(100, 150)
        energies = np.linspace(1e4, 1e5, 100)
        newinput.add_F_tally(
            4,
            ["N", "P"],
            cells,
            energies=energies,
            description="Test F4 tally",
            add_SD=True,
            add_total=True,
            multiplier="1 -52 1",
        )
        assert newinput.other_data["F4"].lines[0] == "F4:N,P\n"
        assert newinput.other_data["FC4"].lines[0] == "FC4 Test F4 tally\n"
        assert len(newinput.other_data["F4"].lines) == 4
        for card in ["F4", "E4"]:
            for line in newinput.other_data[card].lines:
                assert len(line) < 128
        # total adds 1 SD
        assert newinput.other_data["SD4"].lines[-1] == f"SD4 1 {len(cells)}R\n"
        assert newinput.other_data["FM4"].lines[0] == f"FM4 1 -52 1\n"
        cells = ["((1 2 3 4 5 6) < 10)", 12, "(((1 2 3 4 5 6) 18) < 11)"]
        newinput.add_F_tally(
            14,
            ["N"],
            cells,
            energies=energies,
            description="Test F14 tally",
            add_SD=True,
            add_total=True,
            multiplier="1 -52 1",
        )
        assert newinput.other_data["F14"].lines[0] == "F14:N\n"
        assert (
            newinput.other_data["F14"].lines[1]
            == "     ((1 2 3 4 5 6) < 10) 12 (((1 2 3 4 5 6) 18) < 11) T \n"
        )
        assert newinput.other_data["SD14"].lines[-1] == f"SD14 1 {len(cells)}R\n"

    def test_set_cell_void(self):
        newinput = deepcopy(self.testInput)
        Input.set_cell_void(newinput.cells["49"])
        text = newinput.cells["49"].card()
        text = text.replace("\r", "")
        assert text == "49   0     -128 129 48  -49               $imp:n,p=1\n"

    def test_replace_material(self):
        with as_file(resources_inp.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.replace_material(10, "-2", 4)
        assert newinp.cells["21"].get_m() == 10
        assert newinp.cells["21"].get_d() == -2

        newinp.replace_material(0, "10", 10, u_list=[125])
        assert newinp.cells["21"].get_m() == 0

        with as_file(resources_inp.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.replace_material(10, "10", 0, u_list=[125])
        assert newinp.cells["22"].get_m() == 10
        assert newinp.cells["299"].get_m() == 10
        assert newinp.cells["1"].get_m() == 0

    def test_cells_union(self):
        with as_file(resources_inp.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)

        newinp_2 = deepcopy(newinp)

        newinp.cells_union(["1", "22", "299"], None)
        assert "1" in newinp.cells
        assert not "22" in newinp.cells
        assert not "299" in newinp.cells
        text = newinp.cells["1"].card()
        text = text.replace("\r", "")
        assert text == "1 0 ( ( -1 ) : ( -22 ) )  : ( #21 #22    ) imp:n=1 fill=125\n"

        newinp_2.cells_union(["22", "299", "1"], 635)
        assert not "1" in newinp_2.cells
        assert not "22" in newinp_2.cells
        assert not "299" in newinp_2.cells
        assert "635" in newinp_2.cells
        text = newinp_2.cells["635"].card()
        text = text.replace("\r", "")
        assert (
            text == "635 0 ( ( -22 ) : ( #21 #22 ) )  : ( -1 ) imp:n=1\n        U=125\n"
        )

    def test_delete_fill_cards(self):
        with as_file(resources_inp.joinpath("test_universe2.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.delete_fill_cards()
        assert newinp.cells["1"].get_f() is None
        assert (
            newinp.cells["1"].card().replace("\r", "")
            == """1 0 -1 
      imp:n=1               
C a breaking comment
                         $ some dollar comment
                          VOL=1
"""
        )

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
        assert new_cell.card(wrap=False, comment=False).split()[3] == "("
        assert new_cell.card(wrap=False, comment=False).split()[4] == "-128"
        assert new_cell.card(wrap=False, comment=False).split()[-1] == str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[-2] == ")"

        new_cell = Input.add_surface(
            newinput.cells["27"], -sur, None, "intersect", True
        )
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in newinput.cells["27"].values
        )
        assert newinput.cells["27"].input[0].split()[-1] == r"-{:<3}"
        assert new_cell.card(wrap=False, comment=False).split()[3] == "("
        assert new_cell.card(wrap=False, comment=False).split()[4] == "-128"
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
        assert new_cell.card(wrap=False, comment=False).split()[2] == "("
        assert new_cell.card(wrap=False, comment=False).split()[3] == "-1"

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
        assert new_cell.card(wrap=False, comment=False).split()[2] == "("
        assert new_cell.card(wrap=False, comment=False).split()[3] == "-1"
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
        assert mcnp_input.cells["22"].input[0].split()[5] == r":-{:<4}"
        assert any(
            tup and isinstance(tup, tuple) and len(tup) > 0 and tup[0] == sur
            for tup in new_cell.values
        )
        assert new_cell.input[0].split()[5] == r":-{:<4}"
        assert new_cell.card(wrap=False, comment=False).split()[5] == ":" + str(-sur)
        assert new_cell.card(wrap=False, comment=False).split()[4] == ")"
        assert (
            mcnp_input.cells["22"].card(wrap=False, comment=False).split()[-1]
            == "U=125"
        )
        assert new_cell.card(wrap=False, comment=False).split()[2] == "("
        assert new_cell.card(wrap=False, comment=False).split()[3] == "-22"
        assert new_cell.card(wrap=False, comment=False).split()[0] == "50"
        assert mcnp_input.cells["22"].card(wrap=False, comment=False).split()[0] == "50"

    def test_get_density_range(self):
        with as_file(resources_inp.joinpath("test_rho_range.i")) as FILE1:
            inp = Input.from_input(FILE1)
        d_range = inp.get_densities_range()
        assert d_range.loc[30]["Min density [g/cc]"] == 0.945
        assert d_range.loc[30]["Max density [g/cc]"] == 0.946
        assert (
            d_range.loc[7]["Min density [g/cc]"]
            == d_range.loc[7]["Max density [g/cc]"]
            == pytest.approx(0.999978)
        )


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
        assert ["24050", "78195"] == reacfile.get_parents()
        reacfile.reactions[1].daughter == "78195900"

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
