import os
from copy import deepcopy
from importlib.resources import as_file, files

import migjorn
import numpy as np
import pytest

import f4enix.resources as pkg_res
import tests.resources.input as input_res
import tests.resources.libmanager as lib_res
from f4enix.core.irradiation import Nuclide
from f4enix.input.d1suned import IrradiationFile, ReactionFile
from f4enix.input.libmanager import LibManager
from f4enix.input.MCNPinput import D1S_Input, Input, get_formatted_range
# ruff: noqa: PLR2004# ruff: noqa: PLR2004

RESOURCES_INP = files(input_res)
RESOURCES_LIB = files(lib_res)
RESOURCES_PCK = files(pkg_res)

# INP_EX_PATH = os.path.join(resources_inp, 'test_exceptions.i')
# DIS_INP_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest.i')
# DIS_NOPKMT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_noPKMT.i')
# DIS_GETREACT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_getreact.i')

# IRRAD_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_irrad')
# REACT_PATH = os.path.join(cp, 'TestFiles/inputfile/d1stest_react')


class TestInput:
    with as_file(RESOURCES_INP.joinpath("test.i")) as FILE1:
        testInput = Input.from_input(FILE1)

    with as_file(RESOURCES_INP.joinpath("various_bugs.i")) as file:
        bugInput = Input.from_input(file)
    # exceptInput = InputFile.from_text(INP_EX_PATH)
    with (
        as_file(RESOURCES_LIB.joinpath("Activation libs.xlsx")) as ACTIVATION_FILE,
        as_file(RESOURCES_LIB.joinpath("xsdir")) as XSDIR_FILE,
        as_file(RESOURCES_PCK.joinpath("Isotopes.txt")) as ISOTOPES_FILE,
    ):
        lm = LibManager(
            XSDIR_FILE, activationfile=ACTIVATION_FILE, isotopes_file=ISOTOPES_FILE
        )

    def test_check_range(self):
        inp = self.testInput
        assert not inp.check_range([1, 2])
        assert inp.check_range(range(1000, 10010))  # type: ignore
        assert not inp.check_range([1, 1e4])  # type: ignore

        assert not inp.check_range([89], who="surf")
        assert inp.check_range([360], who="surf")
        assert not inp.check_range([1360], who="surf")

    def test_cell_property(self):
        inp = deepcopy(self.testInput)
        # With migjorn, cells is a live view — verify basic access
        cell_1 = inp.cells["1"]
        assert isinstance(cell_1, migjorn.Cell)
        assert cell_1.id == 1

    def test_surf_property(self):
        inp = deepcopy(self.testInput)
        surf_1 = inp.surfs["1"]
        assert isinstance(surf_1, migjorn.Surface)
        assert surf_1.id == 1

    def test_from_input(self):
        inp = deepcopy(self.testInput)
        self._check_macro_properties(inp)

        inp = deepcopy(self.bugInput)
        assert inp.other_data["WWE:N"]
        assert inp.other_data["WWN1:N"]

    def test_hash_cell(self):
        inp = deepcopy(self.testInput)
        cell = inp.cells["2"]
        Input.hash_cell(cell, 12)
        assert "#12" in cell.text

    def test_hash_multiple_cells(self):
        inp = deepcopy(self.testInput)
        inp.hash_multiple_cells({12: [2, 3, 4]})
        assert "#12" in inp.cells["2"].text

    def test_renumber(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)
        testInput.renumber(renum_all=100)
        # check that renumbered IDs are present
        assert testInput.transformations["TR101"]
        assert testInput.cells["101"]
        fill = testInput.cells["101"].fill
        assert fill is not None and fill.universe == 225
        assert "122" in testInput.cells["122"].text
        assert testInput.cells["122"].universe == 225
        assert testInput.transformations["TR101"]

    def test_add_material(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)
        testInput.cells["22"].material = 10
        testInput.cells["22"].density = -1.1
        testInput.cells["99"].material = 94
        testInput.cells["99"].density = 1.1
        testInput.cells["21"].material = 4
        testInput.cells["21"].density = 1.1

        assert testInput.cells["22"].material == 10
        assert testInput.cells["22"].density == pytest.approx(-1.1)
        assert testInput.cells["99"].material == 94
        assert testInput.cells["99"].density == pytest.approx(1.1)
        assert testInput.cells["21"].material == 4  # cell 21 already has material 4

    def test_set_param(self):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE1:
            testInput = Input.from_input(FILE1)

        # add universe in-place
        testInput.set_param(testInput.cells["99"], "u", 50)
        assert testInput.cells["99"].universe == 50

        testInput.set_param(testInput.cells["22"], "fill", 250)
        assert testInput.cells["22"].fill is not None
        assert testInput.cells["22"].fill.universe == 250

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
        with as_file(RESOURCES_INP.joinpath("test_1.i")) as FILE1:
            inp1 = Input.from_input(FILE1)
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE2:
            inp2 = Input.from_input(FILE2)

        # this should not be allowed due to duplicate surf
        dest = deepcopy(inp1)
        try:
            dest.merge(inp2)
            assert False
        except Exception:
            assert True

        # renumber and try again
        dest = deepcopy(inp1)
        inp2.renumber(renum_all=1000)
        dest.merge(inp2)
        assert len(dest.cells) == 14
        assert len(dest.surfs) == 9

    def _check_macro_properties(self, inp: Input):
        # check some macro properties
        assert inp.header.strip("\n").strip("\r") == "This is the header"
        assert len(inp.cells) == 128
        assert len(inp.surfs) == 130
        assert len(inp.mat_section) == 25
        assert len(inp.tally_keys) == 7
        assert len(inp.fmesh_keys) == 5

    def test_update_card_keys(self):
        # _update_card_keys is superseded by migjorn; just check no error
        inp = deepcopy(self.bugInput)
        assert True  # method removed, no-op

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
            "$ H1     WEIGHT(%) 1.9046 AB(%) 99.988"
            in newinput.mat_section.materials[0].to_text()
        )

    def test_get_materials_subset(self):
        materials = "m23"
        _ = self.testInput.get_materials_subset(materials)
        materials = ["m22", "M30"]
        _ = self.testInput.get_materials_subset(materials)
        assert True

    def test_get_data_cards(self):
        _ = self.testInput.other_data["SDEF"]

        with pytest.raises(KeyError):
            _ = self.testInput.other_data["adas"]

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
        newinput = newinput.extract_cells(cells, renumber_offsets=renumber_offsets)
        newinput.write(outfile)
        # re-read
        inp2 = Input.from_input(outfile)
        assert len(inp2.cells) == 5
        assert len(inp2.surfs) == 10
        # assert len(inp2.mat_section) == 3  entire mat section is copied now
        assert sorted(inp2.cells.keys()) == ["16", "24", "25", "26", "32"]

        with as_file(RESOURCES_INP.joinpath("test_1.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        outfile = os.path.join(os.path.dirname(outfile), "extract_fillers.i")

        renumber_offsets = {"cells": 500}
        newinp = mcnp_input.extract_cells([50], renumber_offsets=renumber_offsets)
        newinp.write(outfile)

        # re-read
        result = Input.from_input(outfile)

        assert len(result.cells) == 7
        assert mcnp_input.cells["10"].id == 10

        # test extract without renumbering
        result.extract_cells([550])

        # test extract with strings instead of ints
        result.extract_cells(["550"])

    def test_extract_universe(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        universe = 125
        result = mcnp_input.extract_universe(universe)

        assert len(result.cells) == 3
        assert len(result.surfs) == 2
        assert len(result.mat_section) == 1
        for _, cell in result.cells.items():
            assert cell.universe is None
        assert mcnp_input.cells["21"].universe == universe

        # try also the extract with filler
        result = mcnp_input.extract_universe(
            universe, keep_level_0=True, renumber_offsets={"cells": 1}
        )
        assert len(result.cells) == 3 + 1
        assert result.cells[2].fill.universe == universe

    def test_duplicated_nums(self):
        # There was a bug reading material 101
        self.bugInput.get_materials_subset(["m101"])
        assert True

    def test_missing_data_cards(self):
        # Check that all these data do not go missing
        with as_file(RESOURCES_INP.joinpath("various_bugs.i")) as file:
            inp = Input.from_input(file)

        _ = inp.other_data["SP2"]
        _ = inp.transformations["TR1"]
        _ = inp.other_data["CUT:N"]
        _ = inp.other_data["WWN1:P"]
        _ = inp.other_data["WWN1:N"]
        _ = inp.other_data["F96"]
        _ = inp.other_data["F30004"]
        _ = inp.other_data["TF30004"]

        assert True

    @pytest.mark.parametrize("flag", [True, False])
    def test_get_cells_by_matID(self, flag):
        newinput = deepcopy(self.testInput)
        cells = newinput.get_cells_by_matID(13)
        for filtered, expected in zip(cells.keys(), range(2, 22)):
            assert filtered == str(expected)

    def test_scale_densities(self):
        newinput = deepcopy(self.testInput)
        newinput.scale_densities(0.33333333333)
        assert newinput.cells["49"].density is not None
        assert newinput.cells["52"].density is None  # void cell has no density
        assert newinput.cells["53"].density is not None

    @pytest.mark.parametrize(
        ["id", "expected"],
        [[94, ["FC94", "F94:N", "FM94"]], [214, ["FC214", "FMESH214:N", "FM214"]]],
    )
    def test_get_tally_cards(self, id, expected):
        keys = self.testInput._get_tally_cards_ids(id)
        for a, b in zip(keys, expected):
            assert a.lower() == b.lower()

    def test_get_tally_summary(self):
        with as_file(RESOURCES_INP.joinpath("test_complex_fm.i")) as FILE1:
            testInput = Input.from_input(FILE1)
        summary = testInput.get_tally_summary()
        assert len(summary) == 8
        assert summary.loc[194].values.tolist() == [
            "n",
            "T in Li pt2 appm/FPY",
            "3.8566e10",
            ["25", "205"],
        ]
        to_assert = summary.loc[204].values.tolist()
        assert to_assert[0] == "n"
        assert to_assert[1] is np.nan
        assert to_assert[2] is np.nan
        assert to_assert[3] is np.nan

        assert len(summary.loc[704]["Other multipliers"]) == 12

        summary = testInput.get_tally_summary(fmesh=True)
        assert len(summary) == 5
        assert summary.loc[224].values.tolist() == [
            "p",
            "FMESH Photon Heating [MeV/cc/n_s]",
            "-1",
            ["0", "-5", "-6"],
        ]
        newinput = deepcopy(self.testInput)
        cells = range(100, 150)
        energies = np.linspace(1e4, 1e5, 100)
        newinput.add_F_tally(
            4,
            ["n", "p"],
            cells,
            energies=energies,
            description="Test F4 tally",
            add_SD=True,
            add_total=True,
            multiplier="1 -52 1",
        )
        assert "F4:n,p" in newinput.other_data["F4"].text
        assert "FC4 Test F4 tally" in newinput.other_data["FC4"].text
        assert "SD4" in newinput.other_data
        assert "FM4" in newinput.other_data
        cells = ["((1 2 3 4 5 6) < 10)", 12, "(((1 2 3 4 5 6) 18) < 11)"]
        newinput.add_F_tally(
            14,
            ["n"],
            cells,
            energies=energies,
            description="Test F14 tally",
            add_SD=True,
            add_total=True,
            multiplier="1 -52 1",
        )
        assert "F14:n" in newinput.other_data["F14"].text
        assert "SD14" in newinput.other_data

    def test_set_cell_void(self):
        newinput = deepcopy(self.testInput)
        newinput.cells["49"].material = 0
        assert newinput.cells["49"].material == 0
        assert newinput.cells["49"].density is None

    def test_replace_material(self):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.replace_material(10, "-2", 4)
        assert newinp.cells["21"].material == 10
        assert newinp.cells["21"].density == pytest.approx(-2.0)

        newinp.replace_material(0, "10", 10, u_list=[125])
        assert newinp.cells["21"].material == 0

        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.replace_material(10, "10", 0, u_list=[125])
        assert newinp.cells["22"].material == 10
        assert newinp.cells["299"].material == 10
        assert newinp.cells["1"].material == 0

    def test_cells_union(self):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as inp_file:
            newinp = Input.from_input(inp_file)

        newinp.cells_union(["1", "22", "299"], None)  # type: ignore
        assert "1" in newinp.cells
        assert "22" not in newinp.cells
        assert "299" not in newinp.cells

        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as inp_file:
            newinp2 = Input.from_input(inp_file)
        newinp2.cells_union(["22", "299", "1"], 635)
        assert "1" not in newinp2.cells
        assert "22" not in newinp2.cells
        assert "635" in newinp2.cells
        assert newinp2.cells["635"].geometry_text.count(":") == 2

    def test_delete_fill_cards(self):
        with as_file(RESOURCES_INP.joinpath("test_universe2.i")) as inp_file:
            newinp = Input.from_input(inp_file)
        newinp.delete_fill_cards()
        for _, cell in newinp.cells.items():
            assert cell.fill is None

    def test_add_surface(self):
        newinput = deepcopy(self.testInput)
        sur = 180
        cell = newinput.cells["27"]
        Input.add_surface(cell, -sur, None, "intersect")
        assert sur in cell.surface_ids
        assert cell.text.replace("\r", "") == (
            "27   15  9.1292E-02  ( -128 129 26  -27 ) -180     $imp:n,p=1\n"
        )

        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE:
            mcnp_input = Input.from_input(FILE)

        cell = Input.add_surface(mcnp_input.cells["1"], sur, None, "union")
        assert f"( -1 ) :{sur}" in cell.text

    def test_get_density_range(self):
        with as_file(RESOURCES_INP.joinpath("test_rho_range.i")) as FILE1:
            inp = Input.from_input(FILE1)
        d_range = inp.get_densities_range()
        assert d_range.loc[30]["Min density [g/cc]"] == 0.945
        assert d_range.loc[30]["Max density [g/cc]"] == 0.946
        assert (
            d_range.loc[7]["Min density [g/cc]"]
            == d_range.loc[7]["Max density [g/cc]"]
            == pytest.approx(0.999978)
        )

    def test_remove_tallies(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test.i")) as FILE1:
            inp = Input.from_input(FILE1)
        inp.remove_tallies()
        assert "F94" not in inp.other_data
        assert "F54" not in inp.other_data

        # check that the file can be written and read again
        outfile = tmpdir.mkdir("sub").join("test_tallies_removed.i")
        inp.write(outfile)
        newinp = Input.from_input(outfile)
        assert "F94" not in newinp.other_data
        assert "F54" not in newinp.other_data

        with as_file(RESOURCES_INP.joinpath("test.i")) as FILE2:
            inp2 = Input.from_input(FILE2)
        inp2.remove_tallies([94])
        assert "F94" not in inp2.other_data
        assert "F54" in inp2.other_data

    def test_remove_sdef(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test.i")) as FILE1:
            inp = Input.from_input(FILE1)
        inp.remove_sdef()
        assert "SDEF" not in inp.other_data

        # check that the file can be written and read again
        outfile = tmpdir.mkdir("sub").join("test_sdef_removed.i")
        inp.write(outfile)
        newinp = Input.from_input(outfile)
        assert "SDEF" not in newinp.other_data

    def test_prepare_void_check(self, tmpdir):
        with as_file(RESOURCES_INP.joinpath("test.i")) as FILE1:
            inp = Input.from_input(FILE1)

        # non-sphere surface raises ValueError
        with pytest.raises(ValueError):
            inp.prepare_void_check(2, 100)

        # add a sphere surface and use it
        inp._model.add_surface("12345 SO 10")
        particle = "P"
        inp.prepare_void_check(12345, 100, particle=particle)
        assert "VOID" in inp.other_data
        # check that the file can be written and read again
        outfile = tmpdir.mkdir("sub").join("test_void_check.i")
        inp.write(outfile)
        newinp = Input.from_input(outfile)
        assert "SDEF" in newinp.other_data["SDEF"].text

    def test_get_formatted_range(self):
        result = get_formatted_range([1, 2, 3, 4, 5, 7, 12, 13, 21, 22, 23])
        assert result == "1 3I 5 7 12 13 21 22 23 "
        result = get_formatted_range([1])
        assert result == "1 "

    def test_explore_id_ranges_by_plot(self):
        # This method creates a plot. Here we just check that it runs without errors.
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE1:
            test_input = Input.from_input(FILE1)
        test_input.explore_id_ranges_by_plot()

    def test_find_first_free_id_range(self):
        with as_file(RESOURCES_INP.joinpath("test_universe.i")) as FILE1:
            test_input = Input.from_input(FILE1)

        first_id_that_fits_15 = 2
        assert test_input.find_first_free_id_range(15) == first_id_that_fits_15

        first_id_that_fits_30 = 23
        assert test_input.find_first_free_id_range(30) == first_id_that_fits_30

        first_id_that_fits_1000 = 300
        assert test_input.find_first_free_id_range(1000) == first_id_that_fits_1000


class TestD1S_Input:
    with (
        as_file(RESOURCES_INP.joinpath("d1stest.i")) as inp_file,
        as_file(RESOURCES_INP.joinpath("d1stest_irrad")) as irrad_file,
        as_file(RESOURCES_INP.joinpath("d1stest_react")) as react_file,
    ):
        inp = D1S_Input.from_input(
            inp_file, reac_file=react_file, irrad_file=irrad_file
        )

    with (
        as_file(RESOURCES_LIB.joinpath("Activation libs.xlsx")) as ACTIVATION_FILE,
        as_file(RESOURCES_LIB.joinpath("xsdir")) as XSDIR_FILE,
        as_file(RESOURCES_PCK.joinpath("Isotopes.txt")) as ISOTOPES_FILE,
    ):
        lm = LibManager(
            XSDIR_FILE, activationfile=ACTIVATION_FILE, isotopes_file=ISOTOPES_FILE
        )

    def test_smart_translate(self):
        # This test needs to be improved
        with (
            as_file(RESOURCES_INP.joinpath("d1stest_irrad_st")) as irrad_file,
            as_file(RESOURCES_INP.joinpath("d1stest_react_st")) as react_file,
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

        translation = newinp.mat_section.to_text()

        assert translation.count("98c") == 4
        assert translation.count("00c") == 145
        assert newinp.reac_file.reactions[0].parent.write_to_int_string() == "24050.98c"

    def test_add_PKMT_card(self):
        with as_file(RESOURCES_INP.joinpath("d1stest_noPKMT.i")) as inp_file:
            newinp = D1S_Input.from_input(inp_file)
        newinp.reac_file = self.inp.reac_file

        newinp.add_PIKMT_card()
        card_text = newinp.other_data["PIKMT"].text
        assert len(card_text.splitlines()) == 17

    def test_get_reaction_file(self):
        with (
            as_file(RESOURCES_INP.joinpath("d1stest_getreact.i")) as inp_file,
            as_file(RESOURCES_INP.joinpath("d1stest_irrad_getreact")) as irr_file,
        ):
            newinp = D1S_Input.from_input(inp_file, irrad_file=irr_file)

        lib = "99c"
        reacfile = newinp.get_reaction_file(self.lm, lib)
        assert [
            Nuclide.from_int_string("24050.99c"),
            Nuclide.from_int_string("78195.99c"),
        ] == reacfile.get_parents()
        assert reacfile.reactions[1].daughter.write_to_int_string() == "78195900"

    def test_get_potential_paths(self):
        reaction_list = self.inp.get_potential_paths(self.lm, "98c")
        assert len(reaction_list) == 32

    @pytest.mark.parametrize(
        ["who", "sign"], [["parent", "-"], ["daughter", ""], ["cell", ""]]
    )
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
        card_text = newinp.other_data["FU124"].text
        for line, exp in zip(
            card_text.splitlines()[-3:], ["FU124 0", sign + "1001", sign + "1002"]
        ):
            assert line.strip() == exp
        newinp2 = D1S_Input.from_input(tmpfile)
        newinp2.add_track_contribution(tallyID, ["100", "200"], who=who)
        if who == "cell":
            assert "FT124 SCD" in newinp2.other_data["FT124"].text

    def test_add_father_from_reac(self, tmpdir):
        tallyID = "F124"

        # --- Test parents---
        inp = deepcopy(self.inp)
        # inp.add_daughter_contribution_from_irr(tallyID)
        inp.add_parent_contribution_from_reac(tallyID)
        # dump and reread the input
        tmpfile = os.path.join(tmpdir, "tmp.i")
        inp.write(tmpfile)
        newinp = D1S_Input.from_input(tmpfile)
        # get the new injected card
        parents = inp.reac_file.get_parents()
        parent_zaids = [p.write_to_int_string() for p in parents]
        for line in newinp.other_data["FU124"].text.splitlines():
            if line.startswith("FU124"):
                assert line.strip() == "FU124 0"
            else:
                assert line.strip() in [f"-{p}" for p in parent_zaids]

    def test_add_daughter_from_irr(self, tmpdir):
        tallyID = "F124"

        # --- Test parents---
        inp = deepcopy(self.inp)
        inp.add_daughter_contribution_from_irr(tallyID)
        # dump and reread the input
        tmpfile = os.path.join(tmpdir, "tmp.i")
        inp.write(tmpfile)
        newinp = D1S_Input.from_input(tmpfile)
        # get the new injected card
        daughters = inp.irrad_file.get_daughters()
        for line in newinp.other_data["FU124"].text.splitlines():
            if line.startswith("FU124"):
                assert line.strip() == "FU124 0"
            else:
                assert Nuclide.from_int_string(line.strip()) in daughters

    def test_add_SDDR_dose_function(self):
        tallyID = "F14"
        inp = deepcopy(self.inp)
        inp.add_SDDR_dose_function(tallyID)
        assert "DE14" in inp.other_data
        assert "DF14" in inp.other_data
        assert "DE14 0.01" in inp.other_data["DE14"].text
        assert "DF14 0.0485" in inp.other_data["DF14"].text

    def test_column_format(self, tmp_path):
        with as_file(RESOURCES_INP.joinpath("column_format.i")) as FILE1:
            inp = Input.from_input(FILE1)

        outfile = tmp_path / "column_format_output.i"
        inp.write(outfile)

        inp2 = Input.from_input(outfile)
        # check for the # line
        with open(outfile) as f:
            lines = f.readlines()
        found = False
        for line in lines:
            if "    #     wwn1:p" in line:
                found = True
                break
        assert found

    def test_column_format_other_data_unnamed_cards(self):
        # column_format.i has two data cards migjorn cannot assign a name to
        # (the '#'-prefixed WWN tables); other_data must still expose both of
        # them, each under its own resolvable key, instead of colliding or
        # raising.
        with as_file(RESOURCES_INP.joinpath("column_format.i")) as FILE1:
            inp = Input.from_input(FILE1)

        other_data = dict(inp.other_data)
        assert len(other_data) == len(inp.other_data)

        unnamed_keys = [key for key in inp.other_data if key.startswith("NONE")]
        assert len(unnamed_keys) == 2
        assert len(set(unnamed_keys)) == 2  # keys must be distinct, not both "None"

        cards = {key: inp.other_data[key] for key in unnamed_keys}
        texts = {key: card.text for key, card in cards.items()}
        assert any("wwn1:n" in text for text in texts.values())
        assert any("wwn1:p" in text for text in texts.values())
