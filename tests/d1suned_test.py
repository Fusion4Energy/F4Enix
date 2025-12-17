import os
import pytest
import numpy as np
from importlib.resources import files, as_file

from f4enix.constants import TIME_UNITS
from f4enix.input.d1suned import (
    Reaction,
    ReactionFile,
    Irradiation,
    IrradiationFile,
)
from f4enix.input.libmanager import LibManager
import tests.resources.d1suned as res
import f4enix.resources as pkg_res
from f4enix.input.irradiation import Nuclide
from f4enix.input.irradiation import IrradiationScenario, Pulse

RESOURCES = files(res)
PKG_RESOURCES = files(pkg_res)
# INP = os.path.join(cp, 'TestFiles', 'parserD1S', 'reac_fe')
# Files


class TestIrradiationFile:
    @pytest.mark.parametrize("file", ["irr_test", "irr_test2"])
    def test_fromtext(self, file):
        """
        Test parsing irradiation file
        """
        with as_file(RESOURCES.joinpath(file)) as inp:
            irrfile = IrradiationFile.from_text(inp)
        if file == "irr_test":
            self._assert_file1(irrfile)
        elif file == "irr_test2":
            self._assert_file2(irrfile)

    @staticmethod
    def _assert_file1(irrfile: IrradiationFile):
        assert len(irrfile.irr_schedules) == 6
        TestIrradiation.assert_irr(irrfile.irr_schedules[0])

    @staticmethod
    def _assert_file2(irrfile: IrradiationFile):
        assert len(irrfile.irr_schedules) == 4
        TestIrradiation.assert_irr(irrfile.irr_schedules[0])

    @pytest.mark.parametrize("file", ["irr_test", "irr_test2"])
    def test_write(self, tmpdir, file):
        """
        Test writing irradiation file 1
        """
        with as_file(RESOURCES.joinpath(file)) as inp:
            irrfile = IrradiationFile.from_text(inp)
        irrfile.write(tmpdir)
        outfile = os.path.join(tmpdir, irrfile.name)
        irrfile = IrradiationFile.from_text(outfile)
        if file == "irr_test":
            self._assert_file1(irrfile)
        elif file == "irr_test2":
            self._assert_file2(irrfile)

    def test_from_irrad_schedules(self, tmpdir):
        pulses = [Pulse(10, 5), Pulse(50, 0)] * 2 + [Pulse(100, 10)]
        irr_scenario1 = IrradiationScenario(pulses, name="Scenario 1")

        pulses2 = [Pulse(20, 10), Pulse(40, 0)] * 2 + [Pulse(80, 5)]
        irr_scenario2 = IrradiationScenario(pulses2, name="Scenario 2")

        # multiple scenarios case
        daughter_list = [
            Nuclide.from_formula("Co62"),
            Nuclide.from_formula("Co62m"),
        ]
        irrfile = IrradiationFile.from_irradiation_schedules(
            daughter_list, [irr_scenario1, irr_scenario2], norm=2
        )
        irrfile.write(tmpdir)
        assert len(irrfile.irr_schedules) == 2

        # raise ValueError
        daughter_list = [
            Nuclide.from_formula("Co62"),
            Nuclide.from_formula("Co62m"),
            Nuclide.from_formula("irsCo62m"),
        ]
        with pytest.raises(ValueError):
            irrfile = IrradiationFile.from_irradiation_schedules(
                daughter_list,
                [irr_scenario1, irr_scenario2],
                norm=2,
                scale_IRS={"irsCo62m": [2, 3]},
            )

        # IRS case
        irrfile = IrradiationFile.from_irradiation_schedules(
            daughter_list, [irr_scenario1], norm=2, scale_IRS={"irsCo62m": [1, 2]}
        )
        irrfile.name = "irs_test"
        irrfile.write(tmpdir)
        assert len(irrfile.irr_schedules) == 3
        assert len(irrfile.irr_schedules[0].times) == 2
        assert irrfile.nsc == 2
        assert irrfile.irr_schedules[0].times[0] == irrfile.irr_schedules[0].times[1]
        assert irrfile.irr_schedules[-1].times[0] != irrfile.irr_schedules[-1].times[1]

    def test_get_daughters(self):
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)
        daughters = irrfile.get_daughters()
        expected = ["24051", "25054", "26055", "26059", "27062", "27062900"]
        for one, other in zip(daughters, expected):
            assert one.write_to_int_string() == other

    def test_get_irrad(self):
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)
        # Check the None
        irradiation = irrfile.get_irrad("20051")
        assert irradiation is None
        # Check true value
        irradiation = irrfile.get_irrad("26055")
        assert irradiation.daughter.write_to_int_string() == "26055"

    def test_select_daughters_irradiation_file(self):
        """
        Updates a D1S irradiation file selecting a subset of daughters from a list

        Parameters
        ----------
        daughters : list.
            daughter zaids to be selected

        """
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)
        ans = irrfile.select_daughters_irradiation_file(["24051", "26055"])
        # Keep only useful irradiations
        assert ans is True
        assert len(irrfile.irr_schedules) == 2
        assert irrfile.irr_schedules[0].daughter.write_to_int_string() == "24051"
        assert irrfile.irr_schedules[1].daughter.write_to_int_string() == "26055"

    def test_multiple_irradiations(self, tmpdir):
        """
        Test multiple irradiation lines
        """
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)

        # Test for KeyError when a daughter is missing in times_dict
        missing_key_times = {
            "24051": ["5.982e+00", "5.697e+00"],
            "26055": ["4.487e+00", "6.364e-01"],
            # Missing "25054" and other daughters
        }
        with pytest.raises(KeyError, match="No time correction factors provided"):
            irrfile.add_irradiation_times(missing_key_times)

        # Test for KeyError when an invalid key is provided in times_dict
        invalid_key_times = {
            "24051": ["5.982e+00", "5.697e+00"],
            "25054": ["5.881e+00", "1.829e+00"],
            "99999": ["4.487e+00", "6.364e-01"],  # Invalid key
        }
        with pytest.raises(
            KeyError,
            match="Invalid key '99999' provided. It does not match any daughter.",
        ):
            irrfile.add_irradiation_times(invalid_key_times)

        # Test for ValueError when lists of different lengths are provided
        extra_times = {
            "24051": ["5.982e+00", "5.697e+00"],
            "25054": ["5.881e+00", "1.829e+00"],
            "26055": ["4.487e+00", "6.364e-01"],
            "26059": ["6.645e+00", "5.651e+00"],
            "27062": ["1.336e+00", "4.151e-01"],
            "27062900": ["4.151e-01", "4.151e-01", "1.0e+00"],  # Extra value
        }
        with pytest.raises(
            ValueError,
            match="All input time correction factor lists in `times_dict` must have the same length.",
        ):
            irrfile.add_irradiation_times(extra_times)

        new_times = {
            "24051": ["5.982e+00", "5.697e+00"],
            "25054": ["5.881e+00", "1.829e+00"],
            "26055": ["4.487e+00", "6.364e-01"],
            "26059": ["6.645e+00", "5.651e+00"],
            "27062": ["1.336e+00", "4.151e-01"],
            "27062900": ["4.151e-01", "4.151e-01"],
        }
        irrfile.add_irradiation_times(new_times)

        assert irrfile.nsc == 4
        irrfile.write(tmpdir)

        new_irrfile_path = os.path.join(tmpdir, "irrad")
        new_irrfile = IrradiationFile.from_text(new_irrfile_path)
        assert new_irrfile.nsc == 4
        assert new_irrfile.irr_schedules[0].times[-1] == "5.697e+00"
        new_irrfile.remove_irradiation_time(3)
        assert new_irrfile.nsc == 3
        assert new_irrfile.irr_schedules[0].times[-1] == "5.982e+00"
        irrfile.irr_schedules[0].modify_time_val(3, 4.56)
        assert float(irrfile.irr_schedules[0].times[3]) == 4.56

    def test_rescale_dose(self):
        """
        Test rescaling of dose with different irradiation scenarios
        """
        with as_file(RESOURCES.joinpath("irr_test_rescale")) as inp:
            irrfile_1 = IrradiationFile.from_text(inp)

        df_scaling = irrfile_1.get_scaling_factors_cooling_time(
            1, (9.91e6, TIME_UNITS.SECOND)
        )
        assert pytest.approx(df_scaling.loc[73182, "9910000.0s"], rel=1e-2) == 0.5

        irr_scenario = IrradiationScenario(pulses=[Pulse(10, 1e10)])
        irr_scenario.set_cooling_times(
            [(3600, TIME_UNITS.SECOND), (100, TIME_UNITS.DAY)]
        )
        df_scaling_2 = irrfile_1.get_scaling_factors_new_scenario(1, irr_scenario, 1e10)
        assert pytest.approx(df_scaling_2.loc[73182, "3600s"], rel=1e-2) == 6.704e-3


class TestIrradiation:
    def test_reading(self):
        """
        Test the reading of irradiation line
        """
        text = "   24051     2.896e-07    5.982e+00    5.697e+00     Cr51"
        irr = Irradiation.from_text(text, 2)
        self.assert_irr(irr)

    @staticmethod
    def assert_irr(irr: Irradiation):
        """
        Assert irradiation
        """
        assert irr.daughter.write_to_int_string() == "24051"
        assert irr.lambd == "2.896e-07"
        assert irr.times[0] == "5.982e+00"
        assert irr.times[1] == "5.697e+00"
        assert irr.comment == "Cr51"

    def test_equivalence(self):
        # Equivalent
        text = "   24051     2.896e-07    5.982e+00    5.697e+00     Cr51"
        irr1 = Irradiation.from_text(text, 2)
        text = "   24051     2.896e-07    5.982e+00    5.697     "
        irr2 = Irradiation.from_text(text, 2)
        assert irr1 == irr2

        # Not equal
        text = "   24051     2.896e-07    5.697e+00    5.982e+00     Cr51"
        irr3 = Irradiation.from_text(text, 2)
        text = "   24051     2.896e-07    5.697e+00    Cr51"
        irr4 = Irradiation.from_text(text, 1)
        assert irr1 != irr3
        assert irr1 != {}
        assert irr1 != irr4


class TestReaction:
    def test_fromtext1(self):
        """
        Test different formatting possibilities
        """
        text = "   26054.99c  102  26055     Fe55"
        reaction = Reaction.from_text(text)
        assert reaction.parent.write_to_int_string() == "26054.99c"
        assert reaction.MT == "102"
        assert reaction.daughter.write_to_int_string() == "26055"
        assert reaction.comment == "Fe55"

    def test_fromtext2(self):
        """
        Test different formatting possibilities
        """
        text = "26054.99c 102   26055 Fe55  and some"
        reaction = Reaction.from_text(text)
        assert reaction.parent.write_to_int_string() == "26054.99c"
        assert reaction.MT == "102"
        assert reaction.daughter.write_to_int_string() == "26055"
        assert reaction.comment == "Fe55 and some"

    def test_changelib(self):
        """
        Test change library tag
        """
        parent = Nuclide.from_int_string("26054.99c")
        daughter = Nuclide.from_int_string("26055")
        rec = Reaction(parent, "102", daughter)
        rec.change_lib("31c")
        assert rec.parent.write_to_int_string() == "26054.31c"

    def test_write(self):
        """
        check writing
        """
        text = "26054.99c  102  26055     Fe55 and  some"
        reaction = Reaction.from_text(text)
        ftext = reaction._get_text()
        comptext = ["26054.99c", "102", "26055", "Fe55 and some"]
        assert comptext == ftext


class TestReactionFile:
    @pytest.fixture
    def lm(self):
        xsdirpath = os.path.join(PKG_RESOURCES, "xsdir.txt")
        # isotopes_file = os.path.join(root, 'jade', 'resources', 'Isotopes.txt')
        return LibManager(xsdir_path=xsdirpath)

    def test_fromtext(self):
        """
        right number of reactions
        """
        with as_file(RESOURCES.joinpath("reac_fe")) as inp:
            reac_file = ReactionFile.from_text(inp)
        print(reac_file.reactions)
        assert len(reac_file.reactions) == 11

    def test_write(self, tmpdir):
        """
        writing works
        """
        with as_file(RESOURCES.joinpath("reac_fe")) as inp:
            reac_file = ReactionFile.from_text(inp)

        reac_file.write(tmpdir)
        outpath = os.path.join(tmpdir, "react")
        newfile = ReactionFile.from_text(outpath)
        # Remove the temporary file
        os.remove(outpath)
        # do some operations
        newfile.change_lib("31c")
        assert len(newfile.reactions) == 11
        # Check also first line
        rx = newfile.reactions[0]
        assert rx.parent.write_to_int_string() == "26054.31c"
        assert rx.MT == "102"
        assert rx.daughter.write_to_int_string() == "26055"
        assert rx.comment == "Fe55"

    def test_translation(self, lm: LibManager):
        """
        test translation with libmanager where parents are available

        """
        newlib = "98c"

        with as_file(RESOURCES.joinpath("reac_fe")) as inp:
            reac_file = ReactionFile.from_text(inp)

        reac_file.change_lib(newlib, libmanager=lm)

        for reaction in reac_file.reactions:
            assert reaction.parent.write_to_int_string()[-3:] == newlib

    def test_translation2(self, lm: LibManager):
        """
        test translation with libmanager where parents are not available

        """

        with as_file(RESOURCES.joinpath("reac2")) as inp:
            reac_file = ReactionFile.from_text(inp)

        newlib = "99c"

        reac_file.change_lib(newlib, libmanager=lm)

        for reaction in reac_file.reactions:
            assert reaction.parent.write_to_int_string()[-3:] != newlib

    def test_get_parents(self):
        with as_file(RESOURCES.joinpath("reac_fe")) as inp:
            reac_file = ReactionFile.from_text(inp)
        parents = reac_file.get_parents()
        expected = ["26054", "26056", "26057", "26058"]
        for one, other in zip(parents, expected):
            assert one.write_to_int_string().split(".")[0] == other
