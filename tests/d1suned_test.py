import sys
import os
import pytest
from importlib.resources import files, as_file

from f4enix.input.d1suned import (
    Reaction,
    ReactionFile,
    Irradiation,
    IrradiationFile,
    REACFORMAT,
)
from f4enix.input.libmanager import LibManager
import tests.resources.d1suned as res
import f4enix.resources as pkg_res

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

    def test_get_daughters(self):
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)
        daughters = irrfile.get_daughters()
        assert daughters == ["24051", "25054", "26055", "26059", "27062", "27062900"]

    def test_get_irrad(self):
        with as_file(RESOURCES.joinpath("irr_test")) as inp:
            irrfile = IrradiationFile.from_text(inp)
        # Check the None
        irradiation = irrfile.get_irrad("20051")
        assert irradiation is None
        # Check true value
        irradiation = irrfile.get_irrad("26055")
        assert irradiation.daughter == "26055"

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
        assert irrfile.irr_schedules[0].daughter == "24051"
        assert irrfile.irr_schedules[1].daughter == "26055"

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
        assert irr.daughter == "24051"
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
        assert reaction.parent == "26054.99c"
        assert reaction.MT == "102"
        assert reaction.daughter == "26055"
        assert reaction.comment == "Fe55"

    def test_fromtext2(self):
        """
        Test different formatting possibilities
        """
        text = "26054.99c 102   26055 Fe55  and some"
        reaction = Reaction.from_text(text)
        assert reaction.parent == "26054.99c"
        assert reaction.MT == "102"
        assert reaction.daughter == "26055"
        assert reaction.comment == "Fe55 and some"

    def test_changelib(self):
        """
        Test change library tag
        """
        rec = Reaction("26054.99c", "102", "26055")
        rec.change_lib("31c")
        assert rec.parent == "26054.31c"

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
        assert rx.parent == "26054.31c"
        assert rx.MT == "102"
        assert rx.daughter == "26055"
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
            assert reaction.parent[-3:] == newlib

    def test_translation2(self, lm: LibManager):
        """
        test translation with libmanager where parents are not available

        """

        with as_file(RESOURCES.joinpath("reac2")) as inp:
            reac_file = ReactionFile.from_text(inp)

        newlib = "99c"

        reac_file.change_lib(newlib, libmanager=lm)

        for reaction in reac_file.reactions:
            assert reaction.parent[-3:] != newlib

    def test_get_parents(self):
        with as_file(RESOURCES.joinpath("reac_fe")) as inp:
            reac_file = ReactionFile.from_text(inp)
        parents = reac_file.get_parents()
        assert parents == ["26054", "26056", "26057", "26058"]
