from importlib.resources import as_file, files

import numpy as np
import pytest

import tests.resources.mctal as mctal_res
from f4enix.output.mctal import Mctal

# ruff: noqa: PLR2004

MCTAL_RESOURCES = files(mctal_res)


class TestMctal:
    def test_general(self):
        with as_file(MCTAL_RESOURCES.joinpath("test_m")) as inp:
            MCTAL1 = Mctal(inp)
        assert len(MCTAL1.tallydata) == 11
        assert np.isclose(MCTAL1.tallydata[6]["Value"], 1.61148)

    def test_get_error_summary(self):
        with as_file(MCTAL_RESOURCES.joinpath("error_summary.m")) as inp:
            mctal = Mctal(inp)
        for option in [True, False]:
            df = mctal.get_error_summary(include_abs_err=option)
            assert df["min rel error"].isna().sum() == 4
            assert df["max rel error"].isna().sum() == 4
            assert len(df) == 36
            assert df["min rel error"].max() == 0.4042
            assert df["max rel error"].max() == 0.4042

    def test_error_cmodel(self):
        with as_file(MCTAL_RESOURCES.joinpath("C_Modelm")) as inp:
            Mctal(inp)
        assert True

    def test_error_detector(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_time")) as inp:
            Mctal(inp)
        assert True

    def test_error_tmesh(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_tmesh")) as inp:
            Mctal(inp)
        assert True

    def test_error_cosbin(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_cosbin")) as inp:
            Mctal(inp)
        assert True

    def test_error_radio(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_radio")) as inp:
            Mctal(inp)
        assert True

    def test_detector(self):
        with as_file(MCTAL_RESOURCES.joinpath("detectors.m")) as inp:
            mctal = Mctal(inp)
        assert "Dir" in mctal.tallydata[15].columns
        assert "Dir" in mctal.tallydata[5].columns
        # one energy total per Dir bin
        assert len(mctal.tallydata[5]) == 462 * 2
        assert len(mctal.tallydata[15]) == 703 * 2
        assert len(mctal.totalbin[5]) == 2

    def test_fm(self):
        """This tally makes sure that the Mctal class can handle tallies with a complex
        tally multiplier card even when several results are zero. Tally 704 has 14
        different tally segments and a FM card with 106 reaction rates. This test checks
        that the resulting tallydata has 14*106 = 1484 rows.
        """
        with as_file(MCTAL_RESOURCES.joinpath("mctal_fm")) as inp:
            mctal = Mctal(inp)
        assert mctal.tallydata[704].shape == (1484, 4)

    def test_multi_total(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_multi_total")) as inp:
            mctal = Mctal(inp)
        df = mctal.tallydata[4]
        assert list(df.columns) == ["Cells", "Segments", "Energy", "Value", "Error"]
        assert len(df) == 2 * 3 * 3
        # values in the fixture encode the 1-based (cell, segment, energy) bin
        expected = [
            100 * f + 10 * s + e
            for f in range(1, 3)
            for s in range(1, 4)
            for e in range(1, 4)
        ]
        assert np.allclose(df["Value"], expected)

        row = df[(df["Cells"] == 20) & (df["Segments"] == "total")]
        assert list(row["Energy"]) == [1.0, 20.0, "total"]
        assert np.allclose(row["Value"], [231, 232, 233])

        totals = mctal.totalbin[4]
        assert len(totals) == 2 * 5
        assert ((totals["Segments"] == "total") | (totals["Energy"] == "total")).all()

        mctal.remove_totals(4)
        assert len(mctal.tallydata[4]) == 2 * 2 * 2

    def test_d1s_relative_contribution(self):
        with as_file(MCTAL_RESOURCES.joinpath("mctal_daughter")) as inp:
            mctal = Mctal(inp)
        mctal.set_d1s_relative_contribution(124)
        assert (
            pytest.approx(mctal.tallydata[124].loc[16, "Normalized Value"])
            == 0.18430563795291874
        )
        # test normalization with many bins
        with as_file(MCTAL_RESOURCES.joinpath("mctal_daughter_2")) as inp:
            mctal = Mctal(inp)

        # with pytest.raises(KeyError):
        #     mctal.set_d1s_relative_contribution(124)
        mctal.set_d1s_relative_contribution(5104)
        df = mctal.tallydata[5104]
        assert (
            pytest.approx(
                df[(df["Time"] == 1.0) & (df["Segments"] == 1)][
                    "Normalized Value"
                ].sum()
            )
            == 1.0
        )
