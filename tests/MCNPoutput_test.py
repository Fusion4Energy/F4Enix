from importlib.resources import files, as_file
import pandas as pd
import os

from f4enix.output.MCNPoutput import Output
from f4enix.input.MCNPinput import Input
import tests.resources.output as resources
import pytest

RESOURCES = files(resources)


class TestOutput:
    def test_print_lp_debug(self, tmpdir):
        name = "LPdebug_out_lp."
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)

        # This has no LP
        outp.print_lp_debug(tmpdir)
        checks = outp.get_statistical_checks_tfc_bins()
        assert len(checks) == 12

        with as_file(RESOURCES.joinpath("out_lp.txt")) as file:
            outp = Output(file)
        # This has LP
        outp.print_lp_debug(tmpdir)
        df = pd.read_csv(os.path.join(tmpdir, name + "csv"))
        assert len(df) == 10
        assert len(df.columns) == 6
        df = pd.read_excel(os.path.join(tmpdir, name + "xlsx"))
        assert df["count"].sum() == 10

        outp.print_lp_debug(tmpdir, get_cosine=False)
        df = pd.read_csv(os.path.join(tmpdir, name + "csv"))
        assert len(df) == 10
        assert len(df.columns) == 3

        # LP debug with universes
        with as_file(RESOURCES.joinpath("test_lp_u.o")) as file:
            outp = Output(file)
        with as_file(RESOURCES.joinpath("test_lp_u.i")) as file:
            inp = Input.from_input(file)
        outp.print_lp_debug(tmpdir, input_mcnp=inp)
        df = pd.read_excel(
            os.path.join(tmpdir, "LPdebug_test_lp_u.xlsx"), sheet_name="by universe"
        )
        assert df["count"].values[-1] == 42

    def test_get_tot_lp(self):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)
        assert outp.get_tot_lp() is None

        with as_file(RESOURCES.joinpath("test_lp_u.o")) as file:
            outp = Output(file)
        assert outp.get_tot_lp() == 677

    def test_get_NPS(self):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)
        assert outp.get_NPS() == int(1e4)

    @pytest.mark.parametrize(["table", "shape"], [[60, (4, 11)], [126, (3, 10)]])
    def test_get_table(self, table, shape):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)

        df = outp.get_table(table)
        assert df.shape == shape

    def test_get_fwf_format_from_string(self):
        stringa = "   sdasdaas     scdcsdc    dscds  csc"
        specs = Output._get_fwf_format_from_string(stringa)
        assert specs == [11, 12, 9, 5]

    def test_get_tally_stat_checks(self):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)

        checks = outp.get_tally_stat_checks(46)
        assert len(checks) == 3
        assert checks.loc["passed?"].value_counts()["no"] == 1
        assert checks.loc["passed?"].value_counts()["yes"] == 9
        assert checks.loc["passed?", "PDF slope"] == "yes"

        for cell in [14, 4]:
            try:
                checks = outp.get_tally_stat_checks(cell)
                assert False
            except ValueError:
                assert True

    def test_get_stat_checks_tfc_bins(self):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)
        values = outp.get_statistical_checks_tfc_bins()
        assert values[24] == "All zeros"
        assert values[32] == "Passed"
        assert values[44] == "Missed"

    def test_get_stat_checks_table(self):
        with as_file(RESOURCES.joinpath("test_o")) as file:
            outp = Output(file)
        df = outp.get_stat_checks_table()
        assert len(df) == 12
        assert df["PDF slope"].value_counts()["yes"] == 6
        assert df["PDF slope"].value_counts()["no"] == 1

    @pytest.mark.parametrize(
        ["rel_filepath", "expected"],
        [
            ["out_lp.txt", "6.2"],
            ["FNG-W_Alo", "d1suned411"],
            ["d1sunedout/example", "d1suned411"],
        ],
    )
    def test_get_code_version(self, rel_filepath: os.PathLike, expected: str):
        with as_file(RESOURCES.joinpath(rel_filepath)) as file:
            outp = Output(file)
        assert outp.get_code_version() == expected
