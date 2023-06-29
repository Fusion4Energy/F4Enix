from importlib.resources import files, as_file

from f4enix.output.MCNPoutput import Output
import tests.resources.output as resources
import pytest

RESOURCES = files(resources)


class TestOutput:
    def test_print_lp_debug(self, tmpdir):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)

        # This has no LP
        outp.print_lp_debug(tmpdir)
        checks = outp.get_statistical_checks()
        assert len(checks) == 12

        with as_file(RESOURCES.joinpath('out_lp.txt')) as file:
            outp = Output(file)
        # This has LP
        outp.print_lp_debug(tmpdir)

    def test_get_NPS(self):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)
        assert outp.get_NPS() == int(1e4)

    @pytest.mark.parametrize(['table', 'shape'],
                             [[60, (4, 11)],
                              [126, (3, 10)]])
    def test_get_table(self, table, shape):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)

        df = outp.get_table(table)
        assert df.shape == shape

    def test_get_fwf_format_from_string(self):
        stringa = '   sdasdaas     scdcsdc    dscds  csc'
        specs = Output._get_fwf_format_from_string(stringa)
        assert specs == [11, 12, 9, 5]

    def test_get_tally_stat_checks(self):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)

        checks = outp.get_tally_stat_checks(46)
        assert len(checks) == 3
        assert checks.loc['passed?'].value_counts()['no'] == 1
        assert checks.loc['passed?'].value_counts()['yes'] == 9
        assert checks.loc['passed?', 'PDF slope'] == 'yes'

        for cell in [14, 4]:
            try:
                checks = outp.get_tally_stat_checks(cell)
                assert False
            except ValueError:
                assert True

    def test_get_stat_checks_tfc_bins(self):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)
        values = outp.get_statistical_checks_tfc_bins()
        assert values[24] == 'All zeros'
        assert values[32] == 'Passed'
        assert values[44] == 'Missed'

    def test_get_stat_checks_table(self):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)
        df = outp.get_stat_checks_table()
        assert len(df) == 12
        assert df['PDF slope'].value_counts()['yes'] == 6
        assert df['PDF slope'].value_counts()['no'] == 1
