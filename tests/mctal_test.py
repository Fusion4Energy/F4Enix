import numpy as np
from copy import deepcopy
from importlib.resources import files, as_file

from f4enix.output.mctal import Mctal
import tests.resources.mctal as mctal_res

mctal_resources = files(mctal_res)


class TestMctal:
    def test_general(self):
        with as_file(mctal_resources.joinpath('test_m')) as inp:
            MCTAL1 = Mctal(inp)
        assert len(MCTAL1.tallydata) == 11
        assert np.isclose(MCTAL1.tallydata[6]['Value'], 1.61148)

    def test_get_error_summary(self):
        with as_file(mctal_resources.joinpath('error_summary.m')) as inp:
            mctal = Mctal(inp)
        for option in [True, False]:
            df = mctal.get_error_summary(include_abs_err=option)
            assert df['min rel error'].isna().sum() == 4
            assert df['max rel error'].isna().sum() == 4
            assert len(df) == 36
            assert df['min rel error'].max() == 0.4042
            assert df['max rel error'].max() == 0.4042

    def test_error_cmodel(self):
        with as_file(mctal_resources.joinpath('C_Modelm')) as inp:
            mctal = Mctal(inp)
        assert True

    def test_error_detector(self):
        with as_file(mctal_resources.joinpath('mctal_time')) as inp:
            mctal = Mctal(inp)
        assert True

    def test_error_tmesh(self):
        with as_file(mctal_resources.joinpath('mctal_tmesh')) as inp:
            mctal = Mctal(inp)
        assert True

    def test_error_cosbin(self):
        with as_file(mctal_resources.joinpath('mctal_cosbin')) as inp:
            mctal = Mctal(inp)
        assert True

    def test_error_radio(self):
        with as_file(mctal_resources.joinpath('mctal_radio')) as inp:
            mctal = Mctal(inp)
        assert True
