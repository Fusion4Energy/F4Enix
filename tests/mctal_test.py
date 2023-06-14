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
        df = mctal.get_error_summary()
        assert df['min error'].isna().sum() == 4
        assert df['max error'].isna().sum() == 4
        assert len(df) == 36
        assert df['min error'].max() == 0.4042
        assert df['max error'].max() == 0.4042
