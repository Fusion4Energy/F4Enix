import pytest
import math
from importlib.resources import as_file, files
from f4enix.input.irradiation import (
    Nuclide,
    IrradiationScenario,
    _process_irr_line,
    TCF_Computer,
)
from f4enix.input.d1suned import IrradiationFile
from tests.resources import irradiation

RES = files(irradiation)


class TestNuclide:
    def test_from_formula_basic(self):
        n = Nuclide.from_formula("Li3")
        assert n.zaid == 3003
        assert not n.metastable
        assert not n.IRS_active
        assert n.lib is None

    def test_from_formula_complex(self):
        n = Nuclide.from_formula("irsLi3m.99c")
        assert n.zaid == 3003
        assert n.metastable
        assert n.IRS_active
        assert n.lib == "99c"

    def test_from_int_string_basic(self):
        n = Nuclide.from_int_string("3003")
        assert n.zaid == 3003
        assert not n.metastable
        assert not n.IRS_active
        assert n.lib is None

    def test_from_int_string_complex(self):
        n = Nuclide.from_int_string("9993003900.99c")
        assert n.zaid == 3003
        assert n.metastable
        assert n.IRS_active
        assert n.lib == "99c"

    def test_write_to_formula(self):
        n = Nuclide(3003, metastable=True, IRS_active=True, lib="99c")
        formula = n.write_to_formula()
        assert formula == "irsLi3m.99c"

    def test_write_to_int_string(self):
        n = Nuclide(3003, metastable=True, IRS_active=True, lib="99c")
        int_string = n.write_to_int_string()
        assert int_string == "9993003900.99c"


class TestIrradiationScenario:
    def test_from_fispact(self):
        with as_file(RES.joinpath("inp_fispact.i")) as fisp_file:
            irr_scenario = IrradiationScenario.from_fispact(fisp_file)

        assert irr_scenario.pulses[-1].time == 400
        assert irr_scenario.pulses[-1].intensity == 0

    def test_from_legacy_d1stime(self):
        with as_file(RES.joinpath("inp_d1stime")) as d1s_file:
            irr_scenario = IrradiationScenario.from_legacy_d1stime(d1s_file)

        assert len(irr_scenario.pulses) == 62


class TestTFC_Computer:
    def test_get_lambda(self):
        tcf_computer = TCF_Computer()
        decay_constant = tcf_computer.get_lambda(Nuclide.from_formula("Co60m"))
        assert pytest.approx(decay_constant) == math.log(2) / 628.02

        assert tcf_computer.get_lambda(Nuclide.from_formula("H1")) == 0.0

    def test_compute_correction_factors(self):
        # test en masse all the correction factors of 93c
        with as_file(RES.joinpath("irrad_93c.txt")) as infile:
            irr_file = IrradiationFile.from_text(infile)

        # read the irr scenario
        with as_file(RES.joinpath("irrad_d1stime.i")) as file:
            irr_scenario = IrradiationScenario.from_legacy_d1stime(file)

        tcf = TCF_Computer()
        for irrad in irr_file.irr_schedules:
            nuclide = irrad.daughter
            factors = tcf.compute_correction_factors(
                irr_scenario, [nuclide], norm=2.0e19
            )
            print(nuclide)
            assert pytest.approx(factors[0][0], rel=5e-2) == float(irrad.times[0])
            # assert np.allclose(factors, expected_factors, rtol=2e-2)


def test_process_irr_line():
    line = "(0   /4381 s/     1.111111e16 /50 w/)*4\n"
    pulses = _process_irr_line(line)
    assert len(pulses) == 8
    assert pytest.approx(pulses[0].time) == 4381
    assert (
        pytest.approx(pulses[-1].time) == 50 * 60 * 60 * 24 * 7
    )  # 50 weeks in seconds
    print(pulses)
