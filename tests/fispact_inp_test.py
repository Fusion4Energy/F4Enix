from __future__ import annotations

from unittest.mock import Mock

import pytest

from f4enix.core.irradiation import IrradiationScenario, Pulse
from f4enix.input.fispact_inp import FispactInp
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import Material


@pytest.fixture
def simple_material() -> Material:
    return Material.from_zaids(
        [(1001, 0.5), (8016, 0.5)],
        libman=LibManager(),
        lib="31c",
        name="testmat",
        mat_id=1,
    )


@pytest.fixture
def simple_irradiation_scenario() -> IrradiationScenario:
    pulse_1 = Pulse(time=10, intensity=1.0)
    pulse_2 = Pulse(time=20, intensity=2.0)
    cooling = Pulse(time=30, intensity=0.0)
    return IrradiationScenario(
        pulses=[pulse_1, pulse_2],
        cooling_times=[cooling],
    )


class TestFispactInp:
    def test_valid_add_material_calls(
        self,
        simple_material: Material,
    ):
        inp = FispactInp()

        with pytest.raises(ValueError):
            inp.add_material(simple_material, style="mass")

        with pytest.raises(ValueError):
            inp.add_material(simple_material, style="fuel")

        with pytest.raises(ValueError):
            inp.add_material(simple_material, style="fuel", density=2.0)

        with pytest.raises(ValueError):
            inp.add_material(simple_material, style="fuel", volume=1.0)

        with pytest.raises(ValueError):
            inp.add_material(simple_material, style="mass", volume=1.0)

        inp.add_material(simple_material, style="mass", mass=12.5)
        inp.add_material(simple_material, style="mass", volume=1.0, density=2.0)
        inp.add_material(simple_material, style="fuel", density=2.0, volume=1.0)
        inp.add_material(simple_material, style="fuel", mass=12.5, volume=1.0)
        inp.add_material(simple_material, style="fuel", mass=12.5, density=2.0)

    def test_add_irradiation_scenario_adds_pulses_and_cooling_times(self):
        inp = FispactInp()
        inp.inp = Mock()

        pulse_1 = Pulse(time=10, intensity=1.0)
        pulse_2 = Pulse(time=20, intensity=2.0)
        cooling = Pulse(time=30, intensity=0.0)
        scenario = IrradiationScenario(
            pulses=[pulse_1, pulse_2],
            cooling_times=[cooling],
        )

        inp.add_irradiation_scenario(scenario, norm=2.0)

        inp.inp.addIrradiation.assert_any_call(pulse_1.time, pulse_1.intensity * 2.0)
        inp.inp.addIrradiation.assert_any_call(pulse_2.time, pulse_2.intensity * 2.0)
        inp.inp.addCooling.assert_called_once_with(cooling.time)

    def test_save(
        self,
        simple_material: Material,
        simple_irradiation_scenario: IrradiationScenario,
        tmpdir,
    ):
        inp = FispactInp()
        inp.add_material(simple_material, style="mass", mass=1)
        inp.add_irradiation_scenario(simple_irradiation_scenario, norm=1.0)
        inp.save(tmpdir.join("test.inp"))
