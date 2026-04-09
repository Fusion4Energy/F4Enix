from f4enix.core.material_library import CONCRETE, MICROTHERM
from f4enix.output.decay_pathways import (
    AVAIL_SPECTRUM,
    AVAIL_IRR_SCENARIO,
    PathwayLibrary,
)


class TestPathwayLibrary:
    def test_get_pathways(self):
        plib = PathwayLibrary()

        df = plib.get_pathways(
            spectrum=[AVAIL_SPECTRUM.FirstWall500MW],
            dose=99,
            materials=[CONCRETE, MICROTHERM],
            irradiation_scenario=[AVAIL_IRR_SCENARIO.DT1],
            cooling_times=["10 y"],
        )
        assert len(df) == 9

        df = plib.get_pathways(dose=99)
        assert len(df) == 247

        # temporary fix since this will be changed soon with new data
        # df3 = plib.get_pathways(dose=99)
        # assert len(df3) == 140

    def test_filter_pathways(self):
        plib = PathwayLibrary()
        df = plib.filter_pathways(
            spectrum=[AVAIL_SPECTRUM.FirstWall500MW],
            dose=99,
            materials=[CONCRETE, MICROTHERM],
            irradiation_scenario=[AVAIL_IRR_SCENARIO.DT1],
            cooling_times=["10 y"],
        )
        assert len(df) == 10
