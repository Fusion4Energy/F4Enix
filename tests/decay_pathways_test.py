from f4enix.material_library import CONCRETE, MICROTHERM
from f4enix.output.decay_pathways import AVAILABLE_SPECTRA, PathwayLibrary


class TestPathwayLibrary:
    def test_get_pathways(self):
        plib = PathwayLibrary()
        df = plib.get_pathways(
            spectrum=AVAILABLE_SPECTRA[:-1], dose=95, materials=[CONCRETE, MICROTHERM]
        )
        assert len(df) == 3
