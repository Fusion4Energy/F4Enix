from f4enix.material_library import CONCRETE, MICROTHERM
from f4enix.output.decay_pathways import AVAILABLE_SPECTRA, PathwayLibrary


class TestPathwayLibrary:
    def test_get_pathways(self):
        plib = PathwayLibrary()
        df = plib.get_pathways(
            spectrum=AVAILABLE_SPECTRA[:-1], dose=95, materials=[CONCRETE, MICROTHERM]
        )
        assert len(df) == 38

        df2 = plib.get_pathways(
            spectrum=AVAILABLE_SPECTRA, dose=99, materials=[CONCRETE]
        )
        assert len(df2) > len(df)

        df3 = plib.get_pathways(dose=99)
        assert len(df3) == 140

    def test_filter_pathways(self):
        plib = PathwayLibrary()
        df = plib.filter_pathways(
            spectrum=AVAILABLE_SPECTRA[:-1], dose=95, materials=[CONCRETE, MICROTHERM]
        )
        assert df.reset_index()["material"].unique().tolist() == [
            "Concrete",
            "Microtherm",
        ]

        df2 = plib.filter_pathways(
            spectrum=AVAILABLE_SPECTRA[:-1],
            dose=95,
            materials=[CONCRETE, MICROTHERM],
            cooling_times=["10y"],
        )

        assert len(df) > len(df2)
