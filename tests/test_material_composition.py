from f4enix.material_library import MaterialComposition


class TestMaterialComposition:
    def test_init(self):
        SILVER = MaterialComposition("Pure Silver", [("Ag", 100)])
        assert SILVER.perc == [100]
        assert str(SILVER) == "MaterialComposition(Pure Silver)"

    def test_eq(self):
        SILVER = MaterialComposition("Pure Silver", [("Ag", 100)])
        SILVER2 = MaterialComposition("Pure Silver", [("Ag", 100)])
        SILVER3 = MaterialComposition("Pure Silver", [("Pb", 100)])
        assert SILVER == SILVER2
        assert SILVER != SILVER3
