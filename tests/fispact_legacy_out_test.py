import pytest
from importlib.resources import files, as_file

import tests.resources.fispact_legacy_out as lib_res
from f4enix.output.fispact_legacy_out import FispactZaid, Pathway, PathwayCollection

lib_resources = files(lib_res)


class TestFispactZaid:
    def test_get_str(self):
        zaid = FispactZaid("U", 235)
        assert zaid.get_str() == "U235"

    def test_get_str_metastable(self):
        zaid = FispactZaid("U", 235, metastable=True)
        assert zaid.get_str() == "U235m"


class TestPathway:
    def test_str(self):
        parent = FispactZaid("U", 235)
        daughter = FispactZaid("U", 236)
        pathway = Pathway(parent, daughter, 1, ["fission"])
        assert str(pathway) == "U235 -fission-> U236"

    def test_str_with_intermediates(self):
        parent = FispactZaid("U", 235)
        daughter = FispactZaid("U", 236)
        intermediates = [FispactZaid("U", 236)]
        reactions = ["fission"]
        with pytest.raises(AssertionError):
            pathway = Pathway(parent, daughter, 1, reactions, intermediates)

        reactions = ["fission", "decay"]
        pathway = Pathway(parent, daughter, 1, reactions, intermediates)
        assert str(pathway) == "U235 -fission-> U236 -decay->  U236"


class TestPathwayCollection:
    def test_from_file(self):
        with as_file(lib_resources.joinpath("testSS.out")) as file:
            collection = PathwayCollection.from_file(file)
        assert len(collection.pathways) == 56
        assert str(collection.pathways[0]) == "Mn55 -(n,g)-> Mn56"
        assert str(collection.pathways[1]) == "Fe56 -(n,p)-> Mn56"
        assert len(collection.pathways[-1].intermediates) == 8

    def test_to_dataframe(self):
        with as_file(lib_resources.joinpath("testSS.out")) as file:
            collection = PathwayCollection.from_file(file)
        df = collection.to_dataframe()
        assert len(df) == 56
