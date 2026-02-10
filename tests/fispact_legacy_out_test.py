import pytest
from importlib.resources import files, as_file

import tests.resources.fispact_legacy_out as lib_res
from f4enix.output.fispact_legacy_out import (
    FispactZaid,
    Pathway,
    PathwayCollection,
    FispactOutput,
)

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
        assert str(pathway) == "U235 -fission-> U236 -decay-> U236"


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

    @pytest.mark.parametrize(
        "file", ["testSS.out", "test_weird_pathways.out", "no_pathways.out"]
    )
    def test_weird_pathways(self, file):
        with as_file(lib_resources.joinpath(file)) as file:
            collection = PathwayCollection.from_file(file)
            assert len(collection.pathways) > 0

    def test_parse_pathway(self):
        text = """
 Target nuclide Bi212     99.837% of inventory given by  1 path
 --------------------

 path  1  99.837% Th 232 ---(B)--- Ra228 ---(d)--- Ac228 ---(d)--- Th228 ---(b)--- Ra224 ---(b)--- Rn220 ---(d)--- Po216 ---(d)--- Pb212 ---(d)--- 
                     99.99%(a)      100.00%(b-)     100.00%(b-)     100.00%(a)      100.00%(a)      100.00%(a)      100.00%(a)      100.00%(b-)    
                      0.01%(n,na)                                     0.00%(n,na)     0.00%(n,na)                                                  

 path continued   Pb 212 ---(d)--- Bi212 ---(S)---
                    100.00%(b-)"""
        pathway = PathwayCollection._parse_pathway(text)
        assert pathway.parent.get_str() == "Th232"
        assert pathway.daughter.get_str() == "Bi212"
        assert pathway.reactions[-1] == "(b-)"

        # normal path
        text = """
 path  1 100.000% Th232 ---(D)--- Ra228 ---(D)--- Ac228 ---(d)--- Th228 ---(d)--- Ra224 ---(S)--- 
                    100.00%(a)      100.00%(b-)     100.00%(b-)     100.00%(a)  """
        pathway = PathwayCollection._parse_pathway(text)
        assert len(pathway.reactions) == 4


class TestFispactOutput:
    @pytest.fixture
    def outp(self) -> FispactOutput:
        with as_file(lib_resources.joinpath("testSS.out")) as file:
            output = FispactOutput(file, cooling_times=["1e2", "1e4"])
        return output

    def test_filter_by_cum_dose(self, outp: FispactOutput):
        df = outp.filter_by_cum_dose(95, "1e2")
        assert len(df) == 3
        assert df.iloc[-1]["Cumulative dose sum"] > 95

        df = outp.filter_by_cum_dose(95, "1e2", add_pathways=True)
        assert len(df) == 8

    # def test_filter_bugged(self):
    #     with as_file(lib_resources.joinpath("no_pathways.out")) as file:
    #         output = FispactOutput(file, cooling_times=["1e2", "1e4"])

    #     df = output.filter_by_cum_dose(95, "1e2", add_pathways=True)
