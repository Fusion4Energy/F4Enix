import pytest
from importlib.resources import files, as_file

import tests.resources.fispact_legacy_out as lib_res
from f4enix.output.fispact_legacy_out import (
    Pathway,
    PathwayCollection,
    FispactOutput,
)
from f4enix.core.irradiation import Nuclide
from f4enix.core.irradiation import TCF_Computer

lib_resources = files(lib_res)


class TestPathway:
    def test_str(self):
        parent = Nuclide(element="U", isotope=235)
        daughter = Nuclide(element="U", isotope=236)
        pathway = Pathway(parent, daughter, 1, ["fission"])
        assert str(pathway) == "U235 -fission-> U236"

    def test_str_with_intermediates(self):
        parent = Nuclide(element="U", isotope=235)
        daughter = Nuclide(element="U", isotope=236)
        intermediates = [Nuclide(element="U", isotope=236)]
        reactions = ["fission"]
        with pytest.raises(AssertionError):
            pathway = Pathway(parent, daughter, 1, reactions, intermediates)

        reactions = ["fission", "decay"]
        pathway = Pathway(parent, daughter, 1, reactions, intermediates)
        assert str(pathway) == "U235 -fission-> U236 -decay-> U236"

    def test_reduce(self):
        parent = Nuclide(element="U", isotope=235)
        daughter = Nuclide(element="U", isotope=236)
        intermediates = [Nuclide(element="U", isotope=236, metastable=True)]
        reactions = ["A", "IT"]
        pathway = Pathway(parent, daughter, 100, reactions, intermediates=intermediates)
        reduced = pathway.reduce()

        assert str(reduced) == "U235 -A-> U236"

        pathway2 = Pathway(parent, daughter, 100, reactions[1:])
        reduced = pathway2.reduce()
        assert pathway2 == reduced

        # Try a cutoff
        parent = Nuclide(element="U", isotope=235)
        daughter = Nuclide(element="Co", isotope=60)
        intermediates = [Nuclide(element="Co", isotope=60, metastable=True)]
        reactions = ["A", "IT"]
        pathway = Pathway(parent, daughter, 100, reactions, intermediates=intermediates)

        tfc_computer = TCF_Computer()
        reduce1 = pathway.reduce(600, tfc_computer=tfc_computer)
        reduce2 = pathway.reduce(650, tfc_computer=tfc_computer)

        assert len(reduce1.intermediates) == 1
        assert reduce2.intermediates is None

    def test_from_string(self):
        string = "Ni58 -(n,p)-> Co58m -(IT)-> Co58"
        path = Pathway.from_string(string)
        assert str(path) == string

    def test_equal(self):
        zaid1 = Nuclide(element="U", isotope=235)
        zaid2 = Nuclide(element="U", isotope=235, metastable=True)
        reactions = ["fission"]
        reactions2 = ["fission1"]

        pathway1 = Pathway(zaid1, zaid2, 10, reactions)
        pathway2 = Pathway(zaid1, zaid2, 50, reactions)
        pathway3 = Pathway(zaid1, zaid1, 1, reactions)
        pathway4 = Pathway(zaid1, zaid1, 1, reactions2)
        pathway5 = Pathway(zaid2, zaid2, 1, reactions)

        assert pathway1 == pathway2
        assert pathway1 != pathway3
        assert pathway3 == pathway4  # different spelling of reaction is fine
        assert pathway1 != pathway5

    @pytest.mark.parametrize(
        ["text", "expected"],
        [
            [
                [
                    "path  1  99.915% Cu 63 ---(R)--- Cu 62 ---(S)--- ",
                    "                    100.00%(n,2n)  ",
                ],
                False,
            ],
            [
                [
                    " path  4  42.748% Ni 58 ---(R)--- Co 58m---(b)--- Co 58 ---(S)--- ",
                    "  100.00%(n,p)    100.00%(IT)",
                ],
                False,
            ],
            [
                [
                    "path  1  96.911% Sn112 ---(R)--- Sn111 ---(b)--- In111 ---(S)--- ",
                    "   100.00%(n,2n)   100.00%(b+)   ",
                ],
                True,
            ],
        ],
    )
    def test_is_multistep(self, text, expected):
        path = PathwayCollection._parse_pathway(text)
        assert path.is_multistep() == expected


class TestPathwayCollection:
    def test_from_file(self):
        with as_file(lib_resources.joinpath("testSS.out")) as file:
            collection = PathwayCollection.from_file(file)
        assert len(collection.pathways) == 56
        assert collection.pathways[0].parent.write_to_formula() == "Mn55"
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
        assert pathway.parent.write_to_formula() == "Th232"
        assert pathway.daughter.write_to_formula() == "Bi212"
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
