from importlib.resources import files

import pytest

import tests.resources.ace as res
from f4enix.input.ace import _xsdir_to_ace_filepaths, get_xs
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import Material

RESOURCES = files(res)
LM = LibManager()


class TestACEplots:
    local_ace = {
        ("1001", "FENDL32"): {"path": RESOURCES.joinpath("01H_001")},
        ("1002", "FENDL32"): {"path": RESOURCES.joinpath("01H_002")},
    }

    hydrogen = Material.from_zaids([(1000, 1)], LM, "31c")

    @pytest.mark.parametrize(
        ["this", "ace_files", "group_structure"],
        [
            ["1001", local_ace, None],
            ["1000", local_ace, None],
            [hydrogen, local_ace, "VITAMIN-J-42"],
        ],
    )
    def test_get_xs(self, this, ace_files, group_structure):
        e_li, xs_li = get_xs(
            this=this,
            types=[1, 102],
            library="FENDL32",
            ace_filepaths=ace_files,
            group_structure=group_structure,
        )
        assert len(e_li) == len(xs_li[0]) == len(xs_li[1])

    def test_get_xs_with_xsdir(self):
        e_li, xs_li = get_xs(
            this=self.hydrogen,
            types=["total", "(n,total)"],
            library="31c",
            ace_filepaths=RESOURCES.joinpath("xsdir_test"),
            temperature=300,
        )
        assert len(e_li) == len(xs_li[0]) == len(xs_li[1])

    def test_xsdir_to_ace_filepaths(self):
        file = RESOURCES.joinpath("xsdir_test")
        res = _xsdir_to_ace_filepaths(file)
        for isotope in ["1001", "1002"]:
            assert (isotope, "31c") in res.keys()
