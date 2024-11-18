from importlib.resources import as_file, files

import pytest

import tests.resources.aceplotter as res
from f4enix.input.aceplotter import (
    get_collapsed_cross_section_nuclide,
    get_energy_cross_section,
)

RESOURCES = files(res)


class TestACEplots:
    local_ace = {
        ("1001", "FENDL32"): {"path": RESOURCES.joinpath("01H_001")},
        ("1002", "FENDL32"): {"path": RESOURCES.joinpath("01H_002")},
    }

    def test_get_collapsed_cross_section_nuclide(self):
        collapsed_total_xs, collapsed_total_elastic = (
            get_collapsed_cross_section_nuclide(
                this="1001",
                types=["total", "elastic"],
                library="FENDL32",
                ace_filepaths=self.local_ace,
                group_structure="VITAMIN-J-42",
            )
        )
        assert len(collapsed_total_xs) == len(collapsed_total_elastic) == 42

    @pytest.mark.parametrize("this", ["1001", "H"])
    def test_get_energy_cross_section(self, this):
        e_li, xs_li = get_energy_cross_section(
            this=this,
            types=[105, 16],
            library="FENDL32",
            ace_filepaths=self.local_ace,
        )
        assert len(e_li) == len(xs_li[0]) == len(xs_li[1])
