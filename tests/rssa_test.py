from __future__ import annotations
from importlib.resources import files, as_file
import pytest

from f4enix.output.rssa import RSSA
import tests.resources.rssa as res

RESOURCES = files(res)


class TestRSSA:
    with as_file(RESOURCES.joinpath("small_cyl.w")) as infile:
        rssa = RSSA(infile)

    def test_print(self):
        print(self.rssa)

    # def test_correct_parsing(self):
    #     #TODO
    #     pass

    @pytest.mark.parametrize(
        [
            "particle",
            "z_int",
            "theta_int",
            "norm",
            "value_range",
            "x_range",
            "z_range",
            "outfolder",
            "filename",
        ],
        #
        [
            ["n", 10, 10, 1, None, None, None, None, None],
            ["n", 20, 20, 100, [1e-5, 1e-3], [0, 50], [0, 20], "tmp", "random"],
        ],
    )
    def test_plot_cyl(
        self,
        particle,
        z_int,
        theta_int,
        norm,
        value_range,
        x_range,
        z_range,
        outfolder,
        filename,
        tmpdir,
    ):
        if outfolder == "tmp":
            outfolder = tmpdir

        self.rssa.plot_cyl(
            particle=particle,
            z_int=z_int,
            x_range=x_range,
            theta_int=theta_int,
            z_range=z_range,
            value_range=value_range,
            norm=norm,
            outfolder=outfolder,
            filename=filename,
        )

    def test_properties(self):
        self.rssa.energies
        self.rssa.histories
        self.rssa.mask_neutron_tracks
        self.rssa.mask_photon_tracks
        self.rssa.wgt
        self.rssa.x
        self.rssa.y
        self.rssa.z
        self.rssa.histories
        self.rssa.type
        assert True
