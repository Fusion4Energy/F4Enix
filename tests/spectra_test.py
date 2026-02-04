from f4enix.core.spectra import Spectra, SpectraParsingError
from importlib.resources import files, as_file
import tests.resources.spectra as res
import pytest
import numpy as np
from f4enix.core.egroups import GROUP_STRUCTURES
from pathlib import Path

RES = files(res)


class TestSpectra:
    @pytest.fixture
    def spectra(self) -> Spectra:
        with as_file(RES.joinpath("FLUXESS1")) as fisp_file:
            spectra = Spectra.from_fispact(GROUP_STRUCTURES["VITAMIN-J-175"], fisp_file)
        return spectra

    @pytest.mark.parametrize("file", ["FLUXESS1"])
    def test_from_fispact(self, file):
        # check that spectra are read correctly from fispact files
        with as_file(RES.joinpath(file)) as fisp_file:
            spectra = Spectra.from_fispact(GROUP_STRUCTURES["VITAMIN-J-175"], fisp_file)

    def test_to_fispact_fluxes(self, tmp_path: Path):
        with as_file(RES.joinpath("FLUXESS1")) as fisp_file:
            spectra = Spectra.from_fispact(GROUP_STRUCTURES["VITAMIN-J-175"], fisp_file)

        spectra.to_fispact_fluxes(tmp_path, format=".2E")

        # Check that the file is the same as the original one
        with as_file(RES.joinpath("FLUXESS1")) as f_original_ctxt:
            with open(tmp_path.joinpath("FLUX1"), "r") as f_new, open(
                f_original_ctxt, "r"
            ) as f_original:
                lines_new = f_new.readlines()
                lines_original = f_original.readlines()

                assert lines_new == lines_original

    def test_print_collapse_inp(self, spectra: Spectra, tmp_path: Path):
        spectra.print_collapse_inp(tmp_path)

    def test_print_condensed_inp(self, spectra: Spectra, tmp_path: Path):
        spectra.print_condensed_inp(tmp_path)

    def test_to_arb_flux_file(self, spectra: Spectra, tmp_path: Path):
        spectra.to_arb_flux_file(tmp_path)

        with open(tmp_path.joinpath("arb_flux_FLUX1.txt"), "r") as f:
            # get all tokens
            lines = f.readlines()
            tokens = []
            for line in lines:
                tokens.extend(line.split())
        # Check that the number of tokens is correct
        assert len(tokens) == 2 * len(spectra.spectra_values) + 3
        assert tokens[-1] == "FLUX1"

    def test_wrong_bins(self):
        with pytest.raises(ValueError):
            values = np.random.rand(len(GROUP_STRUCTURES["VITAMIN-J-175"]) - 1)
            Spectra(
                ebins=GROUP_STRUCTURES["VITAMIN-J-175"][:-1],
                spectra_values=values / np.sum(values),
                name="test_spectra",
            )

    def test_not_normalized(self):
        spectrum = Spectra(
            ebins=GROUP_STRUCTURES["VITAMIN-J-175"],
            spectra_values=np.random.rand(len(GROUP_STRUCTURES["VITAMIN-J-175"]) - 1)
            * 100,
            name="test_spectra",
        )
        assert pytest.approx(spectrum.spectra_values.sum()) == 1.0

    @pytest.mark.parametrize(
        "egroup",
        [
            GROUP_STRUCTURES["CCFE-709"],
            GROUP_STRUCTURES["CASMO-16"],
            GROUP_STRUCTURES["VITAMIN-J-175"][:-1],
            np.array(GROUP_STRUCTURES["VITAMIN-J-175"].tolist() + [20.0]),
        ],
    )
    def test_wrong_binning(self, egroup: np.ndarray):
        with pytest.raises(SpectraParsingError):
            with as_file(RES.joinpath("FLUXESS1")) as fisp_file:
                Spectra.from_fispact(egroup, fisp_file)

    def test_plot(self, spectra: Spectra):
        fig, ax = spectra.plot(lethargy=True, add_spectra=[spectra, spectra])
        assert fig is not None
        assert ax is not None

    def test_spaced_name(self):
        with as_file(RES.joinpath("FLUXESS_space_name")) as fisp_file:
            spectra = Spectra.from_fispact(GROUP_STRUCTURES["VITAMIN-J-175"], fisp_file)
        assert spectra.name == "FLUX1_with_space"

    def test_automatic_normalization(self):
        ebins = GROUP_STRUCTURES["VITAMIN-J-175"]
        values = np.random.rand(len(ebins) - 1) * 100
        spectra = Spectra(ebins=ebins, spectra_values=values, name="test_spectra")
        assert np.isclose(np.sum(spectra.spectra_values), 1.0, atol=5e-3)
