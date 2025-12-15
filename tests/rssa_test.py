from importlib.resources import files
from pathlib import Path

import pytest

import tests.resources.rssa as res
from f4enix.output.rssa import RSSA, PlotParameters

RESOURCES = files(res)


@pytest.fixture
def rssa():
    path = Path(RESOURCES.joinpath("small_cyl.w"))  # type: ignore
    return RSSA(path)


def test_read_rssa_parameters(rssa):
    assert isinstance(rssa, RSSA)
    assert rssa.parameters.np1 == -100_000
    assert rssa.parameters.nrss == 72083
    assert rssa.parameters.nrcd == 11
    assert rssa.parameters.njsw == 1
    assert rssa.parameters.niss == 70797
    assert rssa.parameters.niwr == 0
    assert rssa.parameters.mipts == 3
    assert rssa.parameters.kjaq == 0

    assert rssa.parameters.surfaces[0].id == 1
    assert rssa.parameters.surfaces[0].info == -1
    assert rssa.parameters.surfaces[0].type == 1
    assert rssa.parameters.surfaces[0].num_parameters == 0


def test_read_rssa_tracks(rssa):
    assert rssa.tracks.shape == (72083, 11)


def test_neutron_tracks(rssa):
    neutron_tracks = rssa.neutron_tracks
    assert neutron_tracks.shape[0] > 0


def test_photon_tracks(rssa):
    photon_tracks = rssa.photon_tracks
    assert photon_tracks.shape[0] == 0


def test_properties(rssa):
    rssa.x
    rssa.y
    rssa.z
    rssa.energies
    rssa.histories
    rssa.wgt
    rssa.get_summary()
    str(rssa)
    rssa.__repr__()
    assert True


def test_plot_cyl(rssa, tmp_path):
    (
        rssa.plot_cyl()
        .set_particle("n")
        .set_z_limits(-600, 800)
        .set_perimeter_limits(-500, 1000)
        .calculate_bins(bin_width=10)
        .get_particle_current(1e20)
        .set_plot_parameters(
            PlotParameters(
                vmin=1e6,
                vmax=1e14,
                number_of_colors=16,
                legend_orientation="vertical",
                title="Neutron current through the surface [#/cm2/s]",
            )
        )
        .save_figure(tmp_path / "test_plot_cyl.png")
    )
    # Check if the plot file was created
    assert (tmp_path / "test_plot_cyl.png").exists()

    (rssa.plot_cyl().set_particle("p").set_surface_ids([]))


def test_plot_cyl_errors(rssa, tmp_path):
    (
        rssa.plot_cyl()
        .set_particle("n")
        .set_z_limits(-600, 800)
        .set_perimeter_limits(-500, 1000)
        .calculate_bins(bin_width=10)
        .get_particle_current_errors()
        .set_plot_parameters(
            PlotParameters(
                vmin=1e6,
                vmax=1e14,
                number_of_colors=16,
                legend_orientation="vertical",
                title="Neutron current error through the surface [#/cm2/s]",
            )
        )
        .save_figure(tmp_path / "test_plot_cyl_error.png")
    )
    # Check if the plot file was created
    assert (tmp_path / "test_plot_cyl_error.png").exists()


def test_get_ratio(rssa, tmp_path):
    (
        rssa.plot_cyl()
        .set_z_limits(-600, 800)
        .set_perimeter_limits(-500, 1000)
        .calculate_bins(bin_width=10)
        .get_particle_current(1)
        .get_ratio_to(
            rssa.plot_cyl()
            .set_z_limits(-600, 800)
            .set_perimeter_limits(-500, 1000)
            .calculate_bins(bin_width=10)
            .get_particle_current(2)
        )
        .set_plot_parameters(
            PlotParameters(
                vmin=0.1,
                vmax=10,
                number_of_colors=16,
                legend_orientation="vertical",
                title="Neutron to proton current ratio through the surface",
            )
        )
        .save_figure(tmp_path / "test_get_ratio.png")
    )
    # Check if the plot file was created
    assert (tmp_path / "test_get_ratio.png").exists()


def test_plot_plane(rssa, tmp_path):
    (
        rssa.plot_plane()
        .set_particle("n")
        .set_z_limits(-600, 800)
        .calculate_bins(bin_width=10)
        .get_particle_current(1e20)
        .set_plot_parameters(
            PlotParameters(
                vmin=1e6,
                vmax=1e14,
                number_of_colors=16,
                legend_orientation="vertical",
                title="Neutron current through the surface [#/cm2/s]",
            )
        )
        .save_figure(tmp_path / "test_plot_plane.png")
    )
    # Check if the plot file was created
    assert (tmp_path / "test_plot_plane.png").exists()
