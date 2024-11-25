# flake8: noqa: PLR2004
import shutil
from pathlib import Path

import numpy as np
import pytest
import pyvista as pv
from numpy.testing import assert_array_almost_equal

from f4enix.input.ww_gvr.models import CoordinateType, ParticleType, Vectors
from f4enix.input.ww_gvr.weight_window import WW
from tests.test_ww_gvr.resources import expected_values_ww_complex_cart


def test_init_ww_from_ww_file_cart_simple():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    assert (
        ww.file_path == Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    assert ww.particles == [ParticleType.NEUTRON]
    assert ww.geometry.coordinate_type == CoordinateType.CARTESIAN
    assert ww.geometry.director_1 is None
    assert ww.geometry.director_2 is None
    assert ww.geometry._coarse_vectors == Vectors(
        vector_i=np.array([-15, 15]),
        vector_j=np.array([-15, 16]),
        vector_k=np.array([-15, 20]),
    )
    assert ww.geometry._fine_vectors == Vectors(
        vector_i=np.array([2]),
        vector_j=np.array([3]),
        vector_k=np.array([1]),
    )
    assert ww.energies[ParticleType.NEUTRON] == [100.0]

    expected_array = np.array(
        [[[0.11576, 0.093197], [0.67316, 0.5], [0.099821, 0.0898]]]
    )
    assert_array_almost_equal(ww.values[ParticleType.NEUTRON][100.0], expected_array)

    expected_ratios = np.array(
        [[[5.81513476, 5.36497956], [6.74367117, 5.56792873], [6.74367117, 5.56792873]]]
    )

    assert_array_almost_equal(ww.ratios[ParticleType.NEUTRON][100.0], expected_ratios)


def test_init_ww_from_ww_file_cyl_simple():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cyl"
    )
    assert ww.geometry.coordinate_type == CoordinateType.CYLINDRICAL
    assert ww.geometry.director_1 == [0.0, 0.0, 11.0]
    assert ww.geometry.director_2 == [15.0, 0.0, -5.0]
    expected_array = np.array(
        [[[0.5, 0.10463], [0.52965, 0.084479], [0.14258, 0.03275]]]
    )
    assert_array_almost_equal(ww.values[ParticleType.NEUTRON][100.0], expected_array)


def test_init_ww_from_ww_complex_cart():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_complex_cart"
    )

    assert ww.particles == [ParticleType.NEUTRON, ParticleType.PHOTON]

    expected_energies = {
        ParticleType.NEUTRON: [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 100.0],
        ParticleType.PHOTON: [1.2, 2.3],
    }
    assert ww.energies == expected_energies

    expected_values = expected_values_ww_complex_cart.expected_values
    for particle in ww.particles:
        for energy in expected_values[particle].keys():
            assert_array_almost_equal(
                expected_values[particle][energy], ww.values[particle][energy]
            )


def test_values_setter_recalculates_ratios():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )

    # Create a grid with 8 cells for convenience, "ww_simple_cart" only has 6 cells
    dummy_grid = pv.StructuredGrid()
    dummy_grid.dimensions = [3, 3, 3]
    ww.geometry._grid = dummy_grid

    ww.values = {
        ParticleType.NEUTRON: {100.0: np.array([[[1, 2], [1, 1]], [[1, 1], [1, 1]]])}
    }
    expected_ratios = np.array([[[2.0, 2.0], [1.0, 2.0]], [[1.0, 2.0], [1.0, 1.0]]])
    assert_array_almost_equal(ww.ratios[ParticleType.NEUTRON][100.0], expected_ratios)


def test_create_gvr_from_meshtally_file_cart():
    flat_values = [
        9.45568e-05,
        9.54871e-05,
        7.67912e-05,
        7.31554e-05,
        7.30042e-05,
        6.20899e-05,
        5.71612e-05,
        5.69892e-05,
        5.00181e-05,
        4.53375e-05,
        4.50316e-05,
        4.08733e-05,
        7.70699e-05,
        7.64715e-05,
        6.40383e-05,
        6.16851e-05,
        6.22624e-05,
        5.33923e-05,
        5.00923e-05,
        5.02796e-05,
        4.46194e-05,
        4.10111e-05,
        4.10669e-05,
        3.71414e-05,
    ]
    flat_values = np.array(flat_values)
    MAXIMUM_SPLITTING_RATIO = 5.0  # Default one
    MAX_FLUX = flat_values.max()
    expected_values = flat_values / MAX_FLUX * 2 / (MAXIMUM_SPLITTING_RATIO + 1)

    ww = WW.create_gvr_from_meshtally_file(
        Path("tests") / "test_ww_gvr" / "resources" / "meshtally_cart"
    )
    _nested_energies, nested_values = ww._format_nested_energies_and_values()
    result_flat_values = np.array(nested_values[0])

    assert_array_almost_equal(expected_values, result_flat_values)

    softened_ww = WW.create_gvr_from_meshtally_file(
        Path("tests") / "test_ww_gvr" / "resources" / "meshtally_cart",
        softening_factor=0.5,
    )

    assert np.allclose(
        ww.values[ParticleType.NEUTRON][100.0] ** 0.5,
        softened_ww.values[ParticleType.NEUTRON][100.0],
    )


def test_gvr_fails_if_multiple_energies():
    with pytest.raises(ValueError):
        WW.create_gvr_from_meshtally_file(
            Path("tests") / "resources" / "meshtal" / "meshtal_rect_VV"
        )


def test_info():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )

    with open(
        Path("tests") / "test_ww_gvr" / "resources" / "expected_info_simple_cart.txt"
    ) as infile:
        expected = infile.read()

    assert ww.info == expected


def test_repr():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    assert repr(ww) == ww.info


@pytest.mark.parametrize(
    "ww_file",
    [
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart",
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cyl",
        Path("tests") / "test_ww_gvr" / "resources" / "ww_complex_cart",
    ],
)
def test_write_to_ww_file(tmp_path, ww_file):
    original_ww = WW.load_from_ww_file(ww_file)
    original_ww.write_to_ww_file(tmp_path / "test.ww")

    result_ww = WW.load_from_ww_file(tmp_path / "test.ww")

    for particle in original_ww.particles:
        for energy in original_ww.energies[particle]:
            assert_array_almost_equal(
                original_ww.values[particle][energy], result_ww.values[particle][energy]
            )


def test_write_to_ww_file_no_path(tmp_path):
    shutil.copy(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart", tmp_path
    )
    tmp_ww_path = Path(tmp_path / "ww_simple_cart")

    original_ww = WW.load_from_ww_file(tmp_ww_path)
    original_ww.write_to_ww_file()

    result_ww = WW.load_from_ww_file(str(tmp_ww_path) + "_written")

    for particle in original_ww.particles:
        for energy in original_ww.energies[particle]:
            assert_array_almost_equal(
                original_ww.values[particle][energy], result_ww.values[particle][energy]
            )


@pytest.mark.parametrize(
    "ww_file",
    [
        Path("tests") / "test_ww_gvr" / "resources" / "meshtally_cart",
        Path("tests") / "test_ww_gvr" / "resources" / "meshtal_cyl",
    ],
)
def test_write_to_ww_file_gvr(tmp_path, ww_file):
    original_ww = WW.create_gvr_from_meshtally_file(ww_file)
    original_ww.write_to_ww_file(tmp_path / "test.ww")

    result_ww = WW.load_from_ww_file(tmp_path / "test.ww")

    for particle in original_ww.particles:
        for energy in original_ww.energies[particle]:
            assert_array_almost_equal(
                original_ww.values[particle][energy],
                result_ww.values[particle][energy],
                decimal=5,  # The writer only writes 5 significant digits to file
            )


@pytest.mark.parametrize(
    "ww_file",
    [
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart",
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cyl",
        Path("tests") / "test_ww_gvr" / "resources" / "ww_complex_cart",
    ],
)
def test_export_as_vtk(tmp_path, ww_file):
    ww = WW.load_from_ww_file(ww_file)
    ww.export_as_vtk(tmp_path / "test.vts")

    written_grid = pv.read(tmp_path / "test.vts")

    assert (written_grid.points == ww.geometry._grid.points).all()


def test_export_as_vtk_no_path(tmp_path):
    shutil.copy(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart", tmp_path
    )
    tmp_ww_path = Path(tmp_path / "ww_simple_cart")

    ww = WW.load_from_ww_file(tmp_ww_path)
    ww.export_as_vtk()

    written_grid = pv.read(str(tmp_ww_path) + ".vts")

    assert (written_grid.points == ww.geometry._grid.points).all()


def test_export_as_vtk_wrong_suffix(tmp_path):
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    ww.export_as_vtk(tmp_path / "test.wrong")

    written_grid = pv.read(tmp_path / "test.vts")

    assert (written_grid.points == ww.geometry._grid.points).all()


def test_vtk_export_from_gvr(tmp_path):
    ww = WW.create_gvr_from_meshtally_file(
        Path("tests") / "test_ww_gvr" / "resources" / "meshtally_cart"
    )
    ww.export_as_vtk(tmp_path / "test.vts")

    written_grid = pv.read(tmp_path / "test.vts")

    assert isinstance(written_grid, pv.StructuredGrid)
    assert written_grid.bounds[0] == ww.geometry.vectors.vector_i[0]
    assert written_grid.bounds[2] == ww.geometry.vectors.vector_j[0]
    assert written_grid.bounds[4] == ww.geometry.vectors.vector_k[0]


def test_soften():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    original_ratios = ww.ratios[ParticleType.NEUTRON][100.0].copy()

    ww.soften(0.5)

    expected_values = (
        np.array([[[0.11576, 0.093197], [0.67316, 0.5], [0.099821, 0.0898]]]) ** 0.5
    )

    assert_array_almost_equal(ww.values[ParticleType.NEUTRON][100.0], expected_values)

    new_ratios = ww.ratios[ParticleType.NEUTRON][100.0].copy()
    assert float(np.average(new_ratios)) < float(np.average(original_ratios))


def test_soften_quick_return():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    original_ratios = ww.ratios[ParticleType.NEUTRON][100.0].copy()

    # Now we modify directly the values and then do a soften(1.0), this should not
    # recalculate the ratios because the soften value is 1 and therefore the function
    # ends early
    ww._values[ParticleType.NEUTRON][100.0] = np.array(
        [[[1, 2], [43, 23]], [[12, 12], [331, 1]]]
    )
    ww.soften(1.0)
    ratios_after_soften = ww.ratios[ParticleType.NEUTRON][100.0].copy()

    assert_array_almost_equal(original_ratios, ratios_after_soften)


def test_multiply():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )

    expected_values = np.array([[[1.1576, 0.93197], [6.7316, 5], [0.99821, 0.898]]])
    ww.multiply(10.0)
    assert_array_almost_equal(ww.values[ParticleType.NEUTRON][100.0], expected_values)


def test_multiply_quick_return():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )

    expected_values = np.array(
        [[[0.11576, 0.093197], [0.67316, 0.5], [0.099821, 0.0898]]]
    )
    ww.multiply(1)
    assert_array_almost_equal(ww.values[ParticleType.NEUTRON][100.0], expected_values)


def test_add_particle_identical():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    ww.add_particle(norm=1, soft=1)

    assert len(ww.particles) == 2
    assert ParticleType.NEUTRON in ww.particles
    assert ParticleType.PHOTON in ww.particles

    for energy in ww.energies[ParticleType.NEUTRON]:
        assert_array_almost_equal(
            ww.values[ParticleType.NEUTRON][energy],
            ww.values[ParticleType.PHOTON][energy],
        )
        assert_array_almost_equal(
            ww.ratios[ParticleType.NEUTRON][energy],
            ww.ratios[ParticleType.PHOTON][energy],
        )


def test_add_particle_already_exists():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    ww.add_particle(norm=1, soft=1)
    with pytest.raises(ValueError):
        ww.add_particle(norm=1, soft=1)


def test_add_particle_different():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    ww.add_particle(norm=2, soft=0.5)

    for energy in ww.energies[ParticleType.NEUTRON]:
        assert_array_almost_equal(
            ww.values[ParticleType.NEUTRON][energy] * 2**0.5,
            ww.values[ParticleType.PHOTON][energy],
        )


def test_remove_particle():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_complex_cart"
    )
    ww.remove_particle()

    assert len(ww.particles) == 1
    assert ParticleType.NEUTRON in ww.particles
    assert ParticleType.PHOTON not in ww.particles

    assert len(ww.values) == 1


def test_remove_particle_only_one():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    with pytest.raises(ValueError):
        ww.remove_particle()


def test_mitigate_long_histories():
    ww = WW.load_from_ww_file(
        Path("tests") / "test_ww_gvr" / "resources" / "ww_simple_cart"
    )
    previous_max_ratio = np.max(ww.ratios[ParticleType.NEUTRON][100.0])

    ww.mitigate_long_histories(max_ratio=6.57)
    current_max_ratio = np.max(ww.ratios[ParticleType.NEUTRON][100.0])

    expected_values = np.array([[[0.11576, 0.093197], [0.0, 0.5], [0.0, 0.0898]]])

    assert_array_almost_equal(expected_values, ww.values[ParticleType.NEUTRON][100.0])
    assert current_max_ratio < previous_max_ratio
