from pathlib import Path

import pytest
import pyvista

from f4enix.input.ww_gvr.weight_window import WW


def test_read_fmesh_and_export_cartesian_has_correct_parameters(tmpdir):
    ww = WW.create_gvr_from_meshtally_file(
        Path("tests") / "test_ww_gvr" / "resources" / "meshtal_complex_cart"
    )
    ww.export_as_vtk(tmpdir / "test.vts")
    mesh = pyvista.read(tmpdir / "test.vts")
    assert mesh is not None
    assert mesh.bounds == pytest.approx((-200.0, 200.0, -250.0, 300.0, -200.0, 100.0))


def test_write_multi_particle_ww_has_no_blank_line_between_particles(tmpdir):
    """Regression test: no blank line should appear between particle sections."""
    gvr = WW.create_gvr_from_meshtally_file(
        Path("tests") / "test_ww_gvr" / "resources" / "meshtal_cyl"
    )
    gvr.add_particle(norm=1.0, soft=1.0)
    ww_path = tmpdir / "test.ww"
    gvr.write_to_ww_file(ww_path)

    content = Path(ww_path).read_text()
    assert (
        "\n\n" not in content
    ), "Blank line found in WW file (likely between particle sections)"
