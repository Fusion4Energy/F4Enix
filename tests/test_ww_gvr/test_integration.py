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
