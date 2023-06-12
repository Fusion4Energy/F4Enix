from importlib.resources import files, as_file
import pytest
import pyvista as pv

from f4enix.output.plotter import Atlas
import tests.resources.plotter as pkg_res

RESOURCES = files(pkg_res)


class TestAtlas:

    @pytest.fixture
    def mesh(self):
        with as_file(RESOURCES.joinpath('test1.vtk')) as infile:
            mesh = pv.read(infile)
        return mesh

    @pytest.fixture
    def stl(self):
        with as_file(RESOURCES.joinpath('stl.stl')) as infile:
            stl = pv.read(infile)
        return stl

    def test_slice_toroidal(self, mesh, stl):
        atlas = Atlas(mesh, stl)
        mesh_slices, stl_slices = atlas.slice_toroidal(20)
        assert True

    def test_slice_on_axis(self, mesh, stl):
        atlas = Atlas(mesh, stl)
        mesh_slices, stl_slices = atlas.slice_on_axis('y', 3)
        assert True

        # p = pv.Plotter()
        # for mesh_slice, stl_slice in (mesh_slices, stl_slices):
        #     p.add_mesh(mesh_slice, scalars='Error')
        #     p.add_mesh(stl_slice)

        # p.show(jupyter_backend='static')