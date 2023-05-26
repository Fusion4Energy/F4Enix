import os
from importlib.resources import files, as_file
import pytest
import numpy as np

from f4enix.output.pyvistawrap import PyVistaWrapper
import tests.resources.pyvistawrapper.tests as res
import tests.resources.pyvistawrapper.expected as res_exp

resources = files(res)
expected = files(res_exp)


class TestPyVistaWrapper:

    @pytest.mark.parametrize(
            ['file_read', 'list_array_names', 'out_format', 'out_name', 'file_exp'],
            [['example.vts', None, 'csv', 'example_csv.csv', "example_['Values']_csv.csv"],
             ['test_VTK_CUBE_SQUARE.vtr', ["Value - Total"], 'csv', 'test_VTK_CUBE_SQUARE_csv.csv', "test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv"],
             ['meshtal_14.vts', ["Value - Total"], 'csv', 'meshtal_14_csv.csv', "meshtal_14_['Value - Total']_csv.csv"],
             ['cuvmsh_44_CuV_CELF10.vtr', ["Value - Total"], 'csv', 'cuvmsh_44_CuV_CELF10_csv.csv', "cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv"],
             ['PS_NHD_DIV_RHC_INBOARD.vtk', ["NHD[W/cm3]"], 'csv', 'PS_NHD_DIV_RHC_INBOARD_csv.csv', "PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv"],
             ['PS_NHD_DIV_RHC_INBOARD.vtk', ["NHD[W/cm3]"], 'ip_fluent', 'PS_NHD_DIV_RHC_INBOARD_ip_fluent.txt', "PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt"],
             ['PS_NHD_DIV_RHC_INBOARD.vtk', ["NHD[W/cm3]"], 'point_cloud', 'PS_NHD_DIV_RHC_INBOARD_point_cloud.txt', "PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_point_cloud.txt"]
            ]
            )
    def test_write(self, file_read, list_array_names, out_format, out_name,
                   file_exp, tmpdir):
        with as_file(resources.joinpath(file_read)) as inp:
            mesh = PyVistaWrapper.from_file(inp)

        outpath = tmpdir.mkdir('sub_csv')
        mesh.write_mesh(outpath, out_format=out_format,
                        list_array_names=list_array_names)
        outfile = os.path.join(outpath, out_name)

        with as_file(expected.joinpath(file_exp)) as exp:
            with open(outfile, "r") as test, open(exp, 'r') as exp_file:
                for line1, line2 in zip(test, exp_file):
                    assert line1 == line2

        # Also always test the .vtk writing
        mesh.write_mesh(outpath)

    def test_print(self):
        with as_file(resources.joinpath('example.vts')) as inp:
            mesh = PyVistaWrapper.from_file(inp)
        print(mesh)
        print(mesh.__repr__())
        assert True

    @pytest.mark.parametrize(['file', 'array'],
                             [['PS_NHD_DIV_RHC_INBOARD.vtk', 'NHD[W/cm3]'],
                              ['cuvmsh_44_CuV_CELF10.vtr', 'Value - Total']])
    def test_print_array_info(self, file, array):
        with as_file(resources.joinpath(file)) as inp:
            mesh = PyVistaWrapper.from_file(inp)
        mesh.print_array_info(array)
        assert True

    def test_translate(self):
        with as_file(resources.joinpath('PS_NHD_DIV_RHC_INBOARD.vtk')) as inp:
            mesh1 = PyVistaWrapper.from_file(inp)
        mesh1.translate(10, 15, 50)
        with as_file(expected.joinpath('PS_NHD_translated.vtu')) as exp:
            mesh2 = PyVistaWrapper.from_file(exp)

        self._assert_equal_mesh(mesh1, mesh2)

    def test_rotate(self):
        with as_file(resources.joinpath('PS_NHD_DIV_RHC_INBOARD.vtk')) as inp:
            mesh1 = PyVistaWrapper.from_file(inp)
        mesh1.rotate(10, 20, 30)
        with as_file(expected.joinpath('PS_NHD_rotated.vtu')) as exp:
            mesh2 = PyVistaWrapper.from_file(exp)

        self._assert_equal_mesh(mesh1, mesh2)

    def test_scale(self):
        with as_file(resources.joinpath('PS_NHD_DIV_RHC_INBOARD.vtk')) as inp:
            mesh1 = PyVistaWrapper.from_file(inp)
        mesh1.scale([10, 20, 10])
        with as_file(expected.joinpath('PS_NHD_scaled.vtu')) as exp:
            mesh2 = PyVistaWrapper.from_file(exp)

        self._assert_equal_mesh(mesh1, mesh2)

    def test_merge(self):
        with as_file(resources.joinpath('test_VTK_CUBE_SQUARE.vtr')) as inp:
            mesh1 = PyVistaWrapper.from_file(inp)
        with as_file(resources.joinpath('test_VTK_CUBE_SQUARE.vtr')) as inp:
            mesh2 = PyVistaWrapper.from_file(inp)

        mesh1.merge(mesh2)

        assert True

    def test_scalar_multiply(self):
        with as_file(resources.joinpath('PS_NHD_DIV_RHC_INBOARD.vtk')) as inp:
            mesh = PyVistaWrapper.from_file(inp)

        old_max = mesh.mesh['NHD[W/cm3]'].max()
        mesh.scalar_multiply(100, 'NHD[W/cm3]')
        assert mesh.mesh['NHD[W/cm3]'].max() == 100*old_max

    def test_linear_combination(self):
        with as_file(resources.joinpath('PS_NHD_DIV_RHC_INBOARD.vtk')) as inp:
            mesh = PyVistaWrapper.from_file(inp)
        mesh.linear_combination('new', ['NHD[W/cm3]', 'NHD[W/cm3]'], [1, 2])
        assert mesh.mesh['new'].max() == mesh.mesh['NHD[W/cm3]'].max()*3

    def _assert_equal_mesh(self, mesh1: PyVistaWrapper, mesh2: PyVistaWrapper):
        # assert mesh1.mesh.origin == mesh2.mesh.origin
        # assert mesh1.mesh.dimensions == mesh2.mesh.dimensions
        # assert mesh1.mesh.spacing == mesh2.mesh.spacing
        assert np.allclose(mesh1.centers, mesh2.centers, rtol=1e-5)
