import os
from importlib.resources import files, as_file

from f4enix.output.pyvistawrap import PyVistaWrapper
import tests.resources.pyvistawrapper.tests as res
import tests.resources.pyvistawrapper.expected as res_exp

resources = files(res)
expected = files(res_exp)


class TestPyVistaWrapper:
    def test_example_csv(self, tmpdir):
        with as_file(resources.joinpath('example.vts')) as inp:
            mesh = PyVistaWrapper.from_file(inp)

        outpath = tmpdir.mkdir('sub_csv')
        mesh.write_mesh(outpath, out_format='csv')
        outfile = os.path.join(outpath, 'example_csv.csv')

        with as_file(expected.joinpath("example_['Values']_csv.csv")) as exp:
            with open(outfile, "r") as test, open(exp, 'r') as exp_file:
                for line1, line2 in zip(test, exp_file):
                    assert line1 == line2

    # def test_vtk_cube_square_csv(self):
    #     functions.open_mesh("tests/data/test_VTK_CUBE_SQUARE.vtr")
    #     functions.write_mesh("tests/data/test_VTK_CUBE_SQUARE.vtr", ["Value - Total"], "csv")
    #     with open(
    #         "tests/data/expected_results/test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv")
    #     self.assertEqual(expected, result)

    # def test_meshtal_14_csv(self):
    #     functions.open_mesh("tests/data/meshtal_14.vts")
    #     functions.write_mesh("tests/data/meshtal_14.vts", ["Value - Total"], "csv")
    #     with open(
    #         "tests/data/expected_results/meshtal_14_['Value - Total']_csv.csv", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/meshtal_14_['Value - Total']_csv.csv", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/meshtal_14_['Value - Total']_csv.csv")
    #     self.assertEqual(expected, result)

    # def test_cuvmsh_44_celf10_csv(self):
    #     functions.open_mesh("tests/data/cuvmsh_44_CuV_CELF10.vtr")
    #     functions.write_mesh("tests/data/cuvmsh_44_CuV_CELF10.vtr", ["Value - Total"], "csv")
    #     with open(
    #         "tests/data/expected_results/cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv")
    #     self.assertEqual(expected, result)

    # def test_rhc_inboard_csv(self):
    #     functions.open_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk")
    #     functions.write_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk", ["NHD[W/cm3]"], "csv")
    #     with open(
    #         "tests/data/expected_results/PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv")
    #     self.assertEqual(expected, result)

    # def test_rhc_inboard_ip_fluent(self):
    #     functions.open_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk")
    #     functions.write_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk", ["NHD[W/cm3]"], "ip_fluent")
    #     with open(
    #         "tests/data/expected_results/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt")
    #     self.assertEqual(expected, result)

    # def test_rhc_inboard_point_cloud(self):
    #     functions.open_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk")
    #     functions.write_mesh("tests/data/PS_NHD_DIV_RHC_INBOARD.vtk", ["NHD[W/cm3]"], "point_cloud")
    #     with open(
    #         "tests/data/expected_results/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_point_cloud.txt", "r"
    #     ) as infile:
    #         expected = infile.read()
    #     with open("tests/data/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_point_cloud.txt", "r") as infile:
    #         result = infile.read()
    #     os.remove("tests/data/PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_point_cloud.txt")
    #     self.assertEqual(expected, result)