import os
import re
import pytest
import pyvista as pv
from importlib.resources import files, as_file

from f4enix.output.meshtal import Meshtal
from f4enix.output.meshtal import identical_mesh
import tests.resources.meshtal as resources
import tests.resources.meshtal.tests as res
import tests.resources.meshtal.expected as res_exp

resources_write = files(res)
expected = files(res_exp)
RESOURCES = files(resources)


class TestMeshtal:
    def test_thetest(self):
        assert True

    @pytest.mark.parametrize(
        "input_meshtal",
        [
            "meshtal_cuv",
            "meshtal_cyl",
            "meshtal_d1s_CSimpactStudy",
            "meshtal_CUBE_SQUARE",
            "meshtal_CUBE_ONES",
            "test_srcimp",
            "assembly_meshtal_test",
        ],
    )
    def test_mesh_print_tally_info(self, input_meshtal):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = Meshtal(inp, filetype)

        for i in meshtally.mesh.items():
            meshtally.mesh[i[0]].print_info()

        assert True

    def test_same_mesh(self):
        with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        assert meshtally.mesh[124].sameMesh(meshtally.mesh[124], checkErg=True)
        assert meshtally.mesh[124].sameMesh(meshtally.mesh[124])

    @pytest.mark.parametrize(
        "input_meshtal",
        [
            "meshtal_cuv",
            "meshtal_cyl",
            "meshtal_d1s_CSimpactStudy",
            "meshtal_CUBE_SQUARE",
            "meshtal_CUBE_ONES",
            "test_srcimp",
            "assembly_meshtal_test",
            "meshtal_bug",
        ],
    )
    def test_read_mesh(self, input_meshtal):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = Meshtal(inp, filetype)

        for i in meshtally.mesh.items():
            meshtally.readMesh()

        assert True

    @pytest.mark.parametrize(
        "norm",
        [
            "vtot",
            "celf",
        ],
    )
    def test_read_mesh_cuv(self, norm):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath("meshtal_cuv")) as inp:
            meshtally = Meshtal(inp, filetype)

        for i in meshtally.mesh.items():
            meshtally.readMesh(norm=norm)
            meshtally.readMesh(cell_filters=[1, 2], norm=norm)

        assert True

    @pytest.mark.parametrize(
        "input_meshtal", ["meshtal_cuv", "meshtal_cyl", "meshtal_d1s_CSimpactStudy"]
    )
    def test_mesh_print_info(self, input_meshtal):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = Meshtal(inp, filetype)

        for i in meshtally.mesh.items():
            meshtally.print_info()

        assert True

    def test_write_cyl(self, tmpdir):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        outpath = tmpdir.mkdir("sub_cyl")
        meshtally.mesh[124].write(outpath)

    @pytest.mark.parametrize(
        ["file_read", "list_array_names", "out_format", "out_name", "file_exp"],
        [
            [
                "example.vts",
                None,
                "csv",
                "meshtal_cyl_124_csv.csv",
                "example_['Values']_csv.csv",
            ],
            [
                "test_VTK_CUBE_SQUARE.vtr",
                ["Value - Total"],
                "csv",
                "meshtal_cyl_124_csv.csv",
                "test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv",
            ],
            [
                "meshtal_14.vts",
                ["Value - Total"],
                "csv",
                "meshtal_cyl_124_csv.csv",
                "meshtal_14_['Value - Total']_csv.csv",
            ],
            [
                "cuvmsh_44_CuV_CELF10.vtr",
                ["Value - Total"],
                "csv",
                "meshtal_cyl_124_csv.csv",
                "cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "csv",
                "meshtal_cyl_124_csv.csv",
                "PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "ip_fluent",
                "meshtal_cyl_124_ip_fluent.txt",
                "PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "point_cloud",
                "meshtal_cyl_124_point_cloud.txt",
                "PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_point_cloud.txt",
            ],
        ],
    )
    def test_write(
        self, file_read, list_array_names, out_format, out_name, file_exp, tmpdir
    ):
        # Whatever mesh can be read here, the grid will be overridden
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        fmesh = meshtally.mesh[124]

        with as_file(resources_write.joinpath(file_read)) as inp:
            fmesh._read_from_vtk(inp)

        outpath = tmpdir.mkdir("sub_csv")
        fmesh.write(outpath, out_format=out_format, list_array_names=list_array_names)
        outfile = os.path.join(outpath, out_name)

        with as_file(expected.joinpath(file_exp)) as exp:
            with open(outfile, "r") as test, open(exp, "r") as exp_file:
                for line1, line2 in zip(test, exp_file):
                    assert line1 == line2

        # Also always test the .vtk writing
        fmesh.write(outpath)

    @pytest.mark.parametrize(
        "input_meshtal",
        [
            "meshtal_cuv",
            "meshtal_cyl",
            "meshtal_d1s_CSimpactStudy",
            "meshtal_d1s_IVVS_FDR",
            "meshtal_rect_VV",
            "test_srcimp",
            "assembly_meshtal_test",
        ],
    )
    def test_reading(self, input_meshtal):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            Meshtal(inp, filetype)

        assert True

    def test_write_all(self, tmpdir):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        meshtally.write_all(tmpdir)

    @pytest.mark.parametrize(
        "input_meshtal",
        [
            "meshtal_cuv",
            "meshtal_cyl",
            "meshtal_d1s_CSimpactStudy",
            "meshtal_CUBE_SQUARE",
            "meshtal_CUBE_ONES",
            "test_srcimp",
            "assembly_meshtal_test",
        ],
    )
    def test_identical_mesh(self, input_meshtal):
        # To check if the meshtal can be read without any problem"
        filetype = "MCNP"
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = Meshtal(inp, filetype)

        for m, i in enumerate(list(meshtally.mesh.values())[1:]):
            if m == 0:
                j = list(meshtally.mesh.values())[0]
            a, b, c = identical_mesh(i, j)
            j = i

        assert True

    def test_collapse_grid(self):
        with as_file(RESOURCES.joinpath("meshtal_collapse")) as inp:
            meshtal = Meshtal(inp)

        meshtal.readMesh()
        dict_names = {4: ["A", "B"], 14: ["C", "D"]}
        grid = meshtal.collapse_grids(dict_names)
        assert len(grid.array_names) == 4
