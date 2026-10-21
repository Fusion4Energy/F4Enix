import os

import numpy as np
from importlib.resources import as_file, files

from f4enix.core.constants import TIME_UNITS
import pytest

import tests.resources.meshtal as resources
import tests.resources.meshtal.expected as res_exp
import tests.resources.meshtal.tests as res
import tests.resources.irradiation as res_irr

from f4enix.input.MCNPinput import Input
from f4enix.output.meshtal.aux_meshtal_functions import _to_float
from f4enix.output.meshtal.meshtal import Meshtal
from f4enix.input.d1suned import IrradiationFile

resources_write = files(res)
expected = files(res_exp)
RESOURCES = files(resources)
resources_irr = files(res_irr)


class TestMeshtal:
    @pytest.mark.parametrize(
        ["input_meshtal", "filetype"],
        [
            ("meshtal_cuv", "CUV"),
            ("meshtal_cyl", "MCNP"),
            ("meshtal_d1s_CSimpactStudy", "MCNP"),
            ("meshtal_CUBE_SQUARE", "MCNP"),
            ("meshtal_CUBE_ONES", "MCNP"),
            # ("test_srcimp", "CUV"),
            ("assembly_meshtal_test", "MCNP"),
            ("1D_mesh_no_code", "MCNP"),
            ("cdgs_test", "CDGS"),
        ],
    )
    def test_mesh_print_tally_info(self, input_meshtal, filetype):
        """To check if the meshtal can be read without any problem"""
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = Meshtal(filename=inp, filetype=filetype)
            meshtally.readMesh(norm="ctot")

        for _mesh_id, mesh in meshtally.mesh.items():
            mesh.print_info()

    def test_same_mesh(self):
        with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        assert meshtally.mesh[124].sameMesh(meshtally.mesh[124])

    @pytest.mark.parametrize(
        "norm",
        [
            "ctot",
            "celf",
        ],
    )
    def test_read_mesh_cuv(self, norm):
        # To check if the meshtal can be read without any problem"
        filetype = "CUV"
        with as_file(RESOURCES.joinpath("meshtal_cuv")) as inp:
            meshtally = Meshtal(inp, filetype)

        meshtally.readMesh(norm=norm)
        meshtally.readMesh(cell_filters=[1, 2], norm=norm)
        assert meshtally._meshtal_parser.get_meshlist() == ("44",)
        bound = meshtally.mesh[44].grid.bounds
        expected = 5.0
        assert bound[-1] == expected
        meshtally.readMesh(cell_filters=[1], norm=norm)

    @pytest.mark.parametrize(
        ["input_meshtal", "length"], [("1D_mesh", 79), ("1D_mesh_energyonly", 175)]
    )
    def test_1d_features(self, input_meshtal, length):
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh()
        _, data, _ = meshtal.mesh[214].convert2tally()
        assert len(data) == length
        pass

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
                "Tally_124_csv.csv",
                "example_['Values']_csv.csv",
            ],
            [
                "test_VTK_CUBE_SQUARE.vtr",
                ["Value - Total"],
                "csv",
                "Tally_124_csv.csv",
                "test_VTK_CUBE_SQUARE_['Value - Total']_csv.csv",
            ],
            [
                "meshtal_14.vts",
                ["Value - Total"],
                "csv",
                "Tally_124_csv.csv",
                "meshtal_14_['Value - Total']_csv.csv",
            ],
            [
                "cuvmsh_44_CuV_CELF10.vtr",
                ["Value - Total"],
                "csv",
                "Tally_124_csv.csv",
                "cuvmsh_44_CuV_CELF10_['Value - Total']_csv.csv",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "csv",
                "Tally_124_csv.csv",
                "PS_NHD_DIV_RHC_INBOARD_['NHD[W-cm3]']_csv.csv",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "ip_fluent",
                "Tally_124_ip_fluent_NHDWcm3.txt",
                "PS_NHD_DIV_RHC_INBOARD_NHD[W-cm3]_ip_fluent.txt",
            ],
            [
                "PS_NHD_DIV_RHC_INBOARD.vtk",
                ["NHD[W/cm3]"],
                "point_cloud",
                "Tally_124_point_cloud_NHDWcm3.txt",
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
            with open(outfile) as test, open(exp) as exp_file:
                for line1, line2 in zip(test, exp_file):
                    assert line1 == line2

        # Also always test the .vtk writing
        fmesh.write(outpath)

    def test_write_outfile(self, tmpdir):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh()
        fmesh = meshtally.mesh[124]
        outpath = tmpdir.mkdir("sub_csv")
        fmesh.write(outpath, outfile="test")
        assert os.path.exists(os.path.join(outpath, "test.vtr"))

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

    def test_read_list(self, tmpdir):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            meshtally = Meshtal(inp)
        meshtally.readMesh([4, 14])
        meshtally.write_all(tmpdir)

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
            assert i.sameMesh(j)

    def test_collapse_grid(self):
        with as_file(RESOURCES.joinpath("meshtal_collapse")) as inp:
            meshtal = Meshtal(inp)

        meshtal.readMesh()
        dict_names = {4: ["A", "B"], 14: ["C", "D"]}
        grid = meshtal.collapse_grids(dict_names)
        expected = 4
        assert len(grid.array_names) == expected

    def test_transform_grid(self, tmpdir):
        with as_file(RESOURCES.joinpath("meshtal_transform")) as inp:
            meshtal = Meshtal(inp)
        with as_file(RESOURCES.joinpath("transforms.i")) as mcnp_inp:
            input_file = Input.from_input(mcnp_inp)

        meshtal.readMesh()
        meshtal.transform_fmesh(input_file)
        meshtal.write_all(tmpdir)
        assert tuple(meshtal.mesh[2024].grid.bounds) == (
            567.0,
            573.0,
            47.0,
            53.0,
            387.0,
            393.0,
        )
        assert tuple(meshtal.mesh[2024].grid.center) == (570.0, 50.0, 390.0)
        assert tuple(meshtal.mesh[2024].grid.bounds) == (
            567.0,
            573.0,
            47.0,
            53.0,
            387.0,
            393.0,
        )
        assert tuple(meshtal.mesh[2124].grid.center) == (0.0, 0.0, 0.0)
        expected = -3.845138
        assert pytest.approx(meshtal.mesh[2124].grid.bounds[0], 1e-5) == expected
        assert tuple(meshtal.mesh[2224].grid.center) == (1.0, 1.0, 1.0)
        assert tuple(meshtal.mesh[2224].grid.bounds) == (
            -2.0,
            4.0,
            -2.0,
            4.0,
            -2.0,
            4.0,
        )
        assert tuple(meshtal.mesh[2324].grid.center) == (10.0, 10.0, 10.0)
        expected = 6.52463
        assert pytest.approx(meshtal.mesh[2324].grid.bounds[0], 1e-5) == expected
        meshtal.write_all(tmpdir)

    def test_transform_grid_rotation_direction(self):
        # the mesh is defined in the auxiliary system of TR4 and must be
        # moved to the main one: main = R @ aux + o, where the columns of R
        # are the auxiliary axes in main coordinates (MCNP convention)
        with as_file(RESOURCES.joinpath("meshtal_transform")) as inp:
            meshtal = Meshtal(inp)
        with as_file(RESOURCES.joinpath("transforms.i")) as mcnp_inp:
            input_file = Input.from_input(mcnp_inp)
        meshtal.readMesh()
        aux_points = meshtal.mesh[2324].grid.cast_to_structured_grid().points

        meshtal.transform_fmesh(input_file)

        # *TR4  10 10 10  10 100 90  80 10 90  90 90 0
        cos = lambda deg: np.cos(np.radians(deg))  # noqa: E731
        x_aux = [cos(10), cos(100), cos(90)]
        y_aux = [cos(80), cos(10), cos(90)]
        z_aux = [cos(90), cos(90), cos(0)]
        rotation = np.array([x_aux, y_aux, z_aux]).T
        expected = aux_points @ rotation.T + np.array([10, 10, 10])
        assert np.allclose(meshtal.mesh[2324].grid.points, expected)

    def test_time_energy_bins(self):
        # To check if the meshtal can be read without any problem"
        with as_file(RESOURCES.joinpath("meshtal_time_energy_bins")) as inp:
            meshtal = Meshtal(inp)

        meshtal.readMesh()
        assert (
            pytest.approx(meshtal.mesh[54].grid.cell_data["Value - Total_t001_e001"][0])
            == 1.0
        )
        filtered_mesh = meshtal.create_filtered_mesh(
            54,
            binlabels=("Value - Total", "Error - Total"),  # type: ignore
            ebin=0,
            tbin=0,
        )
        assert pytest.approx(filtered_mesh.grid.cell_data["Value - Total"][0]) == 1.0

        filtered_mesh = meshtal.create_filtered_mesh(
            64,
            binlabels=("Value - Total", "Error - Total"),  # type: ignore
            ebin=None,
            tbin=0,
        )
        assert (
            pytest.approx(filtered_mesh.grid.cell_data["Value - Total_e001"][0]) == 1.0
        )

        filtered_mesh = meshtal.create_filtered_mesh(
            64,
            binlabels=("Value - Total", "Error - Total"),  # type: ignore
            ebin=0,
            tbin=None,
        )
        assert (
            pytest.approx(filtered_mesh.grid.cell_data["Value - Total_t001"][0]) == 1.0
        )

        filtered_mesh = meshtal.create_filtered_mesh(
            64,
            binlabels=("Value - Total", "Error - Total"),  # type: ignore
            ebin=None,
            tbin=None,
        )
        assert (
            pytest.approx(filtered_mesh.grid.cell_data["Value - Total_t001_e001"][0])
            == 1.0
        )

    def test_cdgs(self, tmpdir):
        # To check if the meshtal can be read without any problem"
        with as_file(RESOURCES.joinpath("cdgs_test")) as inp:
            meshtal = Meshtal(inp, filetype="CDGS")

        meshtal.readMesh()
        meshtal.mesh[1].write(tmpdir)

    @pytest.mark.parametrize("out", ["col", "cf"])
    def test_three_digit_exponent(self, tmp_path, out):
        # Fortran writes exponents with three digits without the "E"
        text = RESOURCES.joinpath("meshtal_CUBE_SQUARE").read_text()
        old = "4           0.00000E+00"
        assert text.count(old) == 1
        text = text.replace(old, "1.77613-306 2.00000-101")
        if out == "col":
            # drop the Volume and Rslt * Vol columns
            head, data = text.split("Rslt * Vol\n")
            head = head.replace("     Volume    ", "")
            data = "".join(
                " ".join(row.split()[:5]) + "\n" for row in data.splitlines()
            )
            text = head + "\n" + data
        modified = tmp_path / "meshtal"
        modified.write_text(text)

        meshtal = Meshtal(modified)
        meshtal.readMesh()
        grid = meshtal.mesh[list(meshtal.mesh)[0]].grid
        # the grid is stored in single precision, so these tiny values become 0
        assert list(grid.cell_data["Value - Total"]) == [1, 9, 0, 16]
        assert list(grid.cell_data["Error - Total"]) == [0, 0, 0, 0]

    def test_rescale_dose_mesh(self):
        with as_file(RESOURCES.joinpath("meshtal_time_energy_bins")) as inp:
            meshtal = Meshtal(inp)

        meshtal.readMesh()

        irr_file = IrradiationFile.from_text(resources_irr.joinpath("d1stime_o_30d"))
        scaling_factors = irr_file.get_scaling_factors_cooling_time(
            1, (200, TIME_UNITS.DAY)
        )
        rescaled_mesh = meshtal.mesh[44].rescale_dose_map(
            scaling_factors_df=scaling_factors,
            cooling_time_col="200d",
        )
        assert (
            pytest.approx(rescaled_mesh.cell_data["Value - Total_e002"][0])
            == 0.04435259766561963
        )


@pytest.mark.parametrize(
    ["text", "value"],
    [
        ("1.0E-03", 1e-3),
        ("2", 2.0),
        ("1.77613-306", 1.77613e-306),
        ("-3.5+100", -3.5e100),
        ("4.-120", 4e-120),
    ],
)
def test_to_float(text, value):
    assert _to_float(text) == pytest.approx(value)
