import os
import pickle
from copy import copy, deepcopy
from importlib.resources import as_file, files

import aci_tallies
import numpy as np
import pytest
import pyvista as pv

import tests.resources.irradiation as res_irr
import tests.resources.meshtal as resources
import tests.resources.meshtal.expected as res_exp
import tests.resources.meshtal.tests as res
from f4enix.core.constants import TIME_UNITS
from f4enix.input.d1suned import IrradiationFile
from f4enix.input.MCNPinput import Input
from f4enix.input.ww_gvr.weight_window import WW
from f4enix.output.meshtal.fmesh import Bins, Fmesh
from f4enix.output.meshtal.meshtal import Meshtal

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
        assert meshtally.mesh_ids == [44]
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
                "example.vts",
                ["Values"],
                "ip_fluent",
                "Tally_124_ip_fluent_Values.txt",
                "example_Values_ip_fluent.txt",
            ],
            [
                "example.vts",
                ["Values"],
                "point_cloud",
                "Tally_124_point_cloud_Values.txt",
                "example_Values_point_cloud.txt",
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
            fmesh.grid = pv.read(inp)

        outpath = tmpdir.mkdir("sub_csv")
        fmesh.write(outpath, out_format=out_format, list_array_names=list_array_names)
        outfile = os.path.join(outpath, out_name)

        with as_file(expected.joinpath(file_exp)) as exp:
            with open(outfile) as test, open(exp) as exp_file:
                assert test.readlines() == exp_file.readlines()

        # Also always test the .vtk writing
        fmesh.write(outpath)

    @pytest.mark.parametrize("out_format", ["csv", "point_cloud", "ip_fluent"])
    def test_write_small_values(self, tmp_path, out_format):
        with as_file(RESOURCES.joinpath("cdgs_test")) as inp:
            meshtal = Meshtal(inp, filetype="CDGS")
        meshtal.readMesh()
        fmesh = meshtal.mesh[1]
        fmesh.grid["Value - Total_eTot"][0] = 1.23456e-20
        fmesh.write(tmp_path, ["Value - Total_eTot"], out_format=out_format)
        (file,) = tmp_path.iterdir()
        assert "1.23456e-20" in file.read_text()

    def test_write_errors(self, tmp_path):
        with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh()
        fmesh = meshtal.mesh[124]
        with pytest.raises(KeyError, match="Invalid format"):
            fmesh.write(tmp_path, out_format="xls")
        fmesh.grid.point_data["points"] = np.zeros(fmesh.grid.n_points)
        with pytest.raises(ValueError, match="mix between cell and point data"):
            fmesh.write(tmp_path, out_format="csv")
        # either is fine
        fmesh.write(tmp_path, ["points"], out_format="csv", outfile="points")
        fmesh.write(tmp_path, ["Value - Total"], out_format="csv", outfile="cells")
        assert len((tmp_path / "points.csv").read_text().splitlines()) == 18
        assert len((tmp_path / "cells.csv").read_text().splitlines()) == 4

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

        meshtally.readMesh((4, 14))
        assert list(meshtally.mesh) == [4, 14]
        meshtally.readMesh(14)
        assert list(meshtally.mesh) == [14]
        with pytest.raises(KeyError, match="bad tally entry"):
            meshtally.readMesh(5)

    def test_unknown_filetype(self):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            with pytest.raises(KeyError, match="Unknown filetype"):
                Meshtal(inp, filetype="VTK")

    def test_print_info(self):
        with as_file(RESOURCES.joinpath("meshtal_transform")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh(2024)
        fmesh = meshtal.mesh[2024]
        info = fmesh.print_info()
        for name in ("x1bin", "x2bin", "x3bin", "tally", "ebin"):
            assert info[name] is getattr(fmesh, name)
        assert info["geometry"] == "rec"

    def test_same_mesh_position(self):
        with as_file(RESOURCES.joinpath("meshtal_collapse")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh()
        names = {4: ["A", "B"], 14: ["C", "D"]}
        assert meshtal.mesh[4].sameMesh(meshtal.mesh[14])
        meshtal.mesh[14].grid = meshtal.mesh[14].grid.translate((1000, 0, 0))
        assert not meshtal.mesh[4].sameMesh(meshtal.mesh[14])
        with pytest.raises(RuntimeError, match="not compatible"):
            meshtal.collapse_grids(names)

    def test_collapse_binned(self):
        with as_file(RESOURCES.joinpath("meshtal_time_energy_bins")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh()
        with pytest.raises(RuntimeError, match="no binning"):
            meshtal.collapse_grids({54: ["A", "B"]})

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

    def test_unsupported_sections(self):
        # D1SUNED source importance meshes are listed but cannot be read
        with as_file(RESOURCES.joinpath("test_srcimp")) as inp:
            meshtal = Meshtal(inp)
        assert meshtal.mesh_ids == [114, 5114, 8114]
        with pytest.raises(ValueError, match="cannot be read"):
            meshtal.readMesh(5114)

    def test_cdgs_values(self):
        with as_file(RESOURCES.joinpath("cdgs_test")) as inp:
            meshtal = Meshtal(inp, filetype="CDGS")
        meshtal.readMesh()
        fmesh = meshtal.mesh[1]
        assert fmesh.type == "CDGS" and fmesh.particle is None
        assert fmesh.cooling_time == 0.0
        assert fmesh.strength == pytest.approx(2.0416e-3)
        assert list(fmesh.ebin) == [0, 5, 10]
        # the voxel totals add up to the mesh source
        assert fmesh.data[0, -1, ..., 0].sum() == pytest.approx(2.0416e-3, rel=1e-4)
        # voxel 1 (x, y, z indices 0, 0, 0): its two energy bins and Total
        assert fmesh.data[0, :, 0, 0, 0, 0] == pytest.approx(
            [8.81932e-04, 0.0, 8.81932e-04]
        )

    def test_negative_value_run_into_coordinate(self):
        # "112.500-1.57688E-11": the value fills its field
        with as_file(RESOURCES.joinpath("meshtal_bug")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh(114)
        data = meshtal.mesh[114].data
        assert not np.isnan(data).any()
        assert (data[..., 0] < 0).any()

    def test_d1s_discrete_bins(self):
        with as_file(RESOURCES.joinpath("assembly_meshtal_test")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh([44, 9044, 9944])
        assert meshtal.mesh[44].ebin.type == "nuc"
        assert list(meshtal.mesh[44].ebin) == pytest.approx(
            [0.001, 26059.002, 27060.003]
        )
        assert meshtal.mesh[9044].ebin.type == "nuc"  # parent nuclides
        assert meshtal.mesh[9944].ebin.type == "cel"
        assert list(meshtal.mesh[9944].ebin) == pytest.approx([0.001, 10.002, 20.003])
        # 3 nuclides and their total
        assert meshtal.mesh[44].data.shape == (1, 4, 2, 2, 2, 2)
        assert "Value - Total_eTot" in meshtal.mesh[44].grid.array_names

    def test_cylinder_grid_and_1d(self):
        # a full-revolution theta bin is kept as one bin in the data; the grid
        # is the exact cylinder (theta split in 5 degree cells)
        with as_file(RESOURCES.joinpath("1D_mesh")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh(214)
        fmesh = meshtal.mesh[214]
        assert fmesh.geom == "cyl"
        assert fmesh.data.shape == (1, 1, 1, 1, 79, 2)
        assert list(fmesh.x3bin) == [0, 1]
        assert fmesh.grid.n_cells == 79 * 72
        _, df, _ = fmesh.convert2tally()
        assert list(df.columns) == ["Cor A", "Value", "Error"]
        assert list(df["Value"]) == pytest.approx(fmesh.data[0, 0, 0, 0, :, 0])

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
        # the grid is stored in double precision, so these tiny values are kept
        assert list(grid.cell_data["Value - Total"]) == [1, 9, 1.77613e-306, 16]
        assert list(grid.cell_data["Error - Total"]) == [0, 0, 2e-101, 0]

    def test_decay_time_bins_have_no_total(self):
        with as_file(RESOURCES.joinpath("meshtal_time_energy_bins")) as inp:
            meshtal = Meshtal(inp)
        meshtal.readMesh(64)
        fmesh = meshtal.mesh[64]
        assert fmesh.tbin.type == "dtme" and not fmesh.tbin.totalbin
        names = fmesh.grid.array_names
        assert "Value - Total_t003_eTot" in names
        assert not any("tTot" in name for name in names)
        with pytest.warns(DeprecationWarning, match="select_bins"):
            selected = fmesh.get_etbin_data(ebin=-1, tbin=2)
        np.testing.assert_array_equal(selected.data[0, 0], fmesh.data[2, -1])

    def test_convert2tally_theta(self):
        with as_file(RESOURCES.joinpath("meshtal_cyl")) as inp:
            tally = aci_tallies.meshtal.open(inp)[4]
        # one axis only: all its voxels
        for selection, column, n in (
            ({"i": 0, "j": 0}, "Cor C", 10),
            ({"i": 0, "k": 0}, "Cor B", 10),
            ({"j": 0, "k": 0}, "Cor A", 3),
        ):
            _, df, _ = Fmesh.from_aci(tally.isel(**selection)).convert2tally()
            assert list(df.columns) == [column, "Value", "Error"]
            assert len(df) == n
        # two theta bins are a second axis
        fmesh = Fmesh.from_aci(tally.isel(i=0, k=slice(0, 2)))
        with pytest.raises(RuntimeError, match="only be used for 1D meshes"):
            fmesh.convert2tally()

    def test_gvr_of_cylinder_at_origin(self, tmp_path):
        # as MCNP5 writes it, without the VEC
        text = RESOURCES.joinpath("1D_mesh").read_text()
        text = text.replace(
            "origin at   0.00E+00  0.00E+00 -1.00E+03",
            "origin at   0.00E+00  0.00E+00  0.00E+00",
        )
        text = text.replace(", VEC direction   1.00E+00  0.00E+00  0.00E+00", "")
        file = tmp_path / "meshtal"
        file.write_text(text)
        meshtal = Meshtal(file)
        meshtal.readMesh(214)
        assert meshtal.mesh[214].trsf.origin.tolist() == [0, 0, 0]
        assert meshtal.mesh[214].trsf.axis.tolist() == [0, 0, 1]
        ww = WW.create_gvr_from_meshtally_file(file, tally_num=214)
        assert ww.geometry.director_1[:2] == [0, 0]

    def test_rescale_dose_mesh(self):
        with as_file(RESOURCES.joinpath("meshtal_time_energy_bins")) as inp:
            meshtal = Meshtal(inp)

        meshtal.readMesh()

        irr_file = IrradiationFile.from_text(resources_irr.joinpath("d1stime_o_30d"))
        scaling_factors = irr_file.get_scaling_factors_cooling_time(
            (200, TIME_UNITS.DAY)
        )
        rescaled_mesh = meshtal.mesh[44].rescale_dose_map(
            scaling_factors_df=scaling_factors,
            cooling_time_col="200d",
        )
        assert (
            pytest.approx(rescaled_mesh.cell_data["Value - Total_e002"][0])
            == 0.04435259766561963
        )

    def test_from_d1s_meshtal(self):
        with (
            as_file(RESOURCES.joinpath("meshtal_d1s_overflow")) as meshtal_file,
            as_file(RESOURCES.joinpath("input_d1s_overflow")) as input_file,
        ):
            inp = Input.from_input(input_file)
            meshtal = Meshtal.from_d1s_meshtal(meshtal_file, inp)
        meshtal.readMesh()
        assert meshtal.mesh_ids == [44, 54, 64, 74, 84]
        assert list(meshtal.mesh[44].ebin) == pytest.approx(
            [0.001, 41092900.002, 99973182.003]
        )
        assert list(meshtal.mesh[54].tbin) == pytest.approx(
            [0.001, 26059.002, 41092900.003]
        )
        assert list(meshtal.mesh[64].ebin) == pytest.approx(
            [0.001, 27060.002, 41092900.003]
        )

        # standard neutron flux mesh in the same file
        neutron = meshtal.mesh[84]
        assert neutron.particle == "neutron"
        assert neutron.ebin.type == "erg"
        assert list(neutron.ebin) == pytest.approx([0, 1, 20])
        # value and relative error of the first voxel: two energy bins and the total
        assert neutron.data[0, :, 0, 0, 0, :] == pytest.approx(
            np.array(
                [
                    [7.77445e-04, 1.45738e-01],
                    [1.62023e-02, 5.76506e-02],
                    [1.69798e-02, 5.63547e-02],
                ]
            )
        )


class TestMeshtalDeepcopy:
    """deepcopy copies the meshes and shares the reader of the file"""

    @staticmethod
    def _check_copy(meshtal: Meshtal) -> None:
        new = deepcopy(meshtal)
        assert new.mesh_ids == meshtal.mesh_ids
        assert new.mesh.keys() == meshtal.mesh.keys()
        for key, fmesh in meshtal.mesh.items():
            copied = new.mesh[key]
            assert copied is not fmesh
            assert copied.grid is not fmesh.grid
            assert (copied.trsf is None) == (fmesh.trsf is None)
            if fmesh.trsf is not None:
                np.testing.assert_array_equal(copied.trsf.axis, fmesh.trsf.axis)
            name = fmesh.grid.array_names[0]
            np.testing.assert_array_equal(copied.grid[name], fmesh.grid[name])
            np.testing.assert_array_equal(copied.data, fmesh.data)
            # changing the copy leaves the original as it was
            original = fmesh.grid[name].copy()
            data = fmesh.data
            copied.data = copied.data + 1
            assert not np.array_equal(copied.grid[name], original)
            np.testing.assert_array_equal(fmesh.grid[name], original)
            np.testing.assert_array_equal(fmesh.data, data)
            copied.grid[name][:] = -1
            np.testing.assert_array_equal(fmesh.grid[name], original)
            assert copied.select_bins(ebin=0).tally == fmesh.tally
            # a shallow copy shares the grid
            assert copy(fmesh).grid is fmesh.grid
            # pickled
            loaded = pickle.loads(pickle.dumps(fmesh))
            np.testing.assert_array_equal(loaded.data, fmesh.data)
            assert loaded.ebin.type == fmesh.ebin.type
        # the copy can read the file
        new.readMesh(norm="ctot")
        for key, fmesh in meshtal.mesh.items():
            np.testing.assert_array_equal(new.mesh[key].data, fmesh.data)

    @pytest.mark.parametrize(
        ["input_meshtal", "filetype"],
        [
            ("meshtal_CUBE_SQUARE", "MCNP"),
            ("meshtal_cyl", "MCNP"),
            ("meshtal_cuv", "CUV"),
            ("cdgs_test", "CDGS"),
        ],
    )
    def test_deepcopy(self, input_meshtal, filetype):
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtal = Meshtal(inp, filetype=filetype)
            meshtal.readMesh(norm="ctot")
        self._check_copy(meshtal)

    def test_deepcopy_from_d1s_meshtal(self):
        with (
            as_file(RESOURCES.joinpath("meshtal_d1s_overflow")) as meshtal_file,
            as_file(RESOURCES.joinpath("input_d1s_overflow")) as input_file,
        ):
            inp = Input.from_input(input_file)
            meshtal = Meshtal.from_d1s_meshtal(meshtal_file, inp)
            meshtal.readMesh()
        self._check_copy(meshtal)


class TestFmeshValues:
    """The values of a mesh are the arrays of its grid"""

    FILES = [
        ("meshtal_CUBE_SQUARE", "MCNP"),
        ("meshtal_rect_VV", "MCNP"),
        ("meshtal_cyl", "MCNP"),
        ("1D_mesh", "MCNP"),
        ("meshtal_time_energy_bins", "MCNP"),
        ("meshtal_transform", "MCNP"),
        ("meshtal_cuv", "CUV"),
        ("cdgs_test", "CDGS"),
    ]

    @staticmethod
    def _read(input_meshtal: str, filetype: str = "MCNP") -> Meshtal:
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtal = Meshtal(inp, filetype=filetype)
            meshtal.readMesh(norm="ctot")
        return meshtal

    @staticmethod
    def _new_data(fmesh: Fmesh) -> np.ndarray:
        """Different in every bin and voxel, so that a wrong order shows."""
        return np.arange(np.prod(fmesh.shape), dtype=float).reshape(fmesh.shape)

    @pytest.mark.parametrize(["input_meshtal", "filetype"], FILES)
    def test_data_as_read(self, input_meshtal, filetype):
        meshtal = self._read(input_meshtal, filetype)
        for key, fmesh in meshtal.mesh.items():
            file = meshtal._file
            tally = file.to_mesh(key, norm="ctot") if filetype == "CUV" else file[key]
            # aci-tallies: (energy, time, i, j, k)
            for i, array in enumerate((tally.result, tally.rel_error)):
                np.testing.assert_array_equal(
                    fmesh.data[..., i], array.transpose(1, 0, 4, 3, 2)
                )

    @pytest.mark.parametrize(["input_meshtal", "filetype"], FILES)
    def test_assign_data(self, input_meshtal, filetype):
        for fmesh in self._read(input_meshtal, filetype).mesh.values():
            names = fmesh.grid.array_names
            data = self._new_data(fmesh)
            fmesh.data = data
            np.testing.assert_array_equal(fmesh.data, data)
            assert fmesh.grid.array_names == names
            for it in range(fmesh.nt):
                for ie in range(fmesh.ne):
                    selected = fmesh.select_bins(ebin=ie, tbin=it)
                    np.testing.assert_array_equal(
                        selected.data[0, 0], fmesh.data[it, ie]
                    )

    @pytest.mark.parametrize(["input_meshtal", "filetype"], FILES)
    def test_grid_edits_are_used(self, input_meshtal, filetype):
        for fmesh in self._read(input_meshtal, filetype).mesh.values():
            original = fmesh.data
            for it in range(fmesh.nt):
                for ie in range(fmesh.ne):
                    fmesh.grid[fmesh.array_name(it, ie)] *= 2
            np.testing.assert_array_equal(fmesh.data[..., 0], 2 * original[..., 0])
            np.testing.assert_array_equal(fmesh.data[..., 1], original[..., 1])
            np.testing.assert_array_equal(
                fmesh.select_bins(ebin=-1, tbin=-1).data[0, 0], fmesh.data[-1, -1]
            )

    def test_grid_edits_in_convert2tally(self):
        fmesh = self._read("1D_mesh").mesh[214]
        _, before, _ = fmesh.convert2tally()
        fmesh.grid["Value - Total"] *= 2
        _, after, _ = fmesh.convert2tally()
        np.testing.assert_array_equal(after["Value"], 2 * before["Value"])
        np.testing.assert_array_equal(after["Error"], before["Error"])

    def test_cylinder_cells(self):
        """Voxels split in several cells along theta"""
        fmesh = self._read("meshtal_cyl").mesh[4]
        data = self._new_data(fmesh)
        fmesh.data = data
        assert list(fmesh.theta_cells) == [8] * 10
        assert fmesh.grid.n_cells == 8 * data[..., 0].size
        self._assert_cells_where_voxels_are(fmesh)

        fmesh.grid["Value - Total"][0] = -1
        with pytest.raises(ValueError, match="cells of a same voxel"):
            _ = fmesh.data

    def test_theta_bins_of_different_width(self, tmp_path):
        r, z, theta = [0, 1, 3], [0, 2, 5, 6], [0, 0.05, 0.25, 1]
        lines = []
        for i in range(2):
            for j in range(3):
                for k in range(3):
                    center = [
                        sum(b[n : n + 2]) / 2 for b, n in ((r, i), (z, j), (theta, k))
                    ]
                    row = "".join(f"{c:10.3f}" for c in center)
                    lines.append(f"{row} {100 * k + 10 * j + i:.5E} 1.00000E-01")
        bounds = [" ".join(f"{b:.3f}" for b in bins) for bins in (r, z, theta)]
        file = tmp_path / "meshtal"
        file.write_text(
            "mcnp   version 6     ld=02/20/18  probid =  11/18/19 14:41:20\n"
            " TEST\n"
            " Number of histories used for normalizing tallies =             1.00\n\n"
            " Mesh Tally Number         4\n neutron  mesh tally.\n\n"
            " Tally bin boundaries:\n"
            "               origin at  1.00E+00  0.00E+00  0.00E+00 axis in   0.00E+00  "
            "0.00E+00  1.00E+00 direction, VEC direction   1.00E+00  0.00E+00  0.00E+00\n"
            f"    R direction: {bounds[0]}\n"
            f"    Z direction: {bounds[1]}\n"
            f"    Theta direction (revolutions): {bounds[2]}\n"
            "    Energy bin boundaries: 0.00E+00 1.00E+36\n\n"
            "        R         Z         Th    Result     Rel Error\n"
            + "\n".join(lines)
            + "\n"
        )
        meshtal = Meshtal(file)
        meshtal.readMesh()
        fmesh = meshtal.mesh[4]
        # 5 degrees at most: 18, 72 and 270 degrees
        assert list(fmesh.theta_cells) == [4, 15, 54]
        k, j, i = np.indices((3, 3, 2))
        np.testing.assert_array_equal(fmesh.data[0, 0, ..., 0], 100 * k + 10 * j + i)
        self._assert_cells_where_voxels_are(fmesh)

    @staticmethod
    def _assert_cells_where_voxels_are(fmesh: Fmesh) -> None:
        """Each cell has the value of the voxel where its middle is"""
        points = fmesh.grid.cell_centers().points - fmesh.trsf.origin
        axis, vec = fmesh.trsf.axis, fmesh.trsf.vec
        x, y, z = points @ vec, points @ np.cross(axis, vec), points @ axis
        k = np.searchsorted(fmesh.x3bin, (np.arctan2(y, x) / (2 * np.pi)) % 1) - 1
        j = np.searchsorted(fmesh.x2bin, z) - 1
        i = np.searchsorted(fmesh.x1bin, np.hypot(x, y)) - 1
        assert len(set(zip(k, j, i))) == fmesh.nx1 * fmesh.nx2 * fmesh.nx3
        for error, name in enumerate(("Value - Total", "Error - Total")):
            np.testing.assert_array_equal(
                fmesh.grid[name], fmesh.data[0, 0, ..., error][k, j, i]
            )

    def test_rectangular_cells(self):
        fmesh = self._read("meshtal_transform").mesh[2024]
        fmesh.data = data = self._new_data(fmesh)
        points = fmesh.grid.cell_centers().points
        i, j, k = (
            np.searchsorted(bins, points[:, n]) - 1
            for n, bins in enumerate((fmesh.x1bin, fmesh.x2bin, fmesh.x3bin))
        )
        assert len(set(zip(i, j, k))) == fmesh.grid.n_cells
        np.testing.assert_array_equal(
            fmesh.grid["Value - Total"], data[0, 0, ..., 0][k, j, i]
        )

    def test_select_bins(self):
        fmesh = self._read("meshtal_time_energy_bins").mesh[54]
        selected = fmesh.select_bins(ebin=0, tbin=0, binlabels=("flux", "err"))
        assert selected.grid.array_names == ["flux", "err"]
        # it has its own values and grid
        original = fmesh.data
        selected.data = self._new_data(selected)
        assert selected.grid["flux"].max() == selected.data[..., 0].max()
        selected.grid = selected.grid.cast_to_structured_grid()
        selected.grid.points += 1
        np.testing.assert_array_equal(fmesh.data, original)
        # selecting nothing is a copy
        whole = fmesh.select_bins()
        assert whole.grid.array_names == fmesh.grid.array_names
        np.testing.assert_array_equal(whole.data, original)
        whole.grid["Value - Total_t001_e001"][:] = -1
        np.testing.assert_array_equal(fmesh.data, original)
        with pytest.raises(IndexError):
            fmesh.select_bins(ebin=fmesh.ne)
        with pytest.raises(IndexError):
            fmesh.select_bins(tbin=-fmesh.nt - 1)

    def test_transformed_stays_transformed(self):
        meshtal = self._read("meshtal_transform")
        untouched = deepcopy(meshtal)
        with as_file(RESOURCES.joinpath("transforms.i")) as input_file:
            meshtal.transform_fmesh(Input.from_input(input_file))
        moved = 0
        for key, fmesh in meshtal.mesh.items():
            points = np.array(fmesh.grid.points)
            moved += not fmesh.sameMesh(untouched.mesh[key])
            data = self._new_data(fmesh)
            fmesh.data = data
            np.testing.assert_array_equal(fmesh.data, data)
            np.testing.assert_array_equal(fmesh.grid.points, points)
            selected = fmesh.select_bins(ebin=0, tbin=0)
            np.testing.assert_array_equal(selected.grid.points, points)
        assert moved

    def test_bad_data_changes_nothing(self):
        fmesh = self._read("meshtal_CUBE_SQUARE").mesh[124]
        original = fmesh.data
        for bad in (original[..., 0], original[:, :, :, :1]):
            with pytest.raises(ValueError, match="it must be"):
                fmesh.data = bad
        np.testing.assert_array_equal(fmesh.data, original)

    def test_data_is_read_only(self):
        fmesh = self._read("meshtal_CUBE_SQUARE").mesh[124]
        with pytest.raises(ValueError, match="read-only"):
            fmesh.data[..., 0] = 7

    def test_grid_of_another_mesh(self):
        fmesh = self._read("meshtal_CUBE_SQUARE").mesh[124]
        original = fmesh.grid.copy()
        data = fmesh.data
        fmesh.grid.rename_array("Error - Total", "other")
        for use in (lambda: fmesh.data, lambda: fmesh.select_bins(ebin=0)):
            with pytest.raises(ValueError, match="no cell array 'Error - Total'"):
                use()
        fmesh.grid = original.threshold(5, scalars="Value - Total")
        for use in (
            lambda: fmesh.data,
            lambda: fmesh.select_bins(ebin=0),
            lambda: setattr(fmesh, "data", data),
        ):
            with pytest.raises(ValueError, match="cells instead of"):
                use()

    def test_collapsed_grid_has_its_own_values(self):
        meshtal = self._read("meshtal_collapse")
        grid = meshtal.collapse_grids({4: ["v4", "e4"], 14: ["v14", "e14"]})
        for key in (4, 14):
            before = meshtal.mesh[key].data
            np.testing.assert_array_equal(
                grid[f"v{key}"], meshtal.mesh[key].grid["Value - Total"]
            )
            grid[f"v{key}"][:] = -1
            np.testing.assert_array_equal(meshtal.mesh[key].data, before)


class TestBins:
    @pytest.mark.parametrize(
        ["values", "type", "n", "total", "suffixes"],
        [
            ([0, 100], "erg", 1, False, [""]),
            ([0, 1, 20], "erg", 3, True, ["_e001", "_e002", "_eTot"]),
            ([0, 1e20], "tme", 1, False, [""]),
            ([0.001, 10.002], "cel", 3, True, ["_e001", "_e002", "_eTot"]),
            ([-1], "nuc", 1, False, [""]),
            ([1, 2, 3], "line", 4, True, ["_e001", "_e002", "_e003", "_eTot"]),
            ([1, 2, 3], "dtme", 3, False, ["_e001", "_e002", "_e003"]),
        ],
    )
    def test_bins(self, values, type, n, total, suffixes):
        bins = Bins(np.array(values, dtype=float), type)
        assert (bins.n, bins.totalbin) == (n, total)
        assert [bins.suffix("e", i) for i in range(n)] == suffixes
        # as an array
        assert list(bins) == list(np.asarray(bins)) == values
        assert bins[-1] == values[-1] and len(bins) == len(values)

    @pytest.mark.parametrize(
        ["values", "type", "index", "expected", "indices"],
        [
            ([0, 1, 20], "erg", None, [0, 1, 20], [0, 1, 2]),
            ([0, 1, 20], "erg", 1, [1, 20], [1]),
            ([0, 1, 20], "erg", -1, [0, 20], [2]),
            ([0, 100], "erg", -1, [0, 100], [0]),
            ([1, 2], "nuc", 0, [1], [0]),
            ([1, 2], "nuc", 2, [-1], [2]),
            ([1, 2], "dtme", -1, [2], [1]),
        ],
    )
    def test_select(self, values, type, index, expected, indices):
        bins, selected = Bins(np.array(values, dtype=float), type).select(index)
        assert list(bins) == expected and list(selected) == indices
        assert bins.type == type and bins.n == len(indices)

    @pytest.mark.parametrize("index", [3, -4])
    def test_select_out_of_range(self, index):
        with pytest.raises(IndexError):
            Bins(np.array([0.0, 1.0, 20.0]), "erg").select(index)
