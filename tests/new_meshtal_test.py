from importlib.resources import as_file, files

import pytest

import tests.resources.meshtal as resources
import tests.resources.meshtal.expected as res_exp
import tests.resources.meshtal.tests as res
from f4enix.output.new_meshtal import NewMeshtal

resources_write = files(res)
expected = files(res_exp)
RESOURCES = files(resources)


class TestNewMeshtal:
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
        ],
    )
    def test_mesh_print_tally_info(self, input_meshtal, filetype):
        """To check if the meshtal can be read without any problem"""
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = NewMeshtal(filename=inp, filetype=filetype)
            meshtally.readMesh()

        for _mesh_id, mesh in meshtally.mesh.items():
            mesh.print_info()

    def test_same_mesh(self):
        with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
            meshtally = NewMeshtal(inp)
        meshtally.readMesh()
        assert meshtally.mesh[124].sameMesh(meshtally.mesh[124])

    @pytest.mark.parametrize(
        "norm",
        [
            "vtot",
            "celf",
        ],
    )
    def test_read_mesh_cuv(self, norm):
        # To check if the meshtal can be read without any problem"
        filetype = "CUV"
        with as_file(RESOURCES.joinpath("meshtal_cuv")) as inp:
            meshtally = NewMeshtal(inp, filetype)

        for i in meshtally.mesh.items():
            meshtally.readMesh(norm=norm)
            meshtally.readMesh(cell_filters=[1, 2], norm=norm)

    @pytest.mark.parametrize("input_meshtal", ["1D_mesh", "1D_mesh_energyonly"])
    def test_1d_features(self, input_meshtal):
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtal = NewMeshtal(inp)
        meshtal.readMesh()
        meshtal.mesh[214].convert2tally()
