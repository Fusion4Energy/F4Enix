from importlib.resources import as_file, files

import pytest

import tests.resources.meshtal as resources
import tests.resources.meshtal.expected as res_exp
import tests.resources.meshtal.tests as res
from f4enix.output.new_meshtal import NewMeshtal

resources_write = files(res)
expected = files(res_exp)
RESOURCES = files(resources)


class NewMeshtalTest:
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
            "1D_mesh_no_code",
        ],
    )
    def test_mesh_print_tally_info(self, input_meshtal):
        """To check if the meshtal can be read without any problem"""
        with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
            meshtally = NewMeshtal(filename=inp)
