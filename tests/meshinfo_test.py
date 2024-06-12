import numpy as np
from importlib.resources import files, as_file
import pytest

import tests.resources.meshinfo as res
from f4enix.output.meshinfo import (
    MeshInfoFile,
    CoordinateType,
    MeshInfo,
    MeshInfoCyl,
    COLUMN_KEY_MASS_GRAMS,
)

RESOURCES = files(res)


class TestMeshInfoFile:
    def test_read_file_cart(self):
        with as_file(RESOURCES.joinpath("meshinfo_cart")) as infile:
            meshinfofile = MeshInfoFile.from_file(infile)

        meshinfo = meshinfofile.info[9014]
        assert CoordinateType.CART == meshinfo.coordinates
        assert pytest.approx(meshinfo.vector_j[0]) == 734.00
        assert len(meshinfo.vector_k) == 18
        filtered_dataframe = meshinfo.data_mass.get_filtered_dataframe(
            voxels=[1],
            materials=[110],
            cells=[324896],
        )
        result_mass = filtered_dataframe[COLUMN_KEY_MASS_GRAMS].values[0]
        # voxel_volume * cell_proportion_of_volume * density
        expected_mass = 3.38693e03 * 0.31100 * 8.0300
        assert pytest.approx(result_mass) == expected_mass
        # check that the specific cell 940030 is made of material 940001
        cells = meshinfo.data_mass.get_cells_from_materials(materials=[940001])
        assert 940030 in cells

    def test_read_file_cyl(self):
        with as_file(RESOURCES.joinpath("meshinfo_cyl")) as infile:
            meshinfofile = MeshInfoFile.from_file(infile)

        meshinfo = meshinfofile.info[4]
        assert CoordinateType.CYL == meshinfo.coordinates
        assert pytest.approx(0.0) == meshinfo.origin.sum()
        assert 22 == len(meshinfo.vector_i)
        assert pytest.approx(100) == meshinfo.vector_i[-1]
        assert pytest.approx(550) == meshinfo.vector_j[-1]
        assert pytest.approx(1) == meshinfo.vector_k[-1]
        filtered_dataframe = meshinfo.data_mass.get_filtered_dataframe(
            voxels=[108],
            materials=[1],
            cells=[2],
        )
        result_mass = filtered_dataframe[COLUMN_KEY_MASS_GRAMS].values[0]
        # from the file: voxel_volume * cell_proportion_of_volume * density
        expected_mass = 1.38791e03 * 1.00000 * 7.9300
        # WARNING: voxel_volume given by the file is wrong, it considered the theta
        # units as degrees instead of revolutions, we manually correct the mass
        expected_mass /= np.pi * 2
        assert pytest.approx(expected_mass) == result_mass
        # Check that our manual correction was also right by calculating the volume
        # from zero by hand
        manual_volume = np.pi * (3.75**2) * 5 * 7.9300
        # I reduce the precision, surely a harmless decimals thing
        assert pytest.approx(manual_volume, abs=1e-2) == expected_mass

    @pytest.mark.parametrize(
        ["file", "key", "obj"],
        [["meshinfo_cyl", 4, MeshInfoCyl], ["meshinfo_cart", 9014, MeshInfo]],
    )
    def test_save_load(self, file, key, obj, tmpdir):
        with as_file(RESOURCES.joinpath(file)) as infile:
            meshinfofile = MeshInfoFile.from_file(infile)
        meshinfo = meshinfofile.info[key]
        meshinfo.save(tmpdir)

        newmeshinfo = obj.load(tmpdir)
        assert meshinfo == newmeshinfo

    def test_can_read_multiple_meshes_same_file(self):
        with as_file(RESOURCES.joinpath("meshinfo_two_meshes")) as infile:
            meshinfofile = MeshInfoFile.from_file(infile)

        mesh_1 = meshinfofile.info[4]
        mesh_2 = meshinfofile.info[234]
        assert (
            mesh_1.data_mass.df.values[-1].sum() == mesh_2.data_mass.df.values[-1].sum()
        )
