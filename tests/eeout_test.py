from importlib.resources import files, as_file
import pytest

from f4enix.output.eeout import EEOUT
import tests.resources.eeout as eeout_res

eeout_resources = files(eeout_res)


class TestEEOUT:

    @pytest.mark.parametrize(
        ["file", "exp"],
        [
            [
                "cyl_tetra1.eeout",
                {
                    "NUMBER OF PARTICLES": 1,
                    "NUMBER OF NODES": 391,
                    "NUMBER OF MATERIALS": 2,
                    "NUMBER OF INSTANCES": 2,
                    "NUMBER OF 1st TETS": 1306,
                    "NUMBER OF 1st PENTS": 0,
                    "NUMBER OF 1st HEXS": 0,
                    "NUMBER OF 2nd TETS": 0,
                    "NUMBER OF 2nd PENTS": 0,
                    "NUMBER OF 2nd HEXS": 0,
                    "NUMBER OF HISTORIES": 1000000,
                    "NUMBER OF EDITS": 2,
                },
            ],
            [
                "cyl_tetra2.eeout",
                {
                    "NUMBER OF PARTICLES": 1,
                    "NUMBER OF NODES": 2424,
                    "NUMBER OF MATERIALS": 2,
                    "NUMBER OF INSTANCES": 2,
                    "NUMBER OF 1st TETS": 0,
                    "NUMBER OF 1st PENTS": 0,
                    "NUMBER OF 1st HEXS": 0,
                    "NUMBER OF 2nd TETS": 1327,
                    "NUMBER OF 2nd PENTS": 0,
                    "NUMBER OF 2nd HEXS": 0,
                    "NUMBER OF HISTORIES": 1000000,
                    "NUMBER OF EDITS": 2,
                },
            ],
        ],
    )
    def test_read_info(self, file, exp):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)

        eeout.__repr__()
        eeout.__str__()

        assert eeout.info == exp

    @pytest.mark.parametrize("file", ["cyl_tetra1.eeout", "cyl_tetra2.eeout"])
    def test_read_particle_list(self, file):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        assert eeout._read_particle_list() == [1]

    @pytest.mark.parametrize("file", ["cyl_tetra1.eeout", "cyl_tetra2.eeout"])
    def test_get_materials_name(self, file):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        assert eeout._get_materials_name() == {1: "water_001", 2: "steel_002"}

    @pytest.mark.parametrize(
        ["file", "exp_idx"], [["cyl_tetra1.eeout", 334], ["cyl_tetra2.eeout", 1552]]
    )
    def test_read_nodes_xyz(self, file, exp_idx):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        points, idx_elem_type = eeout._read_nodes_xyz()

        assert points.shape == (eeout.n_nodes, 3)
        assert idx_elem_type == exp_idx

    @pytest.mark.parametrize(
        ["file", "exp_idx"], [["cyl_tetra1.eeout", 535], ["cyl_tetra2.eeout", 1756]]
    )
    def test_read_material(self, file, exp_idx):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        mat_ids, idx = eeout._read_material()
        assert len(mat_ids) == eeout.n_elem
        assert set(mat_ids) == {1, 2}
        assert idx == exp_idx

    @pytest.mark.parametrize(
        ["file", "exp_idx"], [["cyl_tetra1.eeout", 1843], ["cyl_tetra2.eeout", 4412]]
    )
    def test_read_connectivity(self, file, exp_idx):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        cells, idx = eeout._read_connectivity()
        assert len(cells) == eeout.n_elem
        assert idx == exp_idx

    @pytest.mark.parametrize("file", ["cyl_tetra1.eeout", "cyl_tetra2.eeout"])
    def test_read_edits(self, file):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        edits, idx = eeout._read_edits()
        assert len(edits) == 2
        assert len(edits["Tally1_par1_FLUX_14"]["values"]) == eeout.n_elem
        assert len(edits["Tally1_par1_FLUX_14"]["errors"]) == eeout.n_elem

    @pytest.mark.parametrize("file", ["cyl_tetra1.eeout", "cyl_tetra2.eeout"])
    def test_read_rho_vol(self, file):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)
        dens, vols = eeout._read_rho_vol()
        assert len(dens) == eeout.n_elem
        assert len(vols) == eeout.n_elem

    @pytest.mark.parametrize("file", ["cyl_tetra1.eeout", "cyl_tetra2.eeout"])
    def test_export(self, file, tmpdir):
        with as_file(eeout_resources.joinpath(file)) as inp:
            eeout = EEOUT(inp)

        eeout.export(tmpdir)
