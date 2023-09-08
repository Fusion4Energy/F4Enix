from importlib.resources import files, as_file

from f4enix.output.eeout import EEOUT
import tests.resources.eeout as eeout_res

eeout_resources = files(eeout_res)


class TestEEOUT:

    with as_file(eeout_resources.joinpath('cyl_tetra1.eeout')) as inp:
        eeout = EEOUT(inp)

    def test_read_info(self):

        self.eeout.__repr__()
        self.eeout.__str__()

        dic_expected = {
            'NUMBER OF PARTICLES': 1,
            'NUMBER OF NODES': 391,
            'NUMBER OF MATERIALS': 2,
            'NUMBER OF INSTANCES': 2,
            'NUMBER OF 1st TETS': 1306,
            'NUMBER OF 1st PENTS': 0,
            'NUMBER OF 1st HEXS': 0,
            'NUMBER OF 2nd TETS': 0,
            'NUMBER OF 2nd PENTS': 0,
            'NUMBER OF 2nd HEXS': 0,
            'NUMBER OF HISTORIES': 1000000,
            'NUMBER OF EDITS': 2}

        assert self.eeout.info == dic_expected

    def test_read_particle_list(self):
        assert self.eeout._read_particle_list() == [1]

    def test_read_nodes_xyz(self):
        nodes_x, nodes_y, nodes_z, idx_elem_type = self.eeout._read_nodes_xyz()

        assert len(nodes_x) == len(nodes_y) == len(nodes_z) == self.eeout.n_nodes
        assert idx_elem_type == 334

    def test_read_material(self):
        mat_ids, idx = self.eeout._read_material()
        assert len(mat_ids) == self.eeout.n_elem
        assert set(mat_ids) == {1, 2}
        assert idx == 535

    def test_read_connectivity(self):
        cells, idx = self.eeout._read_connectivity()
        assert len(cells) == self.eeout.n_elem
        assert idx == 1843

    def test_read_edits(self):
        edits, idx = self.eeout._read_edits()
        assert len(edits) == 2
        assert len(edits['Tally1_par1_FLUX_14']['values']) == self.eeout.n_elem
        assert len(edits['Tally1_par1_FLUX_14']['errors']) == self.eeout.n_elem
