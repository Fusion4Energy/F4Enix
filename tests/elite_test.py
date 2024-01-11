import os
from copy import deepcopy
from importlib.resources import files, as_file
# import tests.resources.input as input_res
import tests.resources.elite as elite_res

from f4enix.input.MCNPinput import Input
from f4enix.input.elite import Elite_Input

resources_elite = files(elite_res)

class TestElite_Input:
    with as_file(resources_elite.joinpath('E-Lite_dummy.i')) as FILE1:
        testInput = Elite_Input.from_input(FILE1)

    def test_extract_sector(self, tmpdir):
        tol = 1e-4
        excel_name = 'E-Lite_dummy_Block_Structure_Summary.xlsx'
        inp = deepcopy(self.testInput)
        outfile = tmpdir.mkdir('sub').join('sector1.i')

        inp.extract_sector(1, excel_file=resources_elite.joinpath(excel_name),
                           outfile=outfile, tol=tol, check_Elite=True)
        # re-read
        inp_sec1 = Input.from_input(outfile)
        assert ') 427024 ) 427016' in inp_sec1.cells['800'].lines[1]
        assert ') :-427024  ) :-427016' in inp_sec1.cells['801'].lines[17]
        pbc_surfs = ['*7', '*8', '*427016', '*427024']
        assert set(pbc_surfs).issubset(set(inp_sec1.surfs.keys()))
        assert inp_sec1.surfs['*7'].scoefs[3] == inp.surfs['7'].scoefs[3]
        assert inp_sec1.surfs['*8'].scoefs[3] == inp.surfs['8'].scoefs[3]
        assert inp_sec1.surfs['110'].scoefs[3] == 2*tol
        assert inp_sec1.surfs['105'].scoefs[3] == - 2*tol
        assert inp_sec1.other_data['SI70'].lines[0].rstrip() ==  'SI70 L 427001'
        assert inp_sec1.other_data['SP70'].lines[0].rstrip() ==  'SP70 1' 

        outfile = os.path.join(os.path.dirname(outfile), 'NBI.i')

        inp.extract_sector('2 & 3', excel_file=resources_elite.joinpath(excel_name),
                           outfile=outfile, tol=tol, check_Elite=False)

        inp_NBI = Input.from_input(outfile)
        assert ') -437544  ) 437543' in inp_NBI.cells['800'].lines[1]
        assert ') :437544 ) :-437543' in inp_NBI.cells['801'].lines[17]
        pbc_surfs = ['*7', '*18', '*437543', '*437544']
        assert set(pbc_surfs).issubset(set(inp_NBI.surfs.keys()))
        assert inp_NBI.surfs['*7'].scoefs[3] == inp.surfs['7'].scoefs[3]
        assert inp_NBI.surfs['*18'].scoefs[3] == inp.surfs['18'].scoefs[3]
        assert inp_NBI.surfs['210'].scoefs[3] == - 2*tol
        assert inp_NBI.surfs['205'].scoefs[3] == 2*tol
        assert inp_NBI.other_data['SI70'].lines[0].rstrip() ==  'SI70 L 435001'
        assert inp_NBI.other_data['SP70'].lines[0].rstrip() ==  'SP70 2'

        outfile = os.path.join(os.path.dirname(outfile), 'sectors6_7.i')

        inp.extract_sector([6, 7], excel_file=resources_elite.joinpath(excel_name),
                           outfile=outfile, tol=tol, check_Elite=False)

        inp_secs6_7 = Input.from_input(outfile)

        assert ') 467024 ) 475016' in inp_secs6_7.cells['800'].lines[1]
        assert ') :-467024  ) :-475016' in inp_secs6_7.cells['801'].lines[17]
        pbc_surfs = ['*20', '*22', '*467024', '*475016']
        assert '21' in inp_secs6_7.surfs.keys()
        assert inp_secs6_7.surfs['21'].scoefs[3] == inp.surfs['21'].scoefs[3]
        assert set(pbc_surfs).issubset(set(inp_secs6_7.surfs.keys()))
        assert inp_secs6_7.surfs['*20'].scoefs[3] == inp.surfs['20'].scoefs[3]
        assert inp_secs6_7.surfs['*22'].scoefs[3] == inp.surfs['22'].scoefs[3]
        assert inp_secs6_7.surfs['110'].scoefs[3] == inp.surfs['110'].scoefs[3]
        assert inp_secs6_7.surfs['105'].scoefs[3] == inp.surfs['105'].scoefs[3]
        assert inp_secs6_7.other_data['SI70'].lines[0].rstrip() ==  'SI70 L 467001 475001'
        assert inp_secs6_7.other_data['SP70'].lines[0].rstrip() ==  'SP70 1 1'