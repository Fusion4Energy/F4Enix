import os
import pytest
import math
import numpy as np
from copy import deepcopy
from importlib.resources import files, as_file
from numjuggler import parser

# import tests.resources.input as input_res
import tests.resources.elite as elite_res

from f4enix.input.MCNPinput import Input
from f4enix.input.elite import Elite_Input

resources_elite = files(elite_res)


class TestElite_Input:
    with as_file(resources_elite.joinpath("E-Lite_dummy.i")) as FILE1:
        testInput = Elite_Input.from_input(FILE1)

    def test_extract_sector(self, tmpdir):
        tol = 1e-4
        excel_name = "E-Lite_dummy_Block_Structure_Summary.xlsx"
        inp = deepcopy(self.testInput)
        outfile = tmpdir.mkdir("sub").join("sector1.i")

        inp.extract_sector(
            1,
            excel_file=resources_elite.joinpath(excel_name),
            outfile=outfile,
            tol=tol,
            check_Elite=True,
        )
        # re-read
        inp_sec1 = Input.from_input(outfile)
        assert ") 427024 ) 427016" in inp_sec1.cells["800"].lines[1]
        assert ") :-427024 ) :-427016" in inp_sec1.cells["801"].lines[17]
        pbc_surfs = ["*7", "*8", "*427016", "*427024"]
        assert set(pbc_surfs).issubset(set(inp_sec1.surfs.keys()))
        assert inp_sec1.surfs["*7"].scoefs[3] == inp.surfs["7"].scoefs[3]
        assert inp_sec1.surfs["*8"].scoefs[3] == inp.surfs["8"].scoefs[3]
        assert inp_sec1.surfs["110"].scoefs[3] == 2 * tol
        assert inp_sec1.surfs["105"].scoefs[3] == -2 * tol
        assert inp_sec1.other_data["SI70"].lines[0].rstrip() == "SI70 L 427001"
        assert inp_sec1.other_data["SP70"].lines[0].rstrip() == "SP70 1"

        outfile = os.path.join(os.path.dirname(outfile), "NBI.i")

        inp.extract_sector(
            "2 & 3",
            excel_file=resources_elite.joinpath(excel_name),
            outfile=outfile,
            tol=tol,
            check_Elite=False,
        )

        inp_NBI = Input.from_input(outfile)
        assert ") -437544 ) 437543" in inp_NBI.cells["800"].lines[1]
        assert ") :437544 ) :-437543" in inp_NBI.cells["801"].lines[17]
        pbc_surfs = ["*7", "*18", "*437543", "*437544"]
        assert set(pbc_surfs).issubset(set(inp_NBI.surfs.keys()))
        assert inp_NBI.surfs["*7"].scoefs[3] == inp.surfs["7"].scoefs[3]
        assert inp_NBI.surfs["*18"].scoefs[3] == inp.surfs["18"].scoefs[3]
        assert inp_NBI.surfs["210"].scoefs[3] == -2 * tol
        assert inp_NBI.surfs["205"].scoefs[3] == 2 * tol
        assert inp_NBI.other_data["SI70"].lines[0].rstrip() == "SI70 L 435001"
        assert inp_NBI.other_data["SP70"].lines[0].rstrip() == "SP70 2"

        outfile = os.path.join(os.path.dirname(outfile), "sectors6_7.i")

        inp.extract_sector(
            [6, 7],
            excel_file=resources_elite.joinpath(excel_name),
            outfile=outfile,
            tol=tol,
            check_Elite=False,
        )

        inp_secs6_7 = Input.from_input(outfile)

        assert ") 467024 ) 475016" in inp_secs6_7.cells["800"].lines[1]
        assert ") :-467024 ) :-475016" in inp_secs6_7.cells["801"].lines[17]
        pbc_surfs = ["*20", "*22", "*467024", "*475016"]
        assert "21" in inp_secs6_7.surfs.keys()
        assert inp_secs6_7.surfs["21"].scoefs[3] == inp.surfs["21"].scoefs[3]
        assert set(pbc_surfs).issubset(set(inp_secs6_7.surfs.keys()))
        assert inp_secs6_7.surfs["*20"].scoefs[3] == inp.surfs["20"].scoefs[3]
        assert inp_secs6_7.surfs["*22"].scoefs[3] == inp.surfs["22"].scoefs[3]
        assert inp_secs6_7.surfs["110"].scoefs[3] == inp.surfs["110"].scoefs[3]
        assert inp_secs6_7.surfs["105"].scoefs[3] == inp.surfs["105"].scoefs[3]
        assert (
            inp_secs6_7.other_data["SI70"].lines[0].rstrip() == "SI70 L 467001 475001"
        )
        assert inp_secs6_7.other_data["SP70"].lines[0].rstrip() == "SP70 1 1"

    def test_set_sdef(self):
        inp = deepcopy(self.testInput)
        mod_data = deepcopy(inp.other_data)
        Elite_Input._set_sdef([4], mod_data)

        assert mod_data["SI70"].input[0].rstrip() == "SI70 L 451001"
        assert mod_data["SP70"].input[0].rstrip() == "SP70 1"

        Elite_Input._set_sdef(["2 & 3"], mod_data)
        assert mod_data["SI70"].input[0].rstrip() == "SI70 L 435001"
        assert mod_data["SP70"].input[0].rstrip() == "SP70 2"

        Elite_Input._set_sdef([1, "2 & 3"], mod_data)
        assert mod_data["SI70"].input[0].rstrip() == "SI70 L 427001 435001"
        assert mod_data["SP70"].input[0].rstrip() == "SP70 1 2"

        Elite_Input._set_sdef(["2 & 3", 4], mod_data)
        assert mod_data["SI70"].input[0].rstrip() == "SI70 L 435001 451001"
        assert mod_data["SP70"].input[0].rstrip() == "SP70 2 1"

        Elite_Input._set_sdef([8, 9], mod_data)
        assert mod_data["SI70"].input[0].rstrip() == "SI70 L 483001 491001"
        assert mod_data["SP70"].input[0].rstrip() == "SP70 1 1"

    def test_initialize_elite(self):
        inp = deepcopy(self.testInput)
        excel_name = "E-Lite_dummy_Block_Structure_Summary.xlsx"
        assert inp._sectors_L0_cells_names == {}
        assert inp._Elite_Input__initialized == False
        inp._initialize_elite(resources_elite.joinpath(excel_name), check_Elite=True)
        assert inp._Elite_Input__initialized == True
        assert 800 not in inp._sectors_L0_cells_names[1]
        assert 451001 in inp._sectors_L0_cells_names[4]
        assert 435001 in inp._sectors_L0_cells_names["2 & 3"]

        inp = deepcopy(self.testInput)
        excel_name = "E-Lite_dummy_Block_Structure_Summary_wrong.xlsx"

        with pytest.raises(
            RuntimeError,
            match="MCNP input is not an E-Lite file, or the Excel Block"
            + " Structure file is not compatible",
        ):
            inp._initialize_elite(
                resources_elite.joinpath(excel_name), check_Elite=True
            )

    def test_modify_graveyard(self):
        # cut/ union graveyard and outercell with planes
        inp = deepcopy(self.testInput)

        sectors = [1]
        outercell = inp.cells["800"]
        graveyard = inp.cells["801"]
        new_outercell, new_gy = inp._modify_graveyard(sectors, outercell, graveyard)
        print_outercell = new_outercell.card(wrap=True, comment=False).splitlines()
        print_gy = new_gy.card(wrap=True, comment=False).splitlines()
        assert print_outercell[0].split()[2][:2] == "("
        assert print_gy[0].split()[2][:2] == "("
        assert ") 427024 ) 427016" in print_outercell[0]
        assert ") :-427024 ) :-427016" in print_gy[0]

        sectors = ["2 & 3"]
        new_outercell, new_gy = inp._modify_graveyard(sectors, outercell, graveyard)
        print_outercell = new_outercell.card(wrap=True, comment=False).splitlines()
        print_gy = new_gy.card(wrap=True, comment=False).splitlines()
        assert print_outercell[0].split()[2][:2] == "("
        assert print_gy[0].split()[2][:2] == "("
        assert ") -437544 ) 437543" in print_outercell[0]
        assert ") :437544 ) :-437543" in print_gy[0]

        sectors = [1, "2 & 3"]
        new_outercell, new_gy = inp._modify_graveyard(sectors, outercell, graveyard)
        print_outercell = new_outercell.card(wrap=True, comment=False).splitlines()
        print_gy = new_gy.card(wrap=True, comment=False).splitlines()
        assert print_outercell[0].split()[2][:2] == "("
        assert print_gy[0].split()[2][:2] == "("
        assert ") 427024 ) 437543" in print_outercell[0]
        assert ") :-427024 ) :-437543" in print_gy[0]

    def test_get_boundaries_angles(self):
        inp = deepcopy(self.testInput)
        # get the angles of the two boundary surfaces, counterclockwise
        assert inp._get_boundaries_angles([1]) == [10, 50]
        assert inp._get_boundaries_angles(["2 & 3", 4]) == [50, 170]

    def test_rotate_plane(self):
        inp = deepcopy(self.testInput)
        tol = 1e-7

        surf = inp.surfs["427016"]
        trans = inp.transformations["TR440"]
        p_coeffs = inp._rotate_plane(trans, np.array(surf.scoefs[:3]))
        # compute rotation matrix
        assert inp._check_tol(tol, list(zip(p_coeffs, [0.76604444, -0.6427876, 0])))

        trans = inp.transformations["TR440"]
        p_coeffs = inp._rotate_plane(trans, np.array([1, 0, 0]))
        # compute rotation matrix
        assert inp._check_tol(tol, list(zip(p_coeffs, [0.866025404, 0.5, 0])))

        surf = inp.surfs["437543"]
        trans = inp.transformations["TR587"]
        p_coeffs = inp._rotate_plane(trans, np.array(surf.scoefs[:3]))
        # compute rotation matrix
        assert inp._check_tol(tol, list(zip(p_coeffs, [0.98480775, -0.17364818, 0])))

    def test_check_tol(self):
        assert not self.testInput._check_tol(1e-5, [(1e-5, 3e-5), (1e-5, 1e-5 - 1e-6)])
        assert self.testInput._check_tol(
            1e-5, [(1e-5, 1e-5 + 1e-6), (1e-5, 1e-5 - 1e-6)]
        )

    def test_modify_boundary(self):
        inp = deepcopy(self.testInput)
        tol = 1e-6
        Elite_Input._modify_boundary(inp.surfs["427016"], 1, 50, 2, tol)
        assert inp.surfs["427016"].input[0][0] == "*"

        Elite_Input._modify_boundary(inp.surfs["110"], 2, 50, 2, tol)
        assert float(inp.surfs["110"].lines[0].rstrip().split()[-1]) == 2 * tol

        Elite_Input._modify_boundary(inp.surfs["205"], 2, 130, 2, tol)
        assert float(inp.surfs["205"].lines[0].rstrip().split()[-1]) == 2 * tol
