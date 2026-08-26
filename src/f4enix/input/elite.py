import copy
import logging
import math
import os

import migjorn
import numpy as np
import pandas as pd

from f4enix.core.constants import (
    PLASMA_CELLS,
    SECTOR_BOUNDARIES,
    SECTOR_BOUNDARIES_ANGLES,
    SECTOR_NAMES,
)
from f4enix.input.MCNPinput import Input


class Elite_Input(Input):
    def __init__(
        self,
        *args,
    ) -> None:
        """Children of the :py:class:`Input`.

        It features the methods for sectors extraction from ITER E-Lite 360
        degree model

        Parameters
        ----------
        *args
            all arguments of parent class Input

        Examples
        --------
        Read the input and the associated Excel file.

        By calling extract_sector you'll be able to get a working MCNP input
        file of the desired sector.

        >>> from f4enix.input.elite import Elite_Input
        ... elite_inp = Elite_Input.from_input('Elite.i')
        ... elite_inp.extract_sector(1, 'Elite_excel.xlsx', 'sector1.i', tol=1e-5, True)
        """
        super().__init__(*args)
        self._initialized = False
        self._block_structure = None
        self._sectors_L0_cells_names: dict[str | int, list[int]] = {}

    def _initialize_elite(self, excel_file: os.PathLike, check_Elite: bool = True):
        # checks if the input is actually e-lite and initializes some variables
        # check is optional
        self._block_structure = pd.read_excel(excel_file)

        # check the envelopes
        if check_Elite:
            inp_cells = self.get_cells_summary()
            inp_L0_cells = inp_cells[inp_cells["universe"].isna()].index.tolist()
            if set(inp_L0_cells) != set(self._block_structure["Cell"].tolist()):
                raise RuntimeError(
                    "MCNP input is not an E-Lite file, or the Excel Block"
                    + " Structure file is not compatible"
                )

        for sec in SECTOR_NAMES:
            # you will have to modify outer cell to put them in the sectors
            self._sectors_L0_cells_names[sec] = self._block_structure[
                self._block_structure["Sector"] == sec
            ]["Cell"].tolist()

        self._initialized = True

    def extract_sector(
        self,
        sectors: int | list[int],
        excel_file: os.PathLike,
        outfile: os.PathLike = "sector",
        tol: float = 1e-5,
        check_Elite: bool = False,
        fix_source: bool = True,
    ) -> None:
        """Writes a working input of a user-selected E-Lite sector.

        The user can provide a sector number or a list of contiguous sector
        numbers in counterclockwise direction (e.g. 1, "2 & 3", [4,5], [9, 1]).
        Currently there is no check of the correctness of the input, the user
        must be careful in providing the correct sector number(s), in the
        correct order.
        The method will extract a working input of the selected sector(s), by
        replacing the boundaries with reflecting surfaces.
        The method is based on an auxiliary Excel file that lists all E-Lite
        envelope containers and their respective sector. The method follows
        these steps:

        - Excel file initialization and model check. If there's no correspondence
          between the envelope structure of E-Lite model and excel file, the
          method aborts. The method can't check if the correct sector number is
          assigned to the envelope container.

        - Collection of all envelope containers and fillers to be extracted

        - Graveyard and outer cell modification

        - Source term fixing

        - L0 surfaces which are enough close to the boundary will be set as reflective

        - L1 surfaces which are too close to the boundary will be translated
          "outwards" with respect to the sector(s), to avoid the arising of
          fatal errors

        Parameters
        ----------
        sectors : int or list
            sector number, or list of contiguous sector numbers in
            counterclockwise direction that will be extracted
        excel_file : os.PathLike
            path to the Excel file that describes the sector structure of the
            version of E-Lite in use
        outfile : os.PathLike, optional
            name of the input file that will be written, by default 'sector'
        tol : int, optional
            determines the maximum distance for two planes to be considered
            equal, and determines the magnitude of the translation of L1 cells.
            It should be chosen basedon the value used in DBCN card, by default 1e-5
        check_Elite : bool, optional
            Automatically checks if the envelope structure of the model is
            consistent with the one reported in the Excel, by default False
        fix_source : bool, optional
            If True, the source will be modified to be consistent with the
            extracted sector(s), by default True
        """
        if not self._initialized:
            self._initialize_elite(excel_file, check_Elite)

        # build list of sectors to be extracted
        if not isinstance(sectors, list):
            sectors = [sectors]

        logging.info("Collecting the cells, surfaces, materials and transf.")
        # collect L0 cells to be extracted
        cells: list[int] = []
        for sector in sectors:
            cells += self._sectors_L0_cells_names[sector]
        # collect L0 surfaces, needed for later
        L0_cells = set(cells)
        L0_sset = []
        for _, cell in enumerate(self.cells.values()):
            if cell.id in L0_cells:
                L0_sset.extend(cell.surface_ids)
        L0_sset = set(L0_sset)
        # backup copy graveyard and outercell, that will be modified
        # append gy and outercell manually as they don't belong to a sector
        cells.append(800)
        cells.append(801)
        # get cells, surfaces and materials to be extracted
        extracted_input = self.extract_cells(
            cells=cells,
            keep_universe=True,
        )

        if fix_source:
            _set_sdef(sectors, extracted_input)
        # set L0 as periodic and modify L1 planes
        _set_boundaries(
            _get_boundaries_angles(sectors),
            tol,
            L0_sset,
            extracted_input,
        )
        # modify graveyard and outercell
        _modify_graveyard(extracted_input, sectors)
        extracted_input.write(outfile)

        logging.info("input written correctly")


def _set_sdef(sectors: list[str | int], inp: Input) -> None:
    # write new SI and SD cards, directly in input attribute
    # the copies are modified, original input is preserved
    new_si = "SI70 L "
    new_sp = "SP70 "
    for sector in sectors:
        new_si = new_si + str(PLASMA_CELLS[sector]) + " "
        if sector != "2 & 3":
            new_sp = new_sp + "1 "
        else:
            new_sp = new_sp + "2 "
    inp.other_data["SI70"] = new_si + "\n"
    inp.other_data["SP70"] = new_sp + "\n"


def _set_boundaries(
    boundaries_angles: list[float],
    tol: float,
    L0_sset: set[int],
    inp: Input,
):
    # define the angles at which there are the planes cutting the sectors
    boundary_angles = []
    for angle in boundaries_angles:
        if angle > 180:
            rev_angle = angle - 180
        else:
            rev_angle = angle + 180
        boundary_angles.append(angle)
        boundary_angles.append(rev_angle)

    for l, angle in enumerate(boundary_angles):
        coeffs = np.array(
            [-math.sin(math.radians(angle)), math.cos(math.radians(angle)), 0]
        )
        # Check all surfaces for each angle
        for surf in inp._model.surfaces():
            # check if the surface belongs to L0 or L1 of the sector(s)
            if surf.id in L0_sset:
                bound_opt = 1
            else:
                bound_opt = 2
            # check if the surface is a plane
            kind = surf.kind.lower()
            if kind in ["p", "px", "py", "pz"]:
                if kind == "p":
                    p_coeffs = np.array(surf.coeffs[:3])
                    d = surf.coeffs[3]
                elif kind == "px":
                    p_coeffs = np.array([1, 0, 0])
                    d = surf.coeffs[0]
                elif kind == "py":
                    p_coeffs = np.array([0, 1, 0])
                    d = surf.coeffs[0]
                elif kind == "pz":
                    p_coeffs = np.array([0, 0, 1])
                    d = surf.coeffs[0]
                p_coeffs = p_coeffs / np.linalg.norm(p_coeffs)
                # check if the surface is a plane passing by the origin
                if _check_tol(tol, [(d, 0)]):
                    # check if the plane has a transformation
                    tr = surf.transform
                    if tr:
                        # if so, apply it and check if it's a plane parallel to
                        # the boundary
                        tr_name = "TR" + str(tr)
                        trans = inp.transformations[tr_name]
                        # check if it's only a translation
                        # if len(trans.values) < 12:
                        #     continue
                        # if it's a rotation, rotate the coefficients
                        p_coeffs = _rotate_plane(trans, p_coeffs)
                    # if not, just check if it's parallel to the boundary
                    # if the plane is similar to boundary planes, modify it
                    if _check_tol(tol, list(zip(p_coeffs, coeffs))):
                        _modify_boundary(surf, bound_opt, angle, l, tol)


def _get_boundaries_angles(sectors: list[int | str]) -> list[float]:
    # get the angles of the two boundary surfaces, counterclockwise
    boundaries_angles = []

    for k, sec in enumerate(sectors):
        bounds = SECTOR_BOUNDARIES_ANGLES[sec]
        if k == 0:
            # get y- plane for first sector in list
            boundaries_angles.append(bounds[0])
        if k == (len(sectors) - 1):
            # get y+ plane for last sector in list
            boundaries_angles.append(bounds[1])

    return boundaries_angles


def _modify_graveyard(inp: Input, sectors: list[int | str]):
    # cut/ union graveyard and outercell with planes
    for k, sector in enumerate(sectors):
        if k == 0:
            Input.add_surface(
                inp.cells["801"],
                -SECTOR_BOUNDARIES[sector][0],
                mode="union",
            )
            Input.add_surface(
                inp.cells["800"],
                SECTOR_BOUNDARIES[sector][0],
                mode="intersect",
            )
        if k == len(sectors) - 1:
            Input.add_surface(
                inp.cells["801"], -SECTOR_BOUNDARIES[sector][1], mode="union"
            )
            Input.add_surface(
                inp.cells["800"],
                SECTOR_BOUNDARIES[sector][1],
                mode="intersect",
            )
        else:
            continue


def _rotate_plane(trans: migjorn.Transform, p_coeffs: np.ndarray) -> np.ndarray:
    # compute rotation matrix
    coeffs = trans.coeffs
    if trans.degrees:
        coeffs = [math.cos(math.radians(v)) for v in coeffs]
    # define rotation matrix
    rot_matrix = np.array(
        [
            [coeffs[3], coeffs[4], coeffs[5]],
            [coeffs[6], coeffs[7], coeffs[8]],
            [coeffs[9], coeffs[10], coeffs[11]],
        ]
    )
    # compute inverse rotation matrix
    inv_rot = np.linalg.inv(rot_matrix)
    # compute rotated plane's coefficients to be checked
    dp = np.dot(inv_rot, p_coeffs)

    return dp


def _check_tol(tol: float, n_coeff: list[tuple[float, float]]) -> bool:
    # for all tuples in the list, checks if the elements in the tuples
    # are within the tolerance
    in_tol = True
    # loop over the first list
    for coeff in n_coeff:
        if not coeff[1] - tol <= coeff[0] <= coeff[1] + tol:
            in_tol = False
            break
    return in_tol


def _modify_boundary(
    surf: migjorn.Surface, bound_opt: int, angle: float, l: float, tol: float
):
    # if in L0, set periodic
    tol_sign = {True: 1, False: -1}

    if bound_opt == 1:
        surf.reflective = True

    # if L1, translate outwards to avoid fatal errors
    elif bound_opt == 2:
        # put the correct sign to the translation
        # check if y+ or y- boundary
        y_plus = l > 1
        # check in which quadrant the plane is
        angle_quadr = not (90 < angle < 270 or -270 < angle < -90)
        # check the sign of the y coefficient of the plane
        y_coeff = surf.coeffs[1] > 0
        # translate the plane and recompute template and input
        sign = tol_sign[y_plus] * tol_sign[angle_quadr] * tol_sign[y_coeff]
        new_offset = surf.coeffs[-1] + sign * 2 * tol
        surf.set_coeff(len(surf.coeffs) - 1, new_offset)
