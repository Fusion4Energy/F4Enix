"""
This module is related to the parsing and manipulation of MCNP input files.

The parser is built on the migjorn Rust-based lossless MCNP parser.
"""

from __future__ import annotations

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or - 
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at: 
    https://eupl.eu/1.2/en/ 
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""
import json
import logging
import os
import re
from copy import deepcopy
from collections.abc import MutableMapping
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import migjorn
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from f4enix.core.constants import (
    PAT_ALL_TALLY_KEYS,
    PAT_CARD_KEY,
    PAT_COMMENT,
    PAT_F_TR_CARD_KEY,
    PAT_FMESH_KEY,
    PAT_NP,
    UNION_INTERSECT_SYMBOLS,
)
from f4enix.core.irradiation import METASTABLE_TAG, Nuclide
from f4enix.input.d1suned import IrradiationFile, Reaction, ReactionFile
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import MatCardsList, Material
from f4enix.input.migjorn_proxies import (
    _CellsProxy,
    _OtherDataProxy,
    _SurfsProxy,
    _TransformsProxy,
)

PAT_MT = re.compile(r"m[tx]\d+", re.IGNORECASE)
PAT_BLANK_LINE = re.compile(r"\n[\s\t]*\n")
ADD_LINE_FORMAT = "         {}\n"

_PAT_MAT_CARD = re.compile(r"^M[TX]?\d", re.IGNORECASE)
_PAT_TR_CARD = re.compile(r"^\*?TR\d", re.IGNORECASE)
_PAT_CONTINUATION = re.compile(r"^[ \t]{5,}|^\t")


class Input:
    def __init__(
        self,
        model: migjorn.Model,
        mat_section: MatCardsList,
    ) -> None:
        """Class representing an MCNP input file.

        Cells, surfaces, materials and transformations are handled explicitly.
        All other datacards are treated generically.
        The parser is backed by migjorn (lossless Rust-based MCNP parser).

        Parameters
        ----------
        model : migjorn.Model
            parsed MCNP model (lossless)
        mat_section : MatCardsList
            material cards section of the input

        Attributes
        ----------
        # TODO: add attributes

        """
        self._model = model
        self._mat_section = mat_section

        # TODO: migjorn should have a header attribute
        lines = self._model.to_source().splitlines(keepends=True)
        header: list[str] = []
        for line in lines:
            if re.match(r"^\d", line.strip()):
                break
            header.append(line.replace("\r", ""))
        self._header = header

    def __deepcopy__(self, memo: dict) -> "Input":
        # migjorn.Model can't be pickled; rebuild from its text source
        new_obj = self.__class__.__new__(self.__class__)
        memo[id(self)] = new_obj
        new_obj.__init__(
            migjorn.Model(self._model.to_source()), deepcopy(self._mat_section, memo)
        )
        return new_obj

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def cells(self) -> _CellsProxy:
        """Proxy mapping of cells keyed by string cell ID.

        Supports ``inp.cells["800"] = cell`` to replace a cell in the model.
        Assigning a mapping replaces all cells: ``inp.cells = other_inp.cells``.
        """
        return _CellsProxy(self._model)

    @cells.setter
    def cells(self, value: list[migjorn.Cell]) -> None:
        # TODO: this can be improved migjorn side
        # delete all cells from model and add the new ones
        for cell in list(self._model.cells):
            self._model.remove_cell(cell.id)
        for cell in value:
            self._model.add_cell(cell.text)

    @property
    def surfs(self) -> _SurfsProxy:
        """Proxy mapping of surfaces keyed by string surface ID (with * prefix if reflective).

        Supports ``inp.surfs["10"] = surf`` to replace a surface in the model.
        Assigning a mapping replaces all surfaces: ``inp.surfs = other_inp.surfs``.
        """
        return _SurfsProxy(self._model)

    @surfs.setter
    def surfs(self, value: list[migjorn.Surface]) -> None:
        # TODO: this can be improved migjorn side
        for surf in list(self._model.surfaces):
            self._model.remove_surface(surf.id)
        for surf in value:
            self._model.add_surface(surf.text)

    @property
    def transformations(self) -> _TransformsProxy:
        """Proxy mapping of TR cards keyed by 'TRn'.

        Supports ``del inp.transformations['TR5']`` to remove a transform.
        Assigning a mapping removes all current transforms (add is unsupported by migjorn).
        Setting is not supported; modify the live handle in-place instead.
        """
        return _TransformsProxy(self._model)

    @transformations.setter
    def transformations(self, value: list[migjorn.Transform]) -> None:
        # TODO: migjorn has no add_transform; this setter is a no-op for now
        raise NotImplementedError

    @property
    def other_data(self) -> _OtherDataProxy:
        """Proxy mapping of non-material, non-transform data card texts.

        Supports ``inp.other_data['SI70'] = 'SI70 L 1\\n'`` to add or replace a card.
        """
        return _OtherDataProxy(self._model)

    @other_data.setter
    def other_data(self, value: list[migjorn.DataCard]) -> None:
        # TODO: there is no setter in migjorn
        raise NotImplementedError

    @property
    def mat_section(self) -> MatCardsList:
        return self._mat_section

    @mat_section.setter
    def mat_section(self, value: MatCardsList) -> None:
        self._mat_section = value

    @property
    def tally_keys(self) -> list[int]:
        keys = []
        for card in self.other_data:
            name = card.name
            m = PAT_ALL_TALLY_KEYS.match(name)
            if (
                m
                and name[0].upper() == "F"
                and name[:2].upper()
                not in ("FM", "FC", "FN", "FU", "FT", "FQ", "FS", "FP")
            ):
                try:
                    keys.append(int(m.group(2)))
                except (IndexError, ValueError, TypeError):
                    pass
        return keys

    @property
    def fmesh_keys(self) -> list[int]:
        keys = []
        for card in self.other_data:
            if PAT_FMESH_KEY.match(card.name):
                try:
                    keys.append(int(re.search(r"\d+", card.name).group()))
                except (AttributeError, ValueError):
                    pass
        return keys

    @property
    def header(self) -> list[str]:
        """Title and leading comment lines before the first cell."""
        return self._header

    @header.setter
    def header(self, value: list[str]) -> None:
        """Set the title and leading comment lines before the first cell."""
        self._header = value

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_input(cls, inputfile: os.PathLike | str) -> "Input":
        """Parse an MCNP input file using migjorn.

        Parameters
        ----------
        inputfile : os.PathLike | str
            path to the MCNP input file

        Returns
        -------
        Input
        """
        name = os.path.basename(str(inputfile)).split(".")[0]
        logging.info(f"Reading file: {name}")
        model = migjorn.Model.from_file(str(inputfile))
        for d in model.diagnostics:
            logging.warning(f"migjorn [{d.severity}]: {d.message}")
        logging.debug("Reading has finished")
        logging.debug("building material section")
        mat_section = MatCardsList.from_migjorn(model)
        logging.debug("Material section built")
        return cls(model, mat_section)

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def write(self, outfilepath: os.PathLike | str, wrap: bool = False) -> None:
        """Write the input to a file.

        Parameters
        ----------
        outfilepath : os.PathLike | str
            output file path
        wrap : bool
            ignored (kept for API compatibility); migjorn output is lossless
        """
        logging.info(f"Writing to {outfilepath}")
        with open(outfilepath, "w", newline="\n") as f:
            f.writelines(self.header)
            for cell in self._model.cells:
                f.write(cell.text.replace("\r", ""))
            # blank line between cells and surfaces
            f.write("\n")
            for surf in self._model.surfaces:
                f.write(surf.text.replace("\r", ""))
            # blank line between surfaces and materials
            f.write("\n")
            # write materials section
            f.write(self.mat_section.to_text())
            # write transformations and other data cards
            for tr in self._model.transforms:
                f.write(tr.text.replace("\r", ""))
            for card in self._model.data_cards:
                # TODO: missing data card text property
                f.write(card.text.replace("\r", ""))

        logging.info("File was written correctly")

    # ------------------------------------------------------------------
    # Structural operations
    # ------------------------------------------------------------------

    def merge(self, other_inp: "Input") -> None:
        """Merge another Input into this one.

        Parameters
        ----------
        other_inp : Input
            input to merge in
        """
        self._model.merge([other_inp._model])
        # TODO: double check here what happens if there are duplicate materials
        self._mat_section.extend(other_inp._mat_section.materials)

    def renumber(
        self,
        cells: int | None = None,
        surfs: int | None = None,
        universes: int | None = None,
        translations: int | None = None,
        renum_all: int | None = None,
    ) -> None:
        """Renumber IDs in the input by a constant offset.

        Parameters
        ----------
        cells : int, optional
            offset for cell IDs
        surfs : int, optional
            offset for surface IDs
        universes : int, optional
            offset for universe IDs
        translations : int, optional
            offset for transformation IDs
        renum_all : int, optional
            applies the same offset to all of the above
        """
        if renum_all is not None:
            cells = surfs = universes = translations = renum_all
        if cells is not None:
            self._model.offset_cells(int(cells))
        if surfs is not None:
            self._model.offset_surfaces(int(surfs))
        if universes is not None:
            self._model.renumber_universes(lambda u: u + int(universes))
        if translations is not None:
            self._model.renumber_transforms(lambda t: t + int(translations))

    def translate(self, newlib: str | dict, libmanager: LibManager) -> None:
        """
        Translate the input to another library

        Parameters
        ----------
        newlib : dict | str
            There are a few ways that newlib can be provided:

            1) str (e.g. 31c), the new library to translate to will be the
            one indicated;

            2) dic (e.g. {'98c' : '99c', '31c: 32c'}), the new library is
            determined based on the old library of the zaid

            3) dic (e.g. {'98c': [list of zaids], '31c': [list of zaids]}),
            the new library to be used is explicitly stated depending
            on the zaidnum.

        libmanager : libmanager.LibManager
            Library manager for the conversion.

        Returns
        -------
        None.

        """

        try:
            if isinstance(newlib, str) and newlib[0] == "{":
                # covert the dic
                newlib = json.loads(newlib)
        except KeyError:
            # It is already a dict, pass
            pass

        self.mat_section.translate(newlib, libmanager)

    def get_materials_subset(self, ids: list[str] | str) -> MatCardsList:
        """Return a subset of materials by ID."""
        if type(ids) is str:
            mats = [self.mat_section[ids.upper()]]
        else:
            mats = [self.mat_section[mid.upper()] for mid in ids]

        return MatCardsList(mats)

    def extract_cells(
        self,
        cells: list[int],
        outfile: os.PathLike | str,
        renumber_offsets: dict | None = None,
        keep_universe: bool = True,
        extract_fillers: bool = True,
    ):
        """given a list of cells, dumps a minimum MCNP working file.

        The file will includes all the requested cells, defined surfaces,
        materials and translations.

        Parameters
        ----------
        cells : list[int]
            desired list of cells
        outfile : os.PathLike | str
            path to the file where the MCNP input needs to be dumped
        renumber_offsets : dict, optional
            apply the self.renumber() function to the extracted input.
            the dict will be passed as keyargs to the function.
            Default is None.
        keep_universe: bool
            If True keeps the 'U=' key in the cell cards, otherwise that is
            removed. Default is True.
        extract_fillers: bool
            if True extract also the cells belonging to a universe that is
            used in a 'FILL=' keyword. This happens recursively. default is
            True.
        """
        logging.info("write MCNP reduced input")
        cell_ids = [int(c) for c in cells]
        newinput = self._extract_cells_as_input(
            cell_ids, keep_universe, extract_fillers
        )
        if renumber_offsets is not None:
            newinput.renumber(**renumber_offsets)
        newinput.write(outfile)

    def _extract_cells_as_input(
        self,
        cell_ids: list[int],
        keep_universe: bool = True,
        extract_fillers: bool = True,
    ) -> "Input":
        """Build a minimal Input containing only the requested cells and their references."""
        # TODO: In migjorn there is already extract_universe, should be something similar
        cset: set[int] = set(cell_ids)
        self._collect_cell_refs(cset, extract_fillers)

        # Collect referenced surface and material IDs
        sset: set[int] = set()
        mset: set[str] = set()
        for cid in cset:
            cell = self._model.cell(cid)
            if cell is None:
                continue
            sset.update(abs(s) for s in cell.signed_surfaces)
            if cell.material and cell.material != 0:
                mset.add(f"M{cell.material}")

        # Build clean source directly from individual card texts — this avoids
        # the stale-handle issue that arises after bulk remove_cell/remove_surface.
        source = self._model.to_source()
        # Grab cell texts from a fresh model (remove/add may invalidate handles)
        full_model = migjorn.Model(source)

        def _cell_text(cid: int) -> str:
            c = full_model.cell(cid)
            if c is None:
                return ""
            text = c.text.replace("\r", "").rstrip("\n") + "\n"
            if not keep_universe:
                # Strip u= parameter inline rather than mutating the model
                text = re.sub(r"\s+[uU]=\S+", "", text)
            return text

        def _surf_text(sid: int) -> str:
            s = full_model.surface(sid)
            return s.text.replace("\r", "").rstrip("\n") + "\n" if s else ""

        cells_src = "".join(_cell_text(cid) for cid in sorted(cset))
        surfs_src = "".join(_surf_text(sid) for sid in sorted(sset))
        title = source.splitlines()[0].replace("\r", "")
        clean_source = title + "\n" + cells_src + "\n" + surfs_src + "\n"
        new_model = migjorn.Model(clean_source)

        # Materials subset
        mat_ids = [mid.upper() for mid in mset]
        materials = self.get_materials_subset(mat_ids) if mat_ids else MatCardsList([])
        if isinstance(materials, Material):
            materials = MatCardsList([materials])

        return Input(new_model, materials)

    def _collect_cell_refs(self, cset: set[int], extract_fillers: bool) -> None:
        """Expand cset to include all referenced (#n complement and fill) cells."""
        cell_set = set(cset)
        while cell_set:
            new_set: set[int] = set()
            uni_set: set[int] = set()
            for cid in cell_set:
                cell = self._model.cell(cid)
                if cell is None:
                    continue
                # Add cells referenced via #n complements
                new_set.update(cell.cell_refs)
                if extract_fillers and cell.fill is not None:
                    uni_set.add(cell.fill.universe)
            if extract_fillers:
                for c in self._model.cells:
                    if c.universe in uni_set:
                        new_set.add(c.id)
            cell_set = new_set - cset
            cset |= cell_set

    def extract_universe(
        self,
        universe: int,
        renumber_offsets: dict | None = None,
        keep_universe: bool = False,
    ) -> "Input":
        """Dump a minimum MCNP working file for the given universe.

        Parameters
        ----------
        universe : int
            universe id to be extracted
        renumber_offsets : dict, optional
            offsets passed to renumber(), by default None
        keep_universe : bool
            keep u= keyword in cell definitions, by default False
        """
        # extract the universe
        extracted_model = self._model.extract_universe(universe)
        extracted_inp = Input(extracted_model)
        mat_ids = []
        for mat in extracted_model.materials:
            mat_ids.append(f"M{mat.id}")
        mat_subset = self.get_materials_subset(mat_ids)
        extracted_inp.mat_section = mat_subset

        # renumber if requested
        if renumber_offsets is not None:
            extracted_inp.renumber(**renumber_offsets)
        # remove u= keywords if requested
        if not keep_universe:
            for cell in extracted_inp._model.cells:
                cell.remove_param("u")

        return extracted_inp

    def get_cells_by_matID(
        self, matID: int | str, deepcopy_flag: bool = True
    ) -> dict[str, migjorn.Cell]:
        """Return all cells assigned to matID.

        Parameters
        ----------
        matID : int | str
            material ID to filter the cells
        deepcopy_flag: bool
            ignored (kept for API compatibility)

        Returns
        -------
        dict[str, migjorn.Cell]
            cells assigned to that material
        """
        logging.debug(f"get cells for material {matID} requested")
        filtered_cells = {}
        for key, cell in self.cells.items():
            if cell.material == int(matID):
                filtered_cells[key] = cell
        return filtered_cells

    def scale_densities(self, factor: float) -> None:
        """Scale the density values of all cells by the same factor. Void
        cells are ignored. Resulting density will be equal to
        original_density*factor

        Parameters
        ----------
        factor : float
            scaling factors for the densities
        """
        for _, cell in self.cells.items():
            if cell.material and cell.material != 0 and cell.density is not None:
                cell.density = cell.density * factor

    def get_densities_range(self) -> pd.DataFrame:
        """Return a DataFrame listing for all material the minimum and maximum
        density values encountered in the input cells.

        Returns
        -------
        pd.DataFrame
            Range of densities used for each material.
        """
        lm = LibManager()
        materials = []
        for _, cell in self.cells.items():
            mat_id = cell.material
            if mat_id == 0 or mat_id is None:
                continue  # ignore void cells
            dens = float(cell.density or 0.0)

            # if the density is negative (mass) leave it as it is,
            # if atomic fraction is used, convert it to mass first
            if dens < 0:
                dens = abs(dens)
            else:
                dens = self.mat_section[f"M{mat_id}"].get_density(dens, lm)

            materials.append([mat_id, dens])

        df = pd.DataFrame(materials, columns=["Material ID", "Density"])
        min_vals = df.groupby("Material ID")["Density"].min()
        max_vals = df.groupby("Material ID")["Density"].max()
        summary = pd.DataFrame()
        summary["Min density [g/cc]"] = min_vals
        summary["Max density [g/cc]"] = max_vals
        summary["Material ID"] = max_vals.index
        return summary.set_index("Material ID").sort_index()

    def get_cells_summary(self) -> pd.DataFrame:
        """Get a summary of infos for each cell

        A DataFrame is returned where for each cell is listed the material,
        density, universe and filler is present.

        Returns
        -------
        pd.DataFrame
            Summary of cells info
        """
        rows = []
        for key, cell in self.cells.items():
            row = {"cell": int(key)}
            row["material"] = cell.material
            row["density"] = cell.density
            row["universe"] = cell.universe
            row["filler"] = cell.fill.universe if cell.fill else None
            rows.append(row)

        df = pd.DataFrame(rows)
        return df.set_index("cell").sort_index()

    def _get_tally_cards_ids(
        self,
        idx: int,
    ) -> list[str]:
        keys = []
        pat = re.compile(r"F[a-zA-Z]*{}$".format(idx))
        for card in self.other_data:
            if pat.match(card.name) is not None:
                keys.append(card.name)
        return keys

    def _retrieve_input(self, tag: str) -> str:
        # get the card text excluding the card name tag and $ comments
        text = self.other_data[tag].text
        first_line = text.splitlines()[0] if text else ""
        inp = first_line.split("$")[0]  # strip inline $ comment
        inp = inp.replace(tag, "").replace(tag.lower(), "").strip()
        return inp

    def _retrieve_FM(self, tag: str) -> list[str]:
        text = self.other_data[tag].text
        # Strip inline $ comments and blank/comment lines
        lines = []
        for ln in text.splitlines():
            ln_clean = ln.split("$")[0].strip()
            if ln_clean and not re.match(r"^[cC](\s|$)", ln_clean):
                lines.append(ln_clean)
        first_line = self._retrieve_input(tag).split("$")[0].split()
        if len(lines) <= 1:
            return first_line
        multi = ["N.A."]
        if first_line:
            multi.append(str(first_line))
        for line in lines[1:]:
            multi.append(line.strip())
        return multi

    def get_tally_summary(self, fmesh: bool = False) -> pd.DataFrame:
        """Get a summary of the tallies defined in the input

        Both normal tallies and fmeshes can be requested. For each tally the
        number, particle, description and multiplier are listed (if available)

        Parameters
        ----------
        fmesh : bool, optional
            if True produced a summary for the fmehses instead of for the
            normal tallies, by default False

        Returns
        -------
        pd.DataFrame
            summary info on defined tallies

        """

        if fmesh:
            tag_tally = "FMESH"
            tallies = self.fmesh_keys
        else:
            tag_tally = "F"
            tallies = self.tally_keys

        rows = []
        for key in tallies:
            desc = np.nan
            particle = np.nan
            multiplier = None
            card_keys = self._get_tally_cards_ids(key)
            for aux_key in card_keys:
                if aux_key[:2] == "FC":
                    desc = self._retrieve_input(aux_key)
                elif aux_key == tag_tally + str(key):
                    card = self.other_data[aux_key]
                    particle = card.particle
                elif aux_key[:2] == "FM":
                    multiplier = self._retrieve_FM(aux_key)

            row = {"Tally": key, "Particle": particle, "Description": desc}

            if multiplier is not None:
                row["Normalization"] = multiplier[0]
                if len(multiplier) > 1:
                    row["Other multipliers"] = multiplier[1:]
                else:
                    row["Other multipliers"] = np.nan
            else:
                row["Normalization"] = np.nan
                row["Other multipliers"] = np.nan

            rows.append(row)

        return pd.DataFrame(rows).set_index("Tally").sort_index()

    def replace_material(
        self,
        new_mat_id: int,
        new_density: str | float,
        old_mat_id: int,
        u_list: list[int] | None = None,
    ) -> None:
        """Replace a material and density in the input with other values.

        Parameters
        ----------
        new_mat_id : int
            id of the new material (0 for void)
        new_density : str | float
            new value for the density (including sign)
        old_mat_id : int
            id of the material to be replaced
        u_list : list[int]
            change the material only if cells belong to one of the universes
            in the list. By default is None, all cells are affected.
        """

        if new_mat_id < 0 or old_mat_id < 0:
            raise ValueError("Wrong values for the material ids")
        for _, cell in self.cells.items():
            in_universe = False
            if u_list is not None:
                if cell.universe in u_list:
                    in_universe = True
            else:
                in_universe = True

            if cell.material == old_mat_id and in_universe:
                cell.material = new_mat_id
                if not cell.is_void:
                    cell.density = float(new_density)

    @staticmethod
    def set_param(
        cell: migjorn.Cell,
        param: str,
        param_value: int | str,
        inplace: bool = True,
    ) -> migjorn.Cell:
        """Add a u= or fill= parameter to a cell (always in-place with migjorn).

        Parameters
        ----------
        cell : migjorn.Cell
            cell to modify
        param : str
            parameter to set, these are the ones accepted by migjorn
        param_value : int | str
            value to set for the parameter
        inplace : bool
            if True, modifies the cell in-place; if False, returns a modified copy

        Returns
        -------
        migjorn.Cell
        """
        # TODO: understand or import what are allowed parameters for mig
        allowed_params = ["u", "fill"]
        if param not in allowed_params:
            raise ValueError(
                f"Invalid parameter '{param}'. Allowed parameters are: {allowed_params}"
            )

        if not inplace:
            cell = deepcopy(cell)

        replaced = cell.set_param(param, str(param_value))
        if not replaced:
            cell.add_param(f"{param}={param_value}")

        return cell

    @staticmethod
    def add_surface(
        cell: migjorn.Cell,
        add_surface: int,
        new_cell_num: int | None = None,
        mode: str = "intersect",
        inplace: bool = True,
    ) -> migjorn.Cell:
        """Add a surface to a cell's geometry as union or intersection.

        Parameters
        ----------
        cell : migjorn.Cell
            cell to modify
        add_surface : int
            signed surface number to add
        new_cell_num : int, optional
            ignored with migjorn (cells are mutated in-place in the model)
        mode : str, optional
            'intersect' (default) or 'union'.
            NOTE: union mode requires migjorn.Cell.add_surface_union (TODO migjorn)
        inplace: bool
            ignored (migjorn mutations are always in-place)
        """
        if not inplace:
            cell = deepcopy(cell)
        if mode.lower() == "intersect":
            cell.add_surface(add_surface)
        elif mode.lower() == "union":
            # TODO migjorn: add Cell.add_surface_union(surface: int) method
            raise NotImplementedError(
                "Union surface addition requires migjorn.Cell.add_surface_union (not yet available)"
            )
        else:
            raise ValueError(f"Invalid mode {mode}. Use 'union' or 'intersect'.")

        return cell

    @staticmethod
    def hash_cell(
        cell: migjorn.Cell,
        hash_id: int,
        new_cell_num: int | None = None,
        inplace: bool = True,
    ) -> migjorn.Cell:
        """Add a #hash_id complement to a cell's geometry.

        Parameters
        ----------
        cell : migjorn.Cell
            cell to modify
        hash_id : int
            ID of the cell to complement
        new_cell_num : int, optional
            ignored with migjorn
        inplace : bool, optional
            ignored (migjorn mutations are always in-place)
        """
        if not inplace:
            cell = deepcopy(cell)
        if new_cell_num is not None:
            # TODO migjorn: add Cell.id setter to allow renumbering
            cell.id = new_cell_num
        cell.add_complement(hash_id)
        return cell

    def hash_multiple_cells(self, hash_dict: dict[int, list[int]]) -> None:
        """Add #hash_id complements to cells.

        Parameters
        ----------
        hash_dict : dict[int, list[int]]
            {hash_id: [cell_ids to add it to]}
        """
        for hash_id, cell_ids in hash_dict.items():
            for cell_num in cell_ids:
                cell = self._model.cell(int(cell_num))
                if cell is not None:
                    cell.add_complement(hash_id)

    def cells_union(
        self,
        cell_num_list: list[str],
        new_cell_num: int | None = None,
    ) -> None:
        """Create a union of cells, replacing them with a single cell.

        Parameters
        ----------
        cell_num_list : list[str]
            cell numbers to unite. The first one will be used as base
            for material, density and parameters.
        new_cell_num : int, optional
            ID for the resulting cell; defaults to the first cell's ID
        """
        # Ensure they are all in the model
        cells: list[migjorn.Cell] = [self._model.cell(int(n)) for n in cell_num_list]
        if None in cells:
            raise ValueError("One or more cell IDs not found in the model.")

        # Ensure they all have the same material
        mat = cells[0].material
        for c in cells:
            if c.material != mat:
                raise ValueError("All cells must have the same material.")

        if new_cell_num is None:
            new_cell_num = cells[0].id

        # Build union geometry text from each cell's geometry part
        geom_parts = [_extract_cell_geometry(c.text) for c in cells]
        union_geom = " : ".join(f"({g})" for g in geom_parts)

        base = cells[0]
        mat_part = "0 " if base.is_void else f"{base.material} {base.density} "
        params_part = " ".join(
            f"{p.key}{':{}'.format(p.particle) if p.particle else ''}={p.value}"
            for p in base.params
        )
        new_text = f"{new_cell_num} {mat_part}{union_geom} {params_part}\n"

        for c in cells:
            self._model.remove_cell(c.id)
        self._model.add_cell(new_text.strip())

    def add_F_tally(
        self,
        tally_ID: int,
        particles: list[str],
        cells: Sequence[int | str],
        add_total: bool = False,
        energies: Sequence[float] | np.ndarray | None = None,
        description: str | None = None,
        multiplier: str | None = None,
        add_SD: bool = True,
    ):
        """Add a F-tally to the input

        Parameters
        ----------
        tally_ID : int
            tally ID
        particles : list[str]
            particle to tally
        cells : list[int | str]
            list of cells to tally. One can also add a string with a complex cell
            definition, like unions of cells and cells in a universe,
            e.g. "((1 2 3) < 10)"
        add_total : bool, optional
            if True adds the total tally, by default False
        energies : list[float], optional
            list of energies for the tally, by default None
        description : str, optional
            description of the tally, by default None
        multiplier : str, optional
            multiplier of the tally, by default None
        add_SD : bool, optional
            if True adds the SD 1 card, by default True
        """
        # Add the tally card
        tally_ID = int(tally_ID)
        self.tally_keys.append(tally_ID)

        particles_str = particles[0]
        if len(particles) > 1:
            for particle in particles[1:]:
                particles_str += "," + particle
        if description is not None:
            self.other_data[f"FC{tally_ID}"] = f"FC{tally_ID} {description}\n"
        # Build tally card text
        cells_list = list(cells)
        if add_total:
            cells_list.append("T")
        cells_str = " ".join(str(c) for c in cells_list)
        self.other_data[f"F{tally_ID}"] = f"F{tally_ID}:{particles_str} {cells_str}\n"
        if energies is not None:
            energies_str = " ".join(f"{e:.4e}" for e in energies)
            self.other_data[f"E{tally_ID}"] = f"E{tally_ID} {energies_str}\n"
        if add_SD:
            repetitions = len(cells_list) - 1
            rep_str = f"1 {repetitions}R" if repetitions else "1"
            self.other_data[f"SD{tally_ID}"] = f"SD{tally_ID} {rep_str}\n"
        if multiplier is not None:
            self.other_data[f"FM{tally_ID}"] = f"FM{tally_ID} {multiplier}\n"

    def add_stopCard(self, nps: int | float = 1e7):
        """
        This does not append a stop card anymore but simply changes the NPS card. If
        a float is provided, it is converted to an integer.

        Parameters
        ----------
        nps : int | float
            number of particles to simulate. Default 1e7

        Returns
        -------
        None.

        """

        line = "NPS " + str(int(nps)) + " \n"
        self.other_data["NPS"] = line

    def check_range(self, range: list[int], who: str = "cell") -> bool:
        """Check if the provided range is not within the used index, i.e., if the range
        is free. Both cells and surfaces can be checked.

        Parameters
        ----------
        range : list[int]
            list of indices to be checked
        who : str, optional
            either 'surf' or 'cell', by default 'cell'

        Returns
        -------
        bool
            True if the provided range is free, False otherwise

        Raises
        ------
        ValueError
            only surf or cell are accepted as who
        """
        if who == "cell":
            used_index = self.cells.keys()
        elif who == "surf":
            used_index = self.surfs.keys()
        else:
            raise ValueError("who can be only 'cell' or 'surf'")
        # check that the provided range is not within the used index
        # need to do the check for strings since there may be asterisks
        for i in range:
            for j in used_index:
                index = j.replace("*", "")
                index = index.replace("+", "")
                if str(i) == index:
                    return False
        return True

    def delete_fill_cards(self) -> None:
        """Delete all fill cards from the input cells."""
        for _, cell in self.cells.items():
            if cell.fill is not None:
                cell.remove_param("fill")

    def remove_tallies(self, tally_ids: list[int] | None = None) -> None:
        """Remove tallies from the input.

        Parameters
        ----------
        tally_ids : list[int]
            list of tally IDs to be removed. Default is None, which means
            that all tallies will be removed.
        """

        if tally_ids is None:
            # Remove all tally-related cards
            for card in self.other_data:
                if PAT_ALL_TALLY_KEYS.match(card.name):
                    del self.other_data[card.name]
        else:
            # Remove only cards matching the provided tally IDs
            for card in self.other_data:
                m = PAT_ALL_TALLY_KEYS.match(card.name)
                if m:
                    num = int(m.group(2))
                    if num in tally_ids:
                        del self.other_data[card.name]

    def remove_sdef(self) -> None:
        """Remove the SDEF card and related source definition cards from the input."""

        for card in self.other_data:
            key_lower = card.name.lower()
            if key_lower.startswith(("sdef", "kcode", "ssr")) or key_lower[:2] in (
                "si",
                "sd",
                "ds",
                "sp",
            ):
                del self.other_data[card.name]

    def prepare_void_check(
        self, surf: migjorn.Surface, nps: int, particle: str = "N"
    ) -> None:
        """Prepare the input for a void check.

        Parameters
        ----------
        surf : migjorn.Surface
            sphere surface to use as source
        nps : int
            number of particles
        particle : str, optional
            particle type, by default "N"

        Raises
        ------
        ValueError
            if the provided surface is not a sphere
        """
        if surf.kind.lower() not in ["so", "sx", "sy", "sz", "s"]:
            raise ValueError("The provided surface is not a sphere")
        # Add if not already in the model; surface may already be present
        if self._model.surface(surf.id) is None:
            self._model.add_surface(surf.text.strip())
        radius = surf.coeffs[-1]
        weight = np.pi * radius**2
        self.remove_sdef()
        self.remove_tallies(None)
        self.other_data["VOID"] = "VOID\n"
        self.other_data["NPS"] = f"NPS {nps}\n"
        self.other_data["SDEF"] = (
            f"SDEF PAR={particle} NRM=-1 SUR={surf.id} WGT={weight} DIR=d1\n"
        )
        self.other_data["SI1"] = "SI1 0 1\n"
        self.other_data["SP1"] = "SP1 -21 1\n"

    def explore_id_ranges_by_plot(self) -> tuple[Figure, Axes]:
        """
        Returns a Figure and Axes object where the cell and surface IDs used in the
        input are plotted. Useful for interactive exploration of the used IDs.

        Example
        -------
        >>> from f4enix.input.MCNPinput import Input
        ... inp = Input.from_input('input.i')
        ... fig, ax = inp.explore_id_ranges_by_plot()
        ... fig.show()
        """
        cell_ids = {int(x) for x in self.cells}
        surface_ids = {int(x.lstrip("*")) for x in self.surfs}

        fig, ax = plt.subplots(figsize=(12, 3))
        ax.scatter(list(cell_ids), [1] * len(cell_ids), s=10, color="blue")
        ax.scatter(list(surface_ids), [2] * len(surface_ids), s=10, color="red")
        ax.set_xlabel("ID number")
        ax.set_yticks([1, 2])
        ax.set_yticklabels(["Cell IDs", "Surface IDs"])
        ax.grid()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))  # Only integer ticks
        ax.set_title("Occupied IDs")
        return fig, ax

    def find_first_free_id_range(self, required_size: int) -> int:
        """
        Given a required size, it finds the first ID that can accommodate both cell and
        surface ID ranges.
        """
        # Extract combined cell and surface IDs
        cells = {int(x) for x in self.cells}
        surfaces = {int(x.lstrip("*")) for x in self.surfs}
        combined_ids = np.array(list(cells.union(surfaces)))
        combined_ids.sort()

        # Add zero to the beginning in case the first IDs are free
        combined_ids = np.insert(combined_ids, 0, 0)

        # Return the first id that satisfies the required size or the last id + 1
        gaps = np.diff(combined_ids) - 1
        for i, gap in enumerate(gaps):
            if gap >= required_size:
                return combined_ids[i] + 1
        return combined_ids[-1] + 1


class D1S_Input(Input):
    def __init__(
        self,
        *args,
        irrad_file: IrradiationFile | None = None,
        reac_file: ReactionFile | None = None,
    ) -> None:
        """Children of the :py:class:`Input`.

        it includes also the reaction and irradiation files necessary for a
        D1S-UNED run and defines additional methods related to them.

        Parameters
        ----------
        args : list
            list of arguments to be passed to the super constructor.
        irrad_file : IrradiationFile, optional
            irradiation file object, by default None
        reac_file : ReactionFile, optional
            readtion file object, by default None

        Attributes
        ----------
        irrad_file : IrradiationFile
            irradiation file object
        reac_file : ReactionFile
            readtion file object

        Examples
        --------
        translate the input defining an activation and transport library.
        the reaction file will be used to identify to which isotopes the
        activation library has to be assigned.

        >>> from f4enix.input.MCNPinput import D1S_Input
        ... from f4enix.input.libmanager import LibManager
        ... d1s_inp = D1S_Input.from_input('d1stest.i', irrad_file='irr_test',
        ...                                reac_file='reac_fe')
        ... d1s_inp.smart_translate('99c', '00c', LibManager())
        """

        super().__init__(*args)
        self.irrad_file = irrad_file
        self.reac_file = reac_file

    def __deepcopy__(self, memo: dict) -> "D1S_Input":
        new_obj = self.__class__.__new__(self.__class__)
        memo[id(self)] = new_obj
        new_obj.__init__(
            migjorn.Model(self._model.to_source()), deepcopy(self._mat_section, memo)
        )
        new_obj.irrad_file = deepcopy(self.irrad_file, memo)
        new_obj.reac_file = deepcopy(self.reac_file, memo)
        return new_obj

    @classmethod
    def from_input(
        cls,
        inputfile: os.PathLike,
        irrad_file: os.PathLike | None = None,
        reac_file: os.PathLike | None = None,
    ) -> D1S_Input:
        """Generate a D1S-UNED input file.

        this includes also the reaction and irradiation files.

        Parameters
        ----------
        inputfile : os.PathLike
            path to the MCNP input (D1S)
        irrad_file : os.PathLike, optional
            path to the irradiation file, by default None (no file associated)
        reac_file : os.PathLike, optional
            path to the reaction file, by default None (no file associated)

        Returns
        -------
        D1S_Input
            generated D1S_Input object
        """
        name = os.path.basename(str(inputfile)).split(".")[0]
        logging.info(f"Reading file: {name}")
        model = migjorn.Model.from_file(str(inputfile))
        for d in model.diagnostics:
            logging.warning(f"migjorn [{d.severity}]: {d.message}")
        logging.debug("Reading has finished")
        materials = MatCardsList.from_migjorn(model)

        newirrad_file = (
            IrradiationFile.from_text(irrad_file) if irrad_file is not None else None
        )
        newreac_file = (
            ReactionFile.from_text(reac_file) if reac_file is not None else None
        )

        return cls(model, materials, irrad_file=newirrad_file, reac_file=newreac_file)

    def get_potential_paths(self, libmanager: LibManager, lib: str) -> list[Reaction]:
        """Given an activation library, return a list of all possible reactions
        paths foreseen by the libmanager that can originate from the material
        section of the input.

        Parameters
        ----------
        libmanager : LibManager
            Handlers of cross sections operations
        lib : str
            activation library to be used (e.g. 99c)

        Returns
        -------
        list[Reaction]
            list of Reaction objects describing all possible paths included in
            LibManager
        """
        reactions = []
        for material in self.mat_section.materials:
            for zaid in material.zaids:
                parent = str(zaid.nuclide.zaid)
                zaidreactions = libmanager.get_reactions(lib, parent)
                # if len(zaidreactions) > 0:
                #     # it is a parent only if reactions are available
                #     parentlist.append(parent)
                for MT, daughter in zaidreactions:
                    reactions.append((parent, MT, daughter))
                    # daughterlist.append(daughter)

        reactions = list(set(reactions))
        reactions.sort()
        # --- Build the reactions and reaction file ---
        reaction_list = []
        for parent, MT, daughter in reactions:
            parent_nuclide = Nuclide(parent, lib=lib)
            if daughter[-3:] == METASTABLE_TAG:
                daughter_nuclide = Nuclide(daughter[:-3], metastable=True)
            else:
                daughter_nuclide = Nuclide(daughter)
            # Build a comment
            comment = "{} -> {}".format(parent_nuclide, daughter_nuclide)

            rx = Reaction(parent_nuclide, MT, daughter_nuclide, comment=comment)
            reaction_list.append(rx)

        return reaction_list

    def get_reaction_file(
        self, libmanager: LibManager, lib: str, set_as_attribute: bool = True
    ) -> ReactionFile:
        """
        Get a reaction file suitable for the input.

        The reaction file is built selecting from all the possible reaction
        paths that can originate in the model due to its material cards only
        the reactions that lead to a daughter listed in the irradiation file.
        By default this is added as the reac_file for the input.

        Parameters
        ----------
        libmanager : LibManager
            Object handling all cross-sections related operations.
        lib : str
            library suffix to be used.
        set_as_attribute: bool
            if True (default) the reactionfile is saved as the self.reac_file

        Returns
        -------
        ReactionFile
            Object representing the react file for D1S.

        Raises
        ------
        ValueError
            if no irradiation files have been assigned yet to the input

        """
        # irrad file is necessary for this operation
        # recover all available daughters
        if self.irrad_file is None:
            raise ValueError("irrad_file attribute cannot be None for this operation")
        else:
            available_daughters = self.irrad_file.get_daughters()

        # Recover all possible reactions
        reactions = self.get_potential_paths(libmanager, lib)

        # perform the selection
        selected_reactions = []
        for reaction in reactions:
            if reaction.daughter in available_daughters:
                # add the reaction to the one to use
                selected_reactions.append(reaction)

        reac_file = ReactionFile(selected_reactions)

        if set_as_attribute:
            self.reac_file = reac_file

        return reac_file

    def smart_translate(
        self,
        activation_lib: str,
        transport_lib: str,
        libmanager: LibManager,
        fix_natural_zaid: bool = False,
    ) -> None:
        """
        Translate the input to another library without relying on old libs.

        Both the activation and transport libraries are changed. The reaction
        file and PIKMT cards are also translated.
        Differently from self.translate(),
        libraries are not modified based on old ones but based on the reaction
        file. This is, to all parent zaids listed in the reactions will be
        assigned the activation_lib, to all others the transport_lib.

        Parameters
        ----------
        activation_lib : str
            library to be used for activation, e.g., 99c
        activation_lib : dict[str, str]
            library to be used for activation, e.g., 31c
        libmanager : LibManager
            Library manager for the conversion.
        fix_natural_zaid: bool
            if True, and additional initial translation with the transport lib
            is done in order to expand the natural zaids. If the transport lib
            do not expand the natural zaid there may be issues.

        Returns
        -------
        None

        """
        if self.reac_file is None:
            raise ValueError("reac_file cannot be None for this operation")

        if fix_natural_zaid:
            # get a first translation to avoid issues with old natural zaids
            self.translate(transport_lib, libmanager)

        active_zaids = []
        transp_zaids = []

        for reaction in self.reac_file.reactions:
            # strip the lib from the parent
            parent = reaction.parent.zaid
            active_zaids.append(str(parent))
            reaction.change_lib(activation_lib)

        # Now check for the remaing materials in the input to be assigned
        # to transport
        for material in self.mat_section.materials:
            for zaid in material.zaids:
                zaidnum = str(zaid.nuclide.zaid)
                if zaidnum not in active_zaids and zaidnum not in transp_zaids:
                    transp_zaids.append(zaidnum)

        newlib = {activation_lib: active_zaids, transport_lib: transp_zaids}

        # trigger translation of the PIKMT card
        self.add_PIKMT_card()

        # Translate the input with the new lib
        self.mat_section.translate(newlib, libmanager)

    def add_PIKMT_card(self) -> None:
        """
        Add a PIKMT card to the input file. If a PKMT card is already present
        this will be overridden.

        Returns
        -------
        None.

        """
        if self.reac_file is None:
            logging.warning("No reaction file has been assigned to the input")
            return
        key = "PIKMT"
        lines = [key + "\n"]
        for parent in self.reac_file.get_parents():
            lines.append("         {}    {}\n".format(parent.zaid, 0))
        self.other_data[key] = "".join(lines)

    def add_track_contribution(
        self, tallykey: str, bins: list[str], who: str = "parent"
    ):
        """
        Given a list of zaid add the FU bin in the requested tallies in order
        to collect the contribution of them to the tally.

        Parameters
        ----------
        tallykey : str
            ID of the tally onto which to operate (e.g. F4).
        bins : list[str]
            list of zaid/cell numbers of the parents/daughters (e.g. 1001).
        who : str, optional
            either 'parent', 'daughter' or 'cell', specifies the types of bin to
            be tracked. The default is 'parent'.

        Raises
        ------
        ValueError
            check for admissible who parameter.

        """
        existing = self.other_data[tallykey].text
        num = str(_get_num_tally(tallykey))

        existing += "FU" + num + " 0\n"

        if who == "parent":
            for zaid in bins:
                existing += ADD_LINE_FORMAT.format("-" + str(zaid))
        elif who in ["daughter", "cell"]:
            for zaid in bins:
                existing += ADD_LINE_FORMAT.format(zaid)
            if who == "cell":
                self.other_data["FT" + num] = "FT" + num + " SCD\n"
        else:
            raise ValueError(who + ' is not an admissible "who" parameters')
        self.other_data[tallykey] = existing

    def add_daughter_contribution_from_irr(self, tallykey: str):
        """Add the daughter contribution to the tally. All the daughters
        listed in the irradiation file will be added to the tally.

        Parameters
        ----------
        tallykey : str
            ID of the tally onto which to operate (e.g. F4).

        Raises
        ------
        ValueError
            if no irradiation file has been assigned to the input.
        """
        if self.irrad_file is None:
            raise ValueError("No irradiation file has been assigned to the input")
        daughters = self.irrad_file.get_daughters()
        # get the zaid numbers only
        daughters_zaids = []
        for daughter in daughters:
            daughters_zaids.append(daughter.write_to_int_string())
        self.add_track_contribution(tallykey, daughters_zaids, who="daughter")

    def add_parent_contribution_from_reac(self, tallykey: str):
        """Add the parent contribution to the tally. All the parents
        listed in the reaction file will be added to the tally.

        Parameters
        ----------
        tallykey : str
            ID of the tally onto which to operate (e.g. F4).

        Raises
        ------
        ValueError
            if no reaction file has been assigned to the input.
        """

        if self.reac_file is None:
            raise ValueError("No reaction file has been assigned to the input")
        parents = list(self.reac_file.get_parents())
        parents_zaids = []
        for parent in parents:
            parents_zaids.append(parent.write_to_int_string())
        self.add_track_contribution(tallykey, parents_zaids, who="parent")

    def add_SDDR_dose_function(self, tallykey: str) -> None:
        """Add the SDDR dose function to the tally.

        Parameters
        ----------
        tallykey : str
            ID of the tally onto which to operate (e.g. F4).
        """

        existing = self.other_data[tallykey].text
        num = str(_get_num_tally(tallykey))

        existing += "FU" + num + " 0\n"
        self.other_data[tallykey] = existing

        self.other_data[f"DE{num}"] = (
            f"DE{num} 0.01 0.015 0.02 0.03 0.04 0.05\n"
            "        0.06 0.07 0.08 0.10 0.15 0.20\n"
            "        0.3 0.40.5 0.6 0.8 1.0\n"
            "        2.0 4.0 6.0 8.0 10.0\n"
        )
        self.other_data[f"DF{num}"] = (
            f"DF{num} 0.0485 0.1254 0.2050 0.2999 0.3381 0.3572\n"
            "        0.3780 0.4066 0.4399 0.5172 0.7523 1.0041\n"
            "        1.5083 1.9958 2.4657 2.9082 3.7269 4.4834\n"
            "        7.4896 12.0153 15.9873 19.9191 23.7600\n"
        )


def _extract_cell_geometry(cell_text: str) -> str:
    """Extract the geometry portion from a cell card text line."""
    # TODO: migjorn side there is likely a better way to extrac the geometry part
    line = cell_text.splitlines()[0].strip().replace("\r", "")
    tokens = line.split()
    if not tokens:
        return ""
    # Token 0: cell number. Token 1: material. Token 2: density (if non-void).
    start = 1
    if len(tokens) > 1 and re.match(r"^-?\d+$", tokens[1]):
        mat = int(tokens[1])
        if mat == 0:
            start = 2  # void: no density field
        elif len(tokens) > 2 and re.match(r"^-?[\d.eE+]+$", tokens[2]):
            start = 3  # non-void: density at index 2
        else:
            start = 2
    # Collect geometry tokens until a keyword (letters that aren't a surface/complement ref)
    geom_tokens = []
    for tok in tokens[start:]:
        if re.match(r"^[a-zA-Z]{2,}", tok) and not tok.startswith("#"):
            break
        geom_tokens.append(tok)
    return " ".join(geom_tokens)


def _get_num_tally(key: str) -> int:
    patnum = re.compile(r"\d+")
    try:
        num = patnum.search(key).group()
    except AttributeError:
        raise ValueError(key + " is not a valid tally ID")
    return int(num)


def get_formatted_range(numbers: list[int]) -> str:
    """Format a range of integers into a compact string that makes use of the interval
    syntax in MCNP.

    Parameters
    ----------
    current_range : list[int]
        List of cell or surface ids to be formatted.

    Returns
    -------
    str
        Formatted string representing the range of ids.
    """
    numbers = sorted(numbers)
    result = ""
    current_range = [numbers[0]]
    for number in numbers[1:]:
        if number == current_range[-1] + 1:
            current_range.append(number)
        else:
            result += _convert_range_to_string(current_range)
            current_range = [number]
    # Handle the last range
    result += _convert_range_to_string(current_range)
    return result


def _convert_range_to_string(current_range: list[int]) -> str:
    result = ""
    if len(current_range) == 1:
        result += f"{current_range[0]} "
    elif len(current_range) == 2:
        result += f"{current_range[0]} {current_range[1]} "
    elif len(current_range) == 3:
        result += f"{current_range[0]} {current_range[1]} {current_range[2]} "
    else:
        result += f"{current_range[0]} "
        result += f"{current_range[-1] - current_range[0] - 1}I "
        result += f"{current_range[-1]} "
    return result
