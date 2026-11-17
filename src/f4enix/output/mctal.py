"""
This module deals with the parsing of MCNP MCTAL files
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

import logging
import os

import aci_tallies
import numpy as np
import pandas as pd
from aci_tallies import mctal as aci_mctal

TOTAL = "total"
# Dataframe column order (the order of the VALS array is given by Tally.axes)
COLUMNS = [
    "Cells",
    "Dir",
    "User",
    "Segments",
    "Multiplier",
    "Cosine",
    "Energy",
    "Time",
    "Cor C",
    "Cor B",
    "Cor A",
]

# (column name, labels, size); labels is None for axes that are not binned
Axis = tuple[str, "list | None", int]


class Header:
    def __init__(self) -> None:
        """General information on the mctal file.

        Attributes
        ----------
        kod : str
            name of the code that was used
        ver : str
            name of the version that was used
        probid : str
            date and time when the problem was run
        knod : int
            dump number
        nps : int
            number of histories that were run
        rnr : int
            number of pseudorandom numbers that were used
        title : str
            problem identification line
        ntal : int
            number of tallies present in the file
        ntals : np.ndarray
            array of tally numbers
        npert : int
            number of perturbations
        """
        self.kod = ""
        self.ver = ""
        self.probid = ""
        self.knod = 0
        self.nps = 0
        self.rnr = 0
        self.title = ""
        self.ntal = 0
        self.ntals = np.array((), dtype=int)
        self.npert = 0


class Tally:
    def __init__(self, tN: int) -> None:
        """Raw content of a single tally of the mctal file.

        Attributes
        ----------
        tallyNumber : int
            tally number
        typeNumber : int
            particle type. If negative, the particle list is in tallyParticles
        detectorType : int
            0=none, 1=point, 2=ring, 3=pinhole radiograph, 4=rectangular
            radiograph, 5=cylindrical radiograph, negative for mesh tallies
        tallyParticles : np.ndarray
            0/1 flags indicating which particle types are used by the tally
        tallyComment : list[str]
            the FC card lines
        mesh : bool
            True if the tally is a mesh tally
        radiograph : bool
            True if the tally is a radiograph tally
        axes : list[tuple[str, list | None, int]]
            (name, labels, size) of each binning, from the slowest to the
            fastest varying one in the VALS array. labels is None if the
            binning is not used.
        values : np.ndarray
            (n_bins, 2) array of values and relative errors
        """
        self.tallyNumber = tN
        self.typeNumber = 0
        self.detectorType = 0
        self.tallyParticles = np.array((), dtype=int)
        self.tallyComment: list[str] = []
        self.mesh = False
        self.radiograph = False
        self.axes: list[Axis] = []
        self.values = np.empty((0, 2))


class Mctal:
    def __init__(self, filepath: os.PathLike | str) -> None:
        """Object responsible for the parsing of MCNP mctal files.

        Parameters
        ----------
        filepath : os.PathLike | str
            path to the mctal file to be parsed

        Attributes
        ----------
        tallydata : dict[int, pd.DataFrame]
            dictionary that at each tally id associate a pandas dataframe
            containing the results. It supports multi-binning tallies.
        totalbin : dict[int, pd.DataFrame | None]
            for each tally, the subset of tallydata rows that contain at least
            one total bin. None if the tally has no total bins.
        header : Header
            it is the parsed data of the mctal file. See the Header doc to
            understand how to access the data.
        tallies : list[Tally]
            raw parsed tallies.

        Examples
        --------
        Parse the mctal file and access the data
        >>> # Import the mctal module
        ... from f4enix.output.mctal import Mctal
        ... # Parse the Mctal file
        ... file = 'mctal'
        ... mctal = Mctal(file)
        ... # get a summary of the min and max errors across tallies
        ... mctal.get_error_summary().sort_values(by='tally num')
            tally num	min rel error	max rel error
        0	    4	        NaN	            1.0000
        22	    6	        0.0007	        0.0272
        19	    14	        NaN	            1.0000
        1	    16	        0.0008	        0.0381
        """
        self.mctalFileName = os.path.basename(filepath)
        self.header = Header()
        self.thereAreNaNs = False

        logging.info("Parsing file: %s", self.mctalFileName)
        self.tallies = self._read(filepath)

        if self.thereAreNaNs:
            logging.warning(
                "The MCTAL file %s contains tallies with NaN values",
                self.mctalFileName,
            )

        self.tallydata, self.totalbin = self._get_dfs()

    def get_error_summary(self, include_abs_err: bool = False) -> pd.DataFrame:
        """Return a dataframe containing a summary of the min and max errror
        registered in each tally. If both value and error are equal to zero,
        the errors will be set to NaN, since it means that nothing has been
        scored in the tally.

        Parameters
        ----------
        include_abs_err : bool, optional
            if True includes the absolute error in addition to the total one,
            by default is False

        Returns
        -------
        pd.DataFrame
            error summary
        """
        rows = []
        for tally, data in self.tallydata.items():
            min_error = data["Error"].min()
            min_idx = data["Error"].idxmin()
            max_error = data["Error"].max()
            max_idx = data["Error"].idxmax()

            min_val = data["Value"].iloc[min_idx]
            max_val = data["Value"].iloc[max_idx]

            # put a NaN, since no particle was tallied in the cell
            if min_error == 0 and min_val == 0:
                min_error = np.nan
            if max_error == 0 and max_val == 0:
                max_error = np.nan

            if include_abs_err:
                rows.append([tally, min_error, min_val, max_error, max_val])
            else:
                rows.append([tally, min_error, max_error])

        if include_abs_err:
            columns = [
                "tally num",
                "min rel error",
                "min abs err",
                "max rel error",
                "max abs err",
            ]
        else:
            columns = ["tally num", "min rel error", "max rel error"]

        df = pd.DataFrame(rows)
        df.columns = columns

        return df

    def remove_totals(self, tally_number: int) -> None:
        """
        Remove all rows containing 'total' in any column from the tally DataFrame
        for the specified tally number in self.tallydata.

        Parameters
        ----------
        tally_number : int
            The tally number whose DataFrame should be cleaned.
        """
        if tally_number not in self.tallydata:
            raise KeyError(f"Tally number {tally_number} not found in tallydata.")
        df = self.tallydata[tally_number]
        # Remove rows where any column contains 'total'
        mask = ~df.apply(
            lambda row: (
                row.astype(str).str.contains("total", case=False, na=False).any()
            ),
            axis=1,
        )
        self.tallydata[tally_number] = df[mask].reset_index(drop=True)

    def set_d1s_relative_contribution(
        self, tally_number: int, user_label: str = "Daughter"
    ) -> None:
        """
        For the given tally number, set the User column to int, rename it to user_label (Parent/Daughter/Cell),
        and for each unique bin defined by the other columns, add a 'Normalized Value' column representing
        the relative contribution of each user (Parent/Daughter/Cell) in that bin.

        Parameters
        ----------
        tally_number : int
            The tally number whose DataFrame should be processed.
        user_label : str
            The new name for the User column (e.g., 'Parent', 'Daughter', 'Cell').
        """
        if tally_number not in self.tallydata:
            raise KeyError(f"Tally number {tally_number} not found in tallydata.")

        self.remove_totals(tally_number)

        df = self.tallydata[tally_number].copy()
        if "User" not in df.columns:
            raise ValueError("No 'User' column found in tally DataFrame.")
        df["User"] = df["User"].astype(int)
        df = df.rename(columns={"User": user_label})
        df = normalize_tally(df, user_label=user_label, inplace=False)

        self.tallydata[tally_number] = df

    def _read(self, filepath: os.PathLike | str) -> list[Tally]:
        file = aci_mctal.open(filepath)
        if file.failed:
            id, why = file.failed[0]
            raise ValueError(
                f"Tally {id} of the MCTAL file {self.mctalFileName} could not be "
                f"read: {why}"
            )
        self._read_header(file)
        tallies = []
        for aci_tally in file.read_all():
            for warning in aci_tally.warnings:
                logging.warning(
                    "MCTAL file %s, tally %s: %s",
                    self.mctalFileName,
                    aci_tally.id,
                    warning,
                )
            tally = _tally_from_aci(aci_tally)
            if np.isnan(tally.values).any():
                self.thereAreNaNs = True
            tallies.append(tally)
        return tallies

    def _read_header(self, file: aci_mctal.MctalFile) -> None:
        header = self.header
        header.kod, header.ver = file.code
        header.probid = file.probid
        header.knod = file.dump
        header.nps = file.histories
        header.rnr = file.randoms
        header.title = file.title
        header.ntals = np.array(file.declared_ids, dtype=int)
        header.ntal = len(header.ntals)
        header.npert = file.perturbations
        if header.ntal == 0:
            raise ValueError(f"No tallies in the MCTAL file {self.mctalFileName}")
        if header.npert != 0:
            raise NotImplementedError(
                "MCTAL files with perturbation cards are not supported"
            )

    def _get_dfs(
        self,
    ) -> tuple[dict[int, pd.DataFrame], dict[int, pd.DataFrame | None]]:
        tallydata = {}
        totalbin = {}

        for tally in self.tallies:
            n_rows = len(tally.values)
            columns = {}
            total_mask = np.zeros(n_rows, dtype=bool)
            # number of rows spanned by one bin of the current axis
            block = n_rows
            for name, labels, size in tally.axes:
                block //= size
                # constant binnings are kept only for single-value tallies
                if labels is None or (size == 1 and n_rows > 1):
                    continue
                column = np.tile(
                    np.repeat(_to_label_array(labels), block),
                    n_rows // (size * block),
                )
                if column.dtype == object:
                    total_mask |= column == TOTAL
                columns[name] = column

            df = pd.DataFrame(
                {name: columns[name] for name in COLUMNS if name in columns}
            )
            df["Value"] = tally.values[:, 0]
            df["Error"] = tally.values[:, 1]

            tallydata[tally.tallyNumber] = df
            totalbin[tally.tallyNumber] = df[total_mask] if total_mask.any() else None

        return tallydata, totalbin


def _to_label_array(labels: list) -> np.ndarray:
    types = {type(label) for label in labels}
    if types <= {int}:
        return np.array(labels, dtype=np.int64)
    if types <= {int, float}:
        return np.array(labels, dtype=float)
    array = np.empty(len(labels), dtype=object)
    array[:] = labels
    return array


# MCTAL detector codes and MCNP 6 particles, in the order of the flags line
_DETECTOR = {
    None: 0,
    "point": 1,
    "ring": 2,
    "pinhole image": 3,
    "rectangular image": 4,
    "cylindrical image": 5,
}
_PARTICLES = "npe|quvfhl+-xyo!<>g/zk%^b_~cw@dtsa*?#"
# MCNP writes FC comments in lines of 75 columns
_COMMENT_WIDTH = 75


def _tally_from_aci(aci_tally: aci_mctal.Tally) -> Tally:
    tally = Tally(aci_tally.id)
    # MCNP 6 form: minus the number of particles, which the flags line gives
    tally.typeNumber = -len(aci_tally.particles)
    tally.tallyParticles = np.array(
        [int(p in aci_tally.particles) for p in _PARTICLES], dtype=int
    )
    mesh = aci_tally.mesh
    # negative for a TMESH mesh tally: minus its geometry (1 rec, 2 cyl, 3 sph)
    tally.detectorType = -mesh["geometry"] if mesh else _DETECTOR[aci_tally.detector]
    tally.mesh = tally.detectorType <= -1
    tally.radiograph = tally.detectorType >= 3
    comment = aci_tally.comment or ""
    tally.tallyComment = [
        comment[i : i + _COMMENT_WIDTH].rstrip()
        for i in range(0, len(comment), _COMMENT_WIDTH)
    ]
    tally.axes = _build_axes(tally, aci_tally)
    # aci-tallies keeps MCNP's eight axes in the order of the VALS array
    tally.values = np.column_stack(
        [aci_tally.result.reshape(-1), aci_tally.rel_error.reshape(-1)]
    )
    return tally


def _listed(axis: aci_tallies.Axis) -> list[float]:
    """The values listed on a bin line of the MCTAL file."""
    if axis.bounds is not None:
        return axis.bounds.tolist()
    if axis.upper is not None:
        return axis.upper.tolist()
    return []


def _build_axes(tally: Tally, aci_tally: aci_mctal.Tally) -> list[Axis]:
    aci_axes = aci_tally.axes
    axes: list[Axis] = []

    if tally.mesh:
        grid = aci_tally.mesh
        if grid is None:
            raise ValueError(f"Tally {tally.tallyNumber}: unreadable TMESH grid")
        # the lower edges of the three axes are used as labels
        edges = [bounds[:-1].tolist() for bounds in grid["bounds"]]
        for name, labels in zip(("Cor C", "Cor B", "Cor A"), edges[::-1]):
            axes.append((name, labels, len(labels)))
    else:
        size = len(aci_axes["f"])
        numbers = aci_tally.numbers
        if numbers:
            facets = aci_tally.facets or [None] * len(numbers)
            labels = [
                _cell_label(number, facet, idx)
                for idx, (number, facet) in enumerate(zip(numbers, facets))
            ]
        else:
            # e.g. detectors, where the cells are not listed
            labels = list(range(1, size + 1))
        _check_size(tally, "Cells", labels, size)
        axes.append(("Cells", labels, size))

    n = len(aci_axes["d"])
    axes.append(("Dir", list(range(n)), n))

    axes.append(_listed_axis(tally, "User", aci_axes["u"]))

    axis = aci_axes["s"]
    size = len(axis)
    n_real = size - 1 if axis.total else size
    labels = list(range(1, n_real + 1)) + [TOTAL] * (size - n_real)
    axes.append(("Segments", labels, size))

    axis = aci_axes["m"]
    n = len(axis)
    # a single multiplier bin (no FM card) is not binned
    if n > 1 or axis.total:
        n_real = n - 1 if axis.total else n
        labels = list(range(n_real)) + [TOTAL] * (n - n_real)
        axes.append(("Multiplier", labels, n))
    else:
        axes.append(("Multiplier", None, 1))

    axis = aci_axes["c"]
    if tally.radiograph:
        # the listed values are the n+1 boundaries of the t-axis grid
        size = len(axis)
        labels = _listed(axis)[:size] or None
        axes.append(("Cosine", labels, size))
    else:
        axes.append(_listed_axis(tally, "Cosine", axis))

    axes.append(_listed_axis(tally, "Energy", aci_axes["e"]))
    axes.append(_listed_axis(tally, "Time", aci_axes["t"]))

    return axes


def _listed_axis(tally: Tally, name: str, axis: aci_tallies.Axis) -> Axis:
    size = len(axis)
    labels: list = _listed(axis)
    # the total bin has no listed boundary
    if axis.total and len(labels) == size - 1:
        labels.append(TOTAL)
    if not labels:
        return (name, None, 1) if size == 1 else (name, list(range(size)), size)
    _check_size(tally, name, labels, size)
    return (name, labels, size)


def _cell_label(number: int | None, facet: int | None, idx: int) -> int | float | str:
    if number is None:
        # a zero stands for a union of cells/surfaces
        return f"Input {idx + 1}"
    if facet is not None:
        # macrobody facet, e.g. 2.2
        return float(f"{number}.{facet}")
    return number


def _check_size(tally: Tally, name: str, labels: list, size: int) -> None:
    if len(labels) != size:
        raise ValueError(
            f"Tally {tally.tallyNumber}: {len(labels)} {name} bins found, "
            f"{size} expected"
        )


def normalize_tally(
    df: pd.DataFrame, user_label: str = "Daughter", inplace: bool = True
) -> pd.DataFrame:
    """
    Adds a 'Normalized Value' column to the DataFrame, representing the relative
    contribution of each user_label in each unique bin defined by the other columns.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame, must contain a column named as user_label and 'Value'.
    user_label : str, optional
        The name of the user column to normalize over (default is "User").
    inplace : bool, optional
        If True, modifies the DataFrame in place. If False, returns a new DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with an added 'Normalized Value' column.
    """
    if user_label not in df.columns:
        raise ValueError(f"No '{user_label}' column found in DataFrame.")
    if "Value" not in df.columns:
        raise ValueError("No 'Value' column found in DataFrame.")

    # Identify columns to group by (all except user_label, Value, Error, Normalized Value)
    group_cols = [
        col
        for col in df.columns
        if col not in [user_label, "Value", "Error", "Normalized Value"]
    ]

    if not inplace:
        df = df.copy()

    if len(group_cols) == 0:
        # No group keys: normalize the whole column
        total = df["Value"].sum()
        df["Normalized Value"] = df["Value"] / total if total != 0 else 0
        return df
    else:
        norm_factor = df.groupby(group_cols)["Value"].transform("sum")
        df["Normalized Value"] = df["Value"] / norm_factor.where(norm_factor != 0, 0)
        return df
