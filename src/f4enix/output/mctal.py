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
from math import prod

import numpy as np
import pandas as pd

from f4enix.core.constants import PAT_MISSING_EXP

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
_CARD_KEYS = frozenset(
    [f"{letter}{flag}" for letter in "fdusmcet" for flag in ("", "t", "c")] + ["vals"]
)

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
        with open(filepath, "r") as infile:
            lines = infile.read().splitlines()
        self.tallies = self._read(lines)

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

    def _read(self, lines: list[str]) -> list[Tally]:
        i = self._parse_header(lines)
        tallies = []
        n_lines = len(lines)
        while i < n_lines:
            key = lines[i][:5].lower()
            if key == "tally":
                tally, i = self._parse_tally(lines, i)
                tallies.append(tally)
            elif key == "kcode":
                break
            else:
                # TFC block, not stored
                i += 1
        return tallies

    def _parse_header(self, lines: list[str]) -> int:
        header = self.header
        tokens = lines[0].split()
        # kod, ver and probid are optional, knod nps rnr are always present
        if len(tokens) >= 3:
            header.knod, header.nps, header.rnr = (int(t) for t in tokens[-3:])
            head = tokens[:-3]
            if len(head) >= 2:
                header.kod, header.ver = head[0], head[1]
                header.probid = " ".join(head[2:])
        header.title = lines[1].strip()

        tokens = lines[2].split()
        header.ntal = int(tokens[1])
        if len(tokens) >= 4:
            header.npert = int(tokens[3])
        if header.ntal == 0:
            raise ValueError(f"No tallies in the MCTAL file {self.mctalFileName}")
        if header.npert != 0:
            raise NotImplementedError(
                "MCTAL files with perturbation cards are not supported"
            )

        i = 3
        ids = []
        while i < len(lines) and lines[i][:5].lower() != "tally":
            ids.extend(lines[i].split())
            i += 1
        header.ntals = np.array(ids, dtype=int)
        return i

    def _parse_tally(self, lines: list[str], i: int) -> tuple[Tally, int]:
        tokens = lines[i].split()
        tally = Tally(int(tokens[1]))
        tally.typeNumber = int(tokens[2])
        if len(tokens) > 3:
            tally.detectorType = int(tokens[3])
        tally.mesh = tally.detectorType <= -1
        tally.radiograph = tally.detectorType >= 3
        i += 1

        if tally.typeNumber < 0:
            tally.tallyParticles = np.array(lines[i].split(), dtype=int)
            i += 1

        while not _is_card(lines[i]):
            tally.tallyComment.append(lines[i][5:].rstrip())
            i += 1

        # {letter: (card tokens, total/cumulative flag, listed values)}
        cards: dict[str, tuple[list[str], str, list[str]]] = {}
        while True:
            tokens = lines[i].split()
            key = tokens[0].lower()
            i += 1
            if key == "vals":
                break
            values = []
            while not _is_card(lines[i]):
                values.extend(lines[i].split())
                i += 1
            cards[key[0]] = (tokens, key[1:], values)

        tally.axes = _build_axes(tally, cards)
        n_tokens = 2 * prod(size for _, _, size in tally.axes)

        tokens = []
        while len(tokens) < n_tokens:
            tokens.extend(lines[i].split())
            i += 1
        if len(tokens) != n_tokens:
            raise ValueError(
                f"Unexpected number of values in tally {tally.tallyNumber} of "
                f"{self.mctalFileName}: expected {n_tokens}, found {len(tokens)}"
            )
        tally.values = _to_floats(tokens).reshape(-1, 2)
        if np.isnan(tally.values).any():
            self.thereAreNaNs = True

        return tally, i

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


def _is_card(line: str) -> bool:
    # card keywords start at column 1, comments and values are indented
    if not line or line[0] == " ":
        return False
    return line.split(maxsplit=1)[0].lower() in _CARD_KEYS


def _to_floats(tokens: list[str]) -> np.ndarray:
    try:
        return np.array(tokens, dtype=float)
    except ValueError:
        return np.array([PAT_MISSING_EXP.sub(r"E\1", t) for t in tokens], dtype=float)


def _to_label_array(labels: list) -> np.ndarray:
    types = {type(label) for label in labels}
    if types <= {int}:
        return np.array(labels, dtype=np.int64)
    if types <= {int, float}:
        return np.array(labels, dtype=float)
    array = np.empty(len(labels), dtype=object)
    array[:] = labels
    return array


def _build_axes(
    tally: Tally, cards: dict[str, tuple[list[str], str, list[str]]]
) -> list[Axis]:
    def card(letter: str) -> tuple[int, str, list[str], list[str]]:
        tokens, flag, values = cards.get(letter, ([letter, "0"], "", []))
        return int(tokens[1]), flag, values, tokens

    axes: list[Axis] = []

    n, _, values, tokens = card("f")
    if tally.mesh:
        # f card: n_bins, unknown, n_cora, n_corb, n_corc; values are the
        # concatenated bin boundaries of the three axes
        dims = [int(t) for t in tokens[3:6]]
        coords = _to_floats(values).tolist()
        edges = []
        start = 0
        for dim in dims:
            # lower edges are used as labels
            edges.append(coords[start : start + dim])
            start += dim + 1
        for name, labels in zip(("Cor C", "Cor B", "Cor A"), edges[::-1]):
            axes.append((name, labels, len(labels)))
    else:
        size = max(n, 1)
        if values:
            labels = [_cell_label(value, idx) for idx, value in enumerate(values)]
        else:
            # e.g. detectors, where the cells are not listed
            labels = list(range(1, size + 1))
        _check_size(tally, "Cells", labels, size)
        axes.append(("Cells", labels, size))

    n, _, _, _ = card("d")
    axes.append(("Dir", list(range(n)) if n > 0 else None, max(n, 1)))

    n, flag, values, _ = card("u")
    axes.append(_listed_axis(tally, "User", n, flag, values))

    n, flag, _, _ = card("s")
    size = max(n, 1)
    n_real = size - 1 if flag == "t" else size
    labels = list(range(1, n_real + 1)) + [TOTAL] * (size - n_real)
    axes.append(("Segments", labels, size))

    n, flag, _, _ = card("m")
    if n > 0:
        n_real = n - 1 if flag == "t" else n
        labels = list(range(n_real)) + [TOTAL] * (n - n_real)
        axes.append(("Multiplier", labels, n))
    else:
        axes.append(("Multiplier", None, 1))

    n, flag, values, _ = card("c")
    if tally.radiograph:
        # the listed values are the n+1 boundaries of the t-axis grid
        size = max(n, 1)
        labels = _to_floats(values[:size]).tolist() or None
        axes.append(("Cosine", labels, size))
    else:
        axes.append(_listed_axis(tally, "Cosine", n, flag, values))

    n, flag, values, _ = card("e")
    axes.append(_listed_axis(tally, "Energy", n, flag, values))

    n, flag, values, _ = card("t")
    axes.append(_listed_axis(tally, "Time", n, flag, values))

    return axes


def _listed_axis(tally: Tally, name: str, n: int, flag: str, values: list[str]) -> Axis:
    size = max(n, 1)
    labels = _to_floats(values).tolist()
    # the total bin has no listed boundary
    if flag == "t" and len(labels) == size - 1:
        labels.append(TOTAL)
    if not labels:
        return (name, None, 1) if size == 1 else (name, list(range(size)), size)
    _check_size(tally, name, labels, size)
    return (name, labels, size)


def _cell_label(value: str, idx: int) -> int | float | str:
    number = float(value)
    if number == 0:
        # a zero stands for a union of cells/surfaces
        return f"Input {idx + 1}"
    if number.is_integer():
        return int(number)
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
