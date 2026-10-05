"""Regression tests comparing the new MCTAL parser with the legacy one."""

from importlib.resources import as_file, files

import numpy as np
import pandas as pd
import pytest

import tests.resources.mctal as mctal_res
from f4enix.output.mctal import TOTAL, Mctal
from tests.legacy.mctal_legacy import Mctal as LegacyMctal

MCTAL_RESOURCES = files(mctal_res)

FIXTURES = [
    "test_m",
    "error_summary.m",
    "C_Modelm",
    "mctal",
    "mctal_time",
    "mctal_tmesh",
    "mctal_cosbin",
    "mctal_radio",
    "detectors.m",
    "mctal_fm",
    "mctal_daughter",
    "mctal_daughter_2",
]
# legacy parser only kept 10 of the 12 radiograph t-axis bins
INCOMPLETE_LEGACY = {"mctal_radio"}


def _labels(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in ("Value", "Error")]


def _keys(df: pd.DataFrame, label_cols: list[str]) -> pd.Series:
    return df[label_cols].astype(str).agg("|".join, axis=1)


def _normalize_new(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "Cells" in df.columns:
        # legacy truncated non-integer cell labels (e.g. macrobody facets)
        df["Cells"] = [
            int(c) if isinstance(c, float) else c for c in df["Cells"].tolist()
        ]
    return df


def _normalize_legacy(df: pd.DataFrame) -> pd.DataFrame:
    # legacy swapped the Cor A and Cor C column names
    return df.rename(columns={"Cor A": "Cor C", "Cor C": "Cor A"})


@pytest.mark.parametrize("filename", FIXTURES)
def test_regression(filename):
    with as_file(MCTAL_RESOURCES.joinpath(filename)) as inp:
        new = Mctal(inp)
        legacy = LegacyMctal(inp)

    assert list(new.tallydata) == list(legacy.tallydata)

    for tally_number, legacy_df in legacy.tallydata.items():
        raw_new_df = new.tallydata[tally_number]
        assert _keys(raw_new_df, _labels(raw_new_df)).is_unique, tally_number
        new_df = _normalize_new(raw_new_df)
        legacy_df = _normalize_legacy(legacy_df)

        assert set(new_df.columns) == set(legacy_df.columns), tally_number
        label_cols = _labels(new_df)

        new_df = new_df.assign(key=_keys(new_df, label_cols))
        legacy_df = legacy_df.assign(key=_keys(legacy_df, label_cols))
        if not legacy_df["key"].is_unique:
            # truncated legacy labels are ambiguous, compare by position
            assert len(new_df) == len(legacy_df), tally_number
            assert np.allclose(new_df["Value"], legacy_df["Value"]), tally_number
            assert np.allclose(new_df["Error"], legacy_df["Error"]), tally_number
            continue

        merged = legacy_df.merge(
            new_df[["key", "Value", "Error"]],
            on="key",
            how="left",
            suffixes=("_old", "_new"),
            indicator=True,
        )
        assert (merged["_merge"] == "both").all(), tally_number

        # legacy set to zero values missing the exponent 'E' (e.g. 1.6-113)
        missing_exp = (
            (merged["Value_old"] == 0)
            & (merged["Error_old"] == 0)
            & (merged["Value_new"].abs() < 1e-99)
        )
        ok = merged[~missing_exp]
        assert np.allclose(ok["Value_old"], ok["Value_new"]), tally_number
        assert np.allclose(ok["Error_old"], ok["Error_new"]), tally_number

        # the new parser may only add total rows that legacy dropped
        extra = new_df[~new_df["key"].isin(legacy_df["key"])]
        if filename not in INCOMPLETE_LEGACY:
            has_total = (extra[label_cols] == TOTAL).any(axis=1)
            assert has_total.all(), tally_number
