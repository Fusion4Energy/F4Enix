"""Module related to the collection of decay pathways for ITER"""

from __future__ import annotations

from importlib.resources import as_file, files

import pandas as pd

import f4enix.resources as pkg_res
from f4enix.input.libmanager import LibManager
from f4enix.material_library import AVAILABLE_MATERIALS, MaterialComposition


def _sort_df_byzaidnum(df: pd.DataFrame, reset_index=True) -> None:
    lm = LibManager()
    if reset_index:
        # store original index to be reapplied later
        index = df.index.names
        df.reset_index(inplace=True)
    df["zaidnum"] = df.apply(
        lambda row: lm.get_zaidnum(row["element"] + str(row["isotope"])), axis=1
    )
    df.sort_values(by="zaidnum", inplace=True)
    del df["zaidnum"]
    if reset_index:
        df.set_index(index, inplace=True)


FIRST_WALL = "First Wall (315g)"
EQ_PORT = "Eq. port interspace (175g)"
N17 = "N-17 neutron emission (175g)"

AVAILABLE_SPECTRA = [FIRST_WALL, EQ_PORT, N17]
AVAILABLE_COOLING_TIMES = [
    "24h",
    "11.6d",
    "30d",
    "180d",
    "1y",
    "10y",
]


class PathwayLibrary:
    def __init__(self, df: pd.DataFrame | None = None):
        """Create a library of decay pathways for ITER

        The object allows to filter the pathways by any combinations of spectrum,
        dose and materials and to get an ordered dataframe that contains only the
        pathways that contribute to the requested dose, for the requested materials,
        spectra and cooling times.

        Parameters
        ----------
        df : pd.DataFrame | None, optional
            Global library of pathways, by default None
        """
        if df is None:
            resources = files(pkg_res)
            with as_file(resources.joinpath("global_filterable.csv")) as infile:
                df = pd.read_csv(infile)

        df["isotope"] = df["isotope"].astype(str) + df["state"].astype(str)
        df["isotope"] = df["isotope"].str.replace("nan", "")
        del df["state"]

        self.library = df.set_index(
            [
                "dose",
                "spectrum",
                "element",
                "isotope",
                # "state",
                "pathway",
                "material",
            ],
        )

    def get_pathways(
        self,
        spectrum: list[str] | None = None,
        dose: int = 95,
        materials: list[MaterialComposition] | None = None,
        cooling_times: list[str] | None = None,
    ) -> pd.DataFrame:
        """Get an ordered datafrme that contains only the pathways that contribute
        to the requested dose, for the requested materials, spectra and cooling times.
        When a pathway is important for more than one material, the maximum dose
        reported in the table will refer to the material with the highest dose.

        Parameters
        ----------
        spectrum : list[str] | None, optional
            allowed spectra are available at f4enix.decay_pathways.AVAILABLE_SPECTRA,
            by default None. If None, all spectra are selected
        dose : int | None, optional
            select the decay pathways that contribute to either 95 or 99 percent,
            by default 95.
        materials : list[MaterialComposition] | None, optional
            available materials are at f4enix.material_library.AVAILABLE_MATERIALS,
            by default None. If None, all materials are selected
        cooling_times : list[str] | None, optional
            select the decay pathways that contribute to the requested cooling times,
            allowed times are available at f4enix.decay_pathways.AVAILABLE_COOLING_TIMES,
            by default None. If None, all cooling times are selected

        Returns
        -------
        pd.DataFrame
            resulting summary of important decay pathways for the requested dose,
            materials, spectra and cooling times.
        """
        # filter the pathways
        df = self.filter_pathways(
            spectrum=spectrum,
            dose=dose,
            materials=materials,
            cooling_times=cooling_times,
        )
        # retain only the maximum dose for each pathway
        dfmax = (
            df.reset_index()
            .groupby(
                [
                    "element",
                    "isotope",
                    "pathway",
                ]
            )
            .max()
        )
        del dfmax["dose"]
        del dfmax["material"]
        del dfmax["spectrum"]
        _sort_df_byzaidnum(dfmax)

        return dfmax

    def filter_pathways(
        self,
        spectrum: list[str] | None = None,
        dose: int = 95,
        materials: list[MaterialComposition] | None = None,
        cooling_times: list[str] | None = None,
    ) -> pd.DataFrame:
        """Filter the pathways by any combinations of spectrum, dose and materials

        Parameters
        ----------
        spectrum : list[str] | None, optional
            allowed spectra are available at f4enix.decay_pathways.AVAILABLE_SPECTRA,
            by default None. If None, all spectra are selected
        dose : int | None, optional
            select the decay pathways that contribute to either 95 or 99 percent,
            by default 95.
        materials : list[MaterialComposition] | None, optional
            available materials are at f4enix.material_library.AVAILABLE_MATERIALS,
            by default None. If None, all materials are selected
        cooling_times : list[str] | None, optional
            select the decay pathways that contribute to the requested cooling times,
            allowed times are available at f4enix.decay_pathways.AVAILABLE_COOLING_TIMES,
            by default None. If None, all cooling times are selected

        Returns
        -------
        pd.DataFrame
            filtereted dataframe containing only the rows related to the subset of
            materials, spectra and dose requested.

        Raises
        ------
        ValueError
            if either a spectrum, cooling time or material is requested that is not available
        """
        if spectrum is None:
            spectrum = AVAILABLE_SPECTRA
        else:
            # check that requested spectrum is available
            for s in spectrum:
                if s not in AVAILABLE_SPECTRA:
                    raise ValueError(f"Spectrum {s} is not available")

        if dose == 95:
            doselabel = "Dose 95%"
        elif dose == 99:
            doselabel = "Dose 99%"
        else:
            raise ValueError("only 95% and 99% doses are available")

        if materials is None:
            materials = AVAILABLE_MATERIALS
        else:
            # check that requested materials are available
            for m in materials:
                if m not in AVAILABLE_MATERIALS:
                    raise ValueError(f"Material {m} is not available")

        if cooling_times is None:
            cooling_times = AVAILABLE_COOLING_TIMES
        else:
            # check that requested cooling times are available
            for t in cooling_times:
                if t not in AVAILABLE_COOLING_TIMES:
                    raise ValueError(f"Cooling time {t} is not available")

        # I need the material names to check the df columns
        mat_names = [m.name for m in materials]
        cooling_times = [f"max dose % at {t}" for t in cooling_times]

        # Select all the rows that match the requested spectrum, dose and materials
        df = self.library.copy()
        df = df.loc[doselabel, spectrum, :, :, :, mat_names]

        return df[cooling_times].dropna(how="all")
