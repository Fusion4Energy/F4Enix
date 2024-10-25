from importlib.resources import as_file, files

import pandas as pd

import f4enix.resources as pkg_res
from f4enix.material_library import AVAILABLE_MATERIALS, MaterialComposition

FIRST_WALL = "First Wall (315g)"
EQ_PORT = "Eq. port interspace (175g)"
N17 = "N-17 neutron emission (175g)"

AVAILABLE_SPECTRA = [FIRST_WALL, EQ_PORT, N17]


class PathwayLibrary:
    def __init__(self, df: pd.DataFrame | None = None):
        if df is None:
            resources = files(pkg_res)
            with as_file(resources.joinpath("global_filterable.csv")) as infile:
                df = pd.read_csv(infile)

        self.library = df.set_index(
            [
                "dose",
                "spectrum",
                "element",
                "isotope",
                "state",
                "pathway",
                "material",
            ],
        )

    def get_pathways(
        self,
        spectrum: list[str] | None = None,
        dose: int | None = None,
        materials: list[MaterialComposition] | None = None,
    ) -> pd.DataFrame:
        if spectrum is None:
            spectrum = AVAILABLE_SPECTRA
        else:
            # check that requested spectrum is available
            for s in spectrum:
                if s not in AVAILABLE_SPECTRA:
                    raise ValueError(f"Spectrum {s} is not available")

        if dose is None or dose == 95:
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

        # I need the material names to check the df columns
        mat_names = [m.name for m in materials]

        # Select all the rows that match the requested spectrum, dose and materials
        df = self.library.copy()
        df = df.loc[doselabel, spectrum, :, :, :, :, mat_names]

        dfmax = (
            df.reset_index()
            .groupby(
                [
                    "element",
                    "isotope",
                    "state",
                    "pathway",
                ]
            )
            .max()
        )
        del dfmax["dose"]
        del dfmax["material"]
        del dfmax["spectrum"]

        return dfmax
