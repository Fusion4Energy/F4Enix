"""Parse legacy FISPACT-II output files to extract the generic pathways."""

from __future__ import annotations

import re
import os
import pandas as pd
import pypact as pp
from dataclasses import dataclass
from pathlib import Path
from f4enix.input.libmanager import LibManager
from f4enix.core.constants import PathLike


perc_pattern = re.compile(r"\d+\.*\d*%")
target_pathway_zaids = re.compile(r"[A-Z][a-z]*\s*\d+m*")
metastable_pat = re.compile(r"\d+m")
isotope_pat = re.compile(r"\d+")
element_pat = re.compile(r"[a-zA-Z]+")
path_id_pat = re.compile(r"\s*path\s+\d+")
reaction_pat = re.compile(r"\([a-zA-Z\d,+-]+\)")
PATHWAY_END = ["(S)", "(L)"]


@dataclass
class FispactZaid:
    """Class to represent a Zaid in a generic pathway in fispact.

    Attributes
    ----------
    element : str
        The element of the isotope.
    isotope : int
        The isotope number.
    metastable : bool
        Whether the isotope is metastable or not.
    """

    element: str
    isotope: int
    metastable: bool = False

    def get_str(self) -> str:
        """Return the Zaid as a string.

        Returns
        -------
        str
            The Zaid as a string.
        """
        if self.metastable:
            return f"{self.element}{self.isotope}m"
        else:
            return f"{self.element}{self.isotope}"


@dataclass
class Pathway:
    """Class to represent a pathway in fispact.

    Attributes
    ----------
    parent : FispactZaid
        The parent FispactZaid object.
    daughter : FispactZaid
        The daughter FispactZaid object.
    perc : float
        The percentage of the pathway contribution to the daughter isotope.
    reactions : list[str]
        The list of reactions in the pathway.
    intermediates : list[FispactZaid], optional
        The list of intermediate FispactZaid objects in the pathway.

    """

    parent: FispactZaid
    daughter: FispactZaid
    perc: float
    reactions: list[str]
    intermediates: list[FispactZaid] = None

    def __post_init__(self):
        # check that reaction list length is always one more than intermediates
        # when it is not None
        if self.intermediates is not None:
            try:
                assert len(self.reactions) == len(self.intermediates) + 1
            except AssertionError as e:
                raise AssertionError(
                    f"Number of reactions is too low for pathway {self}"
                ) from e

    def __str__(self) -> str:
        if self.intermediates is not None:
            text = f"{self.parent.get_str()} "
            for intermediate, reaction in zip(self.intermediates, self.reactions):
                text += f"-{reaction}-> {intermediate.get_str()} "
            text += f"-{self.reactions[-1]}-> {self.daughter.get_str()}"
            return text
        else:
            return f"{self.parent.get_str()} -{self.reactions[0]}-> {self.daughter.get_str()}"


class PathwayCollection:
    def __init__(self, pathways: list[Pathway]) -> None:
        """Collection of pathways. This can be created from a list of
        Pathway objects or directly from a FISPACT legacy output file.

        Parameters
        ----------
        pathways : list[Pathway]
            list of Pathway objects.

        Attributes
        ----------
        pathways : list[Pathway]
            list of Pathway objects.

        Examples
        --------
        >>> from f4enix.output.fispact_legacy_out import PathwayCollection
        ... collection = PathwayCollection.from_file("path/to/fispact/output")
        ... pathway = collection.pathways[0]
        ... print(pathway)
        ... print(pathway.parent, pathway.daughter, pathway.reactions,
        ...       patway.intermediates, pathway.perc)
        Mn55 -(n,g)-> Mn56
        Mn55 Mn56 ['n,g'] None 100.0

        """
        self.pathways = pathways

    @classmethod
    def from_file(cls, file: os.PathLike) -> PathwayCollection:
        """
        Retrieve pathways from a FISPACT legacy output file and return a list of
        Pathway objects.

        Parameters
        ----------
        file : os.PathLike
            The path to the legacy FISPACT output file.

        Returns
        -------
        list[Pathway]
            A list of Pathway objects representing the pathways found in the file.
        """
        lines = []
        possible_ends = ["G E N E R I C   P A T H W", "1 * * * TIME INTERVAL"]
        end_reactions = None
        start_reactions = None
        look_for_end = False
        with open(file, "r", encoding="utf-8") as infile:
            for i, line in enumerate(infile):
                lines.append(line)
                if "Significant loops" in line:
                    start_reactions = i
                    look_for_end = True
                elif look_for_end and any(end in line for end in possible_ends):
                    end_reactions = i
                    break

        if end_reactions is None or start_reactions is None:
            raise ValueError(f"Could not find pathways section in file {file}.")

        lines = lines[start_reactions + 1 : end_reactions]

        paths = []
        for i, line in enumerate(lines):
            if path_id_pat.match(line) is not None:
                pathway = cls._parse_pathway(lines[i : i + 10])  # parse next 10 lines
                paths.append(pathway)

        return PathwayCollection(paths)

    def to_dataframe(self) -> pd.DataFrame:
        """Return the pathways as a pandas DataFrame.
        columns are "Parent", "% contribution", "Intermediates", "Reactions",
        "Daughter".

        Returns
        -------
        pd.DataFrame
            The pathways as a pandas DataFrame.
        """
        # get a libmanager for the sorting
        lm = LibManager()

        rows = []
        for pathway in self.pathways:
            intermediates = []
            if pathway.intermediates is not None:
                for intermediate in pathway.intermediates:
                    intermediates.append(intermediate.get_str())
            rows.append(
                [
                    pathway.parent.get_str(),
                    intermediates,
                    pathway.reactions,
                    pathway.daughter.get_str(),
                    pathway.perc,
                    # additional columns for sorting to be dropped later
                    lm.get_zaidnum(pathway.daughter.get_str()),
                ]
            )

        df = pd.DataFrame(rows)
        df.columns = [
            "Parent",
            "Intermediates",
            "Reactions",
            "Daughter",
            "% contribution",
            "isotope_sort",
        ]
        df.sort_values(by=["isotope_sort"], inplace=True)
        df.set_index(["Daughter", "% contribution"], inplace=True)
        del df["isotope_sort"]
        return df

    @staticmethod
    def _parse_pathway(text: str | list[str]) -> Pathway:
        if isinstance(text, str):
            lines = text.splitlines()
        else:
            lines = text

        reactions = None
        zaids = None
        perc = None
        for j, line in enumerate(lines):
            # be sure no EOL characters are present in the line
            line = line.strip("\n").strip("\rn").strip()
            if path_id_pat.match(line) is not None:
                # Get the zaids
                perc = float(perc_pattern.search(line).group()[:-1])
                tokens = line.split("---")
                # starting from pos 3, every two tokens should be a zaid and a reaction
                zaids = [_get_zaid_from_str(tokens[0].split("%")[-1])]

                for i in range(2, len(tokens[:-1]), 2):
                    zaids.append(_get_zaid_from_str(tokens[i]))

                # get the reactions
                line = lines[j + 1]
                reactions = reaction_pat.findall(line)

                if tokens[-2] in PATHWAY_END:
                    # then there are no more tokens to parse
                    break

            elif "path continued" in line:
                if zaids is None or reactions is None or perc is None:
                    raise ValueError(
                        f"Pathway continuation found but no pathway was being parsed. Text: {text}"
                    )
                # then the path continues on the next line, we need to parse the next line for more zaids and reactions
                tokens = line.split("---")
                zaids.append(_get_zaid_from_str(tokens[0].split("%")[-1]))
                for i in range(2, len(tokens[:-1]), 2):
                    zaids.append(_get_zaid_from_str(tokens[i]))

                line = lines[j + 1]
                reactions.extend(reaction_pat.findall(line))

                if tokens[-2] in PATHWAY_END:
                    # then there are no more tokens to parse
                    break

        if zaids is None or reactions is None or perc is None:
            raise ValueError(f"Could not parse pathway from text: {text}")

        pathway = Pathway(
            parent=zaids[0],
            daughter=zaids[-1],
            reactions=reactions,
            intermediates=zaids[1:-1],
            perc=perc,
        )
        return pathway


class FispactOutput:
    def __init__(self, filepath: PathLike, cooling_times: list[str]):
        """Store data parsed from FISPACT legacy output files

        Parameters
        ----------
        filepath : PathLike
            path to the fispact output
        cooling_times : list[str]
            list of labels associated to the different cooling times. No check is
            performed on the correct length of the label list, the last len(cooling_times)
            timesteps in the inventory data will be associated to these labels.

        Attributes
        ----------
        filepath : Path
            path to the fispact output
        name : str
            name of the fispact output file without extension
        inventory_data : list[pp.TimeStep]
            inventory data parsed from the fispact output file. This is a pypact object.
        """
        self.filepath = Path(filepath)
        self.name = self.filepath.stem

        with pp.Reader(filepath) as output:
            # store inventory data
            self.inventory_data = output.inventory_data

        self.sddr = self._get_sddr(cooling_times)
        self.pathways_collection = PathwayCollection.from_file(self.filepath)

    def _get_sddr(self, cooling_times: list[str]) -> pd.DataFrame:
        dfs = []
        for i, timestep in enumerate(self.inventory_data[-len(cooling_times) :]):
            cooling_time_label = cooling_times[i]

            doses = []
            for nuclide in timestep.nuclides:
                doses.append(
                    {
                        "element": nuclide.element,
                        "isotope": nuclide.isotope,
                        "state": nuclide.state,
                        "dose": nuclide.dose,
                        "cooling time": cooling_time_label,
                    }
                )

            df = pd.DataFrame(doses)
            df = df[df["dose"] > 0]  # filter zero doses
            df["isotope % dose"] = df["dose"] / df["dose"].sum() * 100
            df.sort_values(by="isotope % dose", ascending=False, inplace=True)
            df["Cumulative dose sum"] = df["isotope % dose"].values.cumsum()
            dfs.append(df)

        return pd.concat(dfs)

    def filter_by_cum_dose(
        self, perc: float, label: str, add_pathways: bool = False
    ) -> pd.DataFrame:
        """Filter the SDDR output (at a certain cooling time) to get at least a certain
        percent of the cumulative dose.

        Parameters
        ----------
        perc : float
            percent of cumulative dose to filter at
        label : str
            cooling time label to filter at
        add_pathways : bool, optional
            whether to add pathways rows, by default False
        Returns
        -------
        pd.DataFrame
            filtered dataframe
        """
        # select only requested label
        df = self.sddr[self.sddr["cooling time"] == label].copy()
        last_index = None
        for i, (_, row) in enumerate(df.iterrows()):
            if row["Cumulative dose sum"] > perc:
                last_index = i + 1
                break
        df = df.iloc[:last_index]

        if add_pathways:
            df = self.add_pathways_rows(df)
        return df

    def add_pathways_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add pathways rows to a SDDR dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            filtered dataframe

        Returns
        -------
        pd.DataFrame
            dataframe with pathways rows added
        """
        newrows = []
        for _, row in df.iterrows():
            isotope = str(row["element"]) + str(row["isotope"]) + row["state"]
            found = False
            for pathway in self.pathways_collection.pathways:
                if pathway.daughter.get_str() == isotope:
                    found = True
                    newrow = row.copy()
                    newrow["pathway"] = str(pathway)
                    newrow["pathway % dose"] = (
                        pathway.perc * row["isotope % dose"] / 100
                    )
                    newrows.append(newrow)
            if not found:
                newrow = row.copy()
                newrow["pathway"] = "N.A."
                newrow["pathway % dose"] = row["isotope % dose"]
                newrows.append(newrow)

        newdf = pd.DataFrame(newrows)

        # it may happen that pathways column is not created
        if "pathway % dose" not in newdf.columns:
            raise ValueError(
                "No pathways found for any of the isotopes in the dataframe."
            )
        return newdf.sort_values(by="pathway % dose", ascending=False)


def _get_zaid_from_str(name: str) -> FispactZaid:
    element = element_pat.search(name).group()
    isotope = isotope_pat.search(name).group()
    if metastable_pat.search(name) is not None:
        metastable = True
    else:
        metastable = False

    return FispactZaid(element=element, isotope=isotope, metastable=metastable)
