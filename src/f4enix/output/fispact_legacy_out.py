"""Parse legacy FISPACT-II output files to extract the generic pathways."""

from __future__ import annotations

import re
import os
import pandas as pd
import pypact as pp
from dataclasses import dataclass
from pathlib import Path
from f4enix.input.libmanager import LibManager
from f4enix.core.constants import PathLike, REVERSED_MT_DICT
from f4enix.core.irradiation import Nuclide, TCF_Computer
from copy import deepcopy


perc_pattern = re.compile(r"\d+\.*\d*%")
target_pathway_zaids = re.compile(r"[A-Z][a-z]*\s*\d+m*")
metastable_pat = re.compile(r"\d+m")
isotope_pat = re.compile(r"\d+")
element_pat = re.compile(r"[a-zA-Z]+")
path_id_pat = re.compile(r"\s*path\s+\d+")
reaction_pat = re.compile(r"\([a-zA-Z\d,+-]+\)")
PATHWAY_END = ["(S)", "(L)"]


@dataclass
class Pathway:
    """Class to represent a pathway in fispact.

    Attributes
    ----------
    parent : Nuclide
        The parent Nuclide object.
    daughter : Nuclide
        The daughter Nuclide object.
    perc : float
        The percentage of the pathway contribution to the daughter isotope.
    reactions : list[str]
        The list of reactions in the pathway.
    intermediates : list[Nuclide], optional
        The list of intermediate Nuclide objects in the pathway.

    """

    parent: Nuclide
    daughter: Nuclide
    perc: float
    reactions: list[str]
    intermediates: list[Nuclide] | None = None

    def __eq__(self, other: object) -> bool:
        # they are pathways
        if not isinstance(other, Pathway):
            return False
        # parent and daughters are the same
        if (
            self.parent != other.parent
            or self.daughter != other.daughter
            or self.intermediates != other.intermediates
            # given same intermediates reactions are automatically the same,
            # not checking helps avoid mistakes for small differences in reaction labels
            # or self.reactions != other.reactions
        ):
            return False
        return True

    def __post_init__(self):
        # if an empty list is provided, change it to None
        if isinstance(self.intermediates, list) and len(self.intermediates) == 0:
            self.intermediates = None
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
            text = f"{self.parent.write_to_formula()} "
            for intermediate, reaction in zip(self.intermediates, self.reactions):
                text += f"-{reaction}-> {intermediate.write_to_formula()} "
            text += f"-{self.reactions[-1]}-> {self.daughter.write_to_formula()}"
            return text
        else:
            return f"{self.parent.write_to_formula()} -{self.reactions[0]}-> {self.daughter.write_to_formula()}"

    def is_multistep(self) -> bool:
        """Return whether the pathway is a multistep pathway or not.
        Isomeric transitions are not considered as steps in the pathway.

        Returns
        -------
        bool
            True if the pathway is a multistep pathway, False otherwise.
        """
        if self.intermediates is not None:
            for intermediate in self.intermediates:
                if (
                    intermediate.isotope == self.daughter.isotope
                    and intermediate.element == self.daughter.element
                ):
                    # then this is an isomeric transition
                    # we do not consider it as a step in the pathway
                    continue
                else:
                    return True
        return False

    @classmethod
    def from_string(cls, string: str, perc: float = 100) -> "Pathway":
        """parse from string like: Ni58 -(n,p)-> Co58m -(IT)-> Co58

        Parameters
        ----------
        str : str
            string to parse

        Returns
        -------
        Pathway
            parsed pathway
        """
        reactions = reaction_pat.findall(string)
        # replace the reactions with nothing
        string = reaction_pat.sub("", string)
        zaids = string.split("-->")
        parent = Nuclide.from_formula(zaids[0].strip())
        daughter = Nuclide.from_formula(zaids[-1].strip())
        intermediates = []
        for intermediate in zaids[1:-1]:
            intermediates.append(Nuclide.from_formula(intermediate.strip()))

        return cls(parent, daughter, perc, reactions, intermediates=intermediates)

    def reduce(
        self,
        half_life_cutoff: float | None = None,
        tfc_computer: TCF_Computer | None = None,
    ) -> Pathway:
        """Return a reduced pathway built from the original where isomeric transitions
        have been removed. It is possible to provide a half life cut-off value above
        which metastable steps will be retained.

        Parameters
        ----------
        half_life_cutoff : float | None, optional
            threshold (in seconds) below which isomeric transitions shall be removed
            by default None. If a value is provided, also the tfc_computer must be
            provided.
        tfc_computer : TCF_Computer | None, optional
            TFC_Computer object used to retrieve half-life values, by default None

        Returns
        -------
        Pathway
            New reduced pathway

        Raises
        ------
        ValueError
            If a half life cut-off is provided but no TFC_Computer.
        """
        if half_life_cutoff is not None:
            if tfc_computer is None:
                raise ValueError("TFC_Computer must be provided to compute half-lives")

        it_strings = ["IT", "(IT)"]
        reactions = [self.reactions[0]]
        intermediates = []
        if self.intermediates is None:
            return deepcopy(self)
        for i, intermediate in enumerate(self.intermediates):
            if self.reactions[i + 1] in it_strings:
                # if the metastable is below cut-off, simplify
                if half_life_cutoff and tfc_computer:
                    nuclide = Nuclide.from_formula(
                        str(intermediate.element) + str(intermediate.isotope) + "m",
                    )
                    half_life = tfc_computer.get_half_life(nuclide)
                    if float(half_life) <= half_life_cutoff:
                        continue
                # if no cut-off provided, always simplify
                else:
                    continue
            reactions.append(self.reactions[i + 1])
            intermediates.append(deepcopy(intermediate))
        return Pathway(
            deepcopy(self.parent),
            deepcopy(self.daughter),
            self.perc,
            reactions,
            intermediates=intermediates,
        )

    def get_MT(self) -> int | None:
        """Return the MT number associated to the first reaction in the pathway."""
        # ensure parenthesis are present and no whitespaces
        key = f"({self.reactions[0].replace(' ', '').strip('(').strip(')')})"
        try:
            return REVERSED_MT_DICT[key]
        except KeyError:
            return None


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
                    intermediates.append(intermediate.write_to_formula())
            rows.append(
                [
                    pathway.parent.write_to_formula(),
                    intermediates,
                    pathway.reactions,
                    pathway.daughter.write_to_formula(),
                    pathway.perc,
                    # additional columns for sorting to be dropped later
                    lm.get_zaidnum(pathway.daughter.write_to_formula()),
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
            line = line.strip()
            if path_id_pat.match(line) is not None:
                # Get the zaids
                perc = float(perc_pattern.search(line).group()[:-1])
                tokens = line.split("---")
                # starting from pos 3, every two tokens should be a zaid and a reaction
                zaids = [Nuclide.from_formula(tokens[0].split("%")[-1])]

                for i in range(2, len(tokens[:-1]), 2):
                    zaids.append(Nuclide.from_formula(tokens[i]))

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
                zaids.append(Nuclide.from_formula(tokens[0].split("continued")[-1]))
                for i in range(2, len(tokens[:-1]), 2):
                    zaids.append(Nuclide.from_formula(tokens[i]))

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
    def __init__(
        self,
        filepath: PathLike,
        cooling_times: list[str],
        parse_sddr: bool = True,
        parse_decay_heat: bool = True,
        parse_pathways: bool = True,
    ) -> None:
        """Store data parsed from FISPACT legacy output files

        Parameters
        ----------
        filepath : PathLike
            path to the fispact output
        cooling_times : list[str]
            list of labels associated to the different cooling times. No check is
            performed on the correct length of the label list, the last len(cooling_times)
            timesteps in the inventory data will be associated to these labels.
        parse_sddr : bool, optional
            whether to parse the SDDR data, by default True
        parse_decay_heat : bool, optional
            whether to parse the decay heat data, by default True
        parse_pathways : bool, optional
            whether to parse the pathways data, by default True

        Attributes
        ----------
        filepath : Path
            path to the fispact output
        name : str
            name of the fispact output file without extension
        inventory_data : list[pp.TimeStep]
            inventory data parsed from the fispact output file. This is a pypact object.
        sddr : pd.DataFrame
            SDDR data extracted from the inventory data.
        decay_heat : pd.DataFrame
            Decay heat data extracted from the inventory data.
        pathways_collection : PathwayCollection
            PathwayCollection object containing the pathways extracted from the fispact
            output file.
        """
        self.filepath = Path(filepath)
        self.name = self.filepath.stem

        with pp.Reader(filepath) as output:
            # store inventory data
            self.inventory_data = output.inventory_data

        self.sddr = self._get_sddr(cooling_times) if parse_sddr else None
        self.decay_heat = (
            self._get_decay_heat(cooling_times) if parse_decay_heat else None
        )
        self.pathways_collection = (
            PathwayCollection.from_file(self.filepath) if parse_pathways else None
        )

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

    def _get_decay_heat(self, cooling_times: list[str]) -> pd.DataFrame:
        dfs = []
        for i, timestep in enumerate(self.inventory_data[-len(cooling_times) :]):
            cooling_time_label = cooling_times[i]

            heat_rows = []
            for nuclide in timestep.nuclides:
                heat_rows.append(
                    {
                        "element": nuclide.element,
                        "isotope": nuclide.isotope,
                        "state": nuclide.state,
                        "heat": nuclide.heat,
                        "alpha_heat": nuclide.alpha_heat,
                        "beta_heat": nuclide.beta_heat,
                        "gamma_heat": nuclide.gamma_heat,
                        "cooling time": cooling_time_label,
                    }
                )

            df = pd.DataFrame(heat_rows)
            df = df[df["heat"] > 0]  # filter zero heat
            df["isotope % heat"] = df["heat"] / df["heat"].sum() * 100
            df.sort_values(by="isotope % heat", ascending=False, inplace=True)
            df["Cumulative heat sum"] = df["isotope % heat"].values.cumsum()
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

    def filter_by_cum_heating(
        self, perc: float, label: str, add_pathways: bool = False
    ) -> pd.DataFrame:
        """Filter the decay heating output (at a certain cooling time)
        to get at least a certain percent of the cumulative heat.

        Parameters
        ----------
        perc : float
            percent of cumulative heat to filter at
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
        df = self.decay_heat[self.decay_heat["cooling time"] == label].copy()
        last_index = None
        for i, (_, row) in enumerate(df.iterrows()):
            if row["Cumulative heat sum"] > perc:
                last_index = i + 1
                break
        df = df.iloc[:last_index]

        if add_pathways:
            df = self.add_pathways_rows(df, who="heat")
        return df

    def add_pathways_rows(self, df: pd.DataFrame, who="dose") -> pd.DataFrame:
        """Add pathways rows to a SDDR dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            filtered dataframe
        who : str, optional
            whether to add pathways rows for 'dose' or 'heat' dataframe, by default 'dose'.

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
                if pathway.daughter.write_to_formula() == isotope:
                    found = True
                    newrow = row.copy()
                    newrow["pathway"] = str(pathway)
                    newrow[f"pathway % {who}"] = (
                        pathway.perc * row[f"isotope % {who}"] / 100
                    )
                    newrows.append(newrow)
            if not found:
                newrow = row.copy()
                newrow["pathway"] = "N.A."
                newrow[f"pathway % {who}"] = row[f"isotope % {who}"]
                newrows.append(newrow)

        newdf = pd.DataFrame(newrows)

        # it may happen that pathways column is not created
        if f"pathway % {who}" not in newdf.columns:
            raise ValueError(
                f"No pathways found for any of the isotopes in the {who} dataframe."
            )
        return newdf.sort_values(by=f"pathway % {who}", ascending=False)
