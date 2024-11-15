"""Parse legacy FISPACT-II output files to extract the generic pathways.

"""

from __future__ import annotations

import re
import os
import pandas as pd
from dataclasses import dataclass
from f4enix.input.libmanager import LibManager


perc_pattern = re.compile(r"\d+\.*\d*%")
target_pathway_zaids = re.compile(r"[A-Z][a-z]*\s*\d+m*")
metastable_pat = re.compile(r"\d+m")
isotope_pat = re.compile(r"\d+")
element_pat = re.compile(r"[a-zA-Z]+")
path_id_pat = re.compile(r"\s+path\s+\d+")
reaction_pat = re.compile(r"\([a-zA-Z\d,+-]+\)")


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
            text += f"-{self.reactions[-1]}->  {self.daughter.get_str()}"
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

    @staticmethod
    def _get_zaid_from_str(name: str) -> FispactZaid:
        element = element_pat.search(name).group()
        isotope = isotope_pat.search(name).group()
        if metastable_pat.search(name) is not None:
            metastable = True
        else:
            metastable = False

        return FispactZaid(element=element, isotope=isotope, metastable=metastable)

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
        with open(file, "r", encoding="utf-8") as infile:
            for i, line in enumerate(infile):
                lines.append(line)
                if "Significant loops" in line:
                    start_reactions = i
                elif "G E N E R I C   P A T H W" in line:
                    end_reactions = i
        lines = lines[start_reactions + 1 : end_reactions]

        paths = []
        for i, line in enumerate(lines):
            if path_id_pat.match(line) is not None:
                perc = float(perc_pattern.search(line).group()[:-1])
                zaids = target_pathway_zaids.findall(line)
                try:
                    parent = cls._get_zaid_from_str(zaids[0])
                except IndexError as e:
                    print(line)
                    raise e
                daughter = cls._get_zaid_from_str(zaids[-1])
                if len(zaids) > 2:
                    intermediates = []
                    for zaid in zaids[1:-1]:
                        intermediates.append(cls._get_zaid_from_str(zaid))
                else:
                    intermediates = None

                # if a path is found, the following line will contain details
                # on the reactions
                line = lines[i + 1]
                reactions = reaction_pat.findall(line)

                try:
                    pathway = Pathway(
                        parent=parent,
                        daughter=daughter,
                        reactions=reactions,
                        intermediates=intermediates,
                        perc=perc,
                    )
                except AssertionError as e:
                    # it may be a very long path that continues
                    if "path continued" in lines[i + 4]:
                        other_zaids = target_pathway_zaids.findall(lines[i + 4])
                        intermediates.append(daughter)
                        for zaid in other_zaids[:-1]:
                            intermediates.append(cls._get_zaid_from_str(zaid))
                        daughter = cls._get_zaid_from_str(other_zaids[-1])
                        reactions.extend(reaction_pat.findall(lines[i + 5]))

                        try:
                            pathway = Pathway(
                                parent=parent,
                                daughter=daughter,
                                reactions=reactions,
                                intermediates=intermediates,
                                perc=perc,
                            )
                        except AssertionError:
                            print(lines[i])
                            print(line)
                            raise e
                    else:
                        print(lines[i])
                        print(line)
                        raise e

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
