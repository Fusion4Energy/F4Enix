import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BYTE = np.byte
CHAR = np.char
INT = np.int32
FLOAT = np.float64
LONG = np.int64
NUMBER_OF_EXPECTED_VALUES = 11  # Number of values recorded for each particle
NEUTRON_INDICATOR = 8  # Neutron indicator in the packed variable


@dataclass
class SurfaceParameters:
    id: int
    info: int
    type: int
    num_parameters: int
    parameters: list[int]


@dataclass
class FileParameters:
    np1: int  # Number of histories of the simulation, given as a negative number
    nrss: int  # Number of tracks recorded
    nrcd: int  # Number of values recorded for each particle, it should be 11
    njsw: int  # Number of surfaces in JASW
    niss: int  # Number of different histories that reached the SSW surfaces
    niwr: int  # Number of cells in RSSA file
    mipts: int  # Source particle type
    kjaq: int  # Flag for macrobodies surfaces
    surfaces: list[SurfaceParameters]  # List with surface ids that appear in the file

    def save_to_json(self, path: Path) -> None:
        """Saves the file parameters to a JSON file."""

        # The attributes should have Python basic types not numpy types for
        # serialization
        data = {
            "np1": int(self.np1),
            "nrss": int(self.nrss),
            "nrcd": int(self.nrcd),
            "njsw": int(self.njsw),
            "niss": int(self.niss),
            "niwr": int(self.niwr),
            "mipts": int(self.mipts),
            "kjaq": int(self.kjaq),
            "surfaces": [
                {
                    "id": int(surface.id),
                    "info": int(surface.info),
                    "type": int(surface.type),
                    "num_parameters": int(surface.num_parameters),
                    "parameters": [int(param) for param in surface.parameters],
                }
                for surface in self.surfaces
            ],
        }

        with open(path, "w") as outfile:
            json.dump(data, outfile, indent=4)

    @staticmethod
    def load_from_json(path: Path) -> "FileParameters":
        """Loads the file parameters from a JSON file."""
        with open(path) as infile:
            data = json.load(infile)
            data["surfaces"] = [
                SurfaceParameters(
                    id=surface["id"],
                    info=surface["info"],
                    type=surface["type"],
                    num_parameters=surface["num_parameters"],
                    parameters=surface["parameters"],
                )
                for surface in data["surfaces"]
            ]
            return FileParameters(**data)
