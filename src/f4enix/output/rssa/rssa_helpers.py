from dataclasses import dataclass

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
