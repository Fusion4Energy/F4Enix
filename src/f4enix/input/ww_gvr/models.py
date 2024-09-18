"""
Classes and types used throughout the package.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
from numpy.typing import NDArray


class CoordinateType(Enum):
    CARTESIAN = "cartesian"
    CYLINDRICAL = "cylindrical"


class ParticleType(Enum):
    NEUTRON = "neutron"
    PHOTON = "photon"


@dataclass()
class Vectors:
    """
    A group of three vectors in the form of a Numpy NDArray. They are used to define 3D
    meshes.
    """

    vector_i: NDArray  # Vector defining i direction
    vector_j: NDArray  # Vector defining j direction
    vector_k: NDArray  # Vector defining k direction

    def __iter__(self):
        return iter([self.vector_i, self.vector_j, self.vector_k])

    def __eq__(self, other):
        try:
            np.testing.assert_array_almost_equal(self.vector_i, other.vector_i)
            np.testing.assert_array_almost_equal(self.vector_j, other.vector_j)
            np.testing.assert_array_almost_equal(self.vector_k, other.vector_k)
            return True
        except AssertionError:
            return False


Pathlike = Union[str, Path]
NestedList = List[List[float]]
ValuesByEnergy = Dict[float, np.ndarray]
ValuesByParticle = Dict[ParticleType, ValuesByEnergy]
EnergiesByParticle = Dict[ParticleType, List[float]]
