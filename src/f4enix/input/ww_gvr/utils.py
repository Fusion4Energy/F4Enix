"""
Some utilities for the ww_gvr package.
"""

from typing import Tuple

import numpy as np

from f4enix.input.ww_gvr.models import Vectors


def decompose_b2_vectors(b2_vectors: Vectors) -> Tuple[Vectors, Vectors]:
    """
    Takes a Vectors object with b2 format and returns two Vectors objects, one with the
    coarse vectors and the other with the fine vectors.

    The b2 format of a WW MCNP file is:
            [mesh_position, fine_ints,
            mesh_position, 1.0000, fine_ints,
            mesh_position, 1.0000, fine_ints,
            ...
            mesh_position, 1.0000]

    Coarse vectors have the format: [p1, p2, p3, ...]
    A fine vector shows how many intervals there are between each step of the coarse
    vector's points: [n1, n2, n3, ...]

    Parameters
    ----------
    b2_vectors : Vectors
        Vectors object with b2 format.

    Returns
    -------
    Tuple[Vectors, Vectors]
        Two Vectors objects, one with the coarse vectors and the other with the
        fine vectors.
    """
    coarse = []
    fine = []

    for b2_vector in b2_vectors:
        vector = np.insert(b2_vector, 1, 1.0000)
        coarse.append([vector[i] for i in range(len(vector)) if i % 3 == 0])
        fine.append(
            [
                int(vector[j + 2])
                for j in range(len(vector))
                if j % 3 == 0 and j + 2 < len(vector)
            ]
        )

    coarse_vectors = Vectors(
        vector_i=np.array(coarse[0]),
        vector_j=np.array(coarse[1]),
        vector_k=np.array(coarse[2]),
    )
    fine_vectors = Vectors(
        vector_i=np.array(fine[0]),
        vector_j=np.array(fine[1]),
        vector_k=np.array(fine[2]),
    )

    return coarse_vectors, fine_vectors


def compose_b2_vectors(coarse_vectors: Vectors, fine_vectors: Vectors) -> Vectors:
    """
    The opposite of decompose_b2_vectors.

    Parameters
    ----------
    coarse_vectors : Vectors
        Vectors object with coarse vectors.
    fine_vectors : Vectors
        Vectors object with fine vectors.

    Returns
    -------
    Vectors
        Vectors object with b2 format.
    """
    b2_vectors = []

    for coarse_vector, fine_ints in zip(coarse_vectors, fine_vectors):
        current_vector = []
        for i in range(len(coarse_vector) - 1):
            current_vector.append(coarse_vector[i])
            current_vector.append(1.0000)
            current_vector.append(fine_ints[i])
        current_vector.append(coarse_vector[-1])
        current_vector.append(1)
        current_vector.remove(1)
        b2_vectors.append(current_vector)
    return Vectors(
        vector_i=np.array(b2_vectors[0]),
        vector_j=np.array(b2_vectors[1]),
        vector_k=np.array(b2_vectors[2]),
    )


def build_1d_vectors(coarse_vectors: Vectors, fine_vectors: Vectors) -> Vectors:
    """
    Flatten coarse and fine vectors into 1D Vectors.

    Takes two Vectors objects, one with the coarse vectors and the other with the
    fine vectors which represent the amount of ints between each step of the coarse
    vectors, and returns some flattened vectors with format:
    [mesh_position, mesh_position, mesh_position, ...]

    Parameters
    ----------
    coarse_vectors : Vectors
        Vectors object with coarse vectors.
    fine_vectors : Vectors
        Vectors object with fine vectors.

    Returns
    -------
    Vectors
        Vectors object with 1D vectors.
    """
    b2_vectors = []

    for coarse_vector, fine_ints in zip(coarse_vectors, fine_vectors):
        current_vector = [coarse_vector[0]]
        for i in range(len(coarse_vector) - 1):
            current_vector += np.linspace(
                coarse_vector[i], coarse_vector[i + 1], int(fine_ints[i] + 1)
            ).tolist()[1:]

        b2_vectors.append(current_vector)

    return Vectors(
        vector_i=np.array(b2_vectors[0]),
        vector_j=np.array(b2_vectors[1]),
        vector_k=np.array(b2_vectors[2]),
    )
