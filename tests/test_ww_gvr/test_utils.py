import numpy as np

from f4enix.input.ww_gvr.utils import (
    Vectors,
    build_1d_vectors,
    compose_b2_vectors,
    decompose_b2_vectors,
)


def test_vectors():
    vectors = Vectors(
        vector_i=np.array([0.0, 2.0, 15.0, 1.0]),
        vector_j=np.array([0.0, 3.0, 16.0, 1.0]),
        vector_k=np.array([0.0, 1.0, 20.0, 1.0]),
    )

    iterated = []
    for vector in vectors:
        iterated.append(vector)
    assert iterated[0].tolist() == vectors.vector_i.tolist()

    different = Vectors(
        vector_i=np.array([0.0, 2.0, 16.0, 1.0]),
        vector_j=np.array([0.0, 3.0, 16.0, 1.0]),
        vector_k=np.array([0.0, 1.0, 20.0, 1.0]),
    )
    assert vectors != different


def test_decompose_b2_vectors():
    b2_vectors = Vectors(
        vector_i=np.array([0.0, 2.0, 15.0, 1.0]),
        vector_j=np.array([0.0, 3.0, 16.0, 1.0, 2.0, 19.0, 1.0]),
        vector_k=np.array([0.0, 1.0, 20.0, 1.0]),
    )
    expected_coarse = Vectors(
        vector_i=np.array([0.0, 15.0]),
        vector_j=np.array([0.0, 16.0, 19.0]),
        vector_k=np.array([0.0, 20.0]),
    )
    expected_fine = Vectors(
        vector_i=np.array([2.0]),
        vector_j=np.array([3.0, 2.0]),
        vector_k=np.array([1.0]),
    )

    result_coarse, result_fine = decompose_b2_vectors(b2_vectors)

    assert expected_coarse == result_coarse
    assert expected_fine == result_fine


def test_compose_b2_vectors():
    expected_b2_vectors = Vectors(
        vector_i=np.array([0.0, 2.0, 15.0, 1.0]),
        vector_j=np.array([0.0, 3.0, 16.0, 1.0, 2.0, 19.0, 1.0]),
        vector_k=np.array([0.0, 1.0, 20.0, 1.0]),
    )

    coarse_vectors = Vectors(
        vector_i=np.array([0.0, 15.0]),
        vector_j=np.array([0.0, 16.0, 19.0]),
        vector_k=np.array([0.0, 20.0]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([2.0]),
        vector_j=np.array([3.0, 2.0]),
        vector_k=np.array([1.0]),
    )
    result_b2_vectors = compose_b2_vectors(coarse_vectors, fine_vectors)

    assert expected_b2_vectors == result_b2_vectors


def test_build_1d_vectors():
    coarse_vectors = Vectors(
        vector_i=np.array([-15.0, 5.1, 12.2, 13.3]),
        vector_j=np.array([0.0, 16.0, 19.0]),
        vector_k=np.array([0.0, 1.0]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([2.0, 3.0, 2.0]),
        vector_j=np.array([4.0, 3.0]),
        vector_k=np.array([1.0]),
    )
    expected_vectors = Vectors(
        vector_i=np.array(
            [-15.0, -4.95, 5.1, 7.46666667, 9.83333333, 12.2, 12.75, 13.3]
        ),
        vector_j=np.array([0.0, 4.0, 8.0, 12.0, 16.0, 17.0, 18.0, 19.0]),
        vector_k=np.array([0.0, 1.0]),
    )

    result_vectors = build_1d_vectors(coarse_vectors, fine_vectors)

    assert expected_vectors == result_vectors
