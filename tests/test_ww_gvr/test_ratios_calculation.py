import numpy as np

from f4enix.input.ww_gvr.ratios_calculation import calculate_max_ratio_array


def test_calculate_max_ratio_array():
    values = np.array([[[1, 2], [1, 1]], [[1, 1], [1, 0]]])
    expected = np.array([[[2, 2], [1, 2]], [[1, 2], [1, 1]]])
    result = calculate_max_ratio_array(values)
    assert np.allclose(expected, result)
