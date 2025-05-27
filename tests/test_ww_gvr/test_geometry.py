import numpy as np
import pyvista as pv

from f4enix.input.ww_gvr.geometry import Geometry
from f4enix.input.ww_gvr.models import ParticleType, Vectors
from f4enix.input.ww_gvr.ww_parser import WWHeader, WWHeaderCyl

# ruff: noqa: PLR2004


def test_fill_grid_values_cart():
    header = WWHeader(
        if_=1,
        iv=1,
        ni=1,
        nr=16,
        probid="06/04/21 17:17:28",
        ne=[1],
        nfx=2,
        nfy=2,
        nfz=2,
        origin=[0.0, 0.0, 0.0],
        ncx=2,
        ncy=2,
        ncz=2,
    )
    coarse_vectors = Vectors(
        vector_i=np.array([0, 10, 20]),
        vector_j=np.array([50, 60, 70]),
        vector_k=np.array([150, 160, 170]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([1, 1]),
        vector_j=np.array([1, 1]),
        vector_k=np.array([1, 1]),
    )
    geometry = Geometry(header, coarse_vectors, fine_vectors)

    values = {
        ParticleType.NEUTRON: {100.0: np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])}
    }

    geometry.fill_grid_values(values)
    grid = geometry._grid
    cell_index = grid.find_closest_cell([5, 65, 155])
    value = grid[grid.array_names[0]][cell_index]  # type: ignore

    assert value == 3


def test_fill_grid_values_cyl():
    header = WWHeaderCyl(
        if_=1,
        iv=1,
        ni=1,
        nr=16,
        probid="06/04/21 17:17:28",
        ne=[1],
        nfx=2,
        nfy=2,
        nfz=3,
        origin=[0.0, 0.0, 0.0],
        ncx=2,
        ncy=2,
        ncz=3,
        director_1=[0.0, 0.0, 1.0],
        director_2=[1.0, 0.0, 0.0],
    )
    coarse_vectors = Vectors(
        vector_i=np.array([0, 10, 20]),
        vector_j=np.array([10, 20, 30]),
        vector_k=np.array([0, 0.25, 0.5, 0.75]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([1, 1]),
        vector_j=np.array([1, 1]),
        vector_k=np.array([1, 1, 1]),
    )
    geometry = Geometry(header, coarse_vectors, fine_vectors)

    values = {
        ParticleType.NEUTRON: {
            100.0: np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]], [[9, 10], [11, 12]]])
        }
    }

    geometry.fill_grid_values(values)
    grid = geometry._grid

    # Sample the value at an expected point because .find_closest_cell() does not work
    #  well with cylindrical grids
    sphere = pv.Sphere(center=[-3, 3, 15], radius=1)
    sphere = sphere.sample(grid)  # type: ignore
    result = sphere[grid.array_names[0]][0]  # type: ignore

    assert result == 5


def test_fill_grid_values_cyl_two_theta_ints():
    header = WWHeaderCyl(
        if_=1,
        iv=1,
        ni=1,
        nr=16,
        probid="06/04/21 17:17:28",
        ne=[1],
        nfx=2,
        nfy=2,
        nfz=2,
        origin=[0.0, 0.0, 0.0],
        ncx=2,
        ncy=2,
        ncz=2,
        director_1=[0.0, 0.0, 1.0],
        director_2=[1.0, 0.0, 0.0],
    )
    coarse_vectors = Vectors(
        vector_i=np.array([0, 10, 20]),
        vector_j=np.array([10, 20, 30]),
        vector_k=np.array([0, 0.5, 1]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([1, 1]),
        vector_j=np.array([1, 1]),
        vector_k=np.array([1, 1]),
    )
    geometry = Geometry(header, coarse_vectors, fine_vectors)

    values = {
        ParticleType.NEUTRON: {100.0: np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])}
    }

    geometry.fill_grid_values(values)
    grid = geometry._grid

    # Sample the value at an expected point because .find_closest_cell() does not work
    #  well with cylindrical grids
    sphere = pv.Sphere(center=[0, -5, 15], radius=1)
    sphere = sphere.sample(grid)  # type: ignore
    result = sphere[grid.array_names[0]][0]  # type: ignore

    assert result == 5


def test_decide_plot_parameters_linear():
    header = WWHeader(
        if_=1,
        iv=1,
        ni=1,
        nr=16,
        probid="06/04/21 17:17:28",
        ne=[1],
        nfx=2,
        nfy=2,
        nfz=2,
        origin=[0.0, 0.0, 0.0],
        ncx=2,
        ncy=2,
        ncz=2,
    )
    coarse_vectors = Vectors(
        vector_i=np.array([0, 10, 20]),
        vector_j=np.array([50, 60, 70]),
        vector_k=np.array([150, 160, 170]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([1, 1]),
        vector_j=np.array([1, 1]),
        vector_k=np.array([1, 1]),
    )
    geometry = Geometry(header, coarse_vectors, fine_vectors)

    values = {
        ParticleType.NEUTRON: {100.0: np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])}
    }

    geometry.fill_grid_values(values)

    result = geometry._decide_plot_parameters()

    expected = {
        "scalars": "WW ParticleType.NEUTRON, 100.00 MeV",
        "scalar_bar_args": {
            "interactive": True,
            "title_font_size": 20,
            "label_font_size": 16,
            "shadow": True,
            "italic": True,
            "fmt": "%.e",
            "font_family": "arial",
        },
        "show_edges": False,
        "cmap": "jet",
        "log_scale": False,
        "lighting": False,
    }

    assert result == expected


def test_decide_plot_parameters_log():
    header = WWHeader(
        if_=1,
        iv=1,
        ni=1,
        nr=16,
        probid="06/04/21 17:17:28",
        ne=[1],
        nfx=2,
        nfy=2,
        nfz=2,
        origin=[0.0, 0.0, 0.0],
        ncx=2,
        ncy=2,
        ncz=2,
    )
    coarse_vectors = Vectors(
        vector_i=np.array([0, 10, 20]),
        vector_j=np.array([50, 60, 70]),
        vector_k=np.array([150, 160, 170]),
    )
    fine_vectors = Vectors(
        vector_i=np.array([1, 1]),
        vector_j=np.array([1, 1]),
        vector_k=np.array([1, 1]),
    )
    geometry = Geometry(header, coarse_vectors, fine_vectors)

    values = {
        ParticleType.NEUTRON: {
            100.0: np.array([[[1e3, 2e7], [3e5, 4e2]], [[5e2, 6e6], [7, 8]]])
        }
    }

    geometry.fill_grid_values(values)

    result = geometry._decide_plot_parameters()

    expected = {
        "scalars": "WW ParticleType.NEUTRON, 100.00 MeV",
        "scalar_bar_args": {
            "interactive": True,
            "title_font_size": 20,
            "label_font_size": 16,
            "shadow": True,
            "italic": True,
            "fmt": "%.e",
            "font_family": "arial",
            "n_labels": 10,
        },
        "show_edges": False,
        "cmap": "jet",
        "log_scale": True,
        "lighting": False,
        "clim": [0.001, 10000000],
        "n_colors": 10,
    }

    assert result == expected
