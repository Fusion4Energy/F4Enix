import os
import re
import pytest
import pyvista as pv
from importlib.resources import files, as_file

from f4eparser.output.meshtal import (Meshtal, Fmesh, scalemesh, addmesh,
                                      diffmesh, identical_mesh)
import tests.resources.meshtal as resources

RESOURCES = files(resources)


@pytest.mark.parametrize(
    "input_meshtal",
    [
        "meshtal_cuv",
        "meshtal_cyl",
        "meshtal_d1s_CSimpactStudy",
        "meshtal_CUBE_SQUARE",
        "meshtal_CUBE_ONES",
    ],
)
def test_mesh_print_tally_info(input_meshtal):
    # To check if the meshtal can be read without any problem"
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
        meshtally = Meshtal(inp, filetype)

    for i in meshtally.mesh.items():
        meshtally.mesh[i[0]].print_info()

    assert True


@pytest.mark.parametrize(
    "input_meshtal",
    ["meshtal_cuv", "meshtal_cyl", "meshtal_d1s_CSimpactStudy"]
)
def test_mesh_print_info(input_meshtal):
    # To check if the meshtal can be read without any problem"
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
        meshtally = Meshtal(inp, filetype)

    for i in meshtally.mesh.items():
        meshtally.print_info()

    assert True


# ************** STATUS OF TESTING *****************
# Mesh scale               (scale)      ---> DONE
# Mesh sum                 (sum)        ---> DONE
# Mesh difference          (diff)       ---> DONE
# Energy bin sum           (binsum)     ---> PENDING
# Check identical mesh     (identical)  ---> DONE
# Change correlation       (corr)       ---> DONE
# CuV testing                           ---> PENDING

# 'meshtal_CUBE_SQUARE'
#  ^
#  |   9  16
#  |   1  4
#  Y
#     X --->

# 'meshtal_CUBE_ONES'
#  ^
#  |   1  1
#  |   1  1
#  Y
#     X --->


@pytest.mark.parametrize("input_meshtal", ["meshtal_CUBE_SQUARE"])
def test_mesh_VTKwrite(input_meshtal):
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
        meshtally = Meshtal(inp, filetype)

    meshtally.mesh[124].print_info()
    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    name = "test_VTK_CUBE_SQUARE.vtr"

    meshobj.writeVTK(name)
    assert True


@pytest.mark.parametrize("filename", ["test_VTK_CUBE_SQUARE.vtr"])
def test_mesh_VTKcheck(filename):
    with as_file(RESOURCES.joinpath(filename)) as inp:
        mesh = pv.read(inp)

    assert mesh["Value - Total"][0] == 1
    assert mesh["Value - Total"][1] == 9
    assert mesh["Value - Total"][2] == 4
    assert mesh["Value - Total"][3] == 16


@pytest.mark.parametrize("sfactor", [1, 2, 3, 10, 20])
def test_mesh_scale(sfactor):
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
        meshtally = Meshtal(inp, filetype)
    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    smesh = scalemesh(meshobj, sfactor)

    assert smesh.dat[0][0][0][0] == 1 * sfactor
    assert smesh.dat[0][0][0][1] == 9 * sfactor
    assert smesh.dat[0][0][1][0] == 4 * sfactor
    assert smesh.dat[0][0][1][1] == 16 * sfactor


def test_mesh_sum():
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
        meshtally = Meshtal(inp, filetype)

    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    smesh = addmesh(meshobj, meshobj, f1=1.0, f2=1.0, corr=False)

    assert smesh.dat[0][0][0][0] == 1 * 2
    assert smesh.dat[0][0][0][1] == 9 * 2
    assert smesh.dat[0][0][1][0] == 4 * 2
    assert smesh.dat[0][0][1][1] == 16 * 2


def test_mesh_corr():
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath("meshtal_CUBE_ONES")) as inp:
        meshtally = Meshtal(inp, filetype)

    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    smesh = addmesh(meshobj, meshobj, f1=1.0, f2=1.0, corr=True)

    assert smesh.err[0][0][0][0] == 1
    assert smesh.err[0][0][0][1] == 1
    assert smesh.err[0][0][1][0] == 1
    assert smesh.err[0][0][1][1] == 1

    smesh = addmesh(meshobj, meshobj, f1=1.0, f2=1.0, corr=False)

    assert smesh.err[0][0][0][0] == ((1 + 1) ** 0.5) / 2
    assert smesh.err[0][0][0][1] == ((1 + 1) ** 0.5) / 2
    assert smesh.err[0][0][1][0] == ((1 + 1) ** 0.5) / 2
    assert smesh.err[0][0][1][1] == ((1 + 1) ** 0.5) / 2


def test_mesh_diff():

    filetype = "MCNP"
    with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
        meshtally = Meshtal(inp, filetype)

    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    smesh = diffmesh(meshobj, meshobj)

    assert smesh.dat[0][0][0][0] == 0
    assert smesh.dat[0][0][0][1] == 0
    assert smesh.dat[0][0][1][0] == 0
    assert smesh.dat[0][0][1][1] == 0


def test_mesh_identical():
    filetype = "MCNP"
    with as_file(RESOURCES.joinpath("meshtal_CUBE_SQUARE")) as inp:
        meshtally = Meshtal(inp, filetype)

    meshtally.readMesh([124])
    meshobj = meshtally.mesh[124]

    part, mesh, mtype = identical_mesh(meshobj, meshobj)

    assert part is True
    assert mesh is True
    assert mtype is True


@pytest.mark.parametrize("input_meshtal",
                         ['meshtal_cuv', 'meshtal_cyl',
                          'meshtal_d1s_CSimpactStudy',
                          'meshtal_d1s_IVVS_FDR',
                          'meshtal_rect_VV'])
def test_reading(input_meshtal):
    # To check if the meshtal can be read without any problem"
    filetype = 'MCNP'
    with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
        Meshtal(inp, filetype)

    assert True


@pytest.mark.parametrize("input_meshtal",
                         ['meshtal_rect_VV',
                          'meshtal_cyl',
                          'meshtal_d1s_CSimpactStudy'])
def test_mesh_vtkwrite(input_meshtal, tmpdir):
    filetype = 'MCNP'
    with as_file(RESOURCES.joinpath(input_meshtal)) as inp:
        meshtally = Meshtal(inp, filetype)

    for i, fmesh in meshtally.mesh.items():
        fmesh.print_info()
        meshtally.readMesh(i)

        meshobj = meshtally.mesh[i]
        if meshobj.cart:
            name = 'test_.vtr'
        else:
            name = 'test_.vts'

        outfile = tmpdir.mkdir('sub'+str(i)).join(name)

        meshobj.writeVTK(outfile)
        assert True
