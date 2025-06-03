"""
Classes and functions to read WW files.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, TextIO

import numpy as np
from numpy import float64
from numpy.typing import NDArray

from f4enix.input.ww_gvr.models import NestedList, Vectors
from f4enix.input.ww_gvr.utils import compose_b2_vectors
from f4enix.output.meshtal import Fmesh, Meshtal

CARTESIAN_COORDINATES_ID = 10
ALLOWED_NUMBER_OF_ENERGIES_IN_GVR = 1


@dataclass()
class WWHeader:
    if_: int  # File type. Always 1, unused.
    iv: int  # Time-dependent windows flag. 1/2 = no/yes
    ni: int  # Number of particle types
    nr: int  # 10/16/16 = cartesian/cylindrical/spherical coord
    probid: str  # Date and time of run
    ne: List[int]  # Number of energy windows of each particle
    nfx: int  # Number of fine meshses in i
    nfy: int  # Number of fine meshses in j
    nfz: int  # Number of fine meshses in k
    origin: List[float]  #  Bottom left corner for cart, bottom center for cyl
    ncx: int  # Number of coarse meshses in i
    ncy: int  # Number of coarse meshses in j
    ncz: int  # Number of coarse meshses in k


@dataclass()
class WWHeaderCyl(WWHeader):
    director_1: List[float]  # Vector defining the director 1
    director_2: List[float]  # Vector defining the director 2


@dataclass()
class ParseResult:
    header: WWHeader
    b2_vectors: Vectors
    energies: NestedList
    values: NestedList


def parse(path: Path) -> ParseResult:
    """
    Parse a WW file and return a ParseResult object containing all the information of
    the file.

    Parameters
    ----------
    path : Path
        Path to the WW file.

    Returns
    -------
    ParseResult
        A ParseResult object containing all the information of the file.
    """
    with open(path) as infile:
        header = _parse_header(infile)

        b2_vectors = Vectors(
            vector_i=_read_block_2_vector(infile, header.ncx),
            vector_j=_read_block_2_vector(infile, header.ncy),
            vector_k=_read_block_2_vector(infile, header.ncz),
        )

        energies, values = _read_block_3(infile, header)

    return ParseResult(
        header=header, b2_vectors=b2_vectors, energies=energies, values=values
    )


def _parse_header(infile: TextIO) -> WWHeader:
    words = infile.readline().split()
    if_ = int(words[0])
    iv = int(words[1])
    ni = int(words[2])
    nr = int(words[3])
    probid = " ".join(words[4:])

    words = infile.readline().split()
    ne = [int(word) for word in words]  # Number of energy windows of each particle

    words = infile.readline().split()
    nfx = int(float(words[0]))  # Number of fine meshes in i
    nfy = int(float(words[1]))  # Number of fine meshes in j
    nfz = int(float(words[2]))  # Number of fine meshes in k
    origin = [float(word) for word in words[3:]]

    words = infile.readline().split()
    ncx = int(float(words[0]))  # Number of coarse meshes in i
    ncy = int(float(words[1]))  # Number of coarse meshes in j
    ncz = int(float(words[2]))  # Number of coarse meshes in k

    header_args: Dict[str, Any] = {
        "if_": if_,
        "iv": iv,
        "ni": ni,
        "nr": nr,
        "probid": probid,
        "ne": ne,
        "nfx": nfx,
        "nfy": nfy,
        "nfz": nfz,
        "origin": origin,
        "ncx": ncx,
        "ncy": ncy,
        "ncz": ncz,
    }

    if nr == CARTESIAN_COORDINATES_ID:
        header = WWHeader(**header_args)
    else:  # Cylindrical coordinates
        director_1 = [float(word) for word in words[3:]]
        words = infile.readline().split()
        director_2 = [float(word) for word in words[:3]]
        header_args.update(
            {
                "director_1": director_1,
                "director_2": director_2,
            }
        )
        header = WWHeaderCyl(**header_args)

    return header


def _read_block_2_vector(
    infile: TextIO, number_of_coarse_intervals: int
) -> NDArray[float64]:
    """
    Block 2 vectors combine the info of the coarse and fine meshes, they give the
    position of each coarse mesh point and the amount of fine intervals between them as:
    [mesh_position, fine_ints,
     mesh_position, 1.0000, fine_ints,
     mesh_position, 1.0000, fine_ints,
     ...]
    """
    expected_length = 3 * number_of_coarse_intervals + 1
    vector: List[float] = []
    while len(vector) < expected_length:
        words = infile.readline().split()
        vector.extend([float(word) for word in words])

    return np.array(vector)


def _read_block_3(infile: TextIO, header: WWHeader) -> tuple[NestedList, NestedList]:
    energies: NestedList = []
    values: NestedList = []

    for particle_index in range(header.ni):
        # Read energy bins for this particle
        expected_energy_bins = header.ne[particle_index]
        energies_current_particle: List[float] = []
        while len(energies_current_particle) < expected_energy_bins:
            words = infile.readline().split()
            energies_current_particle.extend([float(word) for word in words])

        # Read the values for this particle, all the energy bins at once
        expected_value_bins = (
            header.nfx * header.nfy * header.nfz * expected_energy_bins
        )
        values_current_particle: List[float] = []
        while len(values_current_particle) < expected_value_bins:
            words = infile.readline().split()
            values_current_particle.extend([float(word) for word in words])

        # Populate the result lists
        energies.append(energies_current_particle)
        values.append(values_current_particle)

    return energies, values


def read_meshtally_file(path: Path, tally_id: int = 0) -> ParseResult:
    meshtal = Meshtal(path)

    # If no tally id is provided, select the first meshtally
    if tally_id == 0:
        tally_id = meshtal.mesh_ids[0]

    meshtal.readMesh(mesh=tally_id)
    mesh = meshtal.mesh[tally_id]

    header = _read_header_from_meshtally_file(mesh)

    b2_vectors = _calculate_b2_vectors_from_meshtally(mesh)

    energies = [[100.0]]  # Only one particle and one energy bin
    values = [list(mesh.data[0, 0, :, :, :, 0].flatten())]

    return ParseResult(
        header=header,
        b2_vectors=b2_vectors,
        energies=energies,
        values=values,
    )


def _read_header_from_meshtally_file(mesh: Fmesh) -> WWHeader:
    if_ = 1  # It is always 1
    iv = 1  # It is always 1, time-dependent windows not supported
    ni = 1  # Number of particle types, always 1 when reading a meshtally
    nr = 10 if mesh.geom == "rec" else 16
    if len(mesh.ebin) - 1 > ALLOWED_NUMBER_OF_ENERGIES_IN_GVR:
        raise ValueError("The meshtally file has more than one energy group.")
    ne = [1]  # Number of energy bins of each particle, always 1 for a meshtally

    # When reading a meshtally we consider that there are no regular intervals between
    # coarse vectors. This means that the length of a coarse vector is equal to the
    # total amount of fine ints at that dimension
    nfx = ncx = len(mesh.x1bin) - 1
    nfy = ncy = len(mesh.x2bin) - 1
    nfz = ncz = len(mesh.x3bin) - 1
    if mesh.trsf and any(mesh.trsf.origin):
        origin = [mesh.trsf.origin[2], mesh.trsf.origin[1], mesh.trsf.origin[0]]
    else:
        origin = [0.0, 0.0, 0.0]

    header_args: Dict[str, Any] = {
        "if_": if_,
        "iv": iv,
        "ni": ni,
        "nr": nr,
        "probid": "Generated by f4enix_ww",
        "ne": ne,
        "nfx": nfx,
        "nfy": nfy,
        "nfz": nfz,
        "origin": origin,
        "ncx": ncx,
        "ncy": ncy,
        "ncz": ncz,
    }

    if nr == CARTESIAN_COORDINATES_ID:
        # Origin should be zero for cartesian as the vectors start already at the origin
        header_args["origin"] = [0, 0, 0]
        header = WWHeader(**header_args)
    else:  # Cylindrical coordinates
        director_1 = mesh.trsf.axis.tolist()
        if mesh.trsf.vec is None:  # MCNP5 Meshtally have no info about vec
            if director_1[0] != 1:
                director_2 = [1.0, 0.0, 0.0]  # Default MCNP vec
            else:
                director_2 = [0.0, 1.0, 0.0]  # Cant be the same as the axis
        else:
            director_2 = mesh.trsf.vec.tolist()
        header_args.update(
            {
                "director_1": director_1,
                "director_2": director_2,
            }
        )
        header = WWHeaderCyl(**header_args)

    return header


def _calculate_b2_vectors_from_meshtally(mesh: Fmesh) -> Vectors:
    coarse_vectors = Vectors(
        vector_i=np.array(mesh.x1bin),
        vector_j=np.array(mesh.x2bin),
        vector_k=np.array(mesh.x3bin),
    )

    fine_vectors = Vectors(
        vector_i=np.ones(mesh.nx1),
        vector_j=np.ones(mesh.nx2),
        vector_k=np.ones(mesh.nx3),
    )

    return compose_b2_vectors(coarse_vectors, fine_vectors)


def write(
    file_path: Path,
    header: WWHeader,
    b2_vectors: Vectors,
    energies: NestedList,
    values: NestedList,
) -> None:
    """
    Write a WW file from the given information.

    Parameters
    ----------
    file_path : Path
        Path to the file to be written.
    header : WWHeader
        The header information of the WW file.
    b2_vectors : Vectors
        The block 2 vectors of the WW file.
    energies : NestedList
        The energies of the WW file.
    values : NestedList
        The values of the WW file.

    Returns
    -------
    None
    """
    with open(file_path, "w") as infile:
        _write_header(infile, header)
        _write_block_2(infile, b2_vectors)
        _write_block_3(infile, energies, values)


def _write_header(f: TextIO, header: WWHeader) -> None:
    f.write(f"{header.if_:>10}{header.iv:>10}{header.ni:>10}{header.nr:>10}")
    f.write(f"{header.probid:>38}\n")

    for e_bins in header.ne:
        f.write(f"{e_bins:>10}")
    f.write("\n")

    f.write(f"{header.nfx:>#9.5g}{header.nfy:>#13.5g}{header.nfz:>#13.5g}")
    f.write(f"{header.origin[0]:>#13.5g}")
    f.write(f"{header.origin[1]:>#13.5g}")
    f.write(f"{header.origin[2]:>#13.5g}\n")

    f.write(f"{header.ncx:>#9.5g}{header.ncy:>#13.5g}{header.ncz:>#13.5g}")
    if isinstance(header, WWHeaderCyl):
        f.write(f"{header.director_1[0]:>#13.5g}")
        f.write(f"{header.director_1[1]:>#13.5g}")
        f.write(f"{header.director_1[2]:>#13.5g}\n")
        f.write(f"{header.director_2[0]:>#9.5g}")
        f.write(f"{header.director_2[1]:>#13.5g}")
        f.write(f"{header.director_2[2]:>#13.5g}")
        f.write(f"{2:>#13.5g}")
    else:
        f.write(f"{1:>#13.5g}")
    f.write("\n")


def _write_block_2(f: TextIO, b2_vectors: Vectors) -> None:
    for b2_vector in b2_vectors:
        vector = b2_vector.tolist()
        packs_of_six = [vector[i : i + 6] for i in range(0, len(vector), 6)]

        for pack in packs_of_six:
            first_number = pack.pop(0)
            f.write(f"{first_number:>#9.5g}")
            for number in pack:
                f.write(f"{number:>#13.5g}")
            f.write("\n")


def _write_block_3(f: TextIO, energies: NestedList, values: NestedList) -> None:
    for particle_index in range(len(energies)):
        particle_ergs = energies[particle_index]
        packs_of_6 = [particle_ergs[i : i + 6] for i in range(0, len(particle_ergs), 6)]
        for energy_pack in packs_of_6:
            first_energy = energy_pack.pop(0)
            f.write(f"{first_energy:>#9.5g}")
            for energy in energy_pack:
                f.write(f"{energy:>#13.5g}")
            f.write("\n")

        particle_vals = values[particle_index]
        packs_val = [particle_vals[i : i + 6] for i in range(0, len(particle_vals), 6)]
        for value_pack in packs_val:
            for value in value_pack:
                f.write(f"{value:>#13.4E}")
            f.write("\n")
