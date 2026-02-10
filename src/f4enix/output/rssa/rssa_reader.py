from typing import BinaryIO

import numpy as np
import polars as pl

from f4enix.output.rssa.rssa_helpers import (
    BYTE,
    FLOAT,
    INT,
    LONG,
    NUMBER_OF_EXPECTED_VALUES,
    FileParameters,
    SurfaceParameters,
)

SCHEMA = pl.Schema(
    {
        "a": int,
        "b": int,
        "wgt": float,
        "erg": float,
        "tme": float,
        "x": float,
        "y": float,
        "z": float,
        "u": float,
        "v": float,
        "c": int,
    },
)


def parse_header(infile: BinaryIO) -> FileParameters:
    first_record = _read_fortran_record(infile)
    # The first line of the file with information like the code version, date and title
    formatted_record_id = first_record.tobytes().decode("UTF-8")
    if "d1suned" in formatted_record_id:
        _last_dump = np.frombuffer(first_record[-4:], INT)
    elif "SF_00001" in formatted_record_id:
        _header = _read_fortran_record(infile)  # code version and other info
    else:
        raise NotImplementedError(
            f"The code that generated this RSSA file has not been implemented"
            f" in this parser, see the code here: {formatted_record_id}..."
        )

    second_record = _read_fortran_record(infile)
    np1 = np.frombuffer(second_record, LONG, 1, 0)[0]
    nrss = np.frombuffer(second_record, LONG, 1, 8)[0]
    nrcd = np.frombuffer(second_record, INT, 1, 16)[0]
    njsw = np.frombuffer(second_record, INT, 1, 20)[0]
    niss = np.frombuffer(second_record, LONG, 1, 24)[0]
    if abs(nrcd) != NUMBER_OF_EXPECTED_VALUES:
        raise NotImplementedError(
            "The amount of values recorded for each particle should be 11 instead of"
            f" {nrcd}..."
        )

    if np1 < 0:
        third_record = _read_fortran_record(infile)
        niwr, mipts, kjaq = np.frombuffer(third_record, INT, 3)
    else:
        raise NotImplementedError("The np1 value is not negative...")

    surfaces = []
    for _ in range(njsw):
        data = _read_fortran_record(infile)
        surf_id = np.frombuffer(data, INT, 1, 0)[0]
        surf_info = np.frombuffer(data, INT, 1, 4)[0] if kjaq == 1 else -1
        surf_type = np.frombuffer(data, INT, 1, 8)[0]
        num_parameters = np.frombuffer(data, INT, 1, 12)[0]
        parameters = np.frombuffer(data, INT, offset=16).tolist()
        surfaces.append(
            SurfaceParameters(
                id=surf_id,
                info=surf_info,
                type=surf_type,
                num_parameters=num_parameters,
                parameters=parameters,
            )
        )

    # we read any extra records as determined by njsw+niwr...
    # no known case of their actual utility
    for _j in range(njsw, njsw + niwr):
        _read_fortran_record(infile)
        raise NotImplementedError(
            "njsw + niwr values are bigger than njsw, behavior not explained"
        )

    # Summary record
    _data = _read_fortran_record(infile)
    # Summary record not processed, its information does not interest us for now

    return FileParameters(
        np1=np1,  # Number of histories of the simulation, given as a negative number
        nrss=nrss,  # Number of tracks recorded
        nrcd=nrcd,  # Number of values recorded for each particle, it should be 11
        njsw=njsw,  # Number of surfaces in JASW
        niss=niss,  # Number of different histories that reached the SSW surfaces
        niwr=niwr,  # Number of cells in RSSA file
        mipts=mipts,  # Source particle type
        kjaq=kjaq,  # Flag for macrobodies surfaces
        surfaces=surfaces,
    )


def parse_tracks(file: BinaryIO) -> pl.DataFrame:
    # Read the whole remaining of the file at once, store all the bytes as a 1D np array
    data = np.fromfile(file, BYTE)

    # Reshape the array so each index holds the information of a single particle
    # we can do this because we know that the particle records have always the same
    # length, 96 bytes
    data = data.reshape(-1, 96)

    # Remove the first and last 4 bytes, these are two integers that tell the record is
    # 88 bytes long
    data = data[:, 4:-4]

    # Convert the array into a 1D array of float numbers instead of simply bytes
    data = np.frombuffer(data.flatten(), FLOAT)

    # Reshape the array so each index holds the information of a single particle
    # all the data is already converted from bytes to floats
    data = data.reshape(-1, 11)

    # Build the DataFrame
    df = pl.DataFrame(data, schema=SCHEMA)

    # Modify the value of "b" for fast filtering of neutrons and photons
    df = df.with_columns(
        (pl.col("b").abs() / (10 ** pl.col("b").abs().log10().floor()))
        .cast(int)
        .alias("b")
    )
    return df


def _read_fortran_record(infile: BinaryIO):
    count_1 = np.fromfile(infile, INT, 1)[0]
    data = np.fromfile(infile, np.byte, count_1)
    count_2 = np.fromfile(infile, INT, 1)[0]
    if count_1 != count_2:
        raise ValueError(
            "The integers that go before and after the Fortran record are not equal..."
        )
    return data
