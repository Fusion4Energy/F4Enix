from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO

import numpy as np
import polars as pl
from polars.io.plugins import register_io_source

from f4enix.output.rssa.rssa_helpers import (
    BYTE,
    FLOAT,
    INT,
    LONG,
    NUMBER_OF_EXPECTED_VALUES,
    FileParameters,
    SurfaceParameters,
)

BYTES_PER_TRACK = 96

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


def parse_header(path: Path | str) -> FileParameters:
    """Reads the header of the RSSA file and returns the file parameters as a
    FileParameters object."""
    with open(path, "rb") as infile:
        return _parse_header_binary(infile)


def parse_tracks(path: Path | str) -> pl.DataFrame:
    """Reads the tracks recorded in the RSSA file and returns them as a polars
    DataFrame."""
    with open(path, "rb") as infile:
        # Skip the header
        _ = _parse_header_binary(infile)
        return _parse_tracks_binary(infile)


def scan_rssa_file(
    path_to_file: Path | str, default_batch_size: int = 100_000_000
) -> pl.LazyFrame:
    """Scans the RSSA to return a polars LazyFrame instead of a DataFrame.

    Parameters
    ----------
    path_to_file: Path | str
        The RSSA file to read.
    default_batch_size: int
        The default batch size to use when reading the file in chunks. This is the
        number of tracks to read at once. The optimal value depends on the size of the
        file and the available memory, but a good starting point is around 100 million
        tracks (~8 Gb).

    Examples
    --------
    >>> from f4enix.output.rssa import RSSA
    >>> from f4enix.output.rssa.rssa_reader import parse_header, scan_rssa_file
    >>> parameters = parse_header('small_cyl.w')
    >>> lazy_tracks = scan_rssa_file('small_cyl.w')
    >>> # We will only read the first 20 tracks of the file
    >>> tracks = lazy_tracks.head(20).collect()
    >>> rssa = RSSA(parameters, tracks)
    """

    def source_generator(
        with_columns: list[str] | None,
        predicate: pl.Expr | None,
        n_rows: int | None,
        batch_size: int | None,
    ) -> Iterator[pl.DataFrame]:
        """
        Generator function that creates the source.
        This function will be registered as IO source.
        """
        if batch_size is None:
            batch_size = default_batch_size

        # Initialize the reader
        with open(path_to_file, "rb") as reader:
            # Skip the header
            _ = _parse_header_binary(reader)

            while n_rows is None or n_rows > 0:
                if n_rows is not None:
                    batch_size = min(batch_size, n_rows)

                number_of_bytes_to_read = (
                    batch_size * BYTES_PER_TRACK
                )  # each particle record is 96 bytes long

                df = _parse_tracks_binary(reader, number_of_bytes_to_read)

                if df.shape[0] == 0:
                    break

                if n_rows is not None:
                    n_rows -= df.shape[0]

                # Apply the column selection and predicate filtering if provided
                if with_columns is not None:
                    df = df.select(with_columns)
                if predicate is not None:
                    df = df.filter(predicate)

                yield df

    return register_io_source(io_source=source_generator, schema=SCHEMA)


def _parse_header_binary(infile: BinaryIO) -> FileParameters:
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


def _parse_tracks_binary(
    file: BinaryIO, number_of_bytes_to_read: int = -1
) -> pl.DataFrame:
    # Read the whole remaining of the file at once, store all the bytes as a 1D np array
    data = np.fromfile(file, BYTE, number_of_bytes_to_read)

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
