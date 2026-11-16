"""Readers of the files `Meshtal` reads, backed by aci-tallies: MCNP/D1SUNED
meshtal files, D1SUNED cell-under-voxel (cuvmsh) files and CDGS files."""

from abc import ABC, abstractmethod
from pathlib import Path

from aci_tallies import Input as AciInput
from aci_tallies import cdgs as aci_cdgs
from aci_tallies import cuv as aci_cuv
from aci_tallies import meshtal as aci_meshtal

from f4enix.output.meshtal.fmesh import Fmesh


class Parser(ABC):
    @abstractmethod
    def __init__(self, filename: str | Path): ...

    @abstractmethod
    def get_meshlist(self) -> tuple: ...

    @abstractmethod
    def get_FMesh(
        self, tally: int, norm: str | None = None, filter: list[int] | None = None
    ) -> Fmesh: ...


class AciMeshtalParser(Parser):
    """Read an MCNP/D1SUNED meshtal file (block and column formats, rectangular and
    cylindrical meshes, D1SUNED nuclide, cell and decay time bins)."""

    def __init__(self, filename: str | Path, input: AciInput | None = None):
        self.filename = filename
        self._file = aci_meshtal.open(filename, input=input)
        self._failed = {id: why for id, why in self._file.failed}

    def get_meshlist(self) -> tuple:
        """Return the tally numbers of the meshes in the file, the ones that cannot
        be read included.

        Returns
        -------
        tuple
            A tuple containing the tally numbers (as strings) of the meshes.
        """
        ids = set(self._file.ids) | set(self._failed)
        return tuple(str(id) for id in sorted(ids))

    def get_FMesh(
        self,
        tally: int,
        norm: str | None = None,
        filter: list[int] | None = None,
    ) -> Fmesh:
        """Returns Fmesh object representing a meshtally.

        Parameters
        ----------
        tally : int
            The mesh tally number to retrieve.
        norm : str | None, optional
            Normalization option for CUV meshes, not used.
        filter : list[int] | None, optional
            List of cells with which voxels in CUV meshes are filtered, not used.

        Returns
        -------
        Fmesh
            An Fmesh object containing the mesh data.

        Raises
        ------
        ValueError
            if the mesh tally cannot be read (e.g. D1SUNED source importance
            meshes).
        """
        tally = int(tally)
        if tally in self._failed:
            raise ValueError(
                f"Mesh tally {tally} of {self.filename} cannot be read: "
                f"{self._failed[tally]}"
            )
        if tally not in self._file:
            raise KeyError("bad tally entry")
        fm = Fmesh.from_aci(self._file[tally])
        fm.type = "meshtal"
        return fm


class AciCUVMeshParser(Parser):
    """Read a D1SUNED cell-under-voxel (cuvmsh) file."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self._file = aci_cuv.open(filename)

    def get_meshlist(self) -> tuple:
        """Return list of mesh stored in the CUV file.

        Returns
        -------
        tuple
            A tuple containing the tally numbers (as strings) of the meshes.
        """
        return tuple(str(id) for id in sorted(self._file.ids))

    def get_FMesh(
        self,
        tally: int,
        norm: str | None = None,
        filter: list[int] | None = None,
    ) -> Fmesh:
        """Returns Fmesh object representing a CUV mesh: each voxel's value is the
        sum over its cells of their value times their volume fraction.

        Parameters
        ----------
        tally : int
            The mesh tally number to retrieve.
        norm : str | None, optional
            'ctot' multiplies each voxel by its volume, 'celf' divides it by the
            volume fraction of the cells kept. Default is None.
        filter : list[int] | None, optional
            Only these cells are summed over. Default is None (all of them).

        Returns
        -------
        Fmesh
            An Fmesh object containing the CUV mesh data.
        """
        if norm not in (None, "ctot", "celf"):
            norm = None  # as the Python parser did, ignore unknown options
        mesh = self._file.to_mesh(int(tally), cells=filter, norm=norm)
        fm = Fmesh.from_aci(mesh)
        fm.type = "CUV"
        return fm


class AciCDGSMeshParser(Parser):
    """Read a D1SUNED Common Decay Gamma Source (CDGS) file, rectangular and
    cylindrical meshes."""

    def __init__(self, filename: str | Path):
        self.filename = filename
        self._file = aci_cdgs.open(filename)

    def get_meshlist(self) -> tuple:
        """Return list of mesh stored in the CDGS file.

        Returns
        -------
        tuple
            A tuple containing the CDGS mesh ids (as strings).
        """
        ids = set(self._file.ids) | {id for id, _ in self._file.failed}
        return tuple(str(id) for id in sorted(ids))

    def get_FMesh(
        self, tally: int, norm: str | None = None, filter: list[int] | None = None
    ) -> Fmesh:
        """Return Fmesh object representing a CDGS source: per voxel and energy bin
        (or line), the sum of its cells' source.
        The tally number is the source mesh number.

        Parameters
        ----------
        tally : int
            The source mesh number to retrieve.
        norm : str | None, optional
            Normalization option, not used in CDGS meshes.
        filter : list[int] | None, optional
            Filter options, not used in CDGS meshes.

        Returns
        -------
        Fmesh
            An Fmesh object containing the mesh data."""
        tally = int(tally)
        if tally not in self._file and tally not in {i for i, _ in self._file.failed}:
            raise KeyError("bad tally entry")
        info = self._file.info(tally)
        fm = Fmesh.from_aci(self._file[tally])
        fm.type = "CDGS"
        fm.particle = None
        fm.cooling_time = info.cooling_time
        fm.strength = info.strength
        fm.comments = info.comment
        return fm
