"""
Parsing of MCNP eeout files.

The main feature is the conversion of the results to a .vtk file.
"""

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or - 
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at: 
    https://eupl.eu/1.2/en/
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""

import re
import os
import logging
import pyvista as pv
import numpy as np

from f4enix.constants import PAT_DIGIT, PAT_SPACE, SCIENTIFIC_PAT

# Identifiers
# idNodes = 'NUMBER OF NODES'
# idTets = 'NUMBER OF 1st TETS'
# id2Tets = 'NUMBER OF 2nd TETS'
# idParticlesType = 'PARTICLE LIST'
# idNodeX = 'NODES X'
# idNodeY = 'NODES Y'
# idNodeZ = 'NODES Z'
# idElem = 'ELEMENT TYPE'
# idConn = 'CONNECTIVITY DATA 1ST ORDER TETS ELEMENT ORDERED'
# id2Conn = 'CONNECTIVITY DATA 2ND ORDER TETS ELEMENT ORDERED'
# idNeighbour = 'NEAREST NEIGHBOR DATA 1ST ORDER TETS'
# id2Neighbour = 'NEAREST NEIGHBOR DATA 2ND ORDER TETS'
# Common patterns
PAT_INFO_NAME = re.compile(r"[A-Za-z\s\d]+")
PAT_VALUE = re.compile(r"\s[\s-]\d")
# patNumber = re.compile(r'\d+')
# patNumberSci = re.compile(r'[-+]*\d+.\d+E[+-]\d+')
# patSpace = re.compile(r'\s+')
# patFlux = re.compile(r'd*4$')
# patHeat = re.compile(r'd*6$')
# # special patterns
# patTets = re.compile(r'\d\d+')
TETRA1 = "1st tetra"
TETRA2 = "2nd tetra"


class EEOUT:

    def __init__(self, filepath: os.PathLike) -> None:
        """Representation of an MCNP .eeout output file

        This class is used to parse and manipulate a .eeout file. The main
        feature is the possibility to convert the results to a .vtk format

        Parameters
        ----------
        filepath : os.PathLike
            path to the .eeout file to be parsed.

        Attributes
        ----------
        filename : str
            name of the original .eeout file
        lines : list[str]
            lines of the file stored in memory
        info : dict
            contains general informations about the file
        n_nodes : int
            number of nodes of the mesh
        n_elem : int
            number of elements of the mesh
        grid : pv.UnstructuredGrid
            vtk object containing all the parsed data from eeout

        Examples
        --------
        Parse an eeout file and export it to .vtu format

        >>> from f4enix.output.eeout import EEOUT
        ... eeout = EEOUT(r'path/to/file.eeout')
        ... eeout.export(r'path/to/output/folder')

        Raises
        ------
        NotImplementedError
            Only 1st and 2nd order tetras are supported
        """
        self.filename = os.path.basename(filepath).split(".")[0]
        logging.info(f"Parsing {self.filename}")
        # store file lines
        lines = []
        with open(filepath, "r", errors="surrogateescape") as infile:
            for line in infile:
                lines.append(line)
        self.lines = lines

        # Get infos
        self.info = self._read_info()

        # get materials name
        mat_names = self._get_materials_name()

        # Determine mesh elements used. Only first and second order tetras
        # are supported. No mixed formulation are supported
        tetra1 = self.info["NUMBER OF 1st TETS"]
        tetra2 = self.info["NUMBER OF 2nd TETS"]
        self.n_nodes = self.info["NUMBER OF NODES"]
        if tetra1 > 0 and tetra2 == 0:
            self.elem_type = TETRA1
            self.n_elem = tetra1
        elif tetra2 > 0 and tetra1 == 0:
            self.elem_type = TETRA2
            self.n_elem = tetra2
        else:
            raise NotImplementedError(
                "Only 1st and 2nd order tetra are supported. No mixed formulation"
            )

        # get the particle list
        self.p_list = self._read_particle_list()
        assert len(self.p_list) == self.info["NUMBER OF PARTICLES"]

        points, idx_elem_type = self._read_nodes_xyz()

        # get the materialss
        material_ids, idx_connect = self._read_material(idx_init=idx_elem_type)

        # TODO the element type should be parsed here in order to allow mixed
        # formulations

        # get connectivity
        elem_connectivity, idx_neigh = self._read_connectivity(start_idx=idx_connect)

        # get edits
        edits, idx_centroids = self._read_edits(start_idx=idx_neigh)

        # read density and volumes
        rhos, volumes = self._read_rho_vol(start_index=idx_centroids)

        # build the UM
        if self.elem_type == TETRA1:
            cell_type = np.full(self.n_elem, pv.CellType.TETRA, dtype=np.uint8)
        elif self.elem_type == TETRA2:
            cell_type = np.full(self.n_elem, pv.CellType.TETRA, dtype=np.uint8)
        else:
            raise NotImplementedError(f"{self.elem_type} not implemented")

        grid = pv.UnstructuredGrid(elem_connectivity.ravel(), cell_type, points)
        grid["material"] = np.vectorize(mat_names.get)(material_ids)
        grid["density"] = rhos
        grid["volume"] = volumes
        for key, field in edits.items():
            grid[f"{key} - value"] = field["values"]
            grid[f"{key} - error"] = field["errors"]

        self.grid = grid

        logging.info("Parsing completed correctly")

    def export(
        self, outfolder: os.PathLike, filename: str = None, format: str = "vtu"
    ) -> None:
        """Export the eeout to different format.

        at the moment only unstructured mesh vtk (vtu) is the only supported
        format.

        Parameters
        ----------
        outfolder : os.PathLike
            path to the output folder where to export.
        filename : str, optional
            name of the file (without format extension), by default None,
            meaning the the original name of the eeout file will be used.
        format : str, optional
            export format, by default 'vtu' which is only one that is currently
            supported.

        Raises
        ------
        ValueError
            if the output folder does not exists.
        """
        if not os.path.exists(outfolder):
            raise ValueError(f"{outfolder} does not exist")

        if format not in ["vtu"]:
            raise NotImplementedError(f"{format} is not a valid format")

        if filename is None:
            filename = self.filename

        outpath = os.path.join(outfolder, filename + "." + format)
        self.grid.save(outpath)

    def _read_info(self) -> dict:
        read = False
        infos = {}
        for line in self.lines:
            # read infos between n. of particles and number of com edits
            if line.find("NUMBER OF PARTICLES") != -1:
                read = True

            if read:
                info_name = PAT_INFO_NAME.match(line).group().strip()
                value = PAT_DIGIT.findall(line)[-1]
                infos[info_name] = int(value)

            if line.find("NUMBER OF EDITS") != -1:
                return infos

        raise RuntimeError("No 'NUMBER OF EDITS' tag was found")

    def _get_materials_name(self) -> dict[int, str]:
        n_materials = self.info["NUMBER OF MATERIALS"]
        for idx, line in enumerate(self.lines):
            pat_flag = re.compile(r"\s+MATERIALS")
            if pat_flag.match(line):
                break

        pat_whatever = re.compile(".+")
        materials = {}
        for i in range(n_materials):
            line = self.lines[i + 1 + idx]
            mat = pat_whatever.match(line).group().strip()
            materials[i + 1] = mat

        return materials

    def _read_particle_list(self) -> int:
        readFlag = False
        for line in self.lines:
            if readFlag:
                particleList = PAT_DIGIT.findall(line)
                particleList = [int(x) for x in particleList]  # get ints
                return particleList

            if line.find("PARTICLE LIST") != -1:
                readFlag = True

        raise RuntimeError("Particles list was not found")

    def _read_material(self, idx_init: int = 0) -> tuple[list[int], int]:
        materialFlag = False
        materialsList = []

        for idx, line in enumerate(self.lines[idx_init:]):
            # No more materials to be read
            if line.find("CONNECTIVITY DATA") != -1:
                del materialsList[-6:]
                return materialsList, idx + idx_init

            if materialFlag:
                strippedline = line.strip()
                matLine = PAT_SPACE.split(strippedline)
                for mat in matLine:
                    materialsList.append(int(mat))

            # trigger
            if line.find("ELEMENT MATERIAL") != -1:
                materialFlag = True

    def _read_nodes_xyz(self) -> tuple[np.ndarray, int]:

        # Reading nodes
        readFlagX = False
        readFlagY = False
        readFlagZ = False

        nodesX = []
        nodesY = []
        nodesZ = []

        pat_trigger_end = re.compile(" ELEMENT TYPE")

        for idx, line in enumerate(self.lines):

            # Reading nodes
            if readFlagX:
                split = PAT_SPACE.split(line)
                for string in split:
                    a = SCIENTIFIC_PAT.search(string)
                    if a is not None:
                        nodesX.append(float(a.group()))

            if readFlagY:
                split = PAT_SPACE.split(line)
                for string in split:
                    a = SCIENTIFIC_PAT.search(string)
                    if a is not None:
                        nodesY.append(float(a.group()))

            if readFlagZ:
                split = PAT_SPACE.split(line)
                for string in split:
                    a = SCIENTIFIC_PAT.search(string)
                    if a is not None:
                        nodesZ.append(float(a.group()))

            if line.find("NODES X") != -1:
                readFlagX = True
            if line.find("NODES Y") != -1:
                readFlagY = True
                readFlagX = False
            if line.find("NODES Z") != -1:
                readFlagZ = True
                readFlagY = False

            # trigger exit
            if pat_trigger_end.match(line) is not None:
                points = np.concatenate(
                    [np.array([nodesX]).T, np.array([nodesY]).T, np.array([nodesZ]).T],
                    axis=1,
                )
                return points, idx

        raise RuntimeError("Could not find the element tag 'ELEMENT TYPE'")

    def _read_connectivity(self, start_idx: int = 0) -> tuple[np.ndarray, int]:
        data = []
        # Reading nodes
        readFlagCon = False

        for idx, line in enumerate(self.lines[start_idx:]):
            if readFlagCon:
                data.append(line)
            if line.find("CONNECTIVITY DATA") != -1:
                readFlagCon = True
                if line.find("NODE ORDERED") != -1:
                    raise NotImplementedError("Only element ordered is supported")

            # trigger exit
            if line.find("NEAREST NEIGHBOR DATA") != -1:
                readFlagCon = False
                break

        del data[-1]
        del data[-1]
        # determine number of nodes per element
        if self.elem_type == TETRA1:
            n_nodes = 4
        elif self.elem_type == TETRA2:
            n_nodes = 10

        # first get all values and straigthen the array
        # (solve multiple line issue)
        conn = []
        for string in data:
            newline = string.strip()
            substrings = PAT_SPACE.split(newline)
            for substring in substrings:
                num = int(PAT_DIGIT.search(substring).group()) - 1
                conn.append(num)

        conn = np.reshape(np.array(conn), (-1, n_nodes))
        # add the n_nodes describer
        describer = np.full((len(conn), 1), n_nodes)

        # finally concat
        connectivity = np.concatenate([describer, conn], axis=1)

        return connectivity, idx + start_idx

    def _read_edits(self, start_idx: int = 0) -> tuple[dict[str, dict[str, list]], int]:

        pat_particle = re.compile(r"(?<=DATA OUTPUT PARTICLE :)\s+\d+")
        pat_edit = re.compile(r"(?<=EDIT LIST :)\s+\d+")
        pat_type = re.compile(r"(?<=TYPE :)\s+\w+")

        valuesFlag = False

        edits = {}

        for idx, line in enumerate(self.lines[start_idx:]):
            # a new edit is found, initialize it
            if line.find("DATA OUTPUT PARTICLE") != -1:
                par = int(pat_particle.search(line).group())
                edit_num = int(pat_edit.search(line).group())
                edit_type = pat_type.search(line).group().strip()
                current_tally = f"Tally{edit_num}_par{par}_{edit_type}"
                logging.info(f"parsing {current_tally}")
                edits[current_tally] = {}

            # We need to read data (either values or errors)
            if valuesFlag:
                # valid data line
                if PAT_VALUE.match(line) is not None:
                    strippedline = line.strip()
                    values = PAT_SPACE.split(strippedline)
                    for value in values:
                        values_list.append(float(value))
                # no more values in this block
                else:
                    edits[current_tally][value_type] = values_list[1:]
                    valuesFlag = False

            # triggers for following line
            if (
                line.find("DATA SETS RESULT TIME BIN") != -1
                or line.find("DATA SETS REL ERROR TIME BIN") != -1
            ):
                valuesFlag = True
                values_list = []
                if "RESULT" in line:
                    value_type = "values"
                elif "REL ERROR" in line:
                    value_type = "errors"
                else:
                    value_type = None

            # triggers the end of the data section
            if line.find("CENTROIDS X") != -1:
                return edits, idx + start_idx

        raise RuntimeError("Could not parse edits properly, missing 'CENTROIDS X'")

    def _read_rho_vol(self, start_index: int = 0) -> tuple[list[float], list[float]]:
        densityFlag = False
        volumeFlag = False
        densityList = []
        volumesList = []
        for line in self.lines[start_index:]:
            if densityFlag or volumeFlag:
                # valid data line
                if PAT_VALUE.match(line) is not None:
                    strippedline = line.strip()
                    splitline = PAT_SPACE.split(strippedline)
                    for value in splitline:
                        if densityFlag:
                            densityList.append(float(value))
                        elif volumeFlag:
                            volumesList.append(float(value))

            if line.find("DENSITY") != -1:
                densityFlag = True
            if line.find("VOLUMES") != -1:
                volumeFlag = True
                densityFlag = False

        return densityList, volumesList

    def __repr__(self) -> str:
        return str(self.info)

    def __str__(self) -> str:
        return str(self.info)
