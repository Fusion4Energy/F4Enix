from pathlib import Path
from typing import List, Optional


import numpy as np

from f4enix.input.ww_gvr import ww_parser
from f4enix.input.ww_gvr.geometry import Geometry
from f4enix.input.ww_gvr.models import (
    CoordinateType,
    EnergiesByParticle,
    NestedList,
    ParticleType,
    Pathlike,
    ValuesByParticle,
)
from f4enix.input.ww_gvr.ww_parser import ParseResult, WWHeader, WWHeaderCyl
from f4enix.input.ww_gvr.ratios_calculation import calculate_max_ratio_array
from f4enix.input.ww_gvr.utils import decompose_b2_vectors


class WW:
    def __init__(
        self,
        file_path: Path,
        geometry: Geometry,
        values: ValuesByParticle,
    ):
        self.file_path = file_path
        self.geometry = geometry
        self.values = values

    @classmethod
    def load_from_ww_file(cls, file_path: Pathlike):
        file_path = Path(file_path)
        parse_result = ww_parser.parse(file_path)

        geometry = cls._create_geometry(parse_result)
        values = cls._nested_to_values_by_particle(parse_result)

        return WW(file_path, geometry, values)

    @classmethod
    def create_gvr_from_meshtally_file(
        cls, file_path: Pathlike, maximum_splitting_ratio: float = 5.0
    ):
        file_path = Path(file_path)
        parse_result = ww_parser.read_meshtally_file(file_path)

        geometry = cls._create_geometry(parse_result)
        values = cls._nested_to_values_by_particle(parse_result)
        cls._convert_flux_values_to_gvr(values, maximum_splitting_ratio)

        return WW(file_path, geometry, values)

    @staticmethod
    def _create_geometry(parse_result: ParseResult) -> Geometry:
        coarse_vectors, fine_vectors = decompose_b2_vectors(parse_result.b2_vectors)
        return Geometry(parse_result.header, coarse_vectors, fine_vectors)

    @staticmethod
    def _nested_to_values_by_particle(parse_result: ParseResult) -> ValuesByParticle:
        particles = [ParticleType.NEUTRON]
        if parse_result.header.ni > 1:
            particles.append(ParticleType.PHOTON)

        values_by_particle: ValuesByParticle = {}
        for i, particle in enumerate(particles):
            values_by_particle[particle] = {}
            energies = parse_result.energies[i]
            reshaped_values = np.array(parse_result.values[i]).reshape(
                len(energies),
                parse_result.header.nfz,
                parse_result.header.nfy,
                parse_result.header.nfx,
            )

            for j, energy in enumerate(energies):
                values_by_particle[particle][energy] = reshaped_values[j]

        return values_by_particle

    def _format_nested_energies_and_values(self) -> tuple[NestedList, NestedList]:
        """The opposite of _nested_to_values_by_particle. Returns energies and values
        in the format that the writer of ww_parser expects."""
        energies = []
        values = []

        for particle in self.particles:
            particle_energies = list(self.values[particle].keys())
            particle_energies = np.array(particle_energies).flatten().tolist()
            energies.append(particle_energies)

            particle_values = list(self.values[particle].values())
            particle_values = np.array(particle_values).flatten().tolist()
            values.append(particle_values)

        return energies, values

    def _calculate_ratios(self) -> ValuesByParticle:
        ratios_by_particle = {}

        for particle in self.particles:
            ratios_by_particle[particle] = {}

            for energy in self.energies[particle]:
                ratios_by_particle[particle][energy] = calculate_max_ratio_array(
                    self.values[particle][energy]
                )

        return ratios_by_particle

    @staticmethod
    def _convert_flux_values_to_gvr(
        values: ValuesByParticle, maximum_splitting_ratio: float
    ) -> None:
        # Only one particle is expected when reading a meshtally file
        all_values = np.concatenate(list(values[ParticleType.NEUTRON].values()))
        max_flux = np.max(all_values)

        for energy in values[ParticleType.NEUTRON].keys():
            # Van Vick / Andrew Davis / Magic algorithm
            values[ParticleType.NEUTRON][energy] /= max_flux
            values[ParticleType.NEUTRON][energy] *= 2 / (maximum_splitting_ratio + 1)

    @property
    def values(self) -> ValuesByParticle:
        return self._values

    @values.setter
    def values(self, values: ValuesByParticle):
        self._values = values
        self._ratios = self._calculate_ratios()
        self.geometry.fill_grid_values(self._values)
        self.geometry.fill_grid_ratios(self._ratios)

    @property
    def ratios(self) -> ValuesByParticle:
        return self._ratios

    @property
    def energies(self) -> EnergiesByParticle:
        energies = {}
        for particle in self.particles:
            energies[particle] = list(self._values[particle].keys())
        return energies

    @property
    def particles(self) -> List[ParticleType]:
        particles = [ParticleType.NEUTRON]
        if len(self.values) > 1:
            particles.append(ParticleType.PHOTON)
        return particles

    @property
    def info(self) -> str:
        text = f" {self.file_path.name} weight window\n"

        vi = self.geometry._coarse_vectors.vector_i
        vj = self.geometry._coarse_vectors.vector_j
        vk = self.geometry._coarse_vectors.vector_k
        text += f"       {'From':-^12} {'To':-^12} {'No. Bins':-^12}\n"
        text += f" I --> {vi[0]:^12.1f} {vi[-1]:^12.1f} {self.geometry.i_ints:^12.0f}\n"
        text += f" J --> {vj[0]:^12.1f} {vj[-1]:^12.1f} {self.geometry.j_ints:^12.0f}\n"
        text += f" K --> {vk[0]:^12.1f} {vk[-1]:^12.1f} {self.geometry.k_ints:^12.0f}\n"

        coordinates = self.geometry.coordinate_type.name.capitalize()
        no_voxels = self.geometry.i_ints * self.geometry.j_ints * self.geometry.k_ints
        text += f"\n {coordinates} coordinates, {no_voxels} voxels.\n"

        text += f"\n The weight window contains {len(self.particles)} particle/s"
        for particle in self.particles:
            text += f"\n\n {particle.name.capitalize():-^50}\n"

            all_values = np.concatenate(list(self.values[particle].values()))
            min_value = np.min(all_values)
            max_value = np.max(all_values)
            amount_positive_values = int(np.sum(all_values > 0))
            percent_positive_values = amount_positive_values / all_values.size * 100
            all_ratios = np.concatenate(list(self.ratios[particle].values()))
            average_ratio = np.average(all_ratios)
            max_ratio = np.max(all_ratios)

            text += f" {'Energy bins:':<20} {self.energies[particle]}"
            text += f"\n {'Min value:':<20} {min_value:.2E}"
            text += f"\n {'Max value:':<20} {max_value:.2E}"
            text += f"\n {'No.Bins > 0 [%]:':<20} {percent_positive_values:.1f}%"
            text += f"\n {'Average ratio:':<20} {average_ratio:.2E}"
            text += f"\n {'Max ratio:':<20} {max_ratio:.2E}"

        return text

    def __repr__(self) -> str:
        return self.info

    def multiply(self, factor: float) -> None:
        if factor == 1.0:
            return

        values = self.values
        for particle in self.particles:
            for energy in self.energies[particle]:
                values[particle][energy] *= factor

        # This line will trigger the setter, which will update the ratios and geometry
        self.values = values

    def soften(self, softening_factor: float) -> None:
        if softening_factor == 1.0:
            return

        values = self.values
        for particle in self.particles:
            for energy in self.energies[particle]:
                values[particle][energy] **= softening_factor

        # This line will trigger the setter, which will update the ratios and geometry
        self.values = values

    def add_particle(self, norm: float = 1, soft: float = 1):
        if len(self.particles) == 2:
            raise ValueError("There are already two particles...")

        values_photon = {
            energy: value * norm**soft
            for energy, value in self.values[ParticleType.NEUTRON].items()
        }

        values = self.values
        values.update({ParticleType.PHOTON: values_photon})
        # This line will trigger the setter, which will update the ratios and geometry
        self.values = values

    def remove_particle(self):
        if len(self.particles) == 1:
            raise ValueError("There is only one particle already...")

        self.values = {ParticleType.NEUTRON: self.values[ParticleType.NEUTRON]}

    def mitigate_long_histories(self, max_ratio: float = 10.0) -> None:
        """The weigth window voxels with a ratio higher than max_ratio will be set to
        zero. This stopes the particle from over-splitting."""
        values = self.values
        for particle in self.particles:
            for energy in self.energies[particle]:
                mask_high_ratios = self.ratios[particle][energy] > max_ratio
                values[particle][energy][mask_high_ratios] = 0.0

        # This line will trigger the setter, which will update the ratios and geometry
        self.values = values

    def write_to_ww_file(self, file_path: Optional[Pathlike] = None) -> None:
        if file_path is None:
            file_path = self.file_path.parent / (self.file_path.stem + "_written")
        file_path = Path(file_path)

        nested_energies, nested_values = self._format_nested_energies_and_values()
        ww_parser.write(
            file_path,
            self._build_header(),
            self.geometry.b2_vectors,
            nested_energies,
            nested_values,
        )

    def _build_header(self) -> WWHeader:
        header_args = {
            "if_": 1,
            "iv": 1,
            "ni": len(self.particles),
            "nr": 10
            if self.geometry.coordinate_type == CoordinateType.CARTESIAN
            else 16,
            "probid": "Created by f4enix_ww",
            "ne": [len(energies) for energies in self.energies.values()],
            "nfx": self.geometry.i_ints,
            "nfy": self.geometry.j_ints,
            "nfz": self.geometry.k_ints,
            "origin": self.geometry.origin,
            "ncx": self.geometry.i_coarse_ints,
            "ncy": self.geometry.j_coarse_ints,
            "ncz": self.geometry.k_coarse_ints,
        }
        if self.geometry.coordinate_type == CoordinateType.CARTESIAN:
            header = WWHeader(**header_args)
        else:  # CoordinateType.CYLINDRICAL
            header_args["director_1"] = self.geometry.director_1
            header_args["director_2"] = self.geometry.director_2
            header = WWHeaderCyl(**header_args)

        return header

    def export_as_vtk(self, file_path: Optional[Pathlike] = None) -> None:
        if file_path is None:
            file_path = self.file_path.parent / (self.file_path.stem)
        file_path = Path(file_path)

        self.geometry.export_as_vtk(file_path)

    def plot(self) -> None:
        self.geometry.plot()
