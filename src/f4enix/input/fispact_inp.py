from pypact import InputData, to_file
from pathlib import Path
from f4enix.input.materials import Material
from f4enix.core.irradiation import IrradiationScenario
from f4enix.input.libmanager import LibManager


class FispactInp:
    def __init__(self, pypact_inp: InputData | None = None):
        if pypact_inp is None:
            self.inp = InputData()
        else:
            self.inp = pypact_inp

    def add_material(
        self,
        material: Material,
        style: str = "mass",
        density: float | None = None,
        mass: float | None = None,
        volume: float | None = None,
    ):
        """Add a material to the FISPACT input

        Parameters
        ----------
        material : Material
            f4enix material definition
        style : str, optional
            The style in which the material is specified, either "mass" or "fuel",
            by default "mass". This will determine if the MASS or FUEL keywords will
            be used
        density : float | None, optional
            The density of the material in g/cm^3, by default None
        mass : float | None, optional
            The mass of the material in kg, by default None
        volume : float | None, optional
            The volume of the material in m^3, by default None

        Raises
        ------
        ValueError
            depending on the style chosen, some arguments must be provided.
            - For "mass" style, mass must be provided or both volume and density must
              be provided.
            - For "fuel" style, density and mass (or volume) or both mass and volume
              must be provided.

        """
        # Add logic to incorporate the material into the input data
        lm = LibManager()
        if style == "mass":
            new_mat = material.switch_fraction("mass", lm, inplace=False)
            # fispact natural abundances are used! be sure these are consistent
            # with the material definition
            for element in new_mat.elements:
                self.inp.addElement(element.name, percentage=-element.get_fraction())
            if mass:
                m = mass
            elif volume and density:
                m = volume * density / 1000  # convert g/cm^3 to kg/m^3
            else:
                raise ValueError(
                    "Either mass or both volume and density must be specified for mass style materials."
                )
            self.inp.setMass(m)

        elif style == "fuel":
            new_mat = material.switch_fraction("atom", lm, inplace=False)
            if density:
                rho = density
            elif volume and mass:
                rho = mass / volume * 1000  #  g/cm^3
            else:
                raise ValueError(
                    "Either density or both mass and volume must be specified for fuel style materials."
                )

            if volume:
                vol = volume
            elif mass and density:
                vol = mass / rho
            else:
                raise ValueError(
                    "Either volume or both mass and density must be specified for fuel style materials."
                )

            tad = new_mat.get_tad(rho, lm)
            for zaid in new_mat.zaids:
                n_atoms = zaid.fraction * tad * vol * 1e6 * 1e24
                self.inp.addIsotope(zaid.name, n_atoms)
        else:
            raise ValueError(
                f"Unknown style: {style}, only 'mass' and 'fuel' are supported."
            )

    def add_irradiation_scenario(
        self, irr_scenario: IrradiationScenario, norm: float = 1
    ):
        """Add an irradiation scenario to the fispact input

        Parameters
        ----------
        irr_scenario : IrradiationScenario
            The irradiation scenario to be added to the fispact input.
        norm : float, optional
            Normalization factor for the irradiation intensities, by default 1.
        """
        for pulse in irr_scenario.pulses:
            # Each intensity will need to be adjusted to the flux of the spectra
            intensity = pulse.intensity * norm
            self.inp.addIrradiation(pulse.time, intensity)
        for time in irr_scenario.cooling_times:
            self.inp.addCooling(time.time)

    def save(self, filename: str | Path):
        """Save the fispact input to a file

        Parameters
        ----------
        filename : str | Path
            The name or path of the file to save the fispact input to.
        """
        to_file(self.inp, filename)
