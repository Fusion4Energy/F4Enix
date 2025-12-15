import re
from pathlib import Path
from importlib.resources import files, as_file
from pypact.input.inputdata import InputData
from pypact.input.serialization import from_file
import numpy as np
import json

from f4enix.constants import TIME_UNITS, TIME_UNITS_CONVERSION, PathLike
from f4enix.input.libmanager import LibManager
from f4enix.constants import SCIENTIFIC_PAT, PAT_DIGIT
from f4enix import resources

LM = LibManager()
RES = files(resources)

# D1STIME PATTERNS
PAT_IRRADIATION = re.compile(r"^\s*irradiation", flags=re.IGNORECASE)

METASTABLE_TAG = "900"
IRS_TAG = "999"


class Pulse:
    def __init__(self, time: float, intensity: float, unit=TIME_UNITS.SECOND) -> None:
        """An object representing a single pulse

        Parameters
        ----------
        time : float
            time duration of the pulse
        intensity : float
            intensity (flux) of the pulse
        unit : TIME_UNITS, optional
            preferred time unit for the pulse, by default TIME_UNITS.SECOND.
            time will always be stored in seconds internally.
        """
        self.time = time * TIME_UNITS_CONVERSION[unit]  # Always store in seconds
        self.intensity = intensity
        self.unit = unit  # preferred unit for display

    def get_time(self, unit: TIME_UNITS) -> float:
        """Get the irradiation time in a specific time unit

        Parameters
        ----------
        unit : TIME_UNITS
            desired time unit

        Returns
        -------
        float
            irradiation time in the desired unit
        """
        return self.time / TIME_UNITS_CONVERSION[unit]

    def __repr__(self) -> str:
        time = self.get_time(self.unit)
        return f"Pulse(time={time} {self.unit}, intensity={self.intensity})"

    def __str__(self) -> str:
        return self.__repr__()


class IrradiationScenario:
    def __init__(
        self,
        pulses: list[Pulse],
        name: str | None = None,
        cooling_times: list[Pulse] | None = None,
    ) -> None:
        """Object representing an irradiation scenario which is characterized by
        a sequence of pulses and cooling times.

        Parameters
        ----------
        pulses : list[Pulse]
            list of irradiation pulses
        name : str | None, optional
            irradiation scenario name, by default None
        cooling_times : list[Pulse] | None, optional
            list of cooling time pulses, by default None. If None, a default cooling
            time of 0s is set.

        Attributes
        ----------
        pulses : list[Pulse]
            list of irradiation pulses
        name : str | None, optional
            irradiation scenario name, by default None
        cooling_times : list[Pulse]
            list of cooling time pulses
        cooling_labels : list[str]
            list of cooling time labels
        """
        self.name = name
        self.pulses = pulses
        if cooling_times is not None:
            self._cooling_times = cooling_times
            self._cooling_labels = [
                f"{pulse.get_time(TIME_UNITS.SECOND)}s" for pulse in cooling_times
            ]
        else:
            self._cooling_times = [
                Pulse(time=0.0, intensity=0.0, unit=TIME_UNITS.SECOND)
            ]
            self._cooling_labels = ["0s"]

    @property
    def cooling_times(self) -> list[Pulse]:
        """Get the cooling times as a list of Pulse objects."""
        return self._cooling_times

    @property
    def cooling_labels(self) -> list[str]:
        """Get the cooling time labels as a list of strings."""
        return self._cooling_labels

    def set_cooling_times(
        self, cooling_times: list[tuple[float, TIME_UNITS]], absolute: bool = True
    ) -> None:
        """set a number of cooling times. Time can be expressed relatively to
        the previous time or absolute after shutdown.

        Parameters
        ----------
        cooling_times : list[tuple[float, TIME_UNITS]]
            list of cooling time durations and their units.
        absolute : bool, optional
            if True, cooling times are absolute after shutdown; if False, they are
            relative to the previous time, by default True
        """
        self._cooling_times = []
        self._cooling_labels = []

        last_cumulative_time = 0.0
        for time_val, time_unit in cooling_times:
            time_sec = time_val * TIME_UNITS_CONVERSION[time_unit]
            # convert the absolute in relative
            if absolute:
                self._cooling_labels.append(f"{time_val}{time_unit.value}")
                time = time_sec - last_cumulative_time
                last_cumulative_time = time_sec
            # keep relatives
            else:
                last_cumulative_time = last_cumulative_time + time_sec
                self._cooling_labels.append(f"{last_cumulative_time}s")
                time = time_sec
            self._cooling_times.append(
                Pulse(time=time, intensity=0.0, unit=TIME_UNITS.SECOND)
            )

    @classmethod
    def from_legacy_d1stime(cls, path_to_file: PathLike) -> "IrradiationScenario":
        """Create an Irradiation Scenario object from a legacy d1stime input format.

        Parameters
        ----------
        path_to_file : PathLike
            path to the input file (d1stime legacy format).

        Returns
        -------
        IrradiationScenario
            The irradiation scenario object.
        """
        with open(path_to_file, "r") as f:
            lines = f.readlines()
        name = Path(path_to_file).stem
        flag_irradiation = False
        pulses = []
        for line in lines:
            irr_line = line
            if PAT_IRRADIATION.match(line):
                flag_irradiation = True
                irr_line = line.split(":", 1)[1].strip()
            elif flag_irradiation and line.strip() == "":
                flag_irradiation = False

            # read irradiation lines only when inside the proper block
            if flag_irradiation:
                pulses.extend(_process_irr_line(irr_line))

        return cls(pulses=pulses, name=name)

    @classmethod
    def from_fispact(cls, path_to_file: PathLike) -> "IrradiationScenario":
        """Create an irradiation scenario from a fispact II input file.

        Parameters
        ----------
        path_to_file : PathLike
            path to the fispact II input file.

        Returns
        -------
        IrradiationScenario
            The irradiation scenario object.
        """
        fisp_inp = InputData()
        from_file(fisp_inp, path_to_file)
        name = Path(path_to_file).stem
        pulses = []
        for time, flux in fisp_inp._irradschedule:
            # time is already converted into seconds by pypact
            pulses.append(Pulse(time=time, intensity=flux, unit=TIME_UNITS.SECOND))

        cooling_times = []
        for cool_time in fisp_inp._coolingschedule:  # already relative in fispact
            cooling_times.append(
                Pulse(time=cool_time, intensity=0.0, unit=TIME_UNITS.SECOND)
            )

        return cls(pulses=pulses, cooling_times=cooling_times, name=name)


class Nuclide:
    def __init__(
        self,
        zaid: int | str,
        metastable: bool = False,
        IRS_active: bool = False,
        lib: str | None = None,
    ) -> None:
        """A general nuclide. Supports metastable, libraries and IRS flags.

        Parameters
        ----------
        zaid : int | str
            ZAID number of the nuclide (e.g. 3003 for Li-3).
        metastable : bool, optional
            true if the nuclide is metastable, by default False
        IRS_active : bool, optional
            true if the IRS flag is active, by default False
        lib : str | None, optional
            library identifier, by default None
        """
        self._zaid = int(zaid)
        self.metastable = metastable
        self.IRS_active = IRS_active
        self.lib = lib

    @classmethod
    def from_formula(cls, zaid_str: str) -> "Nuclide":
        """Create a Nuclide object starting from a formula that looks like "Li3".
        Metastables are supported with an 'm' and IRS flags with the 'irs' tag
        at the beginning of the string. MCNP libraries are also supported using '.XXc'.
        A complete example looks like "irsLi3m.99c".

        Parameters
        ----------
        zaid_str : str
            string describing the nuclide as documented.

        Returns
        -------
        Nuclide
            the created Nuclide object.
        """
        # may or may not have a lib in .XXc format
        pieces = zaid_str.split(".")
        zaid_str = pieces[0]

        # check for metastable
        if zaid_str.endswith("m"):
            zaid_str = zaid_str[:-1]
            metastable = True
        else:
            metastable = False

        # check for IRS
        if zaid_str.startswith("irs"):
            zaid_str = zaid_str[3:]
            IRS_active = True
        else:
            IRS_active = False

        zaid = LM.get_zaidnum(zaid_str)
        if len(pieces) > 1:
            lib = pieces[1]
        else:
            lib = None

        return cls(int(zaid), metastable=metastable, IRS_active=IRS_active, lib=lib)

    @classmethod
    def from_int_string(cls, zaid_int_string: str) -> "Nuclide":
        """Create a nuclide object starting from an integer string that looks like "3003".
        Metastables are supported with a '900' suffix and IRS flags with a '999' prefix.
        MCNP libraries are also supported using '.XXc'.
        A complete example looks like "9993003900.99c".

        Parameters
        ----------
        zaid_int_string : str
            string describing the nuclide as documented.

        Returns
        -------
        Nuclide
            the created Nuclide object.
        """
        # may or may not have a lib in .XXc format
        pieces = zaid_int_string.split(".")
        zaid_str = pieces[0]

        # check for special cases

        if len(zaid_str) > 5:
            # starts with IRS?
            if zaid_str.startswith(IRS_TAG):
                zaid_str = zaid_str[3:]
                IRS_active = True
            else:
                IRS_active = False
            # ends with metastable?
            if zaid_str.endswith(METASTABLE_TAG):
                zaid_str = zaid_str[:-3]
                metastable = True
            else:
                metastable = False
        else:
            metastable = False
            IRS_active = False

        zaid = int(zaid_str)

        if len(pieces) > 1:
            lib = pieces[1]
        else:
            lib = None

        return cls(zaid, metastable=metastable, IRS_active=IRS_active, lib=lib)

    def write_to_formula(self) -> str:
        """Return the formula string representation of the nuclide. E.g. "irsLi3m.99c"."""
        result = ""
        if self.IRS_active:
            result += "irs"
        _, formula = LM.get_zaidname(str(self._zaid))
        result += formula.replace("-", "")
        if self.metastable:
            result += "m"
        if self.lib:
            result += f".{self.lib}"
        return result

    def write_to_int_string(self) -> str:
        """Return the integer string representation of the nuclide. E.g. "9993003900.99c"."""
        zaid_str = str(self._zaid)
        if self.IRS_active:
            zaid_str = "999" + zaid_str
        if self.metastable:
            zaid_str += "900"
        result = zaid_str
        if self.lib:
            result += f".{self.lib}"
        return result

    @property
    def zaid(self) -> int:
        return self._zaid

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Nuclide):
            return False
        return (
            self._zaid == value._zaid
            and self.metastable == value.metastable
            and self.IRS_active == value.IRS_active
            and self.lib == value.lib
        )

    def __repr__(self) -> str:
        return self.write_to_formula()

    def __str__(self) -> str:
        return self.__repr__()


class TCF_Computer:
    def __init__(self) -> None:
        """Auxiliary class to be used to compute time correction factors.

        Attributes
        ----------
        half_lives : dict
            Dictionary containing half-life data for nuclides.
        """
        with as_file(RES.joinpath("half_lives_decay2020.json")) as decay_file:
            with open(decay_file, "r") as f:
                self.half_lives = dict(json.load(f))

    def get_half_life(self, nuclide: Nuclide) -> float | str:
        """Get the half-life for a given nuclide.

        Parameters
        ----------
        nuclide : Nuclide
            The nuclide for which to get the half-life.

        Returns
        -------
        float | str
            The half-life in seconds or "STABLE" if the nuclide is stable.
        """
        nuclide_str = nuclide.write_to_formula()
        nuclide_str = nuclide_str.removeprefix("irs")  # remove IRS if present
        half_life_sec = self.half_lives[nuclide_str]
        return half_life_sec

    def get_lambda(self, nuclide: Nuclide) -> float:
        """Get the decay constant (lambda) for a given nuclide.

        Parameters
        ----------
        nuclide : Nuclide
            The nuclide for which to get the decay constant.

        Returns
        -------
        float
            The decay constant in 1/seconds.
        """
        half_life_sec = self.get_half_life(nuclide)

        if half_life_sec == "STABLE":
            return 0.0

        lambda_value = 0.69314718056 / float(half_life_sec)  # ln(2) / half-life
        return lambda_value

    def compute_correction_factors(
        self,
        scenario: IrradiationScenario,
        nuclides: list[Nuclide],
        norm: float = 1,
    ) -> np.ndarray:
        """Compute time correction factors for D1S methodology.
        N[0] = 0
        N[m] = N[m-1]*exp(-lambda*dt) + I/norm * (1-exp(-lambda*dt))

        Parameters
        ----------
        scenario : IrradiationScenario
            The irradiation scenario. It must include also the cooling time
            equivalent pulses.
        nuclides : list[Nuclide]
            list of nuclides for which to compute the correction factors.
        norm : float, optional
            norm to be used to scale the neutron flux intensity, by default 1

        Returns
        -------
        np.ndarray
            Array of correction factors for each nuclide and each cooling time.
        """
        # get lambda vector
        lambda_vector = np.array(
            [self.get_lambda(nuclide) for nuclide in nuclides]
        )  # .reshape((-1, 1))  # column vector
        N = np.zeros_like(lambda_vector, dtype=float)

        # Compute factor up to the end of irradiation
        for pulse in scenario.pulses:
            N = _compute_factor(N, pulse, lambda_vector, norm)

        # Compute factor during cooling times
        factors = []
        for pulse in scenario.cooling_times:
            N = _compute_factor(N, pulse, lambda_vector, norm)
            factors.append(N.copy())
        return np.array(factors)


def _compute_factor(N, pulse: Pulse, lambda_vector, norm) -> np.ndarray:
    I = pulse.intensity / norm
    t = pulse.time
    exp_term = np.exp(-lambda_vector * t)
    N = N * exp_term + I * (1 - exp_term)
    return N


def _process_irr_line(line: str) -> list[Pulse]:
    pulses = []
    # remove spaces and tabs
    inline = line.replace(" ", "").replace("\t", "").replace("\n", "").replace("\r", "")
    # check for parentheses
    if inline.startswith("("):
        multiplier = inline.split(")*")[-1]
        inline = inline[1 : -(len(multiplier) + 2)]  # remove parentheses and multiplier
        pulses.extend(_process_pulses(inline) * int(multiplier))
    else:
        pulses.extend(_process_pulses(inline))

    return pulses


def _process_pulses(pulse_str: str) -> list[Pulse]:
    pieces = pulse_str.split("/")
    pulses = []
    for i in range(len(pieces) // 2):
        flux = float(pieces[2 * i])
        time_str = pieces[2 * i + 1]
        # In time, separate number and unit
        val = SCIENTIFIC_PAT.search(time_str)
        if not val:
            val = PAT_DIGIT.search(time_str)
            if not val:
                raise ValueError(f"Cannot parse time value from string '{time_str}'")

        val = val.group()
        unit = TIME_UNITS(time_str.replace(val, "").strip().lower())
        pulses.append(Pulse(time=float(val), intensity=flux, unit=unit))

    return pulses
