"""Deals with irradiation scenarios and time correction factors for D1S methodology."""

import json
import re
from copy import deepcopy
from importlib.resources import as_file, files
from pathlib import Path

import numpy as np
import pandas as pd
from pypact.input.inputdata import InputData
from pypact.input.serialization import from_file

from f4enix import resources
from f4enix.core.constants import (
    PAT_DIGIT,
    SCIENTIFIC_PAT,
    TIME_UNITS,
    TIME_UNITS_CONVERSION,
    PathLike,
)
from f4enix.input.libmanager import LibManager

LM = LibManager()
RES = files(resources)

# D1STIME PATTERNS
PAT_IRRADIATION = re.compile(r"^\s*irradiation", flags=re.IGNORECASE)

# ASCII SCENARIO FORMAT PATTERNS
_PAT_ASCII_BLOCK_END = re.compile(r"^\)\s*\*\s*(\d+)$")
_PAT_LABEL_TIME = re.compile(r"^(\d+(?:\.\d+)?)(s|min|h|d|w|m|y)$")

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

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Pulse):
            return False
        return (
            self.time == value.time
            and self.intensity == value.intensity
            and self.unit == value.unit
        )


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
        # make sure all pulses are individual copies
        self.pulses = [deepcopy(pulse) for pulse in pulses]
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

    def get_collapsed_table(self) -> pd.DataFrame:
        """Compile a dataframe with the irradiation sceario data. Pulses or sequences
        of pulses that are repeated are collapsed into single entries with a multiplier."""

        def scan_sequence(pulses: list[Pulse]) -> tuple[int, int]:
            "How many repetitions containing the pulse at index zero are there?"
            multiplier = 1
            for len_sequence in range(1, (len(pulses) // 2 + 1)):
                seq = pulses[:len_sequence]
                for check_idx in range(len_sequence, len(pulses), len_sequence):
                    if pulses[check_idx : check_idx + len_sequence] == seq:
                        multiplier += 1
                    else:
                        break
                if multiplier > 1:
                    return multiplier, len_sequence
            return 1, 1

        def get_record(pulse: Pulse, multiplier: int) -> dict:
            record = {
                "Time": f"{pulse.get_time(pulse.unit)} {pulse.unit.value}",
                "Intensity": pulse.intensity,
                "Repetition": multiplier,
            }
            return record

        remaining_pulses = self.pulses
        records = []
        while len(remaining_pulses) > 0:
            multiplier, len_sequence = scan_sequence(remaining_pulses)
            if multiplier > 1:
                sequence = remaining_pulses[:len_sequence]
                remaining_pulses = remaining_pulses[len_sequence * multiplier :]
                for pulse in sequence:
                    records.append(get_record(pulse, multiplier))
            else:
                records.append(get_record(remaining_pulses[0], multiplier=1))
                remaining_pulses = remaining_pulses[1:]

        return pd.DataFrame(records).set_index(["Repetition", "Time"])

    def to_ascii(self, path: PathLike) -> None:
        """Serialize the irradiation scenario to a human-readable ASCII file.

        The file uses three sections separated by keyword headers:

        - ``NAME`` (optional): the scenario name.
        - ``IRRADIATION``: one pulse per line as ``time unit intensity``.
          Repeated sequences are written compactly as ``( ... ) * N`` blocks.
        - ``COOLING``: one line per cooling time as
          ``cumulative_time unit [label]``.  Times are absolute (measured
          from the end of irradiation).

        Lines starting with ``#`` are treated as comments and ignored on
        reading.

        Parameters
        ----------
        path : PathLike
            Destination file path (conventionally ``*.irr``).
        """
        lines = [
            "# F4Enix Irradiation Scenario",
            "# Cooling times are cumulative (absolute) from end of irradiation",
            "",
        ]

        if self.name is not None:
            lines += [f"NAME  {self.name}", ""]

        lines += ["IRRADIATION", "# time  unit  intensity"]

        remaining = list(self.pulses)
        while remaining:
            mult, len_seq = _scan_sequence(remaining)
            if mult > 1:
                seq = remaining[:len_seq]
                remaining = remaining[len_seq * mult :]
                lines.append("  (")
                for pulse in seq:
                    t = pulse.get_time(pulse.unit)
                    lines.append(f"    {t:g}  {pulse.unit.value}  {pulse.intensity:g}")
                lines.append(f"  ) * {mult}")
            else:
                pulse = remaining[0]
                remaining = remaining[1:]
                t = pulse.get_time(pulse.unit)
                lines.append(f"  {t:g}  {pulse.unit.value}  {pulse.intensity:g}")

        lines += ["", "COOLING", "# cumulative_time  unit  [label]"]

        cumulative_sec = 0.0
        for pulse, label in zip(self._cooling_times, self._cooling_labels):
            cumulative_sec += pulse.time
            time_str, unit_str = _format_cooling_line(cumulative_sec, label)
            lines.append(f"  {time_str}  {unit_str}  {label}")

        lines.append("")

        with open(path, "w") as f:
            f.write("\n".join(lines))

    @classmethod
    def from_ascii(cls, path: PathLike) -> "IrradiationScenario":
        """Load an irradiation scenario from a human-readable ASCII file.

        See :meth:`to_ascii` for the file format specification.

        Parameters
        ----------
        path : PathLike
            Path to the ``.irr`` file to read.

        Returns
        -------
        IrradiationScenario
            The loaded irradiation scenario.
        """
        with open(path, "r") as f:
            raw_lines = f.readlines()

        name = None
        pulses: list[Pulse] = []
        cooling_abs: list[tuple[float, TIME_UNITS]] = []
        cooling_labels: list[str] = []

        section = None
        in_block = False
        block_pulses: list[Pulse] = []

        for raw_line in raw_lines:
            line = raw_line.split("#")[0].strip()
            if not line:
                continue

            upper = line.upper()

            if upper.startswith("NAME"):
                tokens = line.split(None, 1)
                if len(tokens) > 1:
                    name = tokens[1].strip()
                continue

            if upper == "IRRADIATION":
                section = "IRRADIATION"
                continue

            if upper == "COOLING":
                section = "COOLING"
                continue

            if section == "IRRADIATION":
                if line == "(":
                    in_block = True
                    block_pulses = []
                    continue

                m = _PAT_ASCII_BLOCK_END.match(line)
                if m and in_block:
                    n = int(m.group(1))
                    pulses.extend(block_pulses * n)
                    in_block = False
                    block_pulses = []
                    continue

                tokens = line.split()
                if len(tokens) >= 3:
                    pulse = _parse_ascii_pulse_line(tokens)
                    if in_block:
                        block_pulses.append(pulse)
                    else:
                        pulses.append(pulse)

            elif section == "COOLING":
                tokens = line.split()
                if len(tokens) >= 2:
                    time_val = float(tokens[0])
                    unit = TIME_UNITS(tokens[1])
                    label = tokens[2] if len(tokens) >= 3 else f"{tokens[0]}{tokens[1]}"
                    cooling_abs.append((time_val, unit))
                    cooling_labels.append(label)

        # Convert absolute cumulative cooling times to relative deltas
        cooling_pulses: list[Pulse] = []
        prev_sec = 0.0
        for time_val, unit in cooling_abs:
            abs_sec = time_val * TIME_UNITS_CONVERSION[unit]
            rel_sec = abs_sec - prev_sec
            prev_sec = abs_sec
            cooling_pulses.append(
                Pulse(time=rel_sec, intensity=0.0, unit=TIME_UNITS.SECOND)
            )

        scenario = cls(
            pulses=pulses,
            name=name,
            cooling_times=cooling_pulses if cooling_pulses else None,
        )

        if cooling_labels:
            scenario._cooling_labels = cooling_labels

        return scenario


class Nuclide:
    def __init__(
        self,
        zaid: int | str | None = None,
        element: str | None = None,
        isotope: int | str | None = None,
        metastable: bool = False,
        IRS_active: bool = False,
        lib: str | None = None,
    ) -> None:
        """A general nuclide. Supports metastable, libraries and IRS flags.
        As a minimum, either the zaid number or element and isotope must be provided.

        Parameters
        ----------
        zaid : int | str | None
            ZAID number of the nuclide (e.g. 3003 for Li-3).
        element : str | None, optional
            Element symbol of the nuclide (e.g. "Li" for Lithium), by default None
        isotope : int | str | None, optional
            Isotope number of the nuclide (e.g. 3 for Li-3), by default None
        metastable : bool, optional
            true if the nuclide is metastable, by default False
        IRS_active : bool, optional
            true if the IRS flag is active, by default False
        lib : str | None, optional
            library identifier, by default None
        """
        if zaid:
            self._zaid = int(zaid)
            _, formula = LM.get_zaidname(str(self._zaid))
            self._element = formula.split("-")[0].strip()
            self._isotope = int(formula.split("-")[1].strip())
        elif element and isotope:
            self._element = element
            self._isotope = int(isotope)
            self._zaid = LM.get_zaidnum(f"{element}{isotope}")

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
        if self.IRS_active:
            irs = "irs"
        else:
            irs = ""

        if self.metastable:
            metastable = "m"
        else:
            metastable = ""
        if self.lib:
            lib = f".{self.lib}"
        else:
            lib = ""

        return f"{irs}{self.element}{self.isotope}{metastable}{lib}"

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

    @property
    def element(self) -> str:
        return self._element

    @property
    def isotope(self) -> int:
        return self._isotope

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


def _scan_sequence(pulses: list[Pulse]) -> tuple[int, int]:
    """Return (multiplier, len_sequence) for the repeating block starting at pulses[0].

    Returns (1, 1) when no repetition is detected.
    """
    multiplier = 1
    for len_sequence in range(1, (len(pulses) // 2 + 1)):
        seq = pulses[:len_sequence]
        for check_idx in range(len_sequence, len(pulses), len_sequence):
            if pulses[check_idx : check_idx + len_sequence] == seq:
                multiplier += 1
            else:
                break
        if multiplier > 1:
            return multiplier, len_sequence
    return 1, 1


def _format_cooling_line(abs_sec: float, label: str) -> tuple[str, str]:
    """Return (time_str, unit_str) for writing a cooling line.

    If *label* encodes a time value and unit (e.g. ``"1y"``, ``"10d"``), those
    are used directly so the file stays human-readable.  Otherwise the
    absolute seconds value is written with unit ``"s"``.
    """
    m = _PAT_LABEL_TIME.match(label)
    if m:
        return m.group(1), m.group(2)
    return f"{abs_sec:g}", "s"


def _parse_ascii_pulse_line(tokens: list[str]) -> Pulse:
    """Parse a ``time unit intensity`` token list into a :class:`Pulse`."""
    return Pulse(
        time=float(tokens[0]),
        intensity=float(tokens[2]),
        unit=TIME_UNITS(tokens[1]),
    )
