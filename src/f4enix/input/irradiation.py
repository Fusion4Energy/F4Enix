from typing import Any
from f4enix.constants import PathLike, TIME_UNITS, TIME_UNITS_CONVERSION
from f4enix.input.libmanager import LibManager

LM = LibManager()


class Pulse:
    def __init__(self, time: float, intensity: float, unit=TIME_UNITS.SECOND) -> None:
        pass


class IrradiationScenario:
    def __init__(self, pulses: list[Pulse], name: str | None = None):
        self.name = name
        self.pulses = pulses

    @classmethod
    def from_legacy_d1stime(cls, path_to_file: PathLike) -> "IrradiationScenario":
        pass

    @classmethod
    def from_fispact(cls, path_to_file: PathLike) -> "IrradiationScenario":
        pass


class Nuclide:
    def __init__(
        self,
        zaid: int,
        metastable: bool = False,
        IRS_active: bool = False,
        lib: str | None = None,
    ) -> None:
        """A general nuclide. Supports metastable, libraries and IRS flags.

        Parameters
        ----------
        zaid : int
            ZAID number of the nuclide (e.g. 3003 for Li-3).
        metastable : bool, optional
            true if the nuclide is metastable, by default False
        IRS_active : bool, optional
            true if the IRS flag is active, by default False
        lib : str | None, optional
            library identifier, by default None
        """
        self._zaid = zaid
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
        metastable_tag = "900"
        irs_tag = "999"

        if len(zaid_str) > 5:
            # starts with IRS?
            if zaid_str.startswith(irs_tag):
                zaid_str = zaid_str[3:]
                IRS_active = True
            else:
                IRS_active = False
            # ends with metastable?
            if zaid_str.endswith(metastable_tag):
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


class TimeCorrectionFactorComputer:
    def __init__(self, decay_lib: PathLike) -> None:
        pass

    def compute_correction_factor(
        self,
        scenario: IrradiationScenario,
        nuclides: list[Nuclide],
        scale_factors: float = 1,
    ) -> float:
        pass

    def get_lambda(nuclide: Nuclide) -> float:
        pass
