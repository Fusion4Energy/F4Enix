"""Contains material compositions"""

from __future__ import annotations

from typing import List, Tuple


class MaterialComposition:
    def __init__(
        self,
        name: str,
        composition: List[Tuple[str, str | float | int]],
        fraction: bool = False,
    ) -> None:
        """Suitable object to simplify the creation of materials

        Parameters
        ----------
        name : str
            name of the material
        composition : List[Tuple[str, float | str]]
            list of elements and their mass percentage. If the mass is provided
            as string it is converted to float
        fraction : bool
            if True percentages are not to be interpreted as such but instead
            are fractions. By default is False, meaning that percentages are
            expected.

        Attributes
        ----------
        perc : list[float]
            percentages of the different elements composing the material
        elem : list[str]
            list of elements symbols (i.e. Ag)
        name : str
            name of the material
        """
        elements = []
        percentages = []
        for elem, perc in composition:
            perc = float(perc)
            elements.append(elem)
            if fraction:
                percentages.append(perc * 100)
            else:
                percentages.append(perc)

        self.perc = percentages
        self.elem = elements
        self.name = name

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        return f"MaterialComposition({self.name})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MaterialComposition):
            return False
        if (
            self.name == other.name
            and self.perc == other.perc
            and self.elem == other.elem
        ):
            return True
        else:
            return False


# Material library
SILVER = MaterialComposition("Pure Silver", [("Ag", 100)])
TUNGSTEN = MaterialComposition(
    "W with impurities",
    [
        ("C", "0.003"),
        ("N", "0.001"),
        ("O", "0.003"),
        ("MG", "0.0005"),
        ("AL", "0.0015"),
        ("SI", "0.002"),
        ("P", "0.005"),
        ("S", "0.0005"),
        ("TI", "0.001"),
        ("CR", "0.001"),
        ("MN", "0.0005"),
        ("FE", "0.003"),
        ("ZN", "0.0005"),
        ("ZR", "0.001"),
        ("MO", "0.01"),
        ("AG", "0.0005"),
        ("CD", "0.001"),
        ("W", "99.9565"),
        ("PB", "0.001"),
        ("NI", "0.003"),
        ("CU", "0.001"),
        ("NB", "0.001"),
        ("CO", "0.001"),
        ("TA", "0.001"),
        ("CA", "0.001"),
        ("H", "0.0005"),
        ("K", "0.001"),
        ("NA", "0.001"),
        ("AS", "0.0005"),
        ("BA", "0.001"),
    ],
)
LEAD = MaterialComposition("Pure Lead", [("Pb", 100)])
AL6061 = MaterialComposition(
    "Al-6061",
    [
        ("MG", "1.00"),
        ("AL", "97.2"),
        ("SI", "0.6"),
        ("TI", "0.088"),
        ("CR", "0.195"),
        ("MN", "0.088"),
        ("FE", "0.409"),
        ("ZN", "0.146"),
        ("CU", "0.275"),
    ],
)
ALBRZ = MaterialComposition(
    "Al-Brz",
    [
        ("AL", "9.75"),
        ("SI", "0.2"),
        ("MN", "1.0"),
        ("FE", "4.0"),
        ("ZN", "0.02"),
        ("CD", "0.02"),
        ("SN", "0.1"),
        ("PB", "0.02"),
        ("NI", "5.0"),
        ("CU", "79.73"),
        ("NB", "0.1"),
        ("CO", "0.05"),
        ("TA", "0.01"),
    ],
)
BHC = MaterialComposition(
    "Borated Heavy Concrete",
    [
        ("B", "0.41823"),
        ("C", "0.44503"),
        ("N", "0.04945"),
        ("O", "31.89362"),
        ("MG", "0.54392"),
        ("AL", "0.49447"),
        ("SI", "2.20536"),
        ("P", "0.20273"),
        ("S", "0.08406"),
        ("TI", "0.17307"),
        ("CR", "0.00119"),
        ("MN", "0.02576"),
        ("FE", "55.81631"),
        ("ZN", "0.00099"),
        ("ZR", "0.00059"),
        ("CO", "0.00519"),
        ("TA", "0.0002"),
        ("CA", "6.58146"),
        ("NA", "0.44503"),
        ("K", "0.11867"),
        ("CL", "0.01978"),
        ("H", "0.30163"),
        ("BA", "0.00326"),
        ("CE", "0.00376"),
        ("CS", "0.00015"),
        ("DY", "0.0001"),
        ("EU", "0.0001"),
        ("GA", "0.00747"),
        ("HF", "0.0001"),
        ("LA", "0.00267"),
        ("LU", "0.00138"),
        ("ND", "0.00099"),
        ("PR", "0.00702"),
        ("RB", "0.00074"),
        ("SB", "0.0003"),
        ("SC", "0.0001"),
        ("SE", "0.00025"),
        ("SM", "0.0004"),
        ("SR", "0.03704"),
        ("TB", "0.00035"),
        ("TH", "0.00049"),
        ("U", "0.00059"),
        ("V", "0.10532"),
        ("Y", "0.00059"),
        ("YB", "0.0001"),
    ],
)
CONCRETE = MaterialComposition(
    "Concrete",
    [
        ("B", "7.4800E-03"),
        ("C", "5.8868E+00"),
        ("N", "4.9890E-2"),
        ("O", "4.9090E+01"),
        ("MG", "5.1634E-1"),
        ("AL", "1.5465E+00"),
        ("SI", "1.5588E+01"),
        ("P", "2.2450E-2"),
        ("S", "3.8913E-1"),
        ("TI", "7.4830E-2"),
        ("CR", "1.7200E-03"),
        ("MN", "4.3700E-03"),
        ("FE", "9.3291E-1"),
        ("ZN", "3.3400E-03"),
        ("ZR", "1.2700E-03"),
        ("CU", "8.7000E-04"),
        ("CO", "2.5000E-04"),
        ("TA", "2.5000E-04"),
        ("H", "3.6418E-1"),
        ("LI", "9.0000E-04"),
        ("NA", "3.0182E-1"),
        ("CL", "1.2470E-2"),
        ("K", "4.2156E-1"),
        ("CA", "2.4695E+01"),
        ("SC", "2.0000E-04"),
        ("V", "2.1000E-03"),
        ("GA", "2.7700E-03"),
        ("GE", "4.0000E-04"),
        ("SE", "2.0000E-04"),
        ("RB", "5.3400E-03"),
        ("SR", "5.1630E-2"),
        ("Y", "1.1000E-03"),
        ("SB", "2.0000E-04"),
        ("CS", "3.0000E-04"),
        ("BA", "7.2800E-03"),
        ("LA", "1.0200E-03"),
        ("CE", "1.6200E-03"),
        ("PR", "3.8200E-03"),
        ("ND", "1.7000E-03"),
        ("SM", "2.0000E-04"),
        ("EU", "3.0000E-04"),
        ("GD", "4.0000E-04"),
        ("TB", "6.2000E-04"),
        ("DY", "1.4000E-04"),
        ("ER", "4.0000E-04"),
        ("YB", "1.0000E-04"),
        ("LU", "5.0000E-05"),
        ("HF", "2.0000E-04"),
        ("TH", "1.3000E-03"),
        ("U", "4.9900E-03"),
    ],
)
CUCRZR = MaterialComposition(
    "Copper Cromium Zirconium (CuCrZr)",
    [
        ("B", "1.0000E-03"),
        ("O", "3.2000E-2"),
        ("MG", "4.0000E-2"),
        ("AL", "3.0000E-02"),
        ("SI", "4.0000E-2"),
        ("P", "1.4000E-2"),
        ("S", "4.0000E-03"),
        ("CR", "7.5000E-2"),
        ("MN", "2.0000E-03"),
        ("FE", "2.0000E-2"),
        ("ZN", "1.0000E-2"),
        ("ZR", "1.1000E-1"),
        ("SN", "1.0000E-2"),
        ("PB", "1.0000E-2"),
        ("BI", "3.0000E-03"),
        ("NI", "6.0000E-2"),
        ("CU", "9.8870E+01"),
        ("NB", "1.0000E-1"),
        ("CO", "5.0000E-2"),
        ("TA", "1.0000E-2"),
        ("AS", "1.0000E-2"),
        ("SB", "1.1000E-2"),
    ],
)
EUROFER97 = MaterialComposition(
    "Eurofer-97",
    [
        ("B", "2.0000E-03"),
        ("C", "1.2000E-1"),
        ("N", "4.5000E-2"),
        ("O", "1.0000E-2"),
        ("AL", "1.0000E-2"),
        ("SI", "5.0000E-2"),
        ("P", "5.0000E-03"),
        ("S", "5.0000E-03"),
        ("TI", "2.0000E-02"),
        ("CR", "9.0000E+00"),
        ("MN", "6.0000E-1"),
        ("FE", "8.8598E+01"),
        ("ZR", "1.2500E-03"),
        ("MO", "5.0000E-03"),
        ("SN", "1.2500E-03"),
        ("W", "1.1000E+00"),
        ("NI", "1.0000E-2"),
        ("CU", "1.0000E-2"),
        ("NB", "5.0000E-03"),
        ("CO", "1.0000E-2"),
        ("TA", "1.4000E-1"),
        ("V", "2.5000E-1"),
        ("SB", "1.2500E-03"),
        ("AS", "1.2500E-03"),
    ],
)
INCONEL718 = MaterialComposition(
    "Inconel-718",
    [
        ("B", "4.0000E-03"),
        ("C", "5.0000E-2"),
        ("AL", "5.0000E-1"),
        ("SI", "3.5000E-1"),
        ("P", "1.5000E-2"),
        ("S", "1.5000E-2"),
        ("TI", "9.0000E-01"),
        ("CR", "1.9000E+01"),
        ("MN", "3.5000E-01"),
        ("FE", "1.7706E+01"),
        ("MO", "3.0500E+00"),
        ("NI", "5.2500E+01"),
        ("CU", "3.0000E-01"),
        ("NB", "5.1500E+00"),
        ("CO", "1.0000E-01"),
        ("TA", "1.0000E-2"),
    ],
)
NB3SN = MaterialComposition("Nb3Sn", [("SN", "2.9000E+01"), ("NB", "7.1000E+01")])
NBTI = MaterialComposition("NbTi", [("TI", "4.6500E+01"), ("NB", "5.3500E+01")])
SS660 = MaterialComposition(
    "SS660",
    [
        ("B", "1.0000E-2"),
        ("C", "5.5000E-2"),
        ("AL", "3.5000E-01"),
        ("SI", "1.0000E+00"),
        ("P", "2.5000E-02"),
        ("S", "1.5000E-2"),
        ("TI", "2.1000E+00"),
        ("CR", "1.4750E+01"),
        ("MN", "1.5000E+00"),
        ("FE", "5.3395E+01"),
        ("MO", "1.2500E+00"),
        ("NI", "2.5500E+01"),
        ("CO", "5.0000E-2"),
        ("TA", "1.0000E-2"),
        ("V", "5.0000E-01"),
    ],
)
SS316LNIG = MaterialComposition(
    "SS316L(N)-IG",
    [
        ("B", "2.0000E-03"),
        ("C", "3.0000E-2"),
        ("N", "7.0000E-2"),
        ("SI", "5.0000E-01"),
        ("P", "2.5000E-02"),
        ("S", "1.0000E-02"),
        ("TI", "1.0000E-01"),
        ("CR", "1.7500E+01"),
        ("MN", "1.8000E+00"),
        ("FE", "6.4843E+01"),
        ("MO", "2.5000E+00"),
        ("NI", "1.2250E+01"),
        ("CU", "3.0000E-01"),
        ("NB", "1.0000E-02"),
        ("CO", "5.0000E-02"),
        ("TA", "1.0000E-02"),
    ],
)
MICROTHERM = MaterialComposition(
    "Microtherm",
    [
        ("O", 41.42),
        ("Na", 0.1),
        ("Mg", 0.41),
        ("Al", 2.24),
        ("Si", 29.49),
        ("P", 0.01),
        ("S", 0.02),
        ("Cl", 0.01),
        ("Ca", 1.83),
        ("Ti", 23.64),
        ("V", 0.07),
        ("Cr", 0.04),
        ("Fe", 0.2),
        ("Zr", 0.15),
    ],
)
TITANIUM5 = MaterialComposition(
    "Titanium Alloy, Grade 5",
    [
        ("H", 0.000110),
        ("C", 0.000570),
        ("N", 0.000210),
        ("O", 0.001410),
        ("Al", 0.61250),
        ("Ti", 0.893630),
        ("V", 0.040),
        ("Fe", 0.002830),
    ],
    fraction=True,
)
BERYLLIUM = MaterialComposition(
    "Beryllium",
    [
        ("Be", 9.92e-01),
        ("O", 3.69e-03),
        ("Al", 6.00e-04),
        ("C", 9.89e-04),
        ("Fe", 8.00e-04),
        ("Mg", 6.00e-04),
        ("Si", 6.00e-04),
        ("Ni", 1.06e-04),
        ("Mn", 1.50e-05),
        ("Cu", 3.10e-05),
        ("Ti", 6.00e-05),
        ("Co", 6.00e-06),
        ("Pb", 1.97e-05),
        ("Ca", 2.00e-05),
        ("W", 9.99e-05),
        ("Mo", 2.00e-05),
        ("Cr", 3.70e-05),
        ("N", 1.14e-04),
        ("Zr", 4.30e-05),
        ("B", 2.00e-06),
        ("Li", 3.00e-06),
        ("Na", 1.50e-05),
        ("S", 1.00e-05),
        ("F", 5.00e-06),
        ("Cl", 1.00e-05),
        ("P", 4.60e-07),
        ("Nb", 1.20e-07),
        ("Ta", 9.70e-06),
    ],
    fraction=True,
)

AVAILABLE_MATERIALS = [
    SILVER,
    TUNGSTEN,
    LEAD,
    AL6061,
    ALBRZ,
    BHC,
    CONCRETE,
    CUCRZR,
    EUROFER97,
    INCONEL718,
    NB3SN,
    NBTI,
    SS660,
    SS316LNIG,
    MICROTHERM,
    TITANIUM5,
    BERYLLIUM,
]
