# This file has been modified to work with ace files instead of hdf5
from __future__ import annotations

import os
import typing
from itertools import chain
from numbers import Integral, Real
from pathlib import Path
from typing import Callable, Iterable, Sequence

import endf
import matplotlib.axes
import matplotlib.figure

# import f4e_ace.checkvalue as cv
# import f4e_ace.data
import numpy as np

from f4enix.egroups import GROUP_STRUCTURES
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import Material, Zaid

# get a default library manager
LM = LibManager()

# Special MT values
UNITY_MT = -1
XI_MT = -2

# MTs to combine to generate associated plot_types
_INELASTIC = [mt for mt in endf.SUM_RULES[3] if mt != 27]
PLOT_TYPES_MT = {
    "total": endf.SUM_RULES[1],
    "scatter": [2] + _INELASTIC,
    "elastic": [2],
    "inelastic": _INELASTIC,
    "fission": [18],
    "absorption": [27],
    "capture": [101],
    "nu-fission": [18],
    "nu-scatter": [2] + _INELASTIC,
    "unity": [UNITY_MT],
    "slowing-down power": [2] + [XI_MT],
    "damage": [444],
}

# Types of plots to plot linearly in y
PLOT_TYPES_LINEAR = {
    "nu-fission / fission",
    "nu-scatter / scatter",
    "nu-fission / absorption",
    "fission / absorption",
}


ELEMENT_NAMES = list(endf.data.ELEMENT_SYMBOL.values())[1:]


def _xsdir_to_ace_filepaths(filename: str | os.PathLike) -> dict[tuple[str, str], dict]:
    xsdir_ace_files = {}
    from f4enix.input.xsdirpyne import Xsdir

    full_xsdir = Xsdir(filename)
    xsdir_tables = full_xsdir.tables
    for table in xsdir_tables:
        suffix = table.name.split(".")[1]
        if table.zaid.isdigit():
            name = str(table.zaid)
            # name = LM.get_zaidname(table.zaid)
            if (name, suffix) not in xsdir_ace_files.keys():
                xsdir_ace_files[(name, suffix)] = {
                    "path": os.path.join(
                        os.path.dirname(filename), Path(table.filename)
                    )
                }
    return xsdir_ace_files


def _continuous_energy_xs_nuclide(
    this: str,
    types: Sequence[str | int],
    library: str,
    ace_filepaths: str | os.PathLike | dict,
    temperature: float = 294.0,
) -> tuple[np.ndarray, Sequence[Callable]]:
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : str
        Nuclide object to source data from
    types : Iterable of str or Integral
        The type of cross sections to calculate; values can either be those
        in PLOT_TYPES or keys from endf.reaction.REACTION_MT which
        correspond to a reaction description e.g '(n,2n)' or integers which
        correspond to reaction channel (MT) numbers.
    library : str
        The library to use for the cross sections
    ace_filepaths : str or os.PathLike | dict
        Path to the ACE file or xsdir file to use for the plot or dictionary of xsdir files
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : Iterable of Callable
        Requested cross section functions

    """

    # Convert temperature to format needed for access in the library
    strT = f"{int(round(temperature))}K"
    T = temperature

    # Now we can create the data sets to be plotted
    energy_grid = []
    xs = []

    if isinstance(ace_filepaths, str) or isinstance(ace_filepaths, os.PathLike):
        cross_section_lookup_dict = _xsdir_to_ace_filepaths(ace_filepaths)
    else:
        cross_section_lookup_dict = ace_filepaths

    lib = cross_section_lookup_dict[this, library]

    if lib is not None:
        nuc = endf.IncidentNeutron.from_ace(lib["path"])
        # Obtain the nearest temperature

        # endf package does not have a temperature or KTs attribute, so finding it through the key values
        # see issue https://github.com/paulromano/endf-python/issues/13
        reactions = nuc.reactions[1]
        nuc.temperatures = list(reactions.xs.keys())
        nuc.kTs = [
            float(temp.replace("K", "")) * endf.data.K_BOLTZMANN
            for temp in nuc.temperatures
        ]
        total_xs = nuc.reactions[1].xs
        all_energy_values = {}
        for temperature in total_xs.keys():
            tabulated_energies = total_xs[temperature].x
            all_energy_values[temperature] = tabulated_energies
        nuc.energy = all_energy_values

        if strT in nuc.temperatures:
            nucT = strT
        else:
            delta_T = np.array(nuc.kTs) - T * endf.data.K_BOLTZMANN
            closest_index = np.argmin(np.abs(delta_T))
            nucT = nuc.temperatures[closest_index]

        energy_grid = nuc.energy[nucT]

        # Parse the types
        mts = []
        ops = []
        yields = []
        for line in types:
            if line in PLOT_TYPES_MT.keys():
                tmp_mts = [
                    mtj
                    for mti in PLOT_TYPES_MT[line]
                    for mtj in nuc._get_reaction_components(mti)
                ]
                mts.append(tmp_mts)
                if line.startswith("nu"):
                    yields.append(True)
                else:
                    yields.append(False)
                if XI_MT in tmp_mts:
                    ops.append((np.add,) * (len(tmp_mts) - 2) + (np.multiply,))
                else:
                    ops.append((np.add,) * (len(tmp_mts) - 1))
            elif line in endf.reaction.REACTION_MT:
                mt_number = endf.reaction.REACTION_MT[line]
                # cv.check_type("MT in types", mt_number, Integral)
                # cv.check_greater_than("MT in types", mt_number, 0)
                tmp_mts = nuc._get_reaction_components(mt_number)
                mts.append(tmp_mts)
                ops.append((np.add,) * (len(tmp_mts) - 1))
                yields.append(False)
            elif isinstance(line, int):
                # Not a built-in type, we have to parse it ourselves
                # cv.check_type("MT in types", line, Integral)
                # cv.check_greater_than("MT in types", line, 0)
                tmp_mts = nuc._get_reaction_components(line)
                mts.append(tmp_mts)
                ops.append((np.add,) * (len(tmp_mts) - 1))
                yields.append(False)
            else:
                raise TypeError("Invalid type", line)

        for i, mt_set in enumerate(mts):
            # Get the reaction xs data from the nuclide
            funcs = []
            op = ops[i]
            for mt in mt_set:
                if mt == 2:
                    funcs.append(nuc[mt].xs[nucT])
                elif mt in nuc:
                    if yields[i]:
                        # Get the total yield first if available. This will be
                        # used primarily for fission.
                        for prod in chain(nuc[mt].products, nuc[mt].derived_products):
                            if (
                                prod.particle == "neutron"
                                and prod.emission_mode == "total"
                            ):
                                func = Combination(
                                    [nuc[mt].xs[nucT], prod.yield_], [np.multiply]
                                )
                                funcs.append(func)
                                break
                        else:
                            # Total doesn't exist so we have to create from
                            # prompt and delayed. This is used for scatter
                            # multiplication.
                            func = None
                            for prod in chain(
                                nuc[mt].products, nuc[mt].derived_products
                            ):
                                if (
                                    prod.particle == "neutron"
                                    and prod.emission_mode != "total"
                                ):
                                    if func:
                                        func = Combination(
                                            [prod.yield_, func], [np.add]
                                        )
                                    else:
                                        func = prod.yield_
                            if func:
                                funcs.append(
                                    Combination([func, nuc[mt].xs[nucT]], [np.multiply])
                                )
                            else:
                                # If func is still None, then there were no
                                # products. In that case, assume the yield is
                                # one as its not provided for some summed
                                # reactions like MT=4
                                funcs.append(nuc[mt].xs[nucT])
                    else:
                        funcs.append(nuc[mt].xs[nucT])
                elif mt == UNITY_MT:
                    funcs.append(lambda x: 1.0)
                elif mt == XI_MT:
                    awr = nuc.atomic_weight_ratio
                    alpha = ((awr - 1.0) / (awr + 1.0)) ** 2
                    xi = 1.0 + alpha * np.log(alpha) / (1.0 - alpha)
                    funcs.append(lambda x: xi)
                else:
                    funcs.append(lambda x: 0.0)
            funcs = funcs if funcs else [lambda x: 0.0]
            xs.append(Combination(funcs, op))
    else:
        raise ValueError(this + " not in library")

    return energy_grid, xs


def _get_xs_zaid_list(
    zaids: list[Zaid],
    types: Sequence[str | int],
    library: str,
    ace_filepaths: str | os.PathLike | dict,
    temp: float,
    group_structure: str | None = None,
    integral_points: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for zaid in zaids:
        nuc = zaid.element + zaid.isotope
        temp_E, temp_xs = _continuous_energy_xs_nuclide(
            nuc, types, library, ace_filepaths, temp
        )
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        # use the log-log interpolation for the cross sections (5). 2 for lin-lin
        xs[nuc] = [
            endf.function.Tabulated1D(temp_E, temp_xs[line](temp_E), interpolation=5)
            for line in range(len(types))
        ]

    if group_structure:
        energy_grid = GROUP_STRUCTURES[group_structure]
    else:
        # Condense the data for every nuclide
        # First create a union energy grid
        energy_grid = E[0]
        for grid in E[1:]:
            energy_grid = np.union1d(energy_grid, grid)

    # Now we can combine all the nuclidic data
    data = np.zeros((len(types), len(energy_grid)))
    for line in range(len(types)):
        if types[line] == "unity":
            data[line, :] = 1.0
        else:
            for zaid in zaids:
                table = xs[nuc][line]
                nuc = zaid.element + zaid.isotope
                if group_structure:
                    # if we are performing a collapse in groups, perform an integral
                    values = []
                    for i, e in enumerate(energy_grid):
                        if i == 0:
                            # get the minimum in the e range
                            e_0 = np.min(table.x)
                        else:
                            e_0 = energy_grid[i - 1]
                        x = np.logspace(np.log10(e_0), np.log10(e), integral_points)
                        y = table(x)
                        values.append(np.trapz(y, x) / (e - e_0))
                    values = np.array(values)

                else:
                    # if continous, the simple interpolation is ok
                    values = table(energy_grid)

                data[line, :] += zaid.fraction * values

    return energy_grid, data


def get_xs(
    this: str | Material,
    types: list[int | str],
    library: str,
    ace_filepaths: str | dict | os.PathLike,
    group_structure: None | str = None,
    temperature: float = 294.0,
    integral_points: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """

    Parameters
    ----------
    this : str
        the nuclide to collapse in string format e.g. 'Li6'.
    types : list[int  |  str]
        The list of cross section reaction types to consider.
    library : str
        The library to use for cross section data.
    ace_filepaths : str | dict | os.PathLike
        The path to the ACE file or xsdir file to use for the plot or a dictionary of
        xsdir files.
    group_structure : str
        The group structure to use for collapsing cross sections.
    temperature : float, optional
        the temperature at which to evaluate cross sections. Defaults to 294.
    integral_points : int, optional
        The number of points to use for the integral calculation during a xs
        energy collapse. Defaults to 100. If the number of groups is low, it is
        recommended to increase this value to get a better integral of the xs in each
        energy bin.

    Returns
    -------
    _type_
        _description_
    """
    if isinstance(this, str):
        # either a nuclide or an element
        mat = Material.from_zaids([(this, 1)], LM, "00c")
    else:
        mat = this

    # Force normalization and atom fractions
    mat.switch_fraction("mass", LM, inplace=True)
    mat.switch_fraction("atom", LM, inplace=True)
    # Collect all zaids
    zaids = []
    for submat in mat.submaterials:
        zaids.extend(submat.zaidList)

    energy_grid, data = _get_xs_zaid_list(
        zaids,
        types,
        library,
        ace_filepaths,
        temp=temperature,
        group_structure=group_structure,
        integral_points=integral_points,
    )

    return energy_grid, data


# def get_exfor_data(isotope: str, mt: str | int) -> list:
#     """
#     Retrieve EXFOR data for a given isotope and reaction.

#     Args:
#         isotope (str): The isotope for which to retrieve the data.
#         mt (Union[str, int]): The reaction type or MT number.

#     Returns:
#         list: A list of dictionaries containing the retrieved data. Each dictionary
#         contains the following keys:
#             - 'energy': An array of energy values in MeV.
#             - 'cross-section': An array of cross-section values in barns.
#             - 'label': A string representing the author and year of the data.

#     Raises:
#         ValueError: If the specified MT value is not found in the MT reactions dictionary.

#     """
#     from x4i3 import exfor_manager

#     # fmt: off
#     ENDF_X4_dict = {
#         1: "N,TOT", 2: "N,EL", 3: "N,NON", 4: "N,INL", 5: "N,X", 10: "N,TOT",
#         11: "N,2N+D", 16: "N,2N", 17: "N,3N", 18: "N,F", 19: "N,F'", 20: "N,N+F",
#         21: "N,2N+F", 22: "N,N+A", 23: "N,N+3A", 24: "N,2N+A", 25: "N,3N+A",
#         27: "N,ABS", 28: "N,N+P", 29: "N,N+2A", 32: "N,N+D", 33: "N,N+T",
#         34: "N,N+HE3", 37: "N,4N", 38: "N,3N+F", 41: "N,2N+P", 42: "N,3N+P",
#         44: "N,N+2P", 45: "N,N+P+A", 51: "N,N'", 89: "N,N'", 90: "N,N'",
#         91: "N,N'", 101: "N,DIS", 102: "N,G", 103: "N,P", 104: "N,D", 105: "N,T",
#         106: "N,HE3", 107: "N,A", 108: "N,2A", 111: "N,2P", 112: "N,P+A",
#         113: "N,T+2A", 115: "N,P+D", 116: "N,P+T", 117: "N,D+A", 151: "N,RES",
#         201: "N,XN", 202: "N,XG", 203: "N,XP", 204: "N,XD", 205: "N,XT",
#         206: "N,XHE3", 207: "N,XA", 208: "N,XPi_pos", 209: "N,XPi_0",
#         210: "N,XPi_neg", 301: "heating", 444: "damage-energy production",
#         452: "N,nu_tot", 454: "N,ind_FY", 455: "N,nu_d", 456: "N,nu_p",
#         458: "N,rel_fis", 459: "FY_cum", 460: "N,g_bdf", 600: "N,P", 601: "N,P'",
#         649: "N,P'", 650: "N,D", 651: "N,D'", 699: "N,D'", 700: "N,T",
#         701: "N,T'", 749: "N,T'", 750: "N,HE3'", 751: "N,HE3'", 799: "N,HE3'",
#         800: "N,A", 801: "N,A'", 849: "N,A'", 875: "N,2N", 876: "N,2N",
#         889: "N,2N", 890: "N,2N",
#     }
#     # fmt: on

#     data_list = []
#     if mt not in ENDF_X4_dict.keys() and mt not in ENDF_X4_dict.values():
#         raise ValueError(
#             f"mt {mt} not found, available mt values are {ENDF_X4_dict.keys()} and {ENDF_X4_dict.values()}"
#         )

#     if isinstance(mt, int):
#         mt_int = mt
#         mt_equation = ENDF_X4_dict[mt]
#     elif isinstance(mt, str):
#         mt_int = [key for key, value in ENDF_X4_dict.items() if value == mt][0]
#         mt_equation = mt

#     if mt_int == 444:
#         raise ValueError("mt 444 is not a valid mt for this function")

#     db = exfor_manager.X4DBManagerDefault()
#     db_entry = db.retrieve(target=isotope, reaction=mt_equation, quantity="SIG")
#     for _, entry in db_entry.items():
#         datasets = entry.getSimplifiedDataSets()
#         for subentry_key, subentry_value in datasets.items():
#             if (
#                 subentry_value.simplified is True
#                 and len(subentry_value.reaction[0].quantity) == 1
#                 and subentry_value.reaction[0].quantity[0] == "SIG"
#                 and isinstance(subentry_value.reaction[0].quantity, list)
#             ):
#                 x_subentry, y_subentry = [], []
#                 energy = datasets[subentry_key].labels.index("Energy")
#                 xs_data = datasets[subentry_key].labels.index("Data")
#                 for row in subentry_value.data:
#                     x_subentry.append(row[energy])
#                     y_subentry.append(row[xs_data])
#                 data_list.append(
#                     {
#                         "energy": np.array(x_subentry) * 1e6,
#                         "cross-section": np.array(y_subentry),
#                         "label": f"{subentry_value.author[0]} {subentry_value.year}",
#                     }
#                 )
#     return data_list


class Combination:
    """Combination of multiple functions with a user-defined operator

    This class allows you to create a callable object which represents the
    combination of other callable objects by way of a series of user-defined
    operators connecting each of the callable objects.

    Parameters
    ----------
    functions : Iterable of Callable
        Functions to combine according to operations
    operations : Iterable of numpy.ufunc
        Operations to perform between functions; note that the standard order
        of operations will not be followed, but can be simulated by
        combinations of Combination objects. The operations parameter must have
        a length one less than the number of functions.


    Attributes
    ----------
    functions : Iterable of Callable
        Functions to combine according to operations
    operations : Iterable of numpy.ufunc
        Operations to perform between functions; note that the standard order
        of operations will not be followed, but can be simulated by
        combinations of Combination objects. The operations parameter must have
        a length one less than the number of functions.

    """

    def __init__(self, functions, operations):
        self.functions = functions
        self.operations = operations

    def __call__(self, x):
        ans = self.functions[0](x)
        for i, operation in enumerate(self.operations):
            ans = operation(ans, self.functions[i + 1](x))
        return ans

    # def __eq__(self, other):
    #     if isinstance(other, type(self)):
    #         for key, value in self.__dict__.items():
    #             if isinstance(value, np.ndarray):
    #                 if not np.array_equal(value, other.__dict__.get(key)):
    #                     return False
    #             else:
    #                 return value == other.__dict__.get(key)
    #     else:
    #         return False

    #     return True

    @property
    def functions(self):
        return self._functions

    @functions.setter
    def functions(self, functions):
        # cv.check_type("functions", functions, Iterable, Callable)
        self._functions = functions

    @property
    def operations(self):
        return self._operations

    @operations.setter
    def operations(self, operations):
        # cv.check_type("operations", operations, Iterable, np.ufunc)
        # length = len(self.functions) - 1
        # cv.check_length("operations", operations, length, length_max=length)
        self._operations = operations
