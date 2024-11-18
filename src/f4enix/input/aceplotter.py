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
from f4enix.input.materials import Material

# get a default library manager
LM = LibManager()

# Supported keywords for continuous-energy cross section plotting
PLOT_TYPES = [
    "total",
    "scatter",
    "elastic",
    "inelastic",
    "fission",
    "absorption",
    "capture",
    "nu-fission",
    "nu-scatter",
    "unity",
    "slowing-down power",
    "damage",
]

# Supported keywords for multi-group cross section plotting
PLOT_TYPES_MGXS = [
    "total",
    "absorption",
    "scatter",
    "fission",
    "kappa-fission",
    "nu-fission",
    "prompt-nu-fission",
    "deleyed-nu-fission",
    "chi",
    "chi-prompt",
    "chi-delayed",
    "inverse-velocity",
    "beta",
    "decay-rate",
    "unity",
]
# Create a dictionary which can be used to convert PLOT_TYPES_MGXS to the
# f4e_ace.XSdata attribute name needed to access the data
_PLOT_MGXS_ATTR = {
    line: line.replace(" ", "_").replace("-", "_") for line in PLOT_TYPES_MGXS
}
_PLOT_MGXS_ATTR["scatter"] = "scatter_matrix"

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

# Minimum and maximum energies for plotting (units of eV)
_MIN_E = 1.0e-5
_MAX_E = 20.0e6


ELEMENT_NAMES = list(endf.data.ELEMENT_SYMBOL.values())[1:]


def _get_legend_label(this, type) -> str:
    """Gets a label for the element or nuclide or material and reaction plotted"""
    # if isinstance(this, str):
    #     if type in f4e_ace.data.DADZ:
    #         z, a, m = endf.data.zam(this)
    #         da, dz = f4e_ace.data.DADZ[type]
    #         gnds_name = endf.data.gnds_name(z + dz, a + da, m)
    #         return f"{this} {type} {gnds_name}"
    #     return f"{this} {type}"
    # elif this.name == "":
    #     return f"Material {this.id} {type}"
    # else:
    #     return f"{this.name} {type}"
    return "TODO put here the appropriate legend"


def _get_yaxis_label(reactions: dict, divisor_types) -> str:
    """Gets a y axis label for the type of data plotted"""

    heat_values = {"heating", "heating-local", "damage-energy"}

    # if all the types are heating a different stem and unit is needed
    if all(set(value).issubset(heat_values) for value in reactions.values()):
        stem = "Heating"
    elif all(isinstance(item, str) for item in reactions.keys()):
        for nuc_reactions in reactions.values():
            for reaction in nuc_reactions:
                if reaction in heat_values:
                    raise TypeError(
                        "Mixture of heating and Microscopic reactions. "
                        "Invalid type for plotting"
                    )
        stem = "Microscopic"
    elif all(isinstance(item, Material) for item in reactions.keys()):
        stem = "Macroscopic"
    else:
        msg = "Mixture of f4e_ace.Material and elements/nuclides. Invalid type for plotting"
        raise TypeError(msg)

    if divisor_types:
        mid, units = "Data", ""
    else:
        mid = "Cross Section"
        units = {
            "Macroscopic": "[1/cm]",
            "Microscopic": "[b]",
            "Heating": "[eV-barn]",
        }[stem]

    return f"{stem} {mid} {units}"


def _get_title(reactions) -> str:
    """Gets a title for the type of data plotted"""
    if len(reactions) == 1:
        (this,) = reactions
        name = this.name if isinstance(this, Material) else this
        return f"Cross Section Plot For {name}"
    else:
        return "Cross Section Plot"


def plot_xs(
    reactions: dict,
    library: str,
    ace_filepaths: str | os.PathLike,
    divisor_types=None,
    temperature=294.0,
    axis: matplotlib.axes.Axes | None = None,
    enrichment: float | None = None,
    plot_CE: bool = True,
    energy_axis_units: str = "eV",
    **kwargs,
) -> matplotlib.figure.Figure:
    """Creates a figure of continuous-energy cross sections for this item.

    Parameters
    ----------
    reactions : dict
        keys can be either a nuclide or element in string form or an
        Material object. Values are the type of cross sections to
        include in the plot.
    library : str
        The library to use for the cross sections
    ace_filepaths : str | os.PathLike
        Path to the ACE file or xsdir file to use for the plot
    divisor_types : Iterable of str, optional
        The type of cross sections to divide by. If not specified, no division
        will be performed. Defaults to None.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    axis : matplotlib.axes.Axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None.
    plot_CE : bool, optional
        Denotes whether or not continuous-energy will be plotted. Defaults to
        plotting the continuous-energy data.
    energy_axis_units : {'eV', 'keV', 'MeV'}
        Units used on the plot energy axis
    **kwargs :
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        If axis is None, then a Matplotlib Figure of the generated
        cross section will be returned. Otherwise, a value of
        None will be returned as the figure and axes have already been
        generated.

    """
    import matplotlib.pyplot as plt

    assert isinstance(plot_CE, bool)
    if energy_axis_units not in {"eV", "keV", "MeV"}:
        raise ValueError(f"Invalid energy axis units: {energy_axis_units}")

    axis_scaling_factor = {"eV": 1.0, "keV": 1e-3, "MeV": 1e-6}

    # Generate the plot
    if axis is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        fig = None
        ax = axis

    all_types = []

    for this, types in reactions.items():
        all_types = all_types + types

        if plot_CE:
            # cv.check_type("this", this, (str, f4e_ace.Material))
            # Calculate for the CE cross sections
            E, data = continuous_energy_xs(
                this, types, library, ace_filepaths, temperature, enrichment
            )
            if divisor_types:
                # cv.check_length("divisor types", divisor_types, len(types))
                Ediv, data_div = continuous_energy_xs(
                    this, divisor_types, library, ace_filepaths, temperature, enrichment
                )

                # Create a new union grid, interpolate data and data_div on to that
                # grid, and then do the actual division
                Enum = E[:]
                E = np.union1d(Enum, Ediv)
                data_new = np.zeros((len(types), len(E)))

                for line in range(len(types)):
                    data_new[line, :] = np.divide(
                        np.interp(E, Enum, data[line, :]),
                        np.interp(E, Ediv, data_div[line, :]),
                    )
                    if divisor_types[line] != "unity":
                        types[line] = types[line] + " / " + divisor_types[line]
                data = data_new

        E *= axis_scaling_factor[energy_axis_units]

        # Plot the data
        for i in range(len(data)):
            data[i, :] = np.nan_to_num(data[i, :])
            if np.sum(data[i, :]) > 0.0:
                ax.plot(E, data[i, :], label=_get_legend_label(this, types[i]))

    # Set to loglog or semilogx depending on if we are plotting a data
    # type which we expect to vary linearly
    if set(all_types).issubset(PLOT_TYPES_LINEAR):
        ax.set_xscale("log")
        ax.set_yscale("linear")
    else:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(f"Energy [{energy_axis_units}]")
    if plot_CE:
        ax.set_xlim(
            _MIN_E * axis_scaling_factor[energy_axis_units],
            _MAX_E * axis_scaling_factor[energy_axis_units],
        )
    else:
        ax.set_xlim(E[-1], E[0])

    ax.set_ylabel(_get_yaxis_label(reactions, divisor_types))
    ax.legend(loc="best")
    ax.set_title(_get_title(reactions))

    return fig


def continuous_energy_xs(
    this: str | Material,
    types: Sequence,
    library: str,
    ace_filepaths: str | os.PathLike,
    temperature: float = 294.0,
    enrichment: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : str or Material
        Object to source data from. Nuclides and elements should be input as a
        str
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    library : str
        The library to use for the cross sections
    ace_filepaths : str or os.PathLike
        Path to the ACE file or xsdir file to use for the plot
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # # Check types
    # cv.check_type("this", this, (str, f4e_ace.Material))
    # cv.check_type("temperature", temperature, Real)
    # if enrichment:
    #     cv.check_type("enrichment", enrichment, Real)

    if isinstance(this, str):
        if this in ELEMENT_NAMES:
            energy_grid, data = continuous_energy_xs_element(
                this, types, library, ace_filepaths, temperature
            )
        else:
            energy_grid, xs = _continuous_energy_xs_nuclide(
                this,
                types,
                library,
                ace_filepaths,
                temperature,
            )

            # Convert xs (Iterable of Callable) to a grid of cross section values
            # calculated on the points in energy_grid for consistency with the
            # element and material functions.
            data = np.zeros((len(types), len(energy_grid)))
            for line in range(len(types)):
                data[line, :] = xs[line](energy_grid)
    else:
        energy_grid, data = continuous_energy_xs_material(
            this, types, library, ace_filepaths, temperature
        )

    return energy_grid, data


def xsdir_to_ace_filepaths(filename: str | os.PathLike) -> dict[tuple[str, str], dict]:
    xsdir_ace_files = {}
    from f4enix.input.xsdirpyne import Xsdir

    full_xsdir = Xsdir(filename)
    xsdir_tables = full_xsdir.tables
    for table in xsdir_tables:
        suffix = table.name.split(".")[1]
        if table.zaid.isdigit():
            name = LM.get_zaidname(table.zaid)
            if (name, suffix) not in xsdir_ace_files.keys():
                xsdir_ace_files[(name, suffix)] = {
                    "path": str(Path(filename).parent / table.filename)
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
        cross_section_lookup_dict = xsdir_to_ace_filepaths(ace_filepaths)
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
            if line in PLOT_TYPES:
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


def continuous_energy_xs_element(
    this: str,
    types: Sequence[int | str],
    library: str,
    ace_filepaths: str | os.PathLike | dict,
    temperature: float = 294.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates continuous-energy cross sections of an element

    Parameters
    ----------
    this : str
        Element name (e.g. 'Li')
    types : Sequence[int  |  str]
        The type of cross sections to calculate
    library : str
        The library to use for the cross sections
    ace_filepaths : str | os.PathLike | dict
        Path to the ACE file or xsdir file to use for the plot
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Energies at which cross sections are calculated, in units of eV
        Cross sections calculated at the energy grid described by energy_grid
    """
    T = temperature

    # # Expand elements in to nuclides with atomic densities
    nuclides = LM.expand_element(this)
    # # For ease of processing split out the nuclide and its fraction
    # nuc_fractions = {nuclide[0]: nuclide[1] for nuclide in nuclides}
    # # Create a dict of [nuclide name] = nuclide object to carry forward
    # # with a common nuclides format between f4e_ace.Material and Elements
    # nuclides = {nuclide[0]: nuclide[0] for nuclide in nuclides}

    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for idx, row in nuclides.iterrows():
        nuc = idx
        temp_E, temp_xs = continuous_energy_xs(nuc, types, library, ace_filepaths, T)
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        xs[nuc] = [
            endf.function.Tabulated1D(temp_E, temp_xs[line])
            for line in range(len(types))
        ]

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
            for idx, row in nuclides.iterrows():
                nuc = idx
                data[line, :] += row["Mean value"] * xs[nuc][line](energy_grid)

    return energy_grid, data


def continuous_energy_xs_material(
    this, types, library, ace_filepaths, temperature=294.0
):
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : f4e_ace.Material or str
        Object to source data from. Element can be input as str
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    if isinstance(this, f4e_ace.Material):
        if this.temperature is not None:
            T = this.temperature
        else:
            T = temperature
    else:
        raise TypeError("Invalid type for this, must be Material()")

    # Expand elements in to nuclides with atomic densities
    nuc_fractions = this.get_nuclide_atom_densities()
    # Create a dict of [nuclide name] = nuclide object to carry forward
    # with a common nuclides format between f4e_ace.Material and Elements
    nuclides = {nuclide: nuclide for nuclide in nuc_fractions}

    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for nuclide in nuclides.items():
        name = nuclide[0]
        nuc = nuclide[1]
        temp_E, temp_xs = continuous_energy_xs(nuc, types, library, ace_filepaths, T)
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        xs[name] = [
            endf.function.Tabulated1D(temp_E, temp_xs[line])
            for line in range(len(types))
        ]

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
            for nuclide in nuclides.items():
                name = nuclide[0]
                data[line, :] += nuc_fractions[name] * xs[name][line](energy_grid)

    return energy_grid, data


def _multi_group_energy_xs_nuclide(
    this,
    types,
    library,
    ace_filepaths,
    group_structure,
    temperature: float = 294,
) -> tuple[np.ndarray, np.ndarray]:
    groups = GROUP_STRUCTURES[group_structure]

    _, cross_section_functions = _continuous_energy_xs_nuclide(
        this=this,
        types=types,
        library=library,
        ace_filepaths=ace_filepaths,
        temperature=temperature,
    )

    mg_cross_sections = []
    for cross_section_function in cross_section_functions:
        # TODO improve this section by using average
        # here we currently sample the continuous energy cross section at the group boundaries
        # we could take the average value over the group which would be more accurate
        # but this is a simple way to get a rough multigroup cross section
        mg_cross_section = cross_section_function(groups)
        mg_cross_sections.append(mg_cross_section)

    return groups, np.ndarray(mg_cross_sections)


def _multi_group_energy_xs_element(
    this,
    types,
    library,
    ace_filepaths,
    group_structure,
    enrichment=None,
    temperature: float = 294,
) -> tuple[np.ndarray, np.ndarray]:
    energy_groups = GROUP_STRUCTURES[group_structure]

    T = temperature

    # Expand elements in to nuclides with atomic densities
    nuclides = LM.expand_element(this)
    # For ease of processing split out the nuclide and its fraction
    # nuc_fractions = {nuclide[0]: nuclide[1] for nuclide in nuclides}
    # # Create a dict of [nuclide name] = nuclide object to carry forward
    # # with a common nuclides format between f4e_ace.Material and Elements
    # nuclides = {nuclide[0]: nuclide[0] for nuclide in nuclides}

    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for idx, row in nuclides.iterrows():
        temp_E, temp_xs = continuous_energy_xs(idx, types, library, ace_filepaths, T)
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        xs[idx] = [
            endf.function.Tabulated1D(temp_E, temp_xs[line])
            for line in range(len(types))
        ]

    # Now we can combine all the nuclidic data
    data = np.zeros((len(types), len(energy_groups)))
    for line in range(len(types)):
        if types[line] == "unity":
            data[line, :] = 1.0
        else:
            for idx, row in nuclides.iterrows():
                data[line, :] += row["Mean Value"] * xs[idx][line](energy_groups)

    return energy_groups, data


def _multi_group_energy_xs_material(
    this: Material,
    types: Sequence,
    library: str,
    ace_filepaths: str | os.PathLike | dict,
    group_structure: str,
    temperature: float = 294.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : f4e_ace.Material
        This can only be a material here
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    library: str
        The library to use for the cross sections
    ace_filepaths : str or os.PathLike | dict
        Path to the ACE file or xsdir file to use for the plot or dictionary of xsdir files
    group_structure : str
        The group structure to use for collapsing cross sections.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    if not isinstance(this, Material):
        raise TypeError("Invalid type for this, must be Material()")
    T = temperature
    # be sure to work in atom fraction
    this.switch_fraction("atom", LM, inplace=True)

    energy_groups = GROUP_STRUCTURES[group_structure]

    # # Expand elements in to nuclides with atomic densities
    # nuc_fractions = this.get_nuclide_atom_densities()
    # # Create a dict of [nuclide name] = nuclide object to carry forward
    # # with a common nuclides format between f4e_ace.Material and Elements
    # nuclides = {nuclide: nuclide for nuclide in nuc_fractions}

    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for nuclide in this.zaids:
        nuc = nuclide.element + nuclide.isotope
        temp_E, temp_xs = continuous_energy_xs(nuc, types, library, ace_filepaths, T)
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        xs[nuc] = [
            endf.function.Tabulated1D(temp_E, temp_xs[line])
            for line in range(len(types))
        ]

    data = np.zeros((len(types), len(energy_groups)))
    for line in range(len(types)):
        if types[line] == "unity":
            data[line, :] = 1.0
        else:
            for nuclide in this.zaids:
                nuc = nuclide.element + nuclide.isotope
                data[line, :] += nuclide.fraction * xs[nuc][line](energy_groups)

    return energy_groups, data


def multi_group_energy_xs(
    this: str | Material,
    types: Sequence,
    library: str,
    ace_filepaths: str | os.PathLike | dict,
    group_structure: str,
    temperature: float = 294.0,
    enrichment: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : str or f4e_ace.Material
        Object to source data from. Nuclides and elements should be input as a
        str
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    library: str
        The library to use for the cross sections
    ace_filepaths : str or os.PathLike | dict
        Path to the ACE file or xsdir file to use for the plot or dictionary of xsdir files
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # Check types
    # cv.check_type("this", this, (str, f4e_ace.Material))
    # cv.check_type("temperature", temperature, Real)
    # if enrichment:
    #     cv.check_type("enrichment", enrichment, Real)

    if isinstance(this, str):
        if this in ELEMENT_NAMES:
            energy_grid, data = _multi_group_energy_xs_element(
                this,
                types,
                library,
                ace_filepaths,
                group_structure,
                enrichment,
                temperature,
            )
        else:
            energy_grid, data = _multi_group_energy_xs_nuclide(
                this,
                types,
                library,
                ace_filepaths,
                group_structure,
                temperature,
            )
    elif isinstance(this, Material):
        energy_grid, data = _multi_group_energy_xs_material(
            this, types, library, ace_filepaths, group_structure, temperature
        )
    else:
        raise ValueError('Invalid type for "this" argument')
    return energy_grid, data


def get_collapsed_cross_section_nuclide(
    this: str,
    types: list[int | str],
    library: str,
    ace_filepaths: str | dict | os.PathLike,
    group_structure: str,
    temperature: float = 294,
) -> list[list[float]]:
    """
    Calculate the collapsed cross sections for the given nuclide.

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

    Returns
    -------
    list[list[float]]
        A 2D list of collapsed cross sections for each reaction type and group.
    """
    groups = GROUP_STRUCTURES[group_structure]

    energy_grid, cross_section_functions = _continuous_energy_xs_nuclide(
        this=this,
        types=types,
        library=library,
        ace_filepaths=ace_filepaths,
        temperature=temperature,
    )

    # we join the group boundaries with the energy grid to make sure we are getting the
    # cross section at the edges of the group and and sudden resonances
    concatenated_array = np.concatenate((energy_grid, groups))

    sorted_array = np.sort(concatenated_array)

    all_collapsed_cross_sections = []
    # looks over the different cross section reactions requested by the types argument
    for cross_section_function in cross_section_functions:
        group_collapsed_cross_sections = []

        for i in range(len(groups) - 1):
            lower_bound = groups[i]
            upper_bound = groups[i + 1]
            bin_width = upper_bound - lower_bound

            bin_values = sorted_array[
                (sorted_array >= lower_bound) & (sorted_array < upper_bound)
            ]

            cross_sections = np.array(
                [cross_section_function(value) for value in bin_values]
            )

            # Apply the Trapezium Rule
            area = np.trapz(cross_sections, bin_values)

            average_cross_section = area / bin_width

            group_collapsed_cross_sections.append(average_cross_section)

        all_collapsed_cross_sections.append(group_collapsed_cross_sections)

    return all_collapsed_cross_sections


def get_energy_cross_section(
    this: str | Material,
    types: list[int | str],
    library: str,
    ace_filepaths: str | dict | os.PathLike,
    group_structure: None | str = None,
    temperature: float = 294.0,
    enrichment: float | None = None,
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
    enrichment : float | None, optional
        the enrichment of the material. Defaults to None.

    Returns
    -------
    _type_
        _description_
    """
    if group_structure is None:
        energy_grid, data = continuous_energy_xs(
            this, types, library, ace_filepaths, temperature, enrichment
        )
    else:
        energy_grid, data = multi_group_energy_xs(
            this,
            types,
            library,
            ace_filepaths,
            group_structure,
            temperature,
            enrichment,
        )
    return energy_grid, data


def get_exfor_data(isotope: str, mt: str | int) -> list:
    """
    Retrieve EXFOR data for a given isotope and reaction.

    Args:
        isotope (str): The isotope for which to retrieve the data.
        mt (Union[str, int]): The reaction type or MT number.

    Returns:
        list: A list of dictionaries containing the retrieved data. Each dictionary
        contains the following keys:
            - 'energy': An array of energy values in MeV.
            - 'cross-section': An array of cross-section values in barns.
            - 'label': A string representing the author and year of the data.

    Raises:
        ValueError: If the specified MT value is not found in the MT reactions dictionary.

    """
    from x4i3 import exfor_manager

    # fmt: off
    ENDF_X4_dict = {
        1: "N,TOT", 2: "N,EL", 3: "N,NON", 4: "N,INL", 5: "N,X", 10: "N,TOT",
        11: "N,2N+D", 16: "N,2N", 17: "N,3N", 18: "N,F", 19: "N,F'", 20: "N,N+F",
        21: "N,2N+F", 22: "N,N+A", 23: "N,N+3A", 24: "N,2N+A", 25: "N,3N+A",
        27: "N,ABS", 28: "N,N+P", 29: "N,N+2A", 32: "N,N+D", 33: "N,N+T",
        34: "N,N+HE3", 37: "N,4N", 38: "N,3N+F", 41: "N,2N+P", 42: "N,3N+P",
        44: "N,N+2P", 45: "N,N+P+A", 51: "N,N'", 89: "N,N'", 90: "N,N'",
        91: "N,N'", 101: "N,DIS", 102: "N,G", 103: "N,P", 104: "N,D", 105: "N,T",
        106: "N,HE3", 107: "N,A", 108: "N,2A", 111: "N,2P", 112: "N,P+A",
        113: "N,T+2A", 115: "N,P+D", 116: "N,P+T", 117: "N,D+A", 151: "N,RES",
        201: "N,XN", 202: "N,XG", 203: "N,XP", 204: "N,XD", 205: "N,XT",
        206: "N,XHE3", 207: "N,XA", 208: "N,XPi_pos", 209: "N,XPi_0",
        210: "N,XPi_neg", 301: "heating", 444: "damage-energy production",
        452: "N,nu_tot", 454: "N,ind_FY", 455: "N,nu_d", 456: "N,nu_p",
        458: "N,rel_fis", 459: "FY_cum", 460: "N,g_bdf", 600: "N,P", 601: "N,P'",
        649: "N,P'", 650: "N,D", 651: "N,D'", 699: "N,D'", 700: "N,T",
        701: "N,T'", 749: "N,T'", 750: "N,HE3'", 751: "N,HE3'", 799: "N,HE3'",
        800: "N,A", 801: "N,A'", 849: "N,A'", 875: "N,2N", 876: "N,2N",
        889: "N,2N", 890: "N,2N",
    }
    # fmt: on

    data_list = []
    if mt not in ENDF_X4_dict.keys() and mt not in ENDF_X4_dict.values():
        raise ValueError(
            f"mt {mt} not found, available mt values are {ENDF_X4_dict.keys()} and {ENDF_X4_dict.values()}"
        )

    if isinstance(mt, int):
        mt_int = mt
        mt_equation = ENDF_X4_dict[mt]
    elif isinstance(mt, str):
        mt_int = [key for key, value in ENDF_X4_dict.items() if value == mt][0]
        mt_equation = mt

    if mt_int == 444:
        raise ValueError("mt 444 is not a valid mt for this function")

    db = exfor_manager.X4DBManagerDefault()
    db_entry = db.retrieve(target=isotope, reaction=mt_equation, quantity="SIG")
    for _, entry in db_entry.items():
        datasets = entry.getSimplifiedDataSets()
        for subentry_key, subentry_value in datasets.items():
            if (
                subentry_value.simplified is True
                and len(subentry_value.reaction[0].quantity) == 1
                and subentry_value.reaction[0].quantity[0] == "SIG"
                and isinstance(subentry_value.reaction[0].quantity, list)
            ):
                x_subentry, y_subentry = [], []
                energy = datasets[subentry_key].labels.index("Energy")
                xs_data = datasets[subentry_key].labels.index("Data")
                for row in subentry_value.data:
                    x_subentry.append(row[energy])
                    y_subentry.append(row[xs_data])
                data_list.append(
                    {
                        "energy": np.array(x_subentry) * 1e6,
                        "cross-section": np.array(y_subentry),
                        "label": f"{subentry_value.author[0]} {subentry_value.year}",
                    }
                )
    return data_list


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

    def __eq__(self, other):
        if isinstance(other, type(self)):
            for key, value in self.__dict__.items():
                if isinstance(value, np.ndarray):
                    if not np.array_equal(value, other.__dict__.get(key)):
                        return False
                else:
                    return value == other.__dict__.get(key)
        else:
            return False

        return True

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
