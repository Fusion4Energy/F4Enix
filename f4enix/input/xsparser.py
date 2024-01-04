from f4enix.input.acepyne import *
from f4enix.input.libmanager import LibManager
import numpy as np
from numpy import logspace

def get_xs(mat_nuc, mt: int, lib: str, lib_manager: LibManager, 
           material: bool = True):

    if material:
        energies, xs, fractions = get_material_xs(mat_nuc, mt, lib, lib_manager=lib_manager)
        return _build_material_xs(energies, xs, fractions)
    else:
        data_t = get_zaid_xstable(str(mat_nuc)+'.'+lib, lib, lib_manager)
        return get_reaction_xs(data_t, mt)

def get_material_xs(mat, mt: int, lib : str, lib_manager: LibManager):

    mat.translate(lib, lib_manager=lib_manager)
    mat.switch_fraction('atom', lib_manager=lib_manager)

    energies = []
    xs = []
    fractions = []
    for submat in mat.submaterials:
        for zaid in submat.zaidList:
            data_t = get_zaid_xstable(zaid.name, lib, lib_manager)
            
            x, y = get_reaction_xs(data_t, mt)

            if y is None:
                continue

            energies.append(x)
            xs.append(y)
            fractions.append(zaid.fraction)

    return energies, xs, fractions

def get_zaid_xstable(zaid_name: str, lib : str, lib_manager: LibManager):
    # dict of ENDF reactions MT number
    
    # Build dict of nuclear data libraries
    for xs_tab in lib_manager.XS.tables:
        if zaid_name == xs_tab.name:
            break
    filepath = lib_manager.XS.directory + xs_tab.filename.replace("/", "\\")

    try:
        data_XS = Library(filepath)
    except FileNotFoundError:
        print('File not Found Error: ACEfile of nuclide \
                ' + zaid_name + ' not found in library ' + lib + '\n')
        exit()

    try:
        data_XS.read()
    except ValueError:
        s = 'Error in reading ACEfile of nuclide ' + zaid_name + '\n'
        print(s)
        exit()

    if zaid_name.split('.')[1] == '00c':
        ace_zaid = zaid_name.split('.')[0] + '.800nc'
    else:
        ace_zaid = zaid_name

    return data_XS.tables[ace_zaid]

def get_reaction_xs(data_t: NeutronTable, mt: int):

    if mt in data_t.reactions.keys():
        ace_mt = data_t.reactions[mt]
        return data_t.energy[ace_mt.IE:], data_t.reactions[mt].sigma
    
    elif mt == 1:
        return data_t.energy, data_t.sigma_t
    
    elif mt == 101:
        return data_t.energy, data_t.sigma_a

    elif mt == 301:
        return data_t.energy, data_t.heating

    else:
        return None, None

def _build_material_xs(energies: list, xs: list, fractions: list):

    if len(energies) == 0:
        return None, None

    combined_energies = energies[0]

    for en_x in energies[1:]:
        combined_energies = np.union1d(combined_energies, en_x)

    for l, xs_ in enumerate(xs):
        if l == 0:
            mat_xs = fractions[l] * np.interp(combined_energies, energies[l], xs_)
        else:
            mat_xs = mat_xs + fractions[l] * np.interp(combined_energies, energies[l], xs_)
    
    return combined_energies, mat_xs / sum(fractions)


def partial_energy_matrix_mono(E_g: np.ndarray, E_n: np.ndarray,
                               slope: int = -1):
    """Generates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    Here, both of the energy arrays must be monotonic. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : 1d numpy float array 
        Lower resolution energy group structure [MeV] that is of length G+1. 
        Ordered based on slope.
    E_n : 1d numpy float array 
        Higher resolution energy group structure [MeV] that is of length N+1. 
        Ordered based on slope.
    slope : int, optional
        Gives the monotonicity of E_g and E_n.  If positive, then they are 
        monotonicly increasing (lowest-to-highest).  If negative, they are
        monotonicly decreasing (highest-to-lowest).

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    # Some convienence parameters

    G = E_g.shape[0] - 1
    N = E_n.shape[0] - 1

    index_E_n = np.arange(N+1)

    # Get the interior points for each gth group in n-space
    if slope < 0:
        inner_mask = np.array([(E_n <= E_g[g]) & (E_g[g+1] <= E_n) for g in range(G)])
    elif 0 < slope:
        inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])
    else:
        raise ValueError("slope must be positive or negative.")

    # Get the upper and lower nth index for every gth group
    lower_index = np.array([index_E_n[inner_mask[g]][0] for g in range(G)])
    upper_index = np.array([index_E_n[inner_mask[g]][-1] for g in range(G)])

    # Convert the mask to initialize the partial enery matrix
    # Hack off the last index of the mask to make the right size
    pem = np.array(inner_mask[:, :-1], dtype=float)

    # Check for partial contibutions at the edges
    for g in range(G):
        # Lower bound
        lig = lower_index[g]
        if lig != 0:
            pem[g,lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g,uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem

def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    The group structures must have the same monotonicity. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : sequence of floats
        Lower resolution energy group structure [MeV] that is of length G+1. 
    E_n : sequence of floats
        Higher resolution energy group structure [MeV] that is of length N+1. 

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """

    E_g = np.asarray(E_g, dtype=float)
    E_n = np.asarray(E_n, dtype=float)

    if (E_g[:-1] >= E_g[1:]).all() and (E_n[:-1] >= E_n[1:]).all():
        # Both energy arrays are monotonically decreasing
        assert E_g[0] <= E_n[0]
        assert E_n[-1] <= E_g[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, -1)
    elif (E_g[:-1] <= E_g[1:]).all() and (E_n[:-1] <= E_n[1:]).all():
        # Both energy arrays are monotonically increasing
        assert E_n[0] <= E_g[0]
        assert E_g[-1] <= E_n[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, 1)
    else:
        raise ValueError("E_g and E_n are not both monotonic in the same direction.")

    return pem

######################
### Group Collapse ###
######################
def phi_g(E_g: np.ndarray, E_n: np.ndarray, phi_n: np.ndarray):
    """Calculates a lower resolution flux, phi_g, from a lower resolution group stucture E_g, 
    a higher resolution groups E_n, and a higher resolution flux phi_n.

    Parameters
    ----------
    E_g : sequence of floats 
        Lower resolution energy group structure [MeV] that is of length G+1. 
    E_n : sequence of floats 
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : sequence of floats
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over (length N).  

    Returns
    -------
    phi_g : numpy array of floats 
        The flux collapsed to G energy groups.
    """
    pem = partial_energy_matrix(E_g, E_n)
    phi_g = np.dot(pem, phi_n)
    return phi_g


def group_collapse(sigma_n: np.ndarray, phi_n: np.ndarray, 
                   phi_g=None, partial_energies=None, E_g=None, E_n=None, 
                   weights=None):
    """Calculates the group cross-sections for a nuclide for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, 
    and G < N.  

    This function has two optional ways of being called.  If the group boundaries
    E_g and E_n are provided, this will collapse the flux automatically.  However, 
    if a partial energy matrix and flux collapse has already been performed you can
    shortcut their recalculation by calling this function with the phi_g and 
    partial_energies keyword arguments.

    Parameters
    ----------
    sigma_n : array-like of floats
        A high-fidelity cross-section.
    phi_n : array-like of floats
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section 
        over (length N).  
    phi_g : array-like of floats, optional
        The low-fidelity flux [n/cm^2/s] to collapse the fission cross-section 
        down to (length G).  If present, partial_energies is needed as well.
    partial_energies : 2D array-like of floats, optional
        A partial energy matrix as provided by a previous call to the function
        partial_energy_matrix().  If present, phi_g is needed as well.
    E_g : array-like of floats, optional
        Lower resolution energy group structure [MeV] that is of length G+1.
        If present, E_n is needed as well.
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
        If present, E_g is needed as well.

    Returns
    -------
    sigma_g : ndarray
        An array of the collapsed fission cross-section.
    """
    if (phi_g is not None) and (partial_energies is not None):
        pem = partial_energies
    elif (phi_g is None) and (partial_energies is not None):
        pem = partial_energies        
        if weights is None:
            phi_g = np.dot(pem, phi_n)
        else:
            phi_g = np.dot(pem, phi_n * weights)
    elif (E_g is not None) and (E_n is not None):
        pem =  partial_energy_matrix(E_g, E_n)
        if weights is None:
            phi_g = np.dot(pem, phi_n)
        else:
            phi_g = np.dot(pem, phi_n * weights)
    else:
        msg = "Either partial_energies or E_g and E_n must both not be None."
        raise ValueError(msg)

    # Calulate partial group collapse
    if weights is None:
        sigma_g = np.dot(pem, sigma_n * phi_n) / phi_g
    else:
        sigma_g = np.dot(pem, sigma_n * phi_n * weights) / phi_g
    sigma_g[np.isnan(sigma_g)] = 0.0  # handle zero flux that causes NaNs later.
    return sigma_g

def ninespace(start, stop, num=50, endpoint=True):
    """Splits the range into one-minus-log-uniform bins defined by num points.
    In the vernacular, the space is 'split in the nines'.  Note that this assumes
    base = 10.0.

    Parameters
    ----------
    start : number
        The starting value of the sequence.
    stop : number
        The final value of the sequence, unless endpoint
        is False.  In that case, num + 1 values are spaced over the
        interval in nines-space, of which all but the last (a sequence of
        length num) are returned.
    num : integer, optional
        Number of samples to generate. See endpoint.
    endpoint : boolean, optional
        If true, stop is the last sample. Otherwise, it is not included.

    Returns
    -------
    samples : ndarray
        num samples, equally spaced in the nines.

    Examples
    --------
    >>> ninespace(0.9, 0.9999, 4)
        array([0.9, 0.99, 0.999, 0.9999])

    """
    log_start = np.log10(1.0 - start)
    log_stop  = np.log10(1.0 - stop)
    samples = 1.0 - logspace(log_start, log_stop, num, endpoint)
    return samples

def stair_step(x, y):
    """Makes a 1d data set of boundaries (x) and cell-centered values (y) 
    into stair step arrays of the same length.  The returned arrays are 
    suitable for graphing.  This is especially useful in energy vs. spectrum
    data where there are G+1 boundaries and G data points.

    Parameters
    ----------
    x : sequence
        Data of length G+1.
    y : sequence
        Data of length G.

    Returns
    -------
    xss : ndarray 
        Stair-step version of x data, length 2G.
    yss : ndarray 
        Stair-step version of y data, length 2G.

    Examples
    --------
    >>> x = [0.1, 1.0, 10.0, 100.0]
    >>> y = [2.0, 3.0, 4.0]
    >>> bins.stair_step(x, y)
    (array([   0.1,    1. ,    1. ,   10. ,   10. ,  100. ]),
    array([ 2.,  2.,  3.,  3.,  4.,  4.]))

    """
    # Grab the number of points, G
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    G = len(y)
    assert G + 1 == len(x)

    # Initilaize data
    xss = np.empty(2*G, dtype=x.dtype)
    yss = np.empty(2*G, dtype=y.dtype)

    # Set values
    xss[:-1:2] = x[:-1]
    xss[1::2] = x[1:]
 
    yss[::2] = y
    yss[1::2] = y

    return xss, yss

def pointwise_linear_collapse(x_g, x, y):
    """Collapses pointwise data to G groups based on a linear interpolation
    between the points. This is useful for collapsing cross section data.

    Parameters
    ----------
    x_g : array-like
        Group boundaries, length G+1 for G groups, must be monotonic in the 
        same direction as x.
    x : array-like
        Pointwise abscissa to be collapsed, must be monotonic in the same direction
        as x_g and have the same length as y.
    y : array-like
        Pointwise data to be interpolated, must have the same length as x.

    Returns
    -------
    y_g : np.ndarray
        The group collapsed data, length G. 
    """
    return pointwise_collapse(x_g, x, y)

def pointwise_collapse(x_g, x, y, logx=False, logy=False, log=False):
    """Collapses pointwise data to G groups based on a interpolation
    between the points. This is useful for collapsing cross section data.

    Parameters
    ----------
    x_g : array-like
        Group boundaries, length G+1 for G groups, must be monotonic in the 
        same direction as x.
    x : array-like
        Pointwise abscissa to be collapsed, must be monotonic in the same direction
        as x_g and have the same length as y.
    y : array-like
        Pointwise data to be interpolated, must have the same length as x.
    logx: bool, optional, default=False
        lin-log interpolation
    logy: bool, optional, default=False
        log-lin interpolation
    log : bool, optional, default=False
        log-log interpolation

    Returns
    -------
    y_g : np.ndarray
        The group collapsed data, length G. 
    """

    G = x_g.shape[0] - 1
    N = x.shape[0] - 1
    y_g = np.empty(G, dtype='float64')

    # Ensure input is monotonically increasing/decreasing
    if not (np.all(np.diff(x) >= 0.) or np.all(np.diff(x) <= 0)):
            raise ValueError("x must be monotonically "  
                             "increasing/decreasing.")
    if not (np.all(np.diff(x_g) >= 0.) or np.all(np.diff(x_g) <= 0)):
            raise ValueError("x_g must be monotonically "  
                             "increasing/decreasing.")

    reversed = False
    if x_g[0] > x_g[-1]:
        # monotonically decreasing, make increasing so logic is simpler
        x_g = x_g[::-1]
        x = x[::-1]
        y = y[::-1]
        reversed = True

    # Handle logrithmic interpolations
    if logx or log:
        if x_g[0] <= 0.0 or x_g[0] <= 0.0:
            raise ValueError("x values must be positive for logrithmic interpolation")
        x_g = np.log(x_g)
        x = np.log(x)
    if logy or log:
        if y[0] <= 0.0:
            raise ValueError("y values must be positive for logrithmic interpolation")
        y = np.log(y)

    n0 = 0
    n1 = 1
    for g0 in range(G):
        g1 = g0 + 1
        val = 0.0
        while x[n1] <= x_g[g1] and n1 < N:
            if x_g[g0] <= x[n0]:
                # interpolation fully in group, can take midpoint
                val += 0.5 * (y[n1] + y[n0]) * (x[n1] - x[n0])
            else: 
                # lower bound intersection
                ylower = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g0] - x[n0]) + y[n0]
                val += 0.5 * (y[n1] + ylower) * (x[n1] - x_g[g0])
            n0 += 1
            n1 += 1
        # compute end point
        if x_g[g0] <= x[n0]:
            yupper = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g1] - x[n0]) + y[n0]
            val += 0.5 * (yupper + y[n0]) * (x_g[g1] - x[n0])
        else:
            yupper = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g1] - x[n0]) + y[n0]
            ylower = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g0] - x[n0]) + y[n0]
            val += 0.5 * (yupper + ylower) * (x_g[g1] - x_g[g0])
        y_g[g0] = val / (x_g[g1] - x_g[g0])
    if reversed:
        y_g = y_g[::-1]

    # Handle logrithmic interpolation
    if logy or log:
        y_g = np.exp(y_g)

    return y_g

def collapse_xs(x_g, x, y):
    # assuming sorted groups from low to high value

    idx = np.searchsorted(x_g, np.min(x), side='right')

    if idx != 0:
        new_x_g = np.insert(x_g[idx:], 0, x_g[idx - 1])
        xs = pointwise_linear_collapse(new_x_g, x, y)
        xs = np.concatenate((np.zeros(len(x_g)-len(new_x_g)), xs))
    else:
        xs = pointwise_linear_collapse(x_g, x, y)

    return xs


def compute_rr(flux, xs, time:float = -1, a_dens:float = -1):

    if len(flux) == len(xs):
        raise ValueError("group structure of flux and xs is different")
    
    rr = np.dot(flux, xs)
    if time > 0:
        rr = rr * time
    if a_dens > 0:
        rr = rr * a_dens

    return rr
