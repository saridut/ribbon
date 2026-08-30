"""
Functions to construct shapes based on plane spirals.

"""
import math
import numpy as np
from .spiralcartesian import SpiralCornu


def write_xyz(x, y, z, atom_symbol='C', filename='out.xyz', title=''):
    """
    Writes atom positions to a XYZ file.

    Parameters
    ----------
    x, y, z : 1D array_like
        *x*, *y*, and *z* coordinates of the atoms. All three must have the
        same length.
    atom_symbol : str
        Atom symbol. If more than two characters, will be trucated to the first
        two characters.
    filename : path_like
        Name of output file (XYZ format) with extension '.xyz'.
    title : str
        A title string for the XYZ file.

    Returns
    -------
    None

    """
    if not (len(x)==len(y)==len(z)):
        raise ValueError("x, y, and z must have the same length."
                f" len(x) = {len(x)}, len(y) = {len(y)}, len(z) = {len(z)}.")
    if len(atom_symbol) > 2:
        print(f"Truncating atom_symbol(= {atom_symbol} to {atom_symbol[:2]}.")
        sym = atom_symbol[0:2]
    else:
        sym = atom_symbol
    num_atoms = len(x)
    with open(filename, 'w') as fh:
        fh.write(f"{num_atoms}\n")
        fh.write(f"{title}\n")
        for i in range(num_atoms):
            fh.write(f"{sym} {x[i]:g} {y[i]:g} {z[i]:g}\n")


def extrude(xcoords, ycoords, zdist, dz):
    """
    Creates a sheet of monoatomic thickness by extruding a spiral and returns the
    positions of the atoms of the sheet.

    Parameters
    ----------
    xcoords, ycoords : 1D array_like
        X & Y coordinates of the points on the planar curve.
    zdist : float
        Extrusion distance (along *z*-direction)
    dz : float
        Spacing along the *z*-direction

    Returns
    -------
    X, Y, Z : tuple of 1D ndarrays
        Coordinates of the atoms of the sheet.

    """
    x = np.asarray(xcoords, dtype=np.float64)
    y = np.asarray(ycoords, dtype=np.float64)
    nx = x.size
    ny = y.size
    if nx < 2:
        raise ValueError(f"Length of `x` must be >= 2. Length = {nx}.")
    if ny < 2:
        raise ValueError(f"Length of `y` must be >= 2. Length = {ny}.")
    if nx != ny:
        raise ValueError(f"`x` & `y` must have same lengths."
                         f" Length of `x` = {nx} & that of `y` = {ny}.")
    if dz < 0.0:
        raise ValueError(f"`dz`(= {dz:g}) must be >= 0.")
    if (zdist>0.0) and (zdist < dz):
        raise ValueError(f"`zdist`(= {zdist:g}) must be >= dz.")
    if zdist==0.0:
        z = np.array([0.0], dtype=np.float64)
    else:
        z = np.linspace(0.0, zdist, round(1+zdist/dz))
    num_atoms = len(x)*len(z)

    #Create atom coordinate arrays
    X = np.zeros((num_atoms,), dtype=np.float64)
    Y = np.zeros_like(X)
    Z = np.zeros_like(X)

    #Fill coordinate arrays and write to file
    nx = x.size; nz = z.size
    for k in range(nz):
        zval = z[k]
        for i in range(nx):
            xval = x[i]
            yval = y[i]
            j = k*nx + i
            X[j] = xval
            Y[j] = yval
            Z[j] = zval
    return (X, Y, Z)


def double_spiral(spiral, L, ds=0.5, same=False, end=1, f=0.0):
    """
    Create a double spiral from a single spiral.

    The two spirals may be joined at either of the two ends. In addition, before
    joining, a the curvature of fraction of the arclength of both spirals may be
    linearized (if `f` > 0). Linearization is done by replacing a length `fL`
    of the curve by a Cornu spiral, whose curvature varies linearly with
    arclength, such that there is no discontinuity in the curvature of the
    resulting double spiral.

    Parameters
    ----------
    spiral : SpiralPolarBase | SpiralCartesianBase
        A spiral object
    L : float
        Arclength of `spiral`. The double spiral will have arclength `2L`.
    ds : float
        Spacing along the arclength.
    same : bool
        Whether the spirals are on the same side or not.
    end : 0 | 1
        If 0, the ``s=0`` end. If 1, the ``s=L`` end.
    f : float
        Fraction of the arclength `L` over which the curvature will be
        linearized.

    Returns
    -------
    s_ds, kappa_ds, x_ds, y_ds : tuple of 1D ndarrays 
        The arclength coordinates (`s_ds`), curvatures (`kappa_ds`),
        *x*-coordinates (`x_ds`), and *y*-coordinates (`y_ds`) of
        the double spiral.

    Warnings
    --------
    * Linearization may not always be necessary, and in some cases may lead to
      significant changes in the *derivative* of the curvature along the
      arclength.

    * Depending on the parameters of the spiral and which end is chosen for
      joining, the double spiral may self-intersect.

    """
    if L <= 0:
        raise ValueError(f"`L`(= {L:g}) must be > 0.")
    if (ds <= 0) or (ds >= L):
        raise ValueError(f"`ds`(= {ds:g}) must be > 0 and < L(= {L:g}).")
    if end not in [0, 1]:
        raise ValueError(f"`end`(= {end:g}) must be 0 or 1.")
    if (f < 0) or (f > 1):
        raise ValueError(f"`f`(= {f:g}) must be >= 0 and <= 1.")
    #Grid points along the arclength
    s = np.linspace(0.0, L, round(1+L/ds))
    if f == 0:
        #Grid points along the arclength
        kappa = spiral.get_curvature(s)
        x, y = spiral.get_cart_coords(s)
        if end == 0:
            t = spiral.get_tangent(s[0])
        else:
            t = -spiral.get_tangent(s[-1])
            x = np.flipud(x); y = np.flipud(y)
            kappa = np.flipud(kappa)
            s = L - np.flipud(s)
        x -= x[0]; y -= y[0]
        coords = np.column_stack((x, y))
        phi = math.atan2(t[1], t[0]) #Angle w.r.t. x-axis
        #rotate by -phi to align t along x-axis
        ct = math.cos(-phi); st = math.sin(-phi)
        rotmat = np.array([[ct, -st],[st, ct]], dtype=np.float64)
        coords = np.dot(coords, rotmat.T)
    else:
        if end == 0:
            s_j = f*L
            ldist = s_j
            s_ss = s[s>=s_j]
            s_lin = s[s<s_j]
            t = spiral.get_tangent(s_j)
            kappa_ss = spiral.get_curvature(s_ss)
            x_ss, y_ss = spiral.get_cart_coords(s_ss)
        else:
            s_j = (1-f)*L
            ldist = L - s_j
            s_ss = s[s<=s_j]
            s_lin = s[s>s_j]
            t = -spiral.get_tangent(s_j)
            kappa_ss = spiral.get_curvature(s_ss)
            x_ss, y_ss = spiral.get_cart_coords(s_ss)
            x_ss = np.flipud(x_ss); y_ss = np.flipud(y_ss)
            kappa_ss = np.flipud(kappa_ss)
            s_ss = L - np.flipud(s_ss)
            s_lin = L - np.flipud(s_lin)
        x_j, y_j = spiral.get_cart_coords(s_j)
        x_ss -= x_j; y_ss -= y_j

        cornu = SpiralCornu(abs(kappa_ss[0])/ldist)
        t_c = cornu.get_tangent(ldist)
        kappa_lin = cornu.get_curvature(s_lin)
        x_lin, y_lin = cornu.get_cart_coords(s_lin)
        x_j_lin, y_j_lin = cornu.get_cart_coords(ldist)
        coords_lin = np.column_stack((x_lin, y_lin))
        if kappa_ss[0] < 0:
            t_c[:] = -t_c
            kappa_lin[:] = -kappa_lin
            x_lin[:] = -x_lin; y_lin[:] = -y_lin
            x_j_lin = -x_j_lin; y_j_lin = -y_j_lin

        cos_phi = np.linalg.vecdot(t, t_c)/(
            np.linalg.vector_norm(t)*np.linalg.vector_norm(t_c))
        if math.isclose(abs(cos_phi), 1.0):
            cos_phi = np.sign(cos_phi)
        phi = np.acos(cos_phi)
        if t_c[0]*t[1] < t_c[1]*t[0]:
            phi = -phi
        ct = math.cos(-phi); st = math.sin(-phi)
        rotmat = np.array([[ct, -st],[st, ct]], dtype=np.float64)
        coords_ss = np.column_stack((x_ss, y_ss))
        coords_ss = np.dot(coords_ss, rotmat.T)
        coords_ss[:,0] += x_j_lin
        coords_ss[:,1] += y_j_lin

        coords = np.concatenate((coords_lin, coords_ss))
        s = np.concatenate((s_lin, s_ss))
        kappa = np.concatenate((kappa_lin, kappa_ss))
    #Make double spiral
    s_ds = np.concatenate((L-s[:0:-1], L+s))
    if same:
        kappa_ds = np.concatenate((kappa[:0:-1], kappa))
        x_ds = np.concatenate( (-coords[:0:-1,0], coords[:,0]) )
        y_ds = np.concatenate( ( coords[:0:-1,1], coords[:,1]) )
    else:
        kappa_ds = np.concatenate((kappa[:0:-1], -kappa))
        x_ds = np.concatenate( (-coords[:0:-1,0], coords[:,0]) )
        y_ds = np.concatenate( (-coords[:0:-1,1], coords[:,1]) )

    return s_ds, kappa_ds, x_ds, y_ds
