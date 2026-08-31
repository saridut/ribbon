import copy
import math
import numbers
import os
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import RectBivariateSpline
import rotlib

#Classes for real numbers.
RealNumber = (numbers.Real, np.number)

class Ribbon(object):
    """
    Attributes
    ----------
    length : float
        Length of the ribbon.
    width : float
        Width of the ribbon.
    thickness : float
        Thickness of the ribbon. If zero, the ribbon is a 2D surface in
        3-space.
    l, m, n : float or callable
        Curvature along length, width, and thickness directions, respectively.
        If any of the curvatures is callable, it should be of the form *y = f
        (x)*, where *x* (float) is the arclength coordinate and *y* (float) is
        the curvature.
    u, v, w : 1d ndarray
        Coordinates of the reference grid along length, width, and thickness
        directions, respectively.
    mline : (n,3) ndarray
        Coordinates of the grid points on the ribbon midline in current
        configuration.
    msurf : (m,n,3) ndarray
        Coordinates of the grid points on the ribbon midsurface in the current
        configuration.
    grid : (m,n,p,3) ndarray
        Coordinates of the grid points of the ribbon in the current configuration.
    atom_refpos : (n,3) ndarray
        Atom coordinates in the reference state.
    atom_pos : (n,3) ndarray
        Atom coordinates in the current configuration.

    Parameters
    ----------
    length : float
        Length of the ribbon.
    width : float
        Width of the ribbon.
    thickness : float
        Thickness of the ribbon. If zero, the ribbon is a 2D surface in
        3-space.
    gspl : float | None
        Grid spacing in the length direction. If `None`, specify the number of
        grid points using the keyword argument `ngpl`.
    gspw : float | None
        Grid spacing in the width direction. If `None`, specify the number of
        grid points using the keyword argument `ngpw`.
    gspt : float | None
        Grid spacing in the thickness direction. Ignored if `thickness = 0`.
        If `thickness != 0` and `gspt = None`, specify the number of grid points
        using the keyword argument `ngpt`.
    ngpl : int | None, optional
        Number of grid points in the length direction. Must be > 4. Used only
        if `gspl = None`.
    ngpw : int | None, optional
        Number of grid points in the width direction. Must be > 4. Used only if
        `gspw = None`.
    ngpt : int | None, optional
        Number of grid points in the thickness direction. Must be > 2. Used
        only if `thickness != 0` and `gspt = None`.

    """
    def __init__(self, length, width, thickness, gspl, gspw, gspt,
                 ngpl=None, ngpw=None, ngpt=None):
        if not ( isinstance(length, RealNumber) and length > 0 ):
            raise ValueError( f"`length`(= {length:g}) must be an instance of"
                " numbers.Real or numpy.number and must be > 0.")
        else:
            self.length = length

        if not ( isinstance(width, RealNumber) and width > 0 ):
            raise ValueError( f"`width`(= {width:g}) must be an instance of"
                " numbers.Real or numpy.number and must be > 0.")
        else:
            self.width = width

        if not ( isinstance(thickness, RealNumber) and thickness >= 0 ):
            raise ValueError( f"`thickness`(= {thickness:g}) must be an instance of"
                " numbers.Real or numpy.number and must be >= 0.")
        else:
            self.thickness = thickness

        if gspl is not None:
            n = math.ceil(self.length/gspl) + 1
        elif ngpl is not None:
            n = ngpl
        else:
            raise ValueError("`gspl` and `ngpl` cannot both be None.") 
        if n < 4:
            raise ValueError(
                f"Number of grid points (= {n}) along the length direction"
                f" must be >= 4. Reduce grid spacing or increase number of"
                f" grid points."
                )
        self.u = np.linspace(0, self.length, n, dtype=np.float64)

        if gspw is not None:
            n = math.ceil(self.width/gspw) + 1
        elif ngpw is not None:
            n = ngpw
        else:
            raise ValueError("`gspw` and `ngpw` cannot both be None.") 
        if n < 4:
            raise ValueError(
                f"Number of grid points (= {n}) along the width direction"
                f" must be >= 4. Reduce grid spacing or increase number of"
                f" grid points."
                )
        self.v = np.linspace(-self.width/2, self.width/2, n, dtype=np.float64)

        if self.thickness > 0:
            if gspt is not None:
                n = math.ceil(self.thickness/gspt) + 1
            elif ngpt is not None:
                n = ngpt
            else:
                raise ValueError("`gspt` and `ngpt` cannot both be None.") 
            if n < 2:
                raise ValueError(
                    f"Number of grid points (= {n}) along the thickness"
                    f" direction must be >= 2. Reduce grid spacing or increase"
                    f" number of grid points."
                    )
            self.w = np.linspace(-self.thickness/2, self.thickness/2, n,
                                 dtype=np.float64)

        self.mline = np.zeros((self.u.size, 3))
        self._d1 = np.zeros_like(self.mline)
        self._d2 = np.zeros_like(self.mline)
        self._d3 = np.zeros_like(self.mline)
        self.msurf = np.zeros((self.u.size, self.v.size, 3))
        if self.thickness > 0:
            self.grid = np.zeros((self.u.size, self.v.size, self.w.size, 3))
            self._msurf_du = np.zeros_like(self.msurf)
            self._msurf_dv = np.zeros_like(self.msurf)
            self._normals = np.zeros_like(self.msurf)
        else:
            self.grid = np.zeros_like(self.msurf)
            self._msurf_du = np.zeros((0,0))
            self._msurf_dv = np.zeros_like(self._msurf_du)
            self._normals = np.zeros_like(self._msurf_du)
        self.atom_refpos = np.zeros((0,0))
        self.atom_pos = np.zeros_like(self.atom_refpos)
        self._ap_ms = np.zeros_like(self.atom_pos)
        self._ap_du = np.zeros_like(self.atom_pos)
        self._ap_dv = np.zeros_like(self.atom_pos)
        self._ap_normals = np.zeros_like(self.atom_pos)

    def set_atom_refpos(self, atom_refpos, copy=True, atom_pos=None):
        """
        Setter for the reference positions of atoms.

        Parameters
        ----------
        atom_refpos : (n,3) ndarray
            Atom positions in the reference state. Must be within the domain
            [0, :attr:`.length`] x [-:attr:`.width`/2, :attr:`.width`/2] 
            x [-:attr:`.thickness`/2, :attr:`.thickness`/2].
        copy : bool
            Whether to make a copy of the input array `atom_refpos`.
        atom_pos : (n,3) ndarray
            External buffer for :attr:`.atom_pos`. If `None`, an array will be
            created.

        Returns
        -------
        None

        """
        self.atom_refpos = np.asarray(atom_refpos, dtype=np.float64, copy=copy)
        if self.atom_refpos.shape == self.atom_pos.shape:
            self.atom_pos[...] = 0.0
            self._ap_ms[...] = 0.0
            self._ap_du[...] = 0.0
            self._ap_dv[...] = 0.0
            self._ap_normals[...] = 0.0
        else:
            if atom_pos is None:
                self.atom_pos = np.zeros_like(self.atom_refpos)
            else:
                self.atom_pos = atom_pos
                self.atom_pos[...] = 0.0
            self._ap_ms = np.zeros_like(self.atom_pos)
            self._ap_du = np.zeros_like(self.atom_pos)
            self._ap_dv = np.zeros_like(self.atom_pos)
            self._ap_normals = np.zeros_like(self.atom_pos)

    def set_curvatures(self, l, m, n, radius=None, pitch=None):
        """
        Setter for the curvatures along the three directions.

        If any of the curvatures is callable, it should be of the form *y = f
        (x)*, where *x* (float) is the arclength coordinate and *y* (float) is
        the curvature (or radius and pitch).

        Parameters
        ----------
        l : None | float | callable
            Curvature along the length direction. If `None`, `m` should also be
            `None` and `radius` and `pitch` must be specified.
        m : None | float | callable
            Twist along the length direction. If `None`, `l` should also be
            `None` and `radius` and `pitch` must be specified.
        n : float | callable
            Curvature along the width direction 
        radius : None | float
            Radius of curvature of the ribbon. If `None`, `pitch` should also
            be `None` and `l` and `m` must be specified.
        pitch : None | float
            Pitch of the ribbon. If `None`, `radius` should also be `None` and
            `l` and `m` must be specified.

        Returns
        -------
        None
        """
        if (l is not None) and (m is not None):
            if not (isinstance(l, RealNumber) or callable(l) ):
                raise ValueError(f"`l`(= {l}) must be a float or callable.")
            else:
                self.l = l
            if not (isinstance(m, RealNumber) or callable(m) ):
                raise ValueError(f"`m`(= {m}) must be a float or callable.")
            else:
                self.m = m
        elif (l is None) and (m is None) and (radius is not None) and \
                (pitch is not None):
            if not isinstance(radius, RealNumber):
                raise ValueError(f"`radius`(= {radius}) must be float.")
            if not isinstance(pitch, RealNumber):
                raise ValueError(f"`pitch`(= {pitch}) must be a float.")
            if np.isposinf(radius) or np.isposinf(pitch):
                self.l = 0.0; self.m = 0.0
            else:
                fpi2 = 4*np.pi**2
                den = fpi2*radius**2 + pitch**2
                self.l = fpi2*radius/den
                self.m = 2*np.pi*pitch/den
        else:
            raise ValueError(
                    f"Valid cases are (1) `l` and `m` are not None, or (2) `l`"
                    f" and `m` are None and `radius` and `pitch` are not None."
                    f" Input given: `l`(= {l}), `m`(= {m}),"
                    f" `radius`(= {radius}), and `pitch`(= {pitch})."
                    )

        if not (isinstance(n, RealNumber) or callable(n) ):
            raise ValueError(f"`n`(= {n}) must be a float or callable.")
        else:
            self.n = n

        if not (callable(self.l) or callable(self.m)):
            self.meth = 'direct'
        else:
            self.meth = 'ode'
        self.mline[...] = 0.0
        self._d1[...] = 0.0
        self._d2[...] = 0.0
        self._d3[...] = 0.0
        self.msurf[...] = 0.0
        self.grid[...] = 0.0
        self._msurf_du[...] = 0.0
        self._msurf_dv[...] = 0.0
        self._normals[...] = 0.0
        self.atom_pos[...] = 0.0
        self._ap_ms[...] = 0.0
        self._ap_du[...] = 0.0
        self._ap_dv[...] = 0.0
        self._ap_normals[...] = 0.0

    def get_radius(self):
        """
        Returns the radius along the ribbon midline.

        Returns
        -------
        float or tuple of 1D ndarrays
            Radius along the ribbon midline. If all curvatures are constants, a
            *float* is returned. If any of the curvatures is *callable*, a
            tuple *(x,y)* is returned, where *x* contains the arclength
            coordinates of the ribbon midline and *y* contains the
            corresponding radii.
        """
        if not (callable(self.l) or callable(self.m)):
            if self.l == self.m == 0:
                R = math.inf
            else:
                R = self.l/(self.l**2+self.m**2)
            return R
        else:
            if callable(self.l):
                l = [self.l(x) for x in self.u]
            else:
                l = [self.l for x in self.u]
            if callable(self.m):
                m = [self.m(x) for x in self.u]
            else:
                m = [self.m for x in self.u]
            R = [math.inf if x==y==0 else x/(x**2+y**2) for (x,y) in zip(l,m)]
            return (self.u, np.array(R, dtype=np.float64, copy=None, ndmin=1))

    def get_pitch(self):
        """
        Returns the pitch along the ribbon midline.

        Returns
        -------
        float or tuple of 1D ndarrays
            Pitch along the the ribbons midline. If all curvatures are
            constants, a *float* is returned. If any of the curvatures is
            *callable*, a tuple *(x,y)* is returned, where *x* contains the
            arclength coordinates of the ribbon midline and *y* contains the
            corresponding pitch values.
        """
        if not (callable(self.l) or callable(self.m)):
            if self.l == self.m == 0:
                P = math.inf
            else:
                P = 2*math.pi*self.m/(self.l**2+self.m**2)
            return P
        else:
            if callable(self.l):
                l = [self.l(x) for x in self.u]
            else:
                l = [self.l for x in self.u]
            if callable(self.m):
                m = [self.m(x) for x in self.u]
            else:
                m = [self.m for x in self.u]
            P = [math.inf if x==y==0 else 2*math.pi*x/(x**2+y**2)
                 for (x,y) in zip(l,m)]
            return (self.u, np.array(P, dtype=np.float64, copy=None, ndmin=1))


    def get_gauss_curvature(self):
        """
        Returns the gaussian curvature along the ribbon midline.

        Returns
        -------
        float or tuple of 1D ndarrays
            Gaussian curvature along the ribbons midline. If all curvatures are
            constants, a *float* is returned. If any of the curvatures is
            *callable*, a tuple *(x,y)* is returned, where *x* contains the
            arclength coordinates of the ribbon midline and *y* contains the
            corresponding gaussian curvatures.
        """
        if not (callable(self.l) or callable(self.m) or callable(self.n)):
            kg = self.l*self.n - self.m**2
            return kg
        else:
            if callable(self.l):
                l = [self.l(x) for x in self.u]
            else:
                l = [self.l for x in self.u]
            if callable(self.m):
                m = [self.m(x) for x in self.u]
            else:
                m = [self.m for x in self.u]
            if callable(self.n):
                n = [self.n(x) for x in self.u]
            else:
                n = [self.n for x in self.u]
            kg = [x*z-y**2 for (x,y,z) in zip(l,m,n)]
            return (self.u, np.asarray(kg))


    def get_mean_curvature(self):
        """
        Returns the mean curvature along the ribbon midline.

        Returns
        -------
        float or tuple of 1D ndarrays
            Mean curvature along the ribbons midline. If all curvatures are
            constants, a *float* is returned. If any of the curvatures is
            *callable*, a tuple *(x,y)* is returned, where *x* contains the
            arclength coordinates of the ribbon midline and *y* contains the
            corresponding mean curvatures.
        """
        if not (callable(self.l) or callable(self.n)):
            km = 0.5*(self.l+self.n)
            return km
        else:
            if callable(self.l):
                l = [self.l(x) for x in self.u]
            else:
                l = [self.l for x in self.u]
            if callable(self.n):
                n = [self.n(x) for x in self.u]
            else:
                n = [self.n for x in self.u]
            km = [0.5*(x+y) for (x,y) in zip(l,n)]
            return (self.u, np.asarray(km))


    def get_theta(self):
        """
        Returns the angle (in radians) between the principal curvature
        direction and the lengthwise direction along the ribbon midline.

        Returns
        -------
        float or tuple of 1D ndarrays
            Angle (in radians). If all curvatures are constants, a *float* is
            returned. If any of the curvatures is *callable*, a tuple *(x,y)*
            is returned, where *x* contains the arclength coordinates of the
            ribbon midline and *y* contains the corresponding angle.
        """
        if not (callable(self.l) or callable(self.m) or callable(self.n)):
            theta = 0.5*math.atan(2*self.m/(self.l-self.n)) \
                    if self.l != self.n else math.pi/4
            return theta
        else:
            if callable(self.l):
                l = [self.l(x) for x in self.u]
            else:
                l = [self.l for x in self.u]
            if callable(self.m):
                m = [self.m(x) for x in self.u]
            else:
                m = [self.m for x in self.u]
            if callable(self.n):
                n = [self.n(x) for x in self.u]
            else:
                n = [self.n for x in self.u]
            theta = [0.5*math.atan(2*y/(x-z)) if x!=z else math.pi/4
                     for (x,y,z) in zip(l,m,n)]
            return (self.u, np.asarray(theta))


    def create(self, orient_along=[1,0,0]):
        """
        Constructs the ribbon and orients it along `orient_along`.

        Parameters
        ----------
        orient_along : (3,) array_like or None
            If not ``None`` and the curvatures are constant, the ribbon axis
            will be oriented along `orient_along`. This does not need to be a
            unit vector.

        Returns
        -------
        None

        """
        if self.meth == 'direct':
            if self.l==self.m==0:
                self._create_direct_zero()
            else:
                self._create_direct(orient_along)
        elif self.meth == 'ode':
            self._create_ode(orient_along)

    def translate_to_center(self):
        """
        Center is the center of the ribbon midline. This will modify current
        atom positions and all grids.

        """
        center = self.mline.mean(axis=0)
        self.mline -= center
        self.msurf -= center
        self.grid -= center
        self.atom_pos -= center

    def _create_direct_zero(self):
        #Create the midline
        for i in range(self.u.size):
            u = self.u[i]
            self.mline[i,:] = np.array([0,0,u])
            self._d1[i,:] = np.array([1,0,0])
            self._d2[i,:] = np.array([0,1,0])
            self._d3[i,:] = np.array([0,0,1])
        self._create_msga()


    def _create_direct(self, orient_along):
        omega = (  self.l, 0, -self.m)
        omega_mag2 = omega[0]*omega[0] + omega[2]*omega[2]
        iomega_mag2 = 1/omega_mag2
        omega_mag = math.sqrt(omega_mag2)
        iomega_mag = 1/omega_mag

        #Create the midline
        for i in range(self.u.size):
            u = self.u[i]
            sn = math.sin(omega_mag*u)
            cs = math.cos(omega_mag*u)

            self.mline[i,0] = -omega[0]*omega[2]*iomega_mag2*(u - sn*iomega_mag)
            self.mline[i,1] = -omega[0]*iomega_mag2*(1 - cs)
            self.mline[i,2] = omega[2]*omega[2]*iomega_mag2 \
                    *(u - sn*iomega_mag) + sn*iomega_mag

            self._d1[i,0] =  omega[0]*omega[0]*iomega_mag2*(1-cs) + cs
            self._d1[i,1] = -omega[2]*iomega_mag*sn
            self._d1[i,2] = -omega[0]*omega[2]*iomega_mag2*(1-cs)

            self._d2[i,0] = omega[2]*iomega_mag*sn
            self._d2[i,1] = cs
            self._d2[i,2] = omega[0]*iomega_mag*sn

            self._d3[i,0] = -omega[0]*omega[2]*iomega_mag2*(1-cs)
            self._d3[i,1] = -omega[0]*iomega_mag*sn
            self._d3[i,2] =  omega[2]*omega[2]*iomega_mag2*(1-cs) + cs

            #Alternative formulation for omega[0] = -l
            #TODO: check for sign consistency. Ideally, omega[0]=l,
            #omega[2] = m. Check if this may be the case.
           #self.mline[i,0] = omega[0]*omega[2]*iomega_mag2*(u - sn*iomega_mag)
           #self.mline[i,1] = omega[0]*iomega_mag2*(1 - cs)
           #self.mline[i,2] = omega[2]*omega[2]*iomega_mag2 \
           #        *(u - sn*iomega_mag) + sn*iomega_mag

           #self._d1[i,0] = omega[0]*omega[0]*iomega_mag2*(1-cs) + cs
           #self._d1[i,1] = -omega[2]*iomega_mag*sn
           #self._d1[i,2] = omega[0]*omega[2]*iomega_mag2*(1-cs)

           #self._d2[i,0] = omega[2]*iomega_mag*sn
           #self._d2[i,1] = cs
           #self._d2[i,2] = -omega[0]*iomega_mag*sn

           #self._d3[i,0] = omega[0]*omega[2]*iomega_mag2*(1-cs)
           #self._d3[i,1] = omega[0]*iomega_mag*sn
           #self._d3[i,2] = omega[2]*omega[2]*iomega_mag2*(1-cs) + cs
        self._create_msga()

        #Determining the axis of the midline helix: First evaluate the helix at
        #a point one turn away from the starting point.
        if orient_along is None:
            return
        gamag = math.hypot(*orient_along[0:3])
        if math.isclose(gamag, 0.0):
            raise ValueError(f"orient_along = {orient_along} is a zero vector.")
        else:
            gaxis = np.asarray(orient_along[0:3], dtype=np.float64,
                               copy=True)/gamag
        
        u = 2*math.pi*iomega_mag
        p = np.array([
            -omega[0]*omega[2]*iomega_mag2*u,
            0,
            omega[2]*omega[2]*iomega_mag2*u
            ])
        pmag = np.linalg.vector_norm(p)
        if math.isclose(pmag, 0.0):
            axis = np.array([-1,0,0])
        else:
            axis = p/pmag
        self.mline[:,:] = rotlib.align(self.mline, axis, gaxis)
        self._d1[:,:] = rotlib.align(self._d1, axis, gaxis)
        self._d2[:,:] = rotlib.align(self._d2, axis, gaxis)

        shp = self.msurf.shape
        tmp = self.msurf.reshape(shp[0]*shp[1],3)
        self.msurf[:,:] = rotlib.align(tmp, axis, gaxis).reshape(shp)

        shp = self.grid.shape
        tmp = self.grid.reshape(math.prod(shp[0:-1]),3)
        self.grid[...] = rotlib.align(tmp, axis, gaxis).reshape(shp)

        if self.atom_refpos.shape[0] > 0:
            self.atom_pos[:,:] = rotlib.align(self.atom_pos, axis, gaxis)

    @staticmethod
    def _rhs(u, y, l, m):
        ydot = np.zeros_like(y)
        rotlib.normalize_quat(y[0:4])
        omega = [l(u), 0, -m(u)]

        ydot[0] = 0.5*( -omega[0]*y[1] - omega[2]*y[3] )
        ydot[1] = 0.5*(  omega[0]*y[0] - omega[2]*y[2] )
        ydot[2] = 0.5*( -omega[2]*y[1] + omega[0]*y[3] )
        ydot[3] = 0.5*(  omega[2]*y[0] - omega[0]*y[2] )

        ydot[4] = 2*(y[1]*y[3] + y[0]*y[2])
        ydot[5] = 2*(y[2]*y[3] - y[0]*y[1])
        ydot[6] = 2*(y[0]*y[0] + y[3]*y[3]) - 1
        return ydot

    def _create_ode(self, orient_along):
        def rhs(u, y, l, m):
            ydot = np.zeros_like(y)
            rotlib.normalize_quat(y[0:4])
            omega = [l(u), 0, -m(u)]

            ydot[0] = 0.5*( -omega[0]*y[1] - omega[2]*y[3] )
            ydot[1] = 0.5*(  omega[0]*y[0] - omega[2]*y[2] )
            ydot[2] = 0.5*( -omega[2]*y[1] + omega[0]*y[3] )
            ydot[3] = 0.5*(  omega[2]*y[0] - omega[0]*y[2] )

            ydot[4] = 2*(y[1]*y[3] + y[0]*y[2])
            ydot[5] = 2*(y[2]*y[3] - y[0]*y[1])
            ydot[6] = 2*(y[0]*y[0] + y[3]*y[3]) - 1
            return ydot

        l = self.l if callable(self.l) else lambda x: self.l
        m = self.m if callable(self.m) else lambda x: self.m

        y0 = np.array([1,0,0,0,0,0,0], dtype=np.float64)

        solver = ode(rhs)
        solver.set_integrator('dop853', rtol=1e-8, atol=1e-12, max_step=0.001)
        solver.set_f_params(l, m)
        solver.set_initial_value(y0, 0)
        
        for i in range(self.u.size):
            u = self.u[i]
            y = solver.integrate(u, step=False)
            if not solver.successful():
                print('Status = ', solver.get_return_code())
                raise RuntimeError('ODE solver failed.')
            rotlib.normalize_quat(y[0:4])
            self.mline[i,:] = y[4:]
            dcm = rotlib.quat_to_dcm(y[0:4])
            self._d1[i,:] = dcm[0,:]
            self._d2[i,:] = dcm[1,:]
            self._d3[i,:] = dcm[2,:]
        self._create_msga()


    def _create_msga(self):
        #Create the midsurface
        if callable(self.n):
            n = np.broadcast_to( self.n(self.u).reshape(self.u.size,1),
                                (self.u.size,3) )
        else:
            n = np.broadcast_to(self.n,(self.u.size,3))

        for j in range(self.v.size):
            v = self.v[j]
            self.msurf[:,j,:] = (self.mline + v*self._d1 
                                 - 0.5*v*v*n*self._d2)

        rbsX = RectBivariateSpline(self.u, self.v, self.msurf[:,:,0])
        rbsY = RectBivariateSpline(self.u, self.v, self.msurf[:,:,1])
        rbsZ = RectBivariateSpline(self.u, self.v, self.msurf[:,:,2])

        #Calculate the current grid
        if self.thickness == 0:
            self.grid[...] = self.msurf[...]
        else:
            self._msurf_du[:,:,0] = rbsX(self.u, self.v, dx=1, dy=0, grid=True)
            self._msurf_du[:,:,1] = rbsY(self.u, self.v, dx=1, dy=0, grid=True)
            self._msurf_du[:,:,2] = rbsZ(self.u, self.v, dx=1, dy=0, grid=True)

            self._msurf_dv[:,:,0] = rbsX(self.u, self.v, dx=0, dy=1, grid=True)
            self._msurf_dv[:,:,1] = rbsY(self.u, self.v, dx=0, dy=1, grid=True)
            self._msurf_dv[:,:,2] = rbsZ(self.u, self.v, dx=0, dy=1, grid=True)

            self._normals[:,:,:] = np.cross(self._msurf_du, self._msurf_dv, axis=2)
            norm = np.linalg.vector_norm(self._normals, axis=2, keepdims=True)
            self._normals /= norm

            for k in range(self.w.size):
                w = self.w[k]
                self.grid[:,:,k,:] = self.msurf + w*self._normals

        #Atom positions
        if self.atom_refpos.shape[0] > 0:
            self._ap_ms[:,0] = rbsX(self.atom_refpos[:,0], 
                                    self.atom_refpos[:,1], grid=False)
            self._ap_ms[:,1] = rbsY(self.atom_refpos[:,0], 
                                    self.atom_refpos[:,1], grid=False)
            self._ap_ms[:,2] = rbsZ(self.atom_refpos[:,0], 
                                    self.atom_refpos[:,1], grid=False)

            if self.thickness == 0:
                self.atom_pos[:,:] = self._ap_ms[:,:]
            else:
                self._ap_du[:,0] = rbsX(self.atom_refpos[:,0],
                                        self.atom_refpos[:,1], dx=1, dy=0,
                                        grid=False)
                self._ap_du[:,1] = rbsY(self.atom_refpos[:,0],
                                        self.atom_refpos[:,1], dx=1, dy=0,
                                        grid=False)
                self._ap_du[:,2] = rbsZ(self.atom_refpos[:,0],
                                        self.atom_refpos[:,1], dx=1, dy=0,
                                        grid=False)

                self._ap_dv[:,0] = rbsX(self.atom_refpos[:,0],
                                        self.atom_refpos[:,1], dx=0, dy=1,
                                        grid=False)
                self._ap_dv[:,1] = rbsY(self.atom_refpos[:,0], 
                                        self.atom_refpos[:,1], dx=0, dy=1,
                                        grid=False)
                self._ap_dv[:,2] = rbsZ(self.atom_refpos[:,0],
                                        self.atom_refpos[:,1], dx=0, dy=1,
                                        grid=False)

                self._ap_normals[:,:] = np.cross(self._ap_du, self._ap_dv, axis=1)
                norm = np.linalg.vector_norm(self._ap_normals, axis=1, keepdims=True)
                self._ap_normals /= norm
                tmp = np.einsum('ij,i->ij',self._ap_normals, self.atom_refpos[:,2])
                self.atom_pos[:,:] = self._ap_ms + tmp
