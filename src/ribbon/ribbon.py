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
    lm, n : float or callable
        Curvature along length, width, and thickness directions, respectively.
        If any of the curvatures is callable, it should be of the form *y = f
        (x)*, where *x* (float) is the arclength coordinate and *y* (float) is
        the curvature.
    profile_curve : dict
        Specification of the profile curve.
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

        self._profile_curve = {'name': '', 'params': {}}
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


    def set_shape(self, shape, **kwargs):
        """
        Setter for the ribbon shape.

        Parameters
        ----------
        shape : str
            Name of the shape. Must be `'Helicoid'` | `'HelicalRibbon'` |
            `'Torus'` | `'Cylinder'` | `'Plane'` | `'General'`
        kwargs : dict
            Shape parameters for each `shape` as key-value pairs.

            HelicalRibbon
                * `'pitch'`. *float*. Must be non-zero. Positive for
                    right-handed, negative for left-handed.
                * `'radius'`. *float*. Must be non-zero.
            Helicoid
                * `'pitch'`. *float*. Must be non-zero. Positive for
                right-handed, negative for left-handed.
            Torus     
                * `'Radius'`. *float*. Must be non-zero.
                * `'radius'`. *float*. Must be non-zero.
            Cylinder     
                * `'radius'`. *float*. Must be non-zero.
            Planar       
                No parameters required.
            General      
                * `'pitch'`. *float*. Must be non-zero. Positive for right-handed,
                    negative for left-handed.
                * `'radius'`. *float*. Must be non-zero.
                * `'profile'`. *str*. Type of the profile curve.
                * `'profile_params'`. *dict*. Parameters of the profile curve.

        Returns
        -------
        None
        """
        shapes = ['HelicalRibbon', 'Helicoid', 'Torus', 'Cylinder', 'Plane',
                  'General']
        if shape not in shapes:
            raise ValueError(f"`shape`(= {shape}) must be one of:"
                             f" {', '.join(x for x in shapes)}."
                             )
        if shape == 'HelicalRibbon':
            radius = kwargs['radius']
            if not isinstance(radius, RealNumber) or radius == 0.0:
                raise ValueError(f"`radius`(= {radius}) must be a non-zero float.")
            pitch = kwargs['pitch']
            if not isinstance(pitch, RealNumber) or pitch == 0.0:
                raise ValueError(f"`pitch`(= {pitch}) must be a non-zero float.")
            delta = 2*np.pi*radius*self.width/math.sqrt(4*np.pi**2*radius**2
                                                        -self.width**2)
            if pitch < delta:
                raise ValueError("The ribbon will self overlap. Try reducing"
                                 " width, or increasin pitch, or increasing"
                                 " radius.")
            l, m = self._rp2lm(radius, pitch)
            self.lm = (l, m)
            self.n = m**2/l
            #l,m for the profile curve
            l, m = self._rp2lm(radius, -4*math.pi**2*radius**2/pitch)
            self._profile_curve = {'type': 'Helix', 'lm': (l, m)}
        elif shape == 'Helicoid':
            pitch = kwargs['pitch']
            if not isinstance(pitch, RealNumber) or pitch == 0.0:
                raise ValueError(f"`pitch`(= {pitch}) must be a non-zero float.")
            self.lm = (0.0, 2*np.pi/pitch)
            self.n = 0.0
            self._profile_curve = {'type': 'Line'}
        elif shape == 'Cylinder':
            radius = kwargs['radius']
            if not isinstance(radius, RealNumber) or radius == 0.0:
                raise ValueError(f"`radius`(= {radius}) must be a non-zero float.")
            self.lm = (1.0/radius, 0.0)
            self.n = 0.0
            self._profile_curve = {'type': 'Line'}
        elif shape == 'Torus':
            Radius = kwargs['Radius']
            if not isinstance(Radius, RealNumber) or Radius == 0.0:
                raise ValueError(f"`Radius`(= {Radius}) must be a non-zero float.")
            radius = kwargs['radius']
            if not isinstance(radius, RealNumber) or radius == 0.0:
                raise ValueError(f"`radius`(= {radius}) must be a non-zero float.")
            self.lm = (1.0/Radius, 0.0)
            self.n = 1.0/radius
            self._profile_curve = {'type': 'Circle', 'radius': radius}
        elif shape == 'Plane':
            self.lm = (0.0, 0.0)
            self.n = 0.0
            self._profile_curve = {'type': 'Line'}
        elif shape == 'General':
            mline_radius = kwargs['radius']
            if not callable(mline_radius) and \
                    not isinstance(mline_radius, RealNumber):
                raise ValueError(f"`radius`(= {mline_radius}) must be a float.")
            mline_pitch = kwargs['pitch']
            if not callable(mline_pitch) and \
                    not isinstance(mline_pitch, RealNumber):
                raise ValueError(f"`pitch`(= {mline_pitch}) must be a float.")
            self.lm = self._rp2lm(mline_radius, mline_pitch)
            ptype = kwargs['profile']
            ppar = kwargs.get('profile_params', None)
            if ptype == 'Line':
                #No parameters necessary
                self.n = 0.0
                self._profile_curve = {'type': 'Line'}
            elif ptype == 'Circle':
                radius = ppar['radius']
                if not isinstance(radius, RealNumber) or radius == 0.0:
                    raise ValueError(f"Profile parameter `radius`(= {radius})"
                                      " must be a non-zero float.")
                self.n = 1.0/radius
                self._profile_curve = {'type': 'Circle', 'radius': radius}
            elif ptype == 'Helix':
                radius = ppar['radius']
                if not isinstance(radius, RealNumber) or radius == 0.0:
                    raise ValueError(f"Profile parameter `radius`(= {radius})"
                                      " must be a non-zero float.")
                pitch = ppar['pitch']
                if not isinstance(pitch, RealNumber) or pitch == 0.0:
                    raise ValueError(f"Profile parameter `pitch`(= {pitch})"
                                      " must be a non-zero float.")
                if math.copysign(1.0, pitch)*math.copysign(1.0, mline_pitch) > 0:
                    raise ValueError(
                            f"Profile parameter `pitch`(= {pitch})"
                            f" and the shape parameter `pitch`(="
                            f" {mline_pitch}) must have opposite signs. Else,"
                            f" the distortion is too large to be physically"
                            f" meaningful for nanoplatelets."
                            )

                lm = self._rp2lm(radius, pitch)
                self.n = lm[1]**2/lm[0]
                self._profile_curve = {'type': 'Helix', 'lm': lm}
            elif ptype == 'General':
                self.n = ppar['n']
                f = ppar['f']
                self._profile_curve = {'type': 'General', 'f': f}
            else:
                raise ValueError(f"`profile`(= {ptype}) must be one of 'Line'"
                                 " 'Circle', 'Helix', or 'General'.")
        #Zero-out arrays to allow calling this method repeatedly
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
        if not callable(self.lm):
            l = self.lm[0]; m = self.lm[1]
            if l == m == 0.0:
                R = math.inf
            else:
                R = l/(l*l+m*m)
            return R
        else:
            lm = [self.lm(x) for x in self.u] #list of tuples (l,m)
            R = [math.inf if x[0]==x[1]==0.0 else x[0]/(x[0]**2+x[1]**2) for x
                 in lm]
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
        if not callable(self.lm):
            l = self.lm[0]; m = self.lm[1]
            if l == m == 0.0:
                P = math.inf
            else:
                P = 2*math.pi*m/(l*l+m*m)
            return P
        else:
            lm = [self.lm(x) for x in self.u] #list of tuples (l,m)
            P = [math.inf if x[0]==x[1]==0.0 else
                 2*math.pi*x[0]/(x[0]**2+x[1]**2) for x in lm]
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
        if not (callable(self.lm) or callable(self.n)):
            l = self.lm[0]; m = self.lm[1]
            kg = l*self.n - m**2
            return kg
        else:
            lm = [self.lm(x) for x in self.u] #list of tuples (l,m)
            kg = [x[0]*self.n-x[1]**2 for x in lm]
            return (self.u, np.array(kg, dtype=np.float64, copy=None, ndmin=1))


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
        if not callable(self.lm):
            l = self.lm[0]; m = self.lm[1]
            km = 0.5*(l+self.n)
            return km
        else:
            lm = [self.lm(x) for x in self.u] #list of tuples (l,m)
            km = [0.5*(x[0]+self.n) for x in lm]
            return (self.u, np.array(km, dtype=np.float64, copy=None, ndmin=1))


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
        if not callable(self.lm):
            l = self.lm[0]; m = self.lm[1]
            theta = 0.5*math.atan(2*m/(l-self.n)) \
                    if l != self.n else math.pi/4
            return theta
        else:
            lm = [self.lm(x) for x in self.u] #list of tuples (l,m)
            theta = [0.5*math.atan(2*x[1]/(x[0]-self.n)) if x[0]!=self.n else
                     math.pi/4 for x in lm]
            return (self.u, np.array(theta, dtype=np.float64, copy=None, ndmin=1))

    def translate_to_center(self):
        """
        Translates grids and current atom positions so as to bring the center
        of the ribbon midline at (0,0,0).

        Returns
        -------
        None

        """
        center = self.mline.mean(axis=0)
        self.mline -= center
        self.msurf -= center
        self.grid -= center
        self.atom_pos -= center


    def create(self, orient_along=None):
        """
        Constructs the ribbon and orients it along `orient_along`.

        Parameters
        ----------
        orient_along : (3,) array_like or None
            If the curvatures of the ribbon midline are constant, the ribbon
            axis will be oriented along `orient_along`. This does not need to
            be a unit vector. If `None`, the orientation will be remain as is.

        Returns
        -------
        None

        """
        if orient_along is None:
            axis = None
        else:
            mag = math.hypot(*orient_along[0:3])
            if math.isclose(mag, 0.0):
                raise ValueError(f"orient_along = {orient_along} is a zero vector.")
            else:
                axis = np.asarray(orient_along[0:3], dtype=np.float64)/mag
        self._create_msurf(axis)
        self._interpolate()


    def _create_msurf(self, axis):
        """
        Creates the mid-surface.
        """
        #Create the midline
        mline_axis = np.zeros((3,), dtype=np.float64)
        if not callable(self.lm):
            if self.lm[0]==self.lm[1]==0.0:
                self.mline[...] = 0.0
                self._d1[...] = 0.0
                self._d2[...] = 0.0
                self._d3[...] = 0.0
                self.mline[:,2] = self.u
                self._d1[:,0] = 1.0
                self._d2[:,1] = 1.0
                self._d3[:,2] = 1.0
                mline_axis[2] = 1.0
            else:
                self.create_helix(self.u, 'lm', (self.lm,), mline_axis,
                                  self.mline, d1=self._d1, d2=self._d2,
                                  d3=self._d3)
        else:
            self.create_helix_var(self.u, 'lm', (self.lm,), mline_axis,
                                   self.mline, d1=self._d1, d2=self._d2,
                                   d3=self._d3)
        #Create the midsurface
        vline = np.zeros((self.v.size,3), dtype=np.float64)
        #2-D profile curve lies on the x-y plane. It is y = f(x), where x
        #spans the d1 direction, and -y the d2 direction. Needs better
        #explanation here.
        if self._profile_curve['type'] == 'Line':
            vline[:,0] = self.v
        elif self._profile_curve['type'] == 'Circle':
            #Note: Positive radius to move away from d2, negative to move
            #toward d2. 
            r = self._profile_curve['radius']
            vline[:,0] = -r*np.sin(self.v/r)
            vline[:,1] = -r*(1.0 - np.cos(self.v/r))
        elif self._profile_curve['type'] == 'Helix':
            lm = self._profile_curve['lm']
            vline_axis = np.zeros((3,), dtype=np.float64)
            self.create_helix(self.v, 'lm', (lm,), vline_axis,  vline)
            vline = rotlib.align(vline, vline_axis, mline_axis )
            #vline_angle = math.acos(vline_axis[0])
            #mline_angle = math.acos(mline_axis[0])
            #rotate_by = vline_angle - (math.pi/2 + mline_angle)
            #print(f"angv = {math.degrees(vline_angle)}")
            #print(f"angm = {math.degrees(mline_angle)}")
            #print(f"angr = {math.degrees(rotate_by)}")
            #vline = rotlib.aa_rotate_vectors(vline, np.array([0.0, 1.0, 0.0]),
            #                         -rotate_by)
            #print('angle = ', np.degrees(np.acos(np.dot(vline_axis, mline_axis))) )
        elif self._profile_curve['type'] == 'General':
            vline_coords = [self._profile_curve['f'](x) for x in self.v]
            if len(vline_coords[0]) == 2:
                vline[:,:2] = np.asarray(vline_coords, dtype=np.float64)
            elif len(vline_coords[0]) == 3:
                vline[:,:] = np.asarray(vline_coords, dtype=np.float64)
        # Move the profile curve to midline frames to create the surface
        dcm = np.zeros((3,3))
        for i in range(self.u.size):
            dcm[0,:] = self._d1[i,:]
            dcm[1,:] = self._d2[i,:]
            dcm[2,:] = self._d3[i,:]
            vline_shifted = rotlib.dcm_rotate_vectors(vline, dcm )
            self.msurf[i,:,:] = vline_shifted + self.mline[i,:]

        if axis is not None:
            self.mline = rotlib.align(self.mline, mline_axis, axis )
            self._d1 = rotlib.align(self._d1, mline_axis, axis)
            self._d2 = rotlib.align(self._d2, mline_axis, axis)
            self._d3 = rotlib.align(self._d3, mline_axis, axis)
            self.msurf = rotlib.align(self.msurf, mline_axis, axis )
        #End of create midsurface


    def _interpolate(self):
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


    @staticmethod
    def create_helix(s, parametrization, params, axis, coords, d1=None,
                      d2=None, d3=None):
        """
        Creates a helix with constant curvatures. 

        Parameters
        ----------
        s : 1d array_like
            Arclength.
        parametrization : `'rp'` | `'lm'`
            If `'rp'`, helix parameters are radius and pitch. If `'lm'`, the
            parameters are the curvatures `l` and `m` (twist). Use `rp` for
            common use-cases.
        params : tuple
            If `parametrization` = `'rp'`, `params` = (radius, pitch).
            If `parametrization` = `'lm'`, `params` = (lm, ), where lm is a
            tuple (l,m). 
        axis : (3,) ndarray
            On return, the contains the unit vector along the helix axis.
        coords : (n, 3) ndarray
            On return, contains the cartesian coordinates of the helix. `n`
            must equal the length of `s`.
        d1, d2, d3 : (n, 3) ndarray | None
            On return, contains the d1 ( equivalently d2 and d3) director
            vector of the helix. `n` must equal the length of `s`. If None, the
            corresponding director is not returned.

        Returns
        -------
        None

        """
        if axis.shape != (3,):
            raise ValueError(f"axis must have shape (3,).")
        if coords.shape != (len(s),3):
            raise ValueError(f"coords must have shape ({len(s)},3).")
        if d1 is not None:
            if d1.shape != (len(s),3):
                raise ValueError(f"d1 must have shape ({len(s)},3).")
        if d2 is not None:
            if d2.shape != (len(s),3):
                raise ValueError(f"d2 must have shape ({len(s)},3).")
        if d3 is not None:
            if d3.shape != (len(s),3):
                raise ValueError(f"d3 must have shape ({len(s)},3).")
        if parametrization == 'rp':
            l,m = self._rp2lm(params[0], params[1])
        elif parametrization == 'lm':
            l,m = params[0]
        else:
            raise ValueError(f"parametrization(= {parametrization} must be"
                             f" 'rp' or 'lm'.")
        omega = (l, 0.0, -m)
        omega_mag2 = omega[0]*omega[0] + omega[2]*omega[2]
        iomega_mag2 = 1/omega_mag2
        omega_mag = math.sqrt(omega_mag2)
        iomega_mag = 1/omega_mag
        for i in range(len(s)):
            u = s[i]
            sn = math.sin(omega_mag*u)
            cs = math.cos(omega_mag*u)
            coords[i,0] = -omega[0]*omega[2]*iomega_mag2*(u - sn*iomega_mag)
            coords[i,1] = -omega[0]*iomega_mag2*(1 - cs)
            coords[i,2] = omega[2]*omega[2]*iomega_mag2 * (u - sn*iomega_mag)\
                            + sn*iomega_mag
            if d1 is not None:
                d1[i,0] =  omega[0]*omega[0]*iomega_mag2*(1-cs) + cs
                d1[i,1] = -omega[2]*iomega_mag*sn
                d1[i,2] = -omega[0]*omega[2]*iomega_mag2*(1-cs)
            if d2 is not None:
                d2[i,0] = omega[2]*iomega_mag*sn
                d2[i,1] = cs
                d2[i,2] = omega[0]*iomega_mag*sn
            if d3 is not None:
                d3[i,0] = -omega[0]*omega[2]*iomega_mag2*(1-cs)
                d3[i,1] = -omega[0]*iomega_mag*sn
                d3[i,2] =  omega[2]*omega[2]*iomega_mag2*(1-cs) + cs

        #Current axis
        if axis is not None:
            u = 2*math.pi*iomega_mag
            p = np.array([
                -omega[0]*omega[2]*iomega_mag2*u,
                0,
                omega[2]*omega[2]*iomega_mag2*u
                ]) # fmt: skip
            pmag = np.linalg.vector_norm(p)
            if math.isclose(pmag, 0.0):
                axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0
            else:
                axis[:] = p/pmag
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

    @staticmethod
    def create_helix_var(s, parametrization, params, axis, coords, d1=None,
                          d2=None, d3=None):
        """
        Creates a helix with non-constant curvatures. 

        Parameters
        ----------
        s : 1d array_like
            Arclength.
        parametrization : `'rp'` | `'lm'`
            If `'rp'`, helix parameters are radius and pitch. If `'lm'`, the
            parameters are the curvatures `l` and `m` (twist). Use `rp` for
            common use-cases.
        params : tuple
            If `parametrization` = `'rp'`, `params` = (radius, pitch), where
            at least one of radius and pitch must be a callable.
            If `parametrization` = `'lm'`, `params` = (lm, ), where lm is a
            callble returning a tuple (l,m).
        axis : (3,) ndarray
            On return, the contains the unit vector along the helix axis.
        coords : (n, 3) ndarray
            On return, contains the cartesian coordinates of the helix. `n`
            must equal the length of `s`.
        d1, d2, d3 : (n, 3) ndarray | None
            On return, contains the d1 ( equivalently d2 and d3) director
            vector of the helix. `n` must equal the length of `s`. If None, the
            corresponding director is not returned.

        Returns
        -------
        None

        """
        if axis.shape != (3,):
            raise ValueError(f"axis must have shape (3,).")
        if coords.shape != (len(s),3):
            raise ValueError(f"coords must have shape ({len(s)},3).")
        if d1 is not None:
            if d1.shape != (len(s),3):
                raise ValueError(f"d1 must have shape ({len(s)},3).")
        if d2 is not None:
            if d2.shape != (len(s),3):
                raise ValueError(f"d2 must have shape ({len(s)},3).")
        if d3 is not None:
            if d3.shape != (len(s),3):
                raise ValueError(f"d3 must have shape ({len(s)},3).")
        if parametrization == 'rp':
            lm = self._rp2lm(params[0], params[1])
        elif parametrization == 'lm':
            lm = params[0]
        else:
            raise ValueError(f"parametrization(= {parametrization} must be"
                             f" 'rp' or 'lm'.")

        def rhs(u, y, lm):
            ydot = np.zeros_like(y)
            rotlib.quat_normalized(y[0:4])
            l, m = lm(u)
            omega = [l, 0, -m]

            ydot[0] = 0.5*( -omega[0]*y[1] - omega[2]*y[3] )
            ydot[1] = 0.5*(  omega[0]*y[0] - omega[2]*y[2] )
            ydot[2] = 0.5*( -omega[2]*y[1] + omega[0]*y[3] )
            ydot[3] = 0.5*(  omega[2]*y[0] - omega[0]*y[2] )

            ydot[4] = 2*(y[1]*y[3] + y[0]*y[2])
            ydot[5] = 2*(y[2]*y[3] - y[0]*y[1])
            ydot[6] = 2*(y[0]*y[0] + y[3]*y[3]) - 1
            return ydot

        y0 = np.array([1,0,0,0,0,0,0], dtype=np.float64)

        solver = ode(rhs)
        solver.set_integrator('dop853', rtol=1e-8, atol=1e-12, max_step=0.001)
        solver.set_f_params(lm)
        solver.set_initial_value(y0, 0)
        
        for i in range(len(s)):
            u = s[i]
            y = solver.integrate(u, step=False)
            if not solver.successful():
                print('Status = ', solver.get_return_code())
                raise RuntimeError('ODE solver failed.')
            rotlib.quat_normalized(y[0:4])
            coords[i,:] = y[4:]
            if not ( d1 is None and d2 is None and d3 is None ):
                dcm = rotlib.quat_to_dcm(y[0:4])
            if d1 is not None:
                d1[i,:] = dcm[0,:]
            if d2 is not None:
                d2[i,:] = dcm[1,:]
            if d3 is not None:
                d3[i,:] = dcm[2,:]

    @staticmethod
    def _rp2lm(radius, pitch):
        """
        Returns the curvatures l and m from radius and pitch.

        Parameters
        ----------
        radius : float or callable
        pitch : float or callable

        `radius` and `pitch` both cannot be equal to zero. For a circle, radius
        > 0, pitch = 0. For a straight line, the limiting case is radius = 0,
        pitch = inf (corresponding to l = m = 0).

        Returns
        -------
        lm : float or callable
            Will be callable if either of `radius` or `pitch` is callable.
        """
        if not (callable(radius) or callable(pitch)):
            if radius == pitch == 0.0:
                raise ValueError(f"Both radius and pitch cannot be zero.")
            elif math.isinf(pitch):
                lm = (0.0, 0.0)
            else:
                tpi = 2.0*np.pi
                fpi2 = tpi*tpi
                den = fpi2*radius*radius + pitch*pitch
                lm = (fpi2*radius/den, tpi*pitch/den)
        else:
            r = radius if callable(radius) else lambda x: radius
            p = pitch if callable(pitch) else lambda x: pitch
            def lm(x, r=r, p=p):
                tpi = 2.0*np.pi
                fpi2 = tpi*tpi
                rx = r(x)
                px = p(x)
                if rx == px == 0.0:
                    raise ValueError(f"Both radius and pitch cannot be zero."
                                     f"Radius = pitch = 0.0 at {x}.")
                elif math.isinf(px):
                    out = (0.0, 0.0)
                else:
                    den = fpi2*rx*rx + px*px
                    out = (fpi2*rx/den, tpi*px/den)
                return out
        return lm
