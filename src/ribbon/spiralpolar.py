"""
Classes for spirals described by a polar equation *r = f(theta)*.

"""
from abc import ABC, abstractmethod
import numpy as np
from scipy.integrate import solve_ivp


class SpiralPolarBase(ABC):
    """
    Abstract base class for spirals with polar equation *r = f(theta)*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`. This
    need not be the origin of the spiral.  The arclength *s* is measured from
    (:attr:`.r0`, :attr:`.t0`) in the direction of increasing radius.

    Attributes
    ----------
    tincr : bool
        Whether radius increases with increasing *theta*.
    r0, t0 : float
        The spiral curve segment begins at (*r* = `r0`, *theta* = `t0`).

    """
    def __init__(self, tincr, r0):
        """
        Parameters
        ----------
        tincr : bool
            Whether radius increases with increasing *theta*.
        r0 : float
            The spiral curve segment begins at (*r* = `r0`). Must be a
            non-negative value.

        """
        if r0 < 0:
            raise ValueError(f"`r0`(= {r0:g}) must be >= 0.")
        self.tincr = tincr
        self.r0 = r0
        self.t0 = self.r_to_theta(r0)

    @abstractmethod
    def theta_to_r(self, theta):
        """
        Returns the radial coordinates corresponding to the given angular
        coordinates.

        Parameters
        ----------
        theta : array_like
            Angular coordinates in radians.

        Returns
        -------
        float or ndarray
            Radial coordinates (*float* if input is scalar, else *ndarray*)

        """
        pass

    @abstractmethod
    def r_to_theta(self, r):
        """
        Returns the angular coordinates corresponding to the given radial
        coordinates.

        Parameters
        ----------
        r : array_like
            Radial coordinates.

        Returns
        -------
        float or ndarray
            Angular coordinates in radians (*float* if input is scalar, else
            *ndarray*).

        """
        pass

    @abstractmethod
    def _func_arclength(self, theta):
        """
        Returns the arclengths corresponding to the given angular coordinates.

        Parameters
        ----------
        theta : array_like
            Angular coordinates in radians.

        Returns
        -------
        float or ndarray
            Arclengths at given values of the angular coordinates (*float* if
            input is scalar, else *ndarray*) measured from the starting point
            along the direction with increasing radius.

        """
        pass

    def _func_arclength_der(self, theta):
        """
        Returns the derivative of arclength with respect to the angle (theta) at
        the given values of the angular coordinates.

        Parameters
        ----------
        theta : array_like
            Angular coordinates in radians.

        Returns
        -------
        float or ndarray
            Derivative of the arclengths at given values of the angular
            coordinates (*float* if input is scalar, else *ndarray*).

        """
        pass

    @abstractmethod
    def _func_curvature(self, theta):
        """
        Returns the curvature at the given values of the angular coordinates.

        Parameters
        ----------
        theta : array_like
            Angular coordinates in radians.

        Returns
        -------
        float or ndarray
            Curvature at given values of the angular coordinates (*float* if
            input is scalar, else *ndarray*).

        """
        pass

    @abstractmethod
    def get_tangent(self, v, var='s'):
        """
        Returns the unit tangent vector at the given values of polar or
        arclength coordinates.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates or arclengths.
        var : {'r', 'theta', 's'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, 
            angular coordinates --- ``'theta'``, and arclengths --- ``'s'``.

        Returns
        -------
        ndarray
            Unit tangent vectors. Shape is *(2,)* if `v` is scalar, else *(n,2)*
            where *n* is the length of `v`.

        """
        pass

    def _check_bounds(self, v, var='s'):
        """
        Validates that the inputs are within allowable range.

        This function should be overridden by subclasses if needed.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates or arclengths.
        var : {'r', 'theta', 's'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, angular
            coordinates --- ``'theta'``, and arclengths --- ``'s'``.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any input is out-of-range.

        """
        eps = np.finfo(np.float64).eps
        if var not in ['r', 'theta', 's']:
            raise ValueError(f"var(= {var}) must be either 'r', 'theta', or 's'.")
        v_ = np.array(v, dtype=np.float64, copy=None, ndmin=1)
        if var == 's':
            is_less_than = v_ < 0.0
            if np.any(is_less_than):
                i = np.nonzero(is_less_than)[0][0]
                raise ValueError(f"Arclengths must be >= 0."
                                 f" v[{i}] = {v[i]:g}.")
        if var == 'r':
            is_close = np.isclose(v_, self.r0, rtol=1e-8, atol=1e-15)
            is_not_close = np.logical_not( is_close )
            is_less_than = np.logical_and(is_not_close, (v_ < self.r0) )
            if np.any(is_less_than):
                i = np.nonzero(is_less_than)[0][0]
                raise ValueError(
                        f"Radii must be >= `self.r0`(= {self.r0:g})."
                        f" v[{i}] = {v[i]:g}."
                        )
        if var == 'theta':
            is_close = np.isclose(v_, self.t0, rtol=1e-8, atol=1e-15)
            is_not_close = np.logical_not( is_close )
            is_less_than = np.logical_and(is_not_close, (v_ < self.t0) )
            is_greater_than = np.logical_and(is_not_close, (v_ > self.t0) )
            if self.tincr and np.any(is_less_than):
                i = np.nonzero(is_less_than)[0][0]
                raise ValueError(
                    f"Angles must be >= `self.t0`(= {self.t0:g})."
                    f" v[{i}] = {v[i]:g}."
                    )
            elif not self.tincr and np.any(is_greater_than):
                i = np.nonzero(is_greater_than)[0][0]
                raise ValueError(
                    f"Angles must be <= `self.t0`(= {self.t0:g})."
                    f" v[{i}] = {v[i]:g}."
                    )

    def arclength_to_r(self, s):
        """
        Returns the radial coordinates corresponding to the arclengths.

        Parameters
        ----------
        s : array_like
            Arclengths measured from the starting point (as set by the instance
            attribute :attr:`.r0`) along the direction of increasing radius.

        Returns
        -------
        float or ndarray
            Radial coordinates (*float* if input is scalar, else *ndarray*).

        """
        t = self.arclength_to_theta(s)
        return self.theta_to_r(t)

    def arclength_to_theta(self, s, s0=0.0, t0=None):
        """
        Returns the angular coordinates corresponding to the arclengths.

        Parameters
        ----------
        s : array_like
            Arclengths measured from the starting point (as set by the instance
            attribute :attr:`.r0`) along the direction of increasing radius.
        s0 : float
            Value of the arclength at angle `t0`.
        t0 : float | None
            Measure arclength from this value of the angular coordinate.
            If ``None``, set to :attr:`.t0`, which corresponds to the radial
            coordinate :attr:`.r0`.

        Returns
        -------
        float or ndarray
            Radial coordinates (*float* if input is scalar, else *ndarray*).

        """
        s_ = np.array(s, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(s_, 's')
        if t0 is None:
            t0 = self.t0
        theta = np.zeros_like(s_)
        us_, uinds, ucounts = np.unique(s_, return_index=True,
                                        return_counts=True, sorted=True)
        #Function yprime = d theta/ d s, needed for numerical integration to
        #obtain theta as a function of s. Initial condition: theta(s=s0) = t0.
        def func(t, y):
            if self.tincr:
                arclength_der = self._func_arclength_der(y)
            else:
                arclength_der = -self._func_arclength_der(y)
            j = np.logical_not(np.isinf(arclength_der))
            yprime = np.zeros_like(y)
            yprime[j] = 1.0/arclength_der
            return yprime
        if (us_.size == 1) and (us_[0] == s0):
            #Trivial case: theta at s0
            soly = np.array([[t0]], dtype=np.float64)
        else:
            sol = solve_ivp(func, [s0, s_.max()], [t0], method='DOP853',
                            t_eval= us_, max_step=0.1, rtol=1e-6, atol=1e-9)
            soly = sol.y
        for i in range(uinds.size):
            ui = uinds[i]; c = ucounts[i]
            theta[ui:ui+c] = soly[0,i]
        return theta[0] if np.isscalar(s) else theta

    def get_arclength(self, v, var='theta', v0=None):
        """
        Returns the arclength at the given values of polar coordinates.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates.
        var : {'r', 'theta'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, angular
            coordinates --- ``'theta'``. 
        v0 : float or None
            Origin of the arclength. If ``None``, will be set to ``self.t0`` if
            ``var = 'theta'`` or to ``self.r0`` if ``var = 'r'``.

        Returns
        -------
        float or ndarray
            Arclength (*float* if input is scalar, else *ndarray*).

        """
        if var not in ['r', 'theta']:
            raise ValueError(f"`var`(= {var}) must be 'r' or 'theta'.")
        if v0 is None:
            v0 = self.t0 if var=='theta' else self.r0
        v_ = np.array(v, dtype=np.float64, copy=None, ndmin=1)
        if var == 'r':
            t0 = self.r_to_theta(v0)
            t = self.r_to_theta(v_)
        if var == 'theta':
            self._check_bounds(v0, 'theta')
            self._check_bounds(v_, 'theta')
            t0 = v0
            t = v_
        s = np.abs(self._func_arclength(t) - self._func_arclength(t0))
        return s[0] if np.isscalar(v) else s

    def get_curvature(self, v, var='s'):
        """
        Returns the curvature at the given values of polar coordinates or
        arclengths.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates or arclengths.
        var : {'r', 'theta', 's'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, angular
            coordinates --- ``'theta'``, and arclengths --- ``'s'``.

        Returns
        -------
        float or ndarray
            Curvatures (*float* if input is scalar, else *ndarray*).

        """
        if var not in ['r', 'theta', 's']:
            raise ValueError(f"`var`(= {var}) must be either"
                             " 'r', 'theta', or 's'.")
        v_ = np.asarray([v], dtype=np.float64) if np.isscalar(v) \
                else np.asarray(v, dtype=np.float64)
        if var == 'r':
            t = self.r_to_theta(v_)
        if var == 's':
            t = self.arclength_to_theta(v_)
        if var == 'theta':
            self._check_bounds(v_, 'theta')
            t = v_
        kappa = self._func_curvature(t)
        return kappa[0] if np.isscalar(v) else kappa

    def get_cart_coords(self, v, var = 's'):
        """
        Returns the cartesian coordinates at the given values of polar
        coordinates or arclengths.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates or arclengths.
        var : {'r', theta', 's'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, angular
            coordinates --- ``'theta'``, and arclengths --- ``'s'``.

        Returns
        -------
        x, y : tuple
            Cartesian coordinates. If `v` is a scalar, `x` and `y` are floats,
            else they are 1D *ndarray*s.

        """
        r, theta = self.get_polar_coords(v, var)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return (x[0], y[0]) if np.isscalar(v) else (x, y)

    def get_polar_coords(self, v, var='s'):
        """
        Returns the polar coordinates at the given values of `theta`, `r`,
        or arclengths.

        Parameters
        ----------
        v : array_like
            Value of the polar coordinates or arclengths.
        var : {'r', 'theta', 's'}
            Descriptor for `v`. Radial coordinates --- ``'r'``, angular
            coordinates --- ``'theta'``, and arclengths --- ``'s'``.

        Returns
        -------
        r, theta : tuple
            Polar coordinates. If `v` is a scalar, `r` and `theta` are 
            floats, else they are 1D *ndarray*s.

        """
        if var not in ['r', 'theta', 's']:
            raise ValueError(f"var(= {var}) must be either "
                             " 'r', 'theta', or 's'.")
        v_ = np.asarray([v], dtype=np.float64) if np.isscalar(v) \
                else np.asarray(v, dtype=np.float64)
        if var == 'r':
            r = v_
            theta = self.r_to_theta(r)
        elif var == 'theta':
            theta = v_
            r = self.theta_to_r(theta)
        elif var == 's':
            theta = self.arclength_to_theta(v_)
            r = self.theta_to_r(theta)
        out = (r[0], theta[0]) if np.isscalar(v) else (r, theta)
        return out



class SpiralCircleInvolute(SpiralPolarBase):
    """
    Spirals with polar equation: *r = b * (1+theta^2)^(1/2)*, where *b >= 0* is
    the radius of the circle. 

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`, where
    *r0* must be > :attr:`.b`.

    The arclength *s* is measured from (:attr:`.r0`, :attr:`.t0`) in the
    direction of increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0`, *theta* must
    be >= :attr:`.t0`, and *s* must be >= 0.

    Parameters
    ----------
    b : float
        The radius of circle, must be >= 0.
    r0 : 
        Radial coordinate of the starting point of the spiral (must be > `b`).

    """
    def __init__(self, b, r0):
        if b <= 0:
            raise ValueError(f"`b`(= {b:g}) must be > 0.")
        if r0 <= b:
            raise ValueError(f"`r0`(= {r0:g}) must be > `b`(= {b:g}).")
        tincr = True
        self.b = b
        super().__init__(tincr, r0)

    def theta_to_r(self, theta):
        theta_ = np.array(theta, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(theta_, 'theta')
        r = self.b*np.sqrt(1.0 + theta_*theta_)
        return r[0] if np.isscalar(theta) else r

    def r_to_theta(self, r):
        r_ = np.array(r, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(r_, 'r')
        theta2 = (r_/self.b)**2 - 1.0
        theta = np.sqrt(theta2)
        return theta[0] if np.isscalar(r) else theta

    def _func_arclength(self, theta):
        return 0.5*self.b*theta*theta  

    def _func_arclength_der(self, theta):
        return self.b*theta

    def _func_curvature(self, theta):
        return 1.0/(self.b*theta)

    def get_tangent(self, v, var='s'):
        if var not in ['r', 'theta', 's']:
            raise ValueError(f"`var`(= {var}) must be either"
                             " 'r', 'theta', or 's'.")
        v_ = np.array(v, dtype=np.float64, copy=None, ndmin=1)
        if var == 'r':
            theta_ = self.r_to_theta(v_)
        if var == 's':
            theta_ = self.arclength_to_theta(v_)
        if var == 'theta':
            self._check_bounds(v_, 'theta')
            theta_ = v_
        ct = np.cos(theta_); st = np.sin(theta_)
        t2 = theta_*theta_
        den = np.sqrt(t2+(1.0+t2)*(1.0+t2))
        tangent = np.zeros((theta_.size,2), dtype=np.float64)
        tangent[:,0] = (theta_*ct - (1.0+t2)*st)/den
        tangent[:,1] = (theta_*st + (1.0+t2)*ct)/den
        return np.squeeze(tangent)
