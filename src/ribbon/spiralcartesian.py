"""
Classes for spirals described by a Césaro equation :math:`\\kappa = f(s)`, where
:math:`\\kappa` is the curvature and :math:`s` is the arclength.

"""
from abc import ABC, abstractmethod
import math
import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad
from scipy.special import fresnel, sici

class SpiralCartesianBase(ABC):
    @abstractmethod
    def get_tangent(self, s):
        """
        Returns the unit tangent vector at the given values of polar or
        arclength coordinates.

        Parameters
        ----------
        s : array_like
            Value of the polar coordinates or arclengths.

        Returns
        -------
        ndarray
            Unit tangent vectors. Shape is *(2,)* if `s` is scalar, else *(n,2)*
            where *n* is the length of `s`.

        """
        pass

    @abstractmethod
    def get_curvature(self, s):
        """
        Returns the curvature at the given values of arclengths.

        Parameters
        ----------
        s : array_like
            Arclengths.

        Returns
        -------
        float or ndarray
            Curvatures (*float* if input is scalar, else *ndarray*).

        """
        pass

    @abstractmethod
    def get_cart_coords(self, s):
        """
        Returns the cartesian coordinates at the given values of arclengths.

        Parameters
        ----------
        s : array_like
            Arclengths.

        Returns
        -------
        x, y : tuple
            Cartesian coordinates. If `s` is a scalar, `x` and `y` are floats,
            else they are 1D *ndarray*.

        """
        pass

    def _check_bounds(self, s):
        """
        Validates that the inputs are within allowable range.

        Parameters
        ----------
        s : array_like
            Arclengths.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any input is out-of-range.

        """
        s_ = np.array(s, dtype=np.float64, copy=None, ndmin=1)
        if np.any(s_ < 0.0):
            raise ValueError(f"Arclengths must be >= 0.\n"
                            f" Arclengths:\n{s}.")


class SpiralPolynomial(SpiralCartesianBase):
    """
    The Césaro equation :math:`\\kappa = f(s)`, where :math:`f(s)` is a
    polynomial in :math:`s`.

    Warnings
    --------
    Polynomials with degree greater than 2 often result in self-intersecting, or
    at least self-touching spirals.

    Parameters
    ----------
    coeffs : array_like
        Coefficients of the polynomial in order of ascending power, i.e.
        `coeffs[i]` is the coefficient of *s^i*.

    """
    def __init__(self, coeffs):
        if np.allclose(coeffs, 0.0, rtol=1e-8, atol=1e-14):
            raise ValueError(f"`coeffs`(= {coeffs}) must not all be zero.")
        self._func_kappa = Polynomial(coeffs)
        self._func_int_kappa = self._func_kappa.integ(m=1)

    def get_tangent(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        psi = self._func_int_kappa(s_)
        tangent = np.zeros((s_.size,2), dtype=np.float64)
        tangent[:,0] = np.cos(psi)
        tangent[:,1] = np.sin(psi)
        return np.squeeze(tangent)

    def get_curvature(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        kappa = self._func_kappa(s_)
        out = kappa[0] if np.isscalar(s) else kappa
        return out

    def get_cart_coords(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        x = np.zeros_like(s_)
        y = np.zeros_like(s_)
        inds = np.argsort(s_)
        fx = lambda x: np.cos(self._func_int_kappa(x))
        fy = lambda x: np.sin(self._func_int_kappa(x))
        s0 = 0.0; x0 = 0.0; y0 = 0.0
        for i in inds:
            s_i = s_[i]
            int_x_i, err = quad(fx, s0, s_i, limit=400)
            int_y_i, err = quad(fy, s0, s_i, limit=400)
            x_i = x0 + int_x_i
            y_i = y0 + int_y_i
            x[i] = x_i; y[i] = y_i
            s0 = s_i; x0 = x_i; y0 = y_i
        out = (x[0],y[0]) if np.isscalar(s) else (x,y)
        return out


class SpiralCornu(SpiralCartesianBase):
    """
    The Césaro equation is :math:`\\kappa = as`, where :math:`a > 0`.

    Parameters
    ----------
    a : float
        Proportionality constant between curvature and arclength. Must be >
        0.

    """
    def __init__(self, a): 
        if a <= 0.0:
            raise ValueError(f"`a` must be >= 0. a = {a}.")
        self.a = a

    def get_tangent(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        psi = 0.5*self.a*s_*s_
        tangent = np.zeros((s_.size,2), dtype=np.float64)
        tangent[:,0] = np.cos(psi)
        tangent[:,1] = np.sin(psi)
        return np.squeeze(tangent)

    def get_curvature(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        kappa = self.a*s_
        out = kappa[0] if np.isscalar(s) else kappa
        return out

    def get_cart_coords(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        m = math.sqrt(math.pi/self.a)
        y, x = fresnel( s_/m )
        x *= m; y*= m
        out = (x[0],y[0]) if np.isscalar(s) else (x,y)
        return out


class SpiralNielsen(SpiralCartesianBase):
    '''
    The Césaro equation is :math:`\\kappa = a e^{bs}`, where
    :math:`a > 0, b \\neq 0`.

    Parameters
    ----------
    a : float
        Parameter in the Césaro equation. Must be > 0.
    b : float
        Parameter in the Césaro equation. Must be non-zero.

    '''
    def __init__(self, a, b):
        if a <= 0.0:
            raise ValueError(f"`a` must be >= 0. a = {a}.")
        if b == 0.0:
            raise ValueError(f"`b` must be non-zero. b = {b}.")
        self.a = a
        self.b = b

    def get_tangent(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        kappa = self.a*np.exp(self.b*s_)
        m = self.a/self.b
        psi = m*np.exp(self.b*s_) + math.pi/2
        tangent = np.zeros((s_.size,2), dtype=np.float64)
        tangent[:,0] = np.cos(psi)
        tangent[:,1] = np.sin(psi)
        return np.squeeze(tangent)

    def get_curvature(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        kappa = self.a*np.exp(self.b*s_)
        out = kappa[0] if np.isscalar(s) else kappa
        return out

    def get_cart_coords(self, s):
        self._check_bounds(s)
        s_ = np.array(s, dtype=np.float64, ndmin=1)
        m = self.a/self.b
        si0, ci0 = sici(m)
        si, ci = sici(m*np.exp(self.b*s_))
        x = -(si - si0)/self.b
        y =  (ci - ci0)/self.b
        out = (x[0],y[0]) if np.isscalar(s) else (x,y)
        return out
