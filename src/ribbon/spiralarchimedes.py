"""
Classes for spirals described by the polar equation *r = b * theta^n*.

"""
from abc import ABCMeta
import numpy as np
from scipy.special import hyp2f1
from .spiralpolar import SpiralPolarBase

class SpiralArchimedesBase(SpiralPolarBase, metaclass=ABCMeta):
    """
    Base class for archimedean spirals. Polar equation: *r = b * theta^n*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`. The
    arclength *s* is measured from (:attr:`.r0`, :attr:`.t0`) in the direction
    of increasing radius.

    Attributes
    ----------
    tincr : bool
        Whether radius increases with increasing *theta*.
    r0, t0 : float
        The spiral curve segment begins at (*r* = `r0`, *theta* = `t0`). `r0`
        may not be negative.
    b : float
        Parameter in the polar equation. Must be > 0.
    n : float
        Exponent in the polar equation. Must be non-zero. If negative, then
        :attr:`.r0` must be > 0.

    Parameters
    ----------
    b : float
        Parameter in the polar equation. Must be > 0.
    n : float
        Exponent in the polar equation. Must be non-zero. If negative, then
        `r0` must be > 0.
    r0 : float
        The spiral curve segment begins at (*r* = `r0`). Must be a
        non-negative value.

    """

    def __init__(self, b, n, r0):
        if b <= 0:
            raise ValueError(f"`b`(= {b:g}) must be > 0.")
        if n == 0:
            raise ValueError(f"`n`(= {n:g}) must be non-zero.")
        elif n < 0:
            if r0 == 0.0:
                raise ValueError(f"`r0`(= {r0:g}) must be > 0.")
        tincr = True if n > 0 else False
        self.b = b; self.n = n
        super().__init__(tincr, r0)

    def theta_to_r(self, theta):
        theta_ = np.array(theta, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(theta_, 'theta')
        r = np.where(theta_==0, 0.0, self.b*theta_**self.n)
        return r[0] if np.isscalar(theta) else r

    def r_to_theta(self, r):
        r_ = np.array(r, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(r_, 'r') 
        if self.n > 0:
            theta = np.where(r_==0.0, 0.0, (r_/self.b)**(1.0/self.n))
        else:
            theta = (r_/self.b)**(1.0/self.n)
        return theta[0] if np.isscalar(r) else theta

    def arclength_to_theta(self, s, s0=0.0, t0=None):
        s_ = np.array(s, dtype=np.float64, copy=None, ndmin=1)
        self._check_bounds(s_) 
        if (self.r0 == 0) and ( (0 < self.n < 1) or (self.n > 1) ):
            theta = np.zeros_like(s_)
            inz = np.nonzero(s_)
            t0 = 1e-9
            s0 = self.get_arclength(t0, 'theta')
            assert s_[inz].min() >= s0
            snz = super().arclength_to_theta(s_[inz], s0, t0)
            theta[inz] = snz
            out = theta
        else:
            out = super().arclength_to_theta(s_)
        return out[0] if np.isscalar(s) else out

    def get_tangent(self, v, var='s'):
        if var not in ['r', 'theta', 's']:
            raise ValueError(f"`var`(= {var}) must be either 'r', 'theta', or 's'.")
        v_ = np.array(v, dtype=np.float64, copy=None, ndmin=1)
        if var == 'r':
            theta_ = self.r_to_theta(v_)
        if var == 's':
            theta_ = self.arclength_to_theta(v_)
        if var == 'theta':
            self._check_bounds(v_, 'theta')
            theta_ = v_
        ct = np.cos(theta_); st = np.sin(theta_)
        den = np.sqrt(self.n*self.n+theta_*theta_)
        tangent = np.zeros((theta_.size,2), dtype=np.float64)
        tangent[:,0] = (self.n*ct - theta_*st)/den
        tangent[:,1] = (self.n*st + theta_*ct)/den
        if not self.tincr:
            tangent[:,:] = -tangent
        return np.squeeze(tangent)


class SpiralArchimedesGeneral(SpiralArchimedesBase):
    """
    Class implementing a general archimedean spiral with polar equation 
    *r = b * theta^n*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`.  If
    :attr:`.n` > 0, :attr:`.r0` >= 0, otherwise :attr:`.r0` > 0. The arclength
    *s* is measured from (:attr:`.r0`, :attr:`.t0`) in the direction of
    increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0` and *s* must be >=
    0.  If :attr:`.n` > 0, *theta* must be >= :attr:`.t0`, otherwise *theta* <=
    :attr:`.t0`.

    """
    def __init__(self, b, n, r0):
        super().__init__(b, n, r0)

    def _func_arclength(self, theta):
        n = self.n
        tn2 = theta*theta/(n*n)
        if hasattr(theta, '__iter__'):
            tn = np.asarray([x**n if (x != 0) else 0 for x in theta],
                            dtype=np.float64)
        else:
            #theta is a scalar
            tn = theta**n if (theta != 0.0) else 0.0
        num = np.abs(n)*tn*hyp2f1(-0.5, n/2, 1+n/2, -tn2)
        return self.b*num/n

    def _func_arclength_der(self, theta):
        n = self.n
        theta_ = np.array(theta, dtype=np.float64, copy=None, ndmin=1)
        res = np.zeros_like(theta_)
        for i in range(res.size):
            if (theta_[i]==0) and (n==1):
                res[i] = 1.0
            elif (theta_[i]==0) and (0 < n < 1):
                res[i] = np.inf
            elif (theta_[i]==0) and (n > 1):
                res[i] = 0.0
            else:
                t2 = theta_[i]*theta_[i]
                t2n = t2**n
                res[i] = np.sqrt(t2n + n*n*t2n/t2)
        out = self.b*res[0] if np.isscalar(theta) else self.b*res
        return out

    def _func_curvature(self, theta):
        n = self.n; b = self.b
        t2 = theta*theta
        den = (t2+n*n)**1.5
        num = (t2+n*n+n)
        if hasattr(theta, '__iter__'):
            pf = np.asarray([x**(1-n) if (x != 0) else 0 for x in theta],
                            dtype=np.float64)
        else:
            #theta is a scalar
            pf = theta**(1-n)  
        return (pf/b)*(num/den)


class SpiralArchimedes(SpiralArchimedesBase):
    '''
    Class implementing an Archimedean spiral with polar equation 
    *r = b * theta*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`, where
    :attr:`.r0` >= 0. The arclength *s* is measured from (:attr:`.r0`,
    :attr:`.t0`) in the direction of increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0`, *s* must be >= 0,
    and *theta* must be >= :attr:`.t0`.

    '''
    def __init__(self, b, r0):
        super().__init__(b, 1, r0)

    def _func_arclength(self, theta):
        u = np.sqrt(1+theta*theta)
        return 0.5*self.b*( theta*u + np.log(theta+u) )

    def _func_arclength_der(self, theta):
        return self.b*np.sqrt(1+theta*theta)

    def _func_curvature(self, theta):
        tsq = theta*theta
        den = self.b*(tsq + 1.0)**1.5
        return (tsq+2)/den


class SpiralFermat(SpiralArchimedesBase):
    '''
    Class implementing a general archimedean spiral with polar equation 
    *r = b * theta^(1/2)*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`, where
    :attr:`.r0` >= 0. The arclength *s* is measured from (:attr:`.r0`,
    :attr:`.t0`) in the direction of increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0`, *s* must be >= 0,
    and *theta* must be >= :attr:`.t0`.

    Parameters
    ----------
    b : float
        Parameter in the polar equation. Must be > 0.
    r0 : float
        The spiral curve segment begins at (*r* = `r0`). Must be a
        non-negative value.
    '''
    def __init__(self, b, r0):
        super().__init__(b, 1/2, r0)

    def _func_arclength(self, theta):
        return self.b*np.sqrt(theta)*hyp2f1(-0.5, 0.25, 1.25, -4*theta*theta)

    def _func_arclength_der(self, theta):
        theta_ = np.array(theta, dtype=np.float64, copy=None, ndmin=1)
        res = np.full_like(theta_, np.inf)
        inds = np.nonzero(theta_)[0]
        res[inds] = 0.5*self.b*np.sqrt(4*theta_[inds]+1.0/theta_[inds])
        return res[0] if np.isscalar(theta) else res

    def _func_curvature(self, theta):
        r_ = self.theta_to_r(theta)
        r4 = r_**4
        return 2*r_*(4*r4 + 3*self.b**4)/(4*r4 + self.b**4)**1.5


class SpiralHyperbolic(SpiralArchimedesBase):
    '''
    Class implementing a general archimedean spiral with polar equation 
    *r = b * theta^(-1)*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`, where
    :attr:`.r0` > 0. The arclength *s* is measured from (:attr:`.r0`,
    :attr:`.t0`) in the direction of increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0`, *s* must be >= 0, 
    and *theta* must be <= :attr:`.t0`.

    Parameters
    ----------
    b : float
        Parameter in the polar equation. Must be > 0.
    r0 : float
        The spiral curve segment begins at (*r* = `r0`). Must be a
        non-negative value.
    '''
    def __init__(self, b, r0):
        super().__init__(b, -1, r0)

    def _func_arclength(self, theta):
        t2 = theta*theta
        srt = np.sqrt(1+t2)
        return self.b*(-srt/theta + np.log(theta+srt))

    def _func_arclength_der(self, theta):
        t2 = theta*theta
        srt = np.sqrt(1+t2)
        return self.b*srt/t2

    def _func_curvature(self, theta):
        t2 = theta*theta
        t4 = t2*t2
        return t4/(self.b*(t2+1)**1.5)


class SpiralLituus(SpiralArchimedesBase):
    '''
    Class implementing a general archimedean spiral with polar equation 
    *r = b * theta^(-1/2)*.

    The starting point is at *r* = :attr:`.r0` and *theta* = :attr:`.t0`, where
    :attr:`.r0` > 0. The arclength *s* is measured from (:attr:`.r0`,
    :attr:`.t0`) in the direction of increasing radius.

    For all points on this spiral *r* must be >= :attr:`.r0`, *s* must be >= 0, 
    and *theta* must be <= :attr:`.t0`.

    Parameters
    ----------
    b : float
        Parameter in the polar equation. Must be > 0.
    r0 : float
        The spiral curve segment begins at (*r* = `r0`). Must be a
        non-negative value.
    '''
    def __init__(self, b, r0):
        super().__init__(b, -1/2, r0)

    def _func_arclength(self, theta):
        t2 = theta*theta
        return 2*np.sqrt(theta)*hyp2f1(-0.5, -0.25, 0.75, -0.25/t2)

    def _func_arclength_der(self, theta):
        itheta = 1.0/theta
        return 0.5*itheta*np.sqrt(itheta+4*theta)

    def _func_curvature(self, theta):
        t2 = theta*theta
        return (8*t2-2)*(theta/(1+4*t2))**1.5
