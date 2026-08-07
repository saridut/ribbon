#!/usr/bin/env python

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import overload, Literal
import numpy as np
from ._typal import Vecf, Array2f

class SpiralPolarBase(ABC):

    def __init__(self, tincr: bool, r0: float): ...

    @overload
    @abstractmethod
    def theta_to_r(self, theta: float) -> float: ...
    @overload
    @abstractmethod
    def theta_to_r(self, theta: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    @abstractmethod
    def r_to_theta(self, r: float) -> float: ...
    @overload
    @abstractmethod
    def r_to_theta(self, r: Sequence[float] | Vecf) -> Vecf: ...

    @overload 
    @abstractmethod
    def _func_arclength(self, theta: float) -> float: ...
    @overload 
    @abstractmethod
    def _func_arclength(self, theta: Vecf) -> Vecf: ...

    @overload
    def _func_arclength_der(self, theta: float) -> float: ...
    @overload
    def _func_arclength_der(self, theta: Vecf) -> Vecf: ...

    @overload
    @abstractmethod
    def _func_curvature(self, theta: float) -> float: ...
    @overload
    @abstractmethod
    def _func_curvature(self, theta: Vecf) -> Vecf: ...

    @abstractmethod
    def get_tangent(self,
                    v: float | Sequence[float] | Vecf,
                    var: Literal['r', 'theta', 's'] = 's'
                    ) -> Vecf | Array2f: ...

    def _check_bounds(self,
                     v: float | Sequence[float] | Vecf,
                     var: Literal['r', 'theta', 's'] = 's'
                     ) -> None: ...

    @overload
    def arclength_to_r(self,
                       s: float
                       ) -> float: ...
    @overload
    def arclength_to_r(self,
                       s: Sequence[float] | Vecf
                       ) -> Vecf: ...

    @overload
    def arclength_to_theta(self,
                           s: float,
                           s0: float = 0.0,
                           t0: float | None = None,
                           ) -> float: ...
    @overload
    def arclength_to_theta(self,
                           s:  Sequence[float] | Vecf,
                           s0: float = 0.0,
                           t0: float | None = None,
                           ) -> Vecf: ...

    @overload
    def get_arclength(self,
                      v: float,
                      var: Literal['r', 'theta'] = 'theta',
                      v0: float | None = None,
                      ) -> float: ...
    @overload
    def get_arclength(self,
                      v: Sequence[float] | Vecf,
                      var: Literal['r', 'theta'] = 'theta',
                      v0: float | None = None,
                      ) -> Vecf: ...

    @overload
    def get_curvature(self,
                      v: float,
                      var: Literal['r', 'theta', 's'] = 's'
                      ) -> float: ...
    @overload
    def get_curvature(self,
                      v: Sequence[float] | Vecf,
                      var: Literal['r', 'theta', 's'] = 's'
                      ) -> Vecf: ...

    @overload
    def get_cart_coords(self,
                        v: float,
                        var: Literal['r', 'theta', 's'] = 's'
                        ) -> tuple[float,float]: ...
    @overload
    def get_cart_coords(self,
                        v: Sequence[float] | Vecf,
                        var: Literal['r', 'theta', 's'] = 's'
                        ) -> tuple[Vecf,Vecf]: ...

    @overload
    def get_polar_coords(self,
                        v: float,
                        var: Literal['r', 'theta', 's'] = 's'
                        ) -> tuple[float,float]: ...
    @overload
    def get_polar_coords(self,
                        v: Sequence[float] | Vecf,
                        var: Literal['r', 'theta', 's'] = 's'
                        ) -> tuple[Vecf,Vecf]: ...


class SpiralCircleInvolute(SpiralPolarBase):
    def __init__(self, b: float, r0: float): ...

    @overload
    def theta_to_r(self, theta: float) -> float: ...
    @overload
    def theta_to_r(self, theta: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    def r_to_theta(self, r: float) -> float: ...
    @overload
    def r_to_theta(self, r: Sequence[float] | Vecf) -> Vecf: ...

    @overload 
    def _func_arclength(self, theta: float) -> float: ...
    @overload 
    def _func_arclength(self, theta: Vecf) -> Vecf: ...

    @overload
    def _func_arclength_der(self, theta: float) -> float: ...
    @overload
    def _func_arclength_der(self, theta: Vecf) -> Vecf: ...

    @overload
    def _func_curvature(self, theta: float) -> float: ...
    @overload
    def _func_curvature(self, theta: Vecf) -> Vecf: ...

    def get_tangent(self,
                    v: float | Sequence[float] | Vecf,
                    var: Literal['r', 'theta', 's'] = 's'
                    ) -> Vecf | Array2f: ...
