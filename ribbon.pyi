#!/usr/bin/env python

from collections.abc import Callable, Sequence
from _typal import Vecf


class Ribbon(object):
    def __init__(self,
                 length: float,
                 width: float,
                 thickness: float,
                 dl: float,
                 dw: float,
                 dt: float): ...

    def set_curvatures(self,
                       l: float | Callable[[float],float],
                       m: float | Callable[[float],float],
                       n: float | Callable[[float],float]
                       ) -> None: ...

    def get_radius(self) -> float | tuple[Vecf, Vecf]: ...

    def get_pitch(self) -> float | tuple[Vecf, Vecf]: ...

    def get_gauss_curvature(self) -> float | tuple[Vecf, Vecf]: ...

    def get_mean_curvature(self) -> float | tuple[Vecf, Vecf]: ...

    def get_theta(self) -> float | tuple[Vecf, Vecf]: ...

    def create(self, orient_along: Sequence[float] = [1,0,0]) -> None: ...

    def _create_direct_zero(self) -> None: ...

    def _create_direct(self, orient_along: Sequence[float]) -> None: ...

    @staticmethod
    def _rhs(u: float,
             y: Vecf,
             l: Callable[[float],float],
             m: Callable[[float],float]
             ) -> Vecf: ...

    def _create_ode(self, orient_along: Sequence[float]) -> None: ...

    def _create_msga(self) -> None: ...
