#!/usr/bin/env python

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import overload
from ._typal import Vecf, Array2f

class SpiralCartesianBase(ABC):
    @abstractmethod
    def get_tangent(self,
                    s: float | Sequence[float] | Vecf
                    ) -> Vecf | Array2f: ...

    @overload
    @abstractmethod
    def get_curvature(self, s: float) -> float: ...
    @overload
    @abstractmethod
    def get_curvature(self, s: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    @abstractmethod
    def get_cart_coords(self, s: float) -> tuple[float,float]: ...
    @overload
    @abstractmethod
    def get_cart_coords(self,
                        s: Sequence[float] | Vecf
                        ) -> tuple[Vecf,Vecf]: ...

    def _check_bounds(self, s: float | Sequence[float] | Vecf) -> None: ...


class SpiralPolynomial(SpiralCartesianBase):
    def __init__(self, coeffs: Sequence[float]) -> None: ...

    def get_tangent(self,
                    s: float | Sequence[float] | Vecf
                    ) -> Vecf | Array2f: ...

    @overload
    def get_curvature(self, s: float) -> float: ...
    @overload
    def get_curvature(self, s: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    def get_cart_coords(self, s: float) -> tuple[float,float]: ...
    @overload
    def get_cart_coords(self,
                        s: Sequence[float] | Vecf
                        ) -> tuple[Vecf,Vecf]: ...


class SpiralCornu(SpiralCartesianBase):
    def __init__(self, a: float) -> None: ...

    def get_tangent(self,
                    s: float | Sequence[float] | Vecf
                    ) -> Vecf | Array2f: ...

    @overload
    def get_curvature(self, s: float) -> float: ...
    @overload
    def get_curvature(self, s: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    def get_cart_coords(self,
                        s: float
                        ) -> tuple[float,float]: ...
    @overload
    def get_cart_coords(self,
                        s: Sequence[float] | Vecf
                        ) -> tuple[Vecf,Vecf]: ...


class SpiralNielsen(SpiralCartesianBase):
    def __init__(self, a: float, b: float) -> None: ...

    def get_tangent(self,
                    s: float | Sequence[float] | Vecf
                    ) -> Vecf | Array2f: ...

    @overload
    def get_curvature(self, s: float) -> float: ...
    @overload
    def get_curvature(self, s: Sequence[float] | Vecf) -> Vecf: ...

    @overload
    def get_cart_coords(self,
                        s: float
                        ) -> tuple[float,float]: ...
    @overload
    def get_cart_coords(self,
                        s: Sequence[float] | Vecf
                        ) -> tuple[Vecf,Vecf]: ...
