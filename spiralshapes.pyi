#!/usr/bin/env python

from collections.abc import Sequence
from pathlib import Path
from typing import Literal
from .spiralcartesian import SpiralCartesianBase, SpiralCornu
from .spiralpolar import SpiralPolarBase
from ._typal import Vecf


def extrude(xcoords: Sequence[float],
            ycoords: Sequence[float], 
            zdist: float,
            dz: float, 
            atom_symbol: str ='C',
            filename: Path | str ='out.xyz',
            title: str =''
            ) -> tuple[Vecf,Vecf,Vecf]: ...


def double_spiral(spiral: SpiralPolarBase | SpiralCartesianBase,
                  L: float,
                  ds: float = 0.5,
                  same: bool = False,
                  end: Literal[0, 1] = 1,
                  f: float = 0.0
                  ) -> tuple[Vecf,Vecf,Vecf,Vecf]: ...
