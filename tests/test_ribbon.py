#!/usr/bin/env python

from pathlib import Path
import sys
import numpy as np
import pytest
try:
    from ribbon import Ribbon
except:
    srcdir = Path('src').resolve()
    sys.path.insert(0, str(srcdir))
    from ribbon import Ribbon

init_params = [
        (0, 4, 0, 1, 1, 1),
        (2, 4, 0, 1, 1, 1),
        (20, 4, 0, None, 1, 1, {'ngpl': None}),
        (20, 0, 0, 1, 1, 1),
        (20, 2, 0, 1, 1, 1),
        (20, 4, 0, 1, None, 1, {'ngpw': None}),
        (20, 4, -1, 1, 1, 1),
        (20, 4, 1, 1, 1, None, {'ngpt': 1}),
        (20, 4, 1, 1, 1, None, {'ngpt': None}),
        ]

@pytest.fixture(scope='module', params=init_params)
def init_args(request):
    par = request.param
    if isinstance(par[-1], dict):
        args = par[:-1]; kwargs = par[-1]
    else:
        args = par; kwargs = {}
    return args, kwargs


def test_init_error(init_args):
    with pytest.raises(ValueError) as info:
        args = init_args[0]; kwargs = init_args[1]
        ribbon = Ribbon(*args, **kwargs)
    assert info.type is ValueError


def test_radius_pitch():
    l = 0.4; m = 0.4; n = 0.4
    ribbon = Ribbon(20, 4, 0, 0.1, 0.1, 0)
    ribbon.set_curvatures(l, m, n)
    radius = ribbon.get_radius()
    pitch = ribbon.get_pitch()
    ribbon.set_curvatures(None, None, n, radius=radius, pitch=pitch)
    assert np.isclose(ribbon.l, l, 1e-8, 1e-14)
    assert np.isclose(ribbon.m, m, 1e-8, 1e-14)



