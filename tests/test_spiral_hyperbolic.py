#!/usr/bin/env python

from pathlib import Path
import sys
import numpy as np
import pytest
try:
    from ribbon import SpiralHyperbolic
except:
    srcdir = Path('src').resolve()
    sys.path.insert(0, str(srcdir))
    from ribbon import SpiralHyperbolic


@pytest.fixture
def spiral():
    return SpiralHyperbolic(1.0, 0.1)


@pytest.fixture
def spiral_data():
    fn = 'tests/test_data_spiral_hyperbolic.txt'
    data = {}
    with open(fn, 'r') as fh:
        fh.readline() # Throwaway the first line
        head_words = fh.readline().strip(' \n').split()
        for key in head_words:
            data[key] = []
        while line := fh.readline():
            words = line.strip(' \n').split()
            if not words:
                continue
            else:
                vals = [float(x) for x in words]
                for i,v in enumerate(vals):
                    data[head_words[i]].append(v)
    for key,val in data.items():
        data[key] = np.array(val, dtype=np.float64)
    return data
           

def test_init_b():
    with pytest.raises(ValueError) as info:
        spiral = SpiralHyperbolic(0, 1)
    assert info.type is ValueError


def test_init_r0():
    with pytest.raises(ValueError) as info:
        spiral = SpiralHyperbolic(1, -1)
    assert info.type is ValueError


def test_r_to_theta(spiral, spiral_data):
    #Scalar version
    theta = spiral.r_to_theta(spiral_data['r'][0])
    assert np.allclose(theta, spiral_data['theta'][0], 1e-8, 1e-14)
    #Array version
    theta = spiral.r_to_theta(spiral_data['r'])
    assert np.allclose(theta, spiral_data['theta'], 1e-8, 1e-14)


def test_theta_to_r(spiral, spiral_data):
    #Scalar version
    r = spiral.theta_to_r(spiral_data['theta'][0])
    assert np.allclose(r, spiral_data['r'][0], 1e-8, 1e-14)
    #Array version
    r = spiral.theta_to_r(spiral_data['theta'])
    assert np.allclose(r, spiral_data['r'], 1e-8, 1e-14)


def test_arclength_to_theta(spiral, spiral_data):
    #Scalar version
    theta = spiral.arclength_to_theta(spiral_data['s'][0])
    assert np.allclose(theta, spiral_data['theta'][0], 1e-8, 1e-14)
    #Array version
    theta = spiral.arclength_to_theta(spiral_data['s'])
    assert np.allclose(theta, spiral_data['theta'], 1e-8, 1e-14)


def test_get_arclength(spiral, spiral_data):
    #From r
    s = spiral.get_arclength(spiral_data['r'], var='r')
    assert np.allclose(s, spiral_data['s'], 1e-8, 1e-14)
    #From theta
    s = spiral.get_arclength(spiral_data['theta'], var='theta')
    assert np.allclose(s, spiral_data['s'], 1e-8, 1e-14)


def test_get_curvature(spiral, spiral_data):
    #From r
    kappa = spiral.get_curvature(spiral_data['r'], var='r')
    assert np.allclose(kappa, spiral_data['kappa'], 1e-8, 1e-14)
    #From theta
    kappa = spiral.get_curvature(spiral_data['theta'], var='theta')
    assert np.allclose(kappa, spiral_data['kappa'], 1e-8, 1e-14)
    #From s
    kappa = spiral.get_curvature(spiral_data['s'], var='s')
    assert np.allclose(kappa, spiral_data['kappa'], 1e-8, 1e-14)


def test_get_tangent(spiral, spiral_data):
    #From r
    tangent = spiral.get_tangent(spiral_data['r'], var='r')
    assert np.allclose(tangent[:,0], spiral_data['tx'], 1e-8, 1e-14)
    assert np.allclose(tangent[:,1], spiral_data['ty'], 1e-8, 1e-14)
    #From theta
    tangent = spiral.get_tangent(spiral_data['theta'], var='theta')
    assert np.allclose(tangent[:,0], spiral_data['tx'], 1e-8, 1e-14)
    assert np.allclose(tangent[:,1], spiral_data['ty'], 1e-8, 1e-14)
    #From s
    tangent = spiral.get_tangent(spiral_data['s'], var='s')
    assert np.allclose(tangent[:,0], spiral_data['tx'], 1e-8, 1e-14)
    assert np.allclose(tangent[:,1], spiral_data['ty'], 1e-8, 1e-14)


def test_get_polar_coords(spiral, spiral_data):
    #From r
    r, theta = spiral.get_polar_coords(spiral_data['r'], var='r')
    assert np.allclose(r, spiral_data['r'], 1e-8, 1e-14)
    assert np.allclose(theta, spiral_data['theta'], 1e-8, 1e-14)
    #From theta
    r, theta = spiral.get_polar_coords(spiral_data['theta'], var='theta')
    assert np.allclose(r, spiral_data['r'], 1e-8, 1e-14)
    assert np.allclose(theta, spiral_data['theta'], 1e-8, 1e-14)
    #From s
    r, theta = spiral.get_polar_coords(spiral_data['s'], var='s')
    assert np.allclose(r, spiral_data['r'], 1e-8, 1e-14)
    assert np.allclose(theta, spiral_data['theta'], 1e-8, 1e-14)


def test_get_cart_coords(spiral, spiral_data):
    #From r
    x, y = spiral.get_cart_coords(spiral_data['r'], var='r')
    assert np.allclose(x, spiral_data['x'], 1e-8, 1e-14)
    assert np.allclose(y, spiral_data['y'], 1e-8, 1e-14)
    #From theta
    x, y = spiral.get_cart_coords(spiral_data['theta'], var='theta')
    assert np.allclose(x, spiral_data['x'], 1e-8, 1e-14)
    assert np.allclose(y, spiral_data['y'], 1e-8, 1e-14)
    #From s
    x, y = spiral.get_cart_coords(spiral_data['s'], var='s')
    assert np.allclose(x, spiral_data['x'], 1e-8, 1e-14)
    assert np.allclose(y, spiral_data['y'], 1e-8, 1e-14)
