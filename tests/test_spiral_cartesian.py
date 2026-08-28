from pathlib import Path
import sys
import numpy as np
import pytest
try:
    from ribbon import SpiralPolynomial, SpiralCornu, SpiralNielsen
except:
    srcdir = Path('src').resolve()
    sys.path.insert(0, str(srcdir))
    from ribbon import SpiralPolynomial, SpiralCornu, SpiralNielsen

fix_params = [
        ('tests/test_data_spiral_polynomial.txt', SpiralPolynomial, ([1,1,1],)),
        ('tests/test_data_spiral_cornu.txt', SpiralCornu, (1,)),
        ('tests/test_data_spiral_nielsen_1.txt', SpiralNielsen, (1, 0.5)),
        ('tests/test_data_spiral_nielsen_2.txt', SpiralNielsen, (1, -0.5)),
    ]

@pytest.fixture(scope='module', params=fix_params)
def fx_spiral(request):
    fn = request.param[0]
    data = spiral_data(fn)
    spiral = request.param[1]
    args = request.param[2]
    return spiral(*args), data


def spiral_data(fn):
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
           

def test_init_polynomial():
    with pytest.raises(ValueError) as info:
        spiral = SpiralPolynomial([0, 0, 0])
    assert info.type is ValueError


def test_init_cornu():
    with pytest.raises(ValueError) as info:
        spiral = SpiralCornu(-1)
    assert info.type is ValueError

def test_init_nielsen_a():
    with pytest.raises(ValueError) as info:
        spiral = SpiralNielsen(-1, -1)
    assert info.type is ValueError


def test_init_nielsen_b():
    with pytest.raises(ValueError) as info:
        spiral = SpiralNielsen(1, 0)
    assert info.type is ValueError


def test_get_curvature(fx_spiral):
    spiral = fx_spiral[0]
    spiral_data = fx_spiral[1]
    kappa = spiral.get_curvature(spiral_data['s'])
    assert np.allclose(kappa, spiral_data['kappa'], 1e-8, 1e-14)


def test_get_tangent(fx_spiral):
    spiral = fx_spiral[0]
    spiral_data = fx_spiral[1]
    tangent = spiral.get_tangent(spiral_data['s'])
    assert np.allclose(tangent[:,0], spiral_data['tx'], 1e-8, 1e-14)
    assert np.allclose(tangent[:,1], spiral_data['ty'], 1e-8, 1e-14)


def test_get_cart_coords(fx_spiral):
    spiral = fx_spiral[0]
    spiral_data = fx_spiral[1]
    x, y = spiral.get_cart_coords(spiral_data['s'])
    assert np.allclose(x, spiral_data['x'], 1e-8, 1e-14)
    assert np.allclose(y, spiral_data['y'], 1e-8, 1e-14)
