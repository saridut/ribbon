#!/usr/bin/env python

from pathlib import Path
import sys
import numpy as np
try:
    import ribbon
except:
    srcdir = Path('../src').resolve()
    sys.path.insert(0, str(srcdir))
    import ribbon

'''
#Test data for polar spirals
# b = 1.0, r0 = 0.1
fn = 'test_data_circleinvolute.txt'
b = 1.0
r0 = 1.1
spiral = ribbon.SpiralCircleInvolute(b, r0)
dr = 0.2 #Spacing in r
n = 10 #Number of data points
r = r0 + dr*np.arange(0, n, 1, dtype=np.float64)
theta = spiral.r_to_theta(r)
s = spiral.get_arclength(theta, var='theta', v0=None)
kappa = spiral.get_curvature(s, var='s')
tangent = spiral.get_tangent(s, var='s')
x, y = spiral.get_cart_coords(s, var='s')

with open(fn, 'w') as fh:
    fh.write(f"# Circle involute : b = {b}, r0 = {r0}\n")
    header_words = ['s', 'r', 'theta', 'x', 'y', 'kappa', 'tx', 'ty']
    header = ' '.join([f"{x:^21s}" for x in header_words])
    fh.write(header + '\n')
    for i in range(n):
        values = [s[i], r[i], theta[i], x[i], y[i], kappa[i], tangent[i,0],
                  tangent[i,1]]
        buf = ' '.join([f"{x:^ 21.14e}" for x in values])
        fh.write(buf + '\n')
'''

#Test data for cartesian spirals
fn = 'test_data_spiralnielsen_2.txt'
#coeffs = [1.0, 1.0, 1.0]
spiral = ribbon.SpiralNielsen(1.0, -0.5)
ds = 0.5 #Spacing in s
n = 10 #Number of data points
s = ds*np.arange(0, n, 1, dtype=np.float64)
kappa = spiral.get_curvature(s)
tangent = spiral.get_tangent(s)
x, y = spiral.get_cart_coords(s)

with open(fn, 'w') as fh:
    fh.write(f"# Nielsen spiral : a = 1, b = -0.5\n")
    header_words = ['s', 'x', 'y', 'kappa', 'tx', 'ty']
    header = ' '.join([f"{x:^21s}" for x in header_words])
    fh.write(header + '\n')
    for i in range(n):
        values = [s[i], x[i], y[i], kappa[i], tangent[i,0], tangent[i,1]]
        buf = ' '.join([f"{x:^ 21.14e}" for x in values])
        fh.write(buf + '\n')
