#!/usr/bin/env python

from pathlib import Path
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
try:
    import ribbon
except:
    srcdir = Path('../src').resolve()
    sys.path.insert(0, str(srcdir))
    import ribbon

import matplotlib.pyplot as plt


figh, axh = plt.subplots(nrows=1, ncols=2, figsize=(8,3))

L = 20.0
b = 0.2
r0 = 0.0

ds = 0.1
s = np.linspace(0, L, 1+int(L/ds))

spl = ribbon.SpiralArchimedes(b, r0)

ds = ribbon.double_spiral(spl, L, ds=ds, same=False, end=1, f=0.0)
ds_s = ds[0]
ds_kappa = ds[1]
ds_x = ds[2]
ds_y = ds[3]

width = 4.0
dz = 0.1
X, Y, Z = ribbon.extrude(ds_x, ds_y, width, dz)
ribbon.write_xyz(X, Y, Z, filename='dspl.xyz', title='')

axh[0].plot(ds_x, ds_y, ls='-', marker='None', label='_nolegend_')
axh[0].plot(ds_x[0], ds_y[0], ls='None', marker='.', label='_nolegend_')
axh[1].plot(ds_s, ds_kappa, ls='-', marker='None', label='_nolegend_')

#Plot tangent vectors
#axh.arrow(xdata[-1], ydata[-1], tangent[0]/2, tangent[1]/2, width=0.02)
#axh.annotate("", xytext=(0, 0), xy=(xdata[0], ydata[0]),
#             arrowprops=dict(arrowstyle="->"))

axh[0].set_aspect('equal')
plt.show()
#plt.savefig('fn')
