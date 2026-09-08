#!/usr/bin/env python

from pathlib import Path
import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
try:
    from ribbon_ext import Ribbon, write_xyz
except:
    srcdir = Path('../src').resolve()
    sys.path.insert(0, str(srcdir))
    from ribbon import Ribbon, write_xyz

import matplotlib.pyplot as plt


length = 36 #2*np.pi
width = 4.0
thickness = 0.0
#l = 0.5
#l = lambda x: 0.5 *(1-x/length)
#l = lambda x: 0.2 *(x**0.4)
#m = 0.5 #lambda x: 0.5 #*(1-x/length) #+ve right handed, -ve left handed
#n = 0.4 #lambda x: 0.5*np.cos(x) #0.5
#n = lambda x: 0.5*(1-x/length)
    
#Atoms reference position
#atom_coords = []
#for x in np.arange(0, length, 0.4):
#    for y in np.arange(-width/2, width/2, 0.4):
#            atom_coords.append([x,y,0])

ribbon = Ribbon(length, width, thickness, 0.1, 0.1, 0.0)
#ribbon.set_atom_refpos(np.asarray(atom_coords))
ribbon.set_shape('HelicalRibbon', radius=1.5, pitch=10)
#ribbon.set_shape('Helicoid', pitch=16)
#ribbon.set_shape('Cylinder', radius=2)
#ribbon.set_shape('Torus', Radius=10 , radius=1)
#ribbon.set_shape('Plane')
#ribbon.set_shape('General', radius=2.5, pitch=10,
#                   profile='Circle', profile_params={'radius': 20.0})
#ribbon.set_shape('General', radius=1.2, pitch=6, profile='Line')
#ribbon.set_shape('General', radius=1, pitch=6,
#                 profile='Helix', profile_params={'radius': 1.0, 'pitch': -6})
#Catenary
#def f(s, a=2.5): 
#    x = a * math.asinh(s/a)
#    y = a * math.cosh(x/a)
#    return x, y
#n = 1.0/2.5
#ribbon.set_shape('General', radius=1, pitch=6,
#                 profile='General', profile_params={'n': n, 'f': f})

#Varying curvature
#ribbon.set_shape('General', radius=lambda x: 1.0 + 0.05*x,
#                 pitch=0.0,
#                 profile='Circle', profile_params={'radius': -5.0})

#ribbon.create(orient_along=None)
ribbon.create(orient_along=[1,0,0])

#out = ribbon.get_radius()
#if isinstance(out, tuple):
#    for x,y in zip(out[0],out[1]):
#        print(f"{x:g}  {y:g}") 
#else:
#    print(f"R = {out}")
print(f"R = {ribbon.get_radius()}\n"
      f"P = {ribbon.get_pitch()}\n"
      f"kg = {ribbon.get_gauss_curvature()}\n"
      f"km = {ribbon.get_mean_curvature()}\n"
      f"theta = {ribbon.get_theta()}")

#write_xyz(ribbon.atom_pos[:,0], ribbon.atom_pos[:,1], ribbon.atom_pos[:,2]) 

#raise SystemExit()
figh, axh = plt.subplots(nrows=1, ncols=1, figsize=(12,9),
                subplot_kw={'projection':'3d', 'proj_type': 'persp'}
                )

axh.plot(ribbon.mline[:,0], ribbon.mline[:,1], ribbon.mline[:,2], '-k', lw=1.2)
axh.plot(ribbon.mline[0,0], ribbon.mline[0,1], ribbon.mline[0,2], 'ob')
#axh.plot_surface(ribbon.msurf[:,:,0], ribbon.msurf[:,:,1],
#            ribbon.msurf[:,:,2], facecolor='r', edgecolor='0.6', lw=0.7, alpha=0.4,
#            rstride=4, cstride=8)

axh.plot(ribbon.msurf[0,:,0], ribbon.msurf[0,:,1], ribbon.msurf[0,:,2], '-g', lw=1.0)
axh.plot(ribbon.msurf[-1,:,0], ribbon.msurf[-1,:,1], ribbon.msurf[-1,:,2], '-g', lw=1.0)
axh.plot(ribbon.msurf[:,0,0], ribbon.msurf[:,0,1], ribbon.msurf[:,0,2], '-r', lw=1.0)
axh.plot(ribbon.msurf[:,-1,0], ribbon.msurf[:,-1,1], ribbon.msurf[:,-1,2], '-r', lw=1.0)

j = -1
axh.plot_wireframe(ribbon.msurf[:j,:,0], ribbon.msurf[:j,:,1],
            ribbon.msurf[:j,:,2], linestyles='-', linewidths=0.5,
            rstride=2, cstride=8, color='0.7')

frl = np.linspace(0, ribbon.u.size-1, 10, dtype=np.int32)
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d1[frl,0], ribbon._d1[frl,1], ribbon._d1[frl,2],
           length=1.0, color='r')
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d2[frl,0], ribbon._d2[frl,1], ribbon._d2[frl,2],
           length=1.0, color='g')
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d3[frl,0], ribbon._d3[frl,1], ribbon._d3[frl,2],
           length=1.0, color='b')

axh.set_aspect('equal', 'box')
#axh.set_xlabel('x')
#axh.set_ylabel('y')
#axh.set_zlabel('z')
#axh.set_xticks(axh.get_xlim())
#axh.set_yticks([])
#axh.set_zticks([])
axh.set_axis_off()
plt.show()

