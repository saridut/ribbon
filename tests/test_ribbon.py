#!/usr/bin/env python

import math
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
from ribbon import Ribbon

#######################################################################
###FIGURE
fig_prop = {'figsize'        : (4, 3),
            'facecolor'      : '0.75',
            'edgecolor'      : 'white', 
            'subplot.left'   : 0.125,
            'subplot.right'  : 0.9, 
            'subplot.bottom' : 0.1,  
            'subplot.top'    : 0.9, 
            'autolayout'     : True
            }

###AXES
axes_prop = {'axisbelow'     : True, #'line',
             'unicode_minus' : False,
             'facecolor'     : 'white',
             'edgecolor'     : 'black',
             'linewidth'     : 0.75,
             'grid'          : False
             }

###LINES
lines_prop = {'linewidth'       : 0.8,
             'linestyle'       : '-',
             'color'           : 'black',
             'marker'          : 'None',
             'markeredgewidth' : 0.5,
             'markersize'      : 3,
             'dash_joinstyle'  : 'miter',
             'dash_capstyle'   : 'butt',
             'solid_joinstyle' : 'miter',
             'solid_capstyle'  : 'projecting',
             'antialiased'     : True
             }

###PATCH
patch_prop = {'linewidth'    : 1.0,
              'facecolor'    : 'blue',
              'edgecolor'    : 'black',
              'antialiased'  : True
              }

###TICK
xtick_prop = {'major.size' : 8, 
              'minor.size' : 4,
              'major.width': 0.4,
              'minor.width' : 0.4,
              'labelsize'  : 'small',
              'direction': 'in'
              }


ytick_prop = {'major.size' : 8, 
              'minor.size' : 4,
              'major.width': 0.4,
              'minor.width' : 0.4,
              'labelsize'  : 'small',
              'direction': 'in'
              }
###GRIDS
grid_prop = {'color'     : 'b0b0b0',
             'linestyle' : '-',
             'linewidth' : 0.4, # in points
             'alpha'     : 0.6  # transparency, between 0.0 and 1.0
             }

###LEGEND
leg_prop = {'fancybox'     : False,
            'numpoints'    : 1,    
            'fontsize'     : 'x-small',
            'borderpad'    : 0.5, 
            'markerscale'  : 1.0, 
            'labelspacing' : 0.5,
            'handlelength' : 2., 
            'handleheight' : 0.7,
            'handletextpad': 0.8,
            'borderaxespad': 0.5,
            'columnspacing': 2.,
            'frameon'      : False
            }

###FONT
#font_prop = {'family' : 'serif',
#              'serif' : 'Times',
#              'weight': 'medium',
#              'size'  : 10
#              }
              
###TEXT & LATEX
text_prop = {'usetex': True,
            'latex.preamble': 
                r'\usepackage{amssymb}'+
                r'\usepackage{amsmath}'+
                r'\usepackage{sansmathfonts}'+ 
                r'\usepackage[T1]{fontenc}'
                 }

###PS & EPS BACKEND
ps_prop = {'useafm': True}

#mpl.rc('figure', **fig_prop  )
#mpl.rc('axes',   **axes_prop )
#mpl.rc('lines',  **lines_prop)
#mpl.rc('patch',  **patch_prop)
#mpl.rc('xtick',  **xtick_prop)
#mpl.rc('ytick',  **ytick_prop)
#mpl.rc('grid',  **grid_prop)
#mpl.rc('legend', **leg_prop  )
##mpl.rc('font',   **font_prop )
#mpl.rc('text',   **text_prop )
#mpl.rc('ps',     **ps_prop   )

#######################################################################
#mpl.use('PDF')
import matplotlib.pyplot as plt

#Plotting
figh, axh = plt.subplots(nrows=1, ncols=1)

length = 20
width = 2
thickness = 0.0
l = 0.4
#l = lambda x: 0.5 *(1-x/length)
#l = lambda x: 0.2 *(x**0.4)
m = 0.4 #lambda x: 0.5 #*(1-x/length) #+ve right handed, -ve left handed
n = 0.4 #lambda x: 0.5*np.cos(x) #0.5
#n = lambda x: 0.5*(1-x/length)
    
#   atom_coords = []
#   for x in np.arange(0, length+0.025, 0.4):
#       for y in np.arange(-width/2, 0.025+width/2, 0.4):
#           for z in np.arange(-thickness/2, 0.025+thickness/2, 0.3):
#               atom_coords.append([x,y,z])

    #ribbon = Ribbon(length, width, thickness, 0.1, 0.1, 0.08, np.asarray(atom_coords))
ribbon = Ribbon(length, width, thickness, 0.1, 0.1, 0.08)
ribbon.set_curvatures(l, m, n)
#out = ribbon.get_radius()
#if isinstance(out, tuple):
#    for x,y in zip(out[0],out[1]):
#        print(f"{x:g}  {y:g}") 
#else:
#    print(f"R = {out}")
#print(f"R = {ribbon.get_radius()}\n"
#      f"P = {ribbon.get_pitch()}\n"
#      f"kg = {ribbon.get_gauss_curvature()}\n"
#      f"km = {ribbon.get_mean_curvature()}\n"
#      f"theta = {ribbon.get_theta()}")
#raise SystemExit()
ribbon.create(orient_along=[0,0,1])


figh, axh = plt.subplots(nrows=1, ncols=1 ,
        subplot_kw={'projection':'3d', 'proj_type': 'ortho'},
                         figsize=(12,9))

axh.plot(ribbon.mline[:,0], ribbon.mline[:,1], ribbon.mline[:,2], '-k', lw=1.2)
axh.plot(ribbon.mline[0,0], ribbon.mline[0,1], ribbon.mline[0,2], 'or')

#axh.plot_surface(ribbon.msurf[:,:,0], ribbon.msurf[:,:,1],
#            ribbon.msurf[:,:,2], facecolor='r', edgecolor='0.6', lw=0.7, alpha=0.4,
#            rstride=4, cstride=8)

axh.plot_wireframe(ribbon.msurf[:,:,0], ribbon.msurf[:,:,1],
            ribbon.msurf[:,:,2], linestyles='-', linewidths=0.5,
            rstride=2, cstride=8, color='0.5')

frl = np.linspace(0,ribbon.u.size-1,10, dtype=np.int32)
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d1[frl,0], ribbon._d1[frl,1], ribbon._d1[frl,2],
           length=1.2, color='r')
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d2[frl,0], ribbon._d2[frl,1], ribbon._d2[frl,2],
           length=1.2, color='g')
axh.quiver(ribbon.mline[frl,0], ribbon.mline[frl,1], ribbon.mline[frl,2],
           ribbon._d3[frl,0], ribbon._d3[frl,1], ribbon._d3[frl,2],
           length=1.2, color='b')

axh.set_aspect('equal', 'box')
#axh.set_xlabel('x')
#axh.set_ylabel('y')
#axh.set_zlabel('z')
#axh.set_xticks(axh.get_xlim())
#axh.set_yticks([])
#axh.set_zticks([])
axh.set_axis_off()
plt.show()

