#!/usr/bin/env python

import math
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
from spirauliya.spiralpolar import SpiralCircleInvolute as spiral
#from spirals import double_spiral
#from spirals import SpiralNielsen as spiral
#from spirals import SpiralPolynomial as spiral
#from spirals import SpiralCornu as cspiral
#from spirals import extrude


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

L = 15.0
cnst_b = -0.1
cnst_a = 2.0
r0 = cnst_b+1e-9

#sa = spiral([1.0, 0.1, 0.3])
sa = spiral(1, 1.01)
#print(type(sa.super()).__name__)
#print(sa.r_to_theta.__doc__)
#sa = cspiral(1.0)
#print(f"b = {sa.b}\n n = {sa.n}\n r0 = {sa.r0}\n t0 = {sa.t0}\n tincr = {sa.tincr}")
#print(f"b = {sa.b}\n r0 = {sa.r0}\n t0 = {sa.t0}\n tincr = {sa.tincr}")
#print(sa._func_kappa)
#print(sa._func_int_kappa)
#print(sa.get_arclength('r', 10.0))

n = 200
v = np.linspace(0.0, 5, n)
#t = sa.arclength_to_theta(v)
#print('t = ', t)
#al = sa.get_arclength(v, 'theta')
#print('al = ', al)

#raise SystemExit()
#kappa = sa.get_curvature(v, 'theta')
#tangent = sa.get_tangent(5)
#print(tangent, np.linalg.norm(tangent))
#print(kappa)
#axh.plot(v, t, ls='-', marker='None', label='_nolegend_')
#axh.plot(al, v, ls='-', marker='None', label='_nolegend_')

#xdata, ydata = sa.get_cart_coords(v)
#axh.plot(xdata, ydata, ls='-', c='r', marker='None', label='_nolegend_')

ds_data = double_spiral(sa, L, ds=1.0, same=False, end=1, f=0.05)
s = ds_data[0]
kappa = ds_data[1]
x = ds_data[2]
y = ds_data[3]
#for u,v in zip(s,kappa):
#    print(f"{u:g} {v:g}")
extrude(x, y, 0, 1.0, 'C', 'extrusion.xyz', 'This file')

#axh.plot(s, kappa, ls='-', marker='.', label='_nolegend_')
axh.plot(x, y, ls='-', marker='None', label='_nolegend_')

#axh.plot(xdata, ydata, ls='-', marker='o', c='r', label='_nolegend_')
#axh.arrow(xdata[-1], ydata[-1], tangent[0]/2, tangent[1]/2, width=0.02)
#axh.annotate("", xytext=(0, 0), xy=(xdata[0], ydata[0]),
#             arrowprops=dict(arrowstyle="->"))


#axh.plot(xdata[0], ydata[0], ls='None', marker='o', c='r', label='_nolegend_')

#func_kappa = lambda x: np.tanh(x)
#sl = SpiralGeneral(func_kappa)
#s = sl.get_arclength(sf.tmin, sf.tmin+4*math.pi)
#xdata, ydata = sl.get_cart_coords(0, 5, 0.1)
#axh.plot(xdata, ydata, ls=':', marker='.', c = 'k',
#    label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))


axh.set_aspect('equal')
#axh.set_axis_off()
axh.grid(True)
#axh.tick_params(which='both', top=True, right=True, color='0.5')
#axh.set_xscale('linear')
#axh.set_yscale('linear')
#legend = axh.legend(loc='best', ncol=1)

#axh.set_xlim(-2, 2)
#axh.set_ylim(-2, 2)

#axh.xaxis.set_minor_locator(ticr.AutoMinorLocator(n=4))
#axh.yaxis.set_minor_locator(ticr.AutoMinorLocator(n=4))

plt.show()
#plt.close()
