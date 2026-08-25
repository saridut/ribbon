#!/usr/bin/env python

import math
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
from scipy.special import fresnel
from scipy.integrate import romb
from spirals import SpiralFermat


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

mpl.rc('figure', **fig_prop  )
mpl.rc('axes',   **axes_prop )
mpl.rc('lines',  **lines_prop)
mpl.rc('patch',  **patch_prop)
mpl.rc('xtick',  **xtick_prop)
mpl.rc('ytick',  **ytick_prop)
mpl.rc('grid',  **grid_prop)
mpl.rc('legend', **leg_prop  )
#mpl.rc('font',   **font_prop )
mpl.rc('text',   **text_prop )
mpl.rc('ps',     **ps_prop   )

#######################################################################
mpl.use('PDF')
import matplotlib.pyplot as plt

#Plotting
figh, axh = plt.subplots(nrows=1, ncols=1)

L = 100
cnst_a = 1
cnst_b = 1.6
#cnst_b = 0.11
rmin = 0.0

for cnst_b in [1.6]:
    sa = SpiralArchimedes(cnst_b, rmin)
#sa = SpiralFermat(cnst_b, rmin)
#sa = SpiralLogarithmic(cnst_a, cnst_b, rmin)
#sa = SpiralClothoid(cnst_a)

#ds = DoubleSpiral(sa, 2*L, 'D', is_tapered=True, taper_angle=math.radians(35))

    n = 257
    s = np.linspace(0, L, n)
    kappa = sa.get_curvature(s, 's')
    ds = L/(n-1)
    xdata, ydata = sa.get_cart_coords(0, L, ds, 's')
    av_kappa = romb(kappa*s, dx=ds)/romb(s, dx=ds)
    av_rho = 1/av_kappa
    #av_rho = 1/np.mean(kappa)
    print(av_rho)
    #with open('sa_b_1.6_L_100_rmin_0.txt', 'w') as fh:
    #    fh.write('x,y,s,kappa\n')
    #    for i in range(len(s)):
    #        fh.write('%g,%g,%g,%g\n'%(xdata[i], ydata[i], s[i], kappa[i]))
    #axh.plot(xdata, ydata, ls='-', lw=4, marker='None', label='_nolegend_')


#xdata, ydata = sa.get_cart_coords(0, L, 0.1, 's')
#axh.plot(xdata, ydata, ls='-', lw=1.2, marker='None', c = '0.6',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))
#
#axh.plot([0], [0], ls='None', marker='x', c = '0.6',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))
#axh.plot(-xdata, ydata, ls='-', lw=1.2, marker='None', c = 'r',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))
#axh.plot([-cnst_a/2], [cnst_a/2], ls='None', marker='x', c = 'r',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))

#xdata, ydata, n = ds.get_cart_coords(0.1)
#axh.plot(xdata[0:n], ydata[0:n], ls='-', marker='None', c = '0.6',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))
#axh.plot(xdata[n:], ydata[n:], ls='-', marker='None', c = 'r',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))
#axh.plot([-cnst_a/2, cnst_a/2], [-cnst_a/2, cnst_a/2], ls='None', marker='x', c = '0.6',
#     label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))

#func_kappa = lambda x: np.tanh(x)
#sl = SpiralGeneral(func_kappa)
#s = sl.get_arclength(sf.tmin, sf.tmin+4*math.pi)
#xdata, ydata = sl.get_cart_coords(0, 5, 0.1)
#axh.plot(xdata, ydata, ls=':', marker='.', c = 'k',
#    label=r'$b = %g, r_{\mathrm{min}} = %g$'%(cnst_b, rmin))


axh.set_aspect('equal')
axh.set_axis_off()
#axh.grid(True)
#axh.tick_params(which='both', top=True, right=True, color='0.5')
#axh.set_xscale('linear')
#axh.set_yscale('linear')
#legend = axh.legend(loc='best', ncol=1)

#axh.set_xlim(-2, 2)
#axh.set_ylim(-2, 2)

#axh.xaxis.set_minor_locator(ticr.AutoMinorLocator(n=4))
#axh.yaxis.set_minor_locator(ticr.AutoMinorLocator(n=4))

#axh.set_xlabel(r'$s\;(\mathrm{nm})$', labelpad=5.0)
#axh.set_ylabel(r'$\rho\;(\mathrm{nm})$', labelpad=5.0)
#plt.savefig('fig-cspiral.pdf', transparent=True)
plt.close()
