import sys
sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_vbin import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from lev_contour import *
#
def main(fname='../outputFILES/spec_surface.h5', oname='./ps_files/spec_surface',
         windx=1, clim=[-3.,0.], xlim=[-3.,3.],
         ylim=[-3.,3.]):
#+
# NAME:
#	plot_surfb_vbin
#
# PURPOSE:
#	This procedure plots surface brightnes of spec_vbin.eo output
#
#-
#
#------------------read all information from hdf5-file------------------
#

   npoints, points_xcoord, points_ycoord, xic1, xic2, xobs, \
      alpha, gamma, int2d_tot, int2d_abs, int2d_emi=read_surfb_vbin(fname)

   int2d_tot=int2d_tot/xic1

   int2d_plot = int2d_tot*1.
   indx1 = np.where(int2d_tot > 0.)
   indx2 = np.where(int2d_tot <= 0.)   
   int2d_plot[indx1] = np.log10(int2d_tot[indx1])
   int2d_plot[indx2] = -15.0

   clevels, ticks = get_clevels(clim=clim)
#
#-----------------------plotting----------------------------------------
#

   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)


   #
   plt.ion()
   plt.show()
#
   ax = fig1.subplots()
#   ax.set_title(r'(\log (I/I_{c}')
   ax.set_xlabel(r'$x_p$')
   ax.set_ylabel(r'$y_p$')
   ax.axis('equal') #isotropic
   ax.set_xlim(xlim)
   ax.set_ylim(ylim)
   
   
#;plot_triangles, '../outputFILES', obs_indx=1, oname=oname, xlim=xlim, ylim=ylim


#   x_coord=np.zeros(shape=(ntheta,nr))
#   y_coord=np.zeros(shape=(ntheta,nr))
#   for i in range(0,nr):
#       for j in range(0,ntheta):
#          x_coord[j][i] = r[i]*np.sin(theta[j])
#          y_coord[j][i] = r[i]*np.cos(theta[j])
#
   contourplot = ax.tricontourf(points_ycoord, points_xcoord, np.transpose(int2d_plot),
                  levels=clevels,
                  extend='both',
                  cmap='jet')


#ax2.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
#cntr2 = ax2.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

   contourplot.cmap.set_over('red')
   contourplot.cmap.set_under('black')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax,orientation='vertical',ticks=ticks)
   cbar.set_label(r'$\log (I/I_{c})$')
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


#fname='modelw_surfb01/spec_surface.h5'
#oname='ps_files/spec_surface01.ps'  
fname='../outputFILES/spec_surface.h5'

windx = 1
main(fname=fname,windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
