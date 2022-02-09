import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

from const_cgs import *
from read_vbin import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from read_models import *
from lev_contour import *
import imageio
import os
#
#
#
def main(fname='../outputFILES/spec_int2d.h5', oname1='./ps_files/int2d',
         oname2='./ps_files/tau2d', clim_int2d=[-5.,-1.], clim_tau=[-2.,0.],
         windx=1):
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

   xic1, vth_fiducial, xnue0, xobs, alpha, gamma, \
      nx, nz, x_int2d, z_int2d, \
      int2d, iemi2d, iabs2d, tau2d, vn2d = read_int2d(fname)
#
#-----------------------plotting----------------------------------------
#
   fig1 = plt.figure(windx,figsize=(8,10),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ncol=3
   nrow=3
   ax = fig1.subplots()

   clevels, ticks = get_clevels(clim=clim_int2d)
   cmap_lev = get_cmap('jet')
   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs)
   ax.set_title(titlestr)
   ax.set_xlabel(r'$z_p$')
   ax.set_ylabel(r'$x_p$')
   ax.set_xlim(np.min(x_int2d),np.max(x_int2d))
   ax.set_ylim(np.min(z_int2d),np.max(z_int2d))
   
   contourplot = ax.contourf(x_int2d, z_int2d, np.log10(int2d),
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()   
   cbar = fig1.colorbar(contourplot,ax=ax,orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(I)$')

#overplot a circle
   nd=101
   xc=np.sin(np.linspace(0.,2.*np.pi))
   yc=np.cos(np.linspace(0.,2.*np.pi))   
   ax.plot(xc,yc,color='black', zorder=1)
   
   onamea = oname1+'.png'
   onameb = oname1+'.ps'
   fig1.savefig(onamea, bbox_inches='tight')
   fig1.savefig(onameb, bbox_inches='tight')   
#
#-----------------------plotting----------------------------------------
#
   fig2 = plt.figure(windx,figsize=(8,10),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ncol=3
   nrow=3
   ax = fig2.subplots()

   clevels, ticks = get_clevels(clim=clim_tau)
   cmap_lev = get_cmap('jet')
   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs)
   ax.set_title(titlestr)
   ax.set_xlabel(r'$z_p$')
   ax.set_ylabel(r'$x_p$')
   ax.set_xlim(np.min(x_int2d),np.max(x_int2d))
   ax.set_ylim(np.min(z_int2d),np.max(z_int2d))
   contourplot = ax.contourf(x_int2d, z_int2d, np.log10(tau2d),
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()   
   cbar = fig2.colorbar(contourplot,ax=ax,orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\tau_z)$')

#overplot a circle
   nd=101
   xc=np.sin(np.linspace(0.,2.*np.pi))
   yc=np.cos(np.linspace(0.,2.*np.pi))   
   ax.plot(xc,yc,color='black', zorder=1)   

   onamea = oname2+'.png'
   onameb = oname2+'.ps'
   fig2.savefig(onamea, bbox_inches='tight')
   fig2.savefig(onameb, bbox_inches='tight')      
#


#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


fname='../outputFILES/spec_int2d_00001.h5'

windx = 1
main(fname=fname,windx=windx)



sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
