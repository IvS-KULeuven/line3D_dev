import matplotlib.pyplot as plt
import numpy as np
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
def main(fname='../outputFILES/spec_surface.h5', oname2='./ps_files/surfb',
         oname3='./ps_files/surfbc',
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
      np_surfb, nzeta_surfb, p_surfb, zeta_surfb, \
      iem2d_surfb, iemi2d_surfb, iabs2d_surfb, icont2d_surfb = read_surfb(fname)

   xp_surfb = np.zeros(shape=(nzeta_surfb,np_surfb))
   yp_surfb = np.zeros(shape=(nzeta_surfb,np_surfb))                        
   for i in range(0,np_surfb):
      for j in range(0,nzeta_surfb):
         xp_surfb[j][i] = p_surfb[i]*np.sin(zeta_surfb[j])
         yp_surfb[j][i] = p_surfb[i]*np.cos(zeta_surfb[j])            
#
#-----------------------plotting----------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ncol=3
   nrow=3
   ax2 = fig2.subplots()

   iem2d_surfb = iem2d_surfb/icont2d_surfb
   iemi2d_surfb = iemi2d_surfb/icont2d_surfb
   iabs2d_surfb = iabs2d_surfb/icont2d_surfb

   clevels, ticks = get_clevels(clim=[0.,1.000001])
   cmap_lev = get_cmap('jet')
   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs)
   ax2.set_title(titlestr)
   ax2.set_xlabel(r'$x_p$')
   ax2.set_ylabel(r'$y_p$')
   ax2.axis('equal') #isotropic
   ax2.set_xlim(-np.max(p_surfb),np.max(p_surfb))
   ax2.set_ylim(-np.max(p_surfb),np.max(p_surfb))
   contourplot = ax2.contourf(xp_surfb, yp_surfb, iem2d_surfb,
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()   
   cbar = fig2.colorbar(contourplot,ax=ax2,orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$I/I_{cont}$')

#overplot a circle
   nd=101
   xc=np.sin(np.linspace(0.,2.*np.pi))
   yc=np.cos(np.linspace(0.,2.*np.pi))   
   ax2.plot(xc,yc,color='black', zorder=1)
#
#-------------------------continuum surface brightness----------------------------
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1./.4
   ysize=xsize/aspect_ratio
   fig3 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ax3 = fig3.subplots(1,2)

   icont2d_surfb
   icont2d_surfb = icont2d_surfb/xic1

#   indx=np.where(iem2d_surfb > 1.0000001)
#   print(indx)
#   exit()

   clevels, ticks = get_clevels(clim=[0.2,0.7])
   cmap_lev = get_cmap('Blues')   
#   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs_surfb[i])
#   ax3[0].set_title(titlestr)
   ax3[0].set_xlabel(r'$x_p$')
   ax3[0].set_ylabel(r'$y_p$')
   ax3[0].axis('equal') #isotropic
   ax3[0].set_xlim(-np.max(p_surfb),np.max(p_surfb))
   ax3[0].set_ylim(-np.max(p_surfb),np.max(p_surfb))
   contourplot = ax3[0].contourf(xp_surfb, yp_surfb, icont2d_surfb,
                              levels=clevels,
                              extend='both',
                                 cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('white')
   contourplot.changed()   
   cbar = fig3.colorbar(contourplot,ax=ax3[0],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$I/I_{core}$')

#overplot a circle
#   nd=101
#   xc=np.sin(np.linspace(0.,2.*np.pi))
#   yc=np.cos(np.linspace(0.,2.*np.pi))   
#   ax3[0].plot(xc,yc,color='black', zorder=1)
   

   clevels, ticks = get_clevels(clim=[0.,1.000001])
#   clevels, ticks = get_clevels(clim=[0.,0.2])   
   
   cmap_lev = get_cmap('jet')   
   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs)
   ax3[1].set_title(titlestr)
   ax3[1].set_xlabel(r'$x_p$')
   ax3[1].set_ylabel(r'$y_p$')
   ax3[1].axis('equal') #isotropic
   ax3[1].set_xlim(-np.max(p_surfb),np.max(p_surfb))
   ax3[1].set_ylim(-np.max(p_surfb),np.max(p_surfb))
   contourplot = ax3[1].contourf(xp_surfb, yp_surfb, iem2d_surfb,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()
   cbar = fig3.colorbar(contourplot,ax=ax3[1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$I/I_{cont}$')
#
#---------------------------------------------------------------------------------
#
#
   onamea = oname2+'.png'
   onameb = oname2+'.ps'
   fig2.savefig(onamea, bbox_inches='tight')
   fig2.savefig(onameb, bbox_inches='tight')   

   onamea = oname3+'.png'
   onameb = oname3+'.ps'
   fig3.savefig(onamea, bbox_inches='tight')
   fig3.savefig(onameb, bbox_inches='tight')   

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


fname='../outputFILES/nico_wr3d/spec_surface_00001.h5'
fname='../outputFILES/spec_surface_00001.h5'

windx = 1
main(fname=fname,windx=windx)



sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
