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
#from colorspacious import cspace_converter
#
def main(fnames=['../outputFILES/spec_surface.h5'], windx=1, create_animation=False):
#+
# NAME:
#	plot_surfb_vbin
#
# PURPOSE:
#	This procedure plots surface brightnes of spec_vbin.eo output
#
#-
#

   animation_files = []
   animation_images = []   

   clevels, ticks = get_clevels(clim=[0.,1.000001])
   cmap_lev = get_cmap('jet')
#
#-----------------------------------------------------------------------
#
   indx=0
   for fname in fnames:
      print('reading file', fname)
#
#------------------read all information from hdf5-file------------------      
#
      xic1, vth_fiducial, xnue0, xobs, alpha, gamma, \
         np_surfb, nzeta_surfb, p_surfb, zeta_surfb, \
         iem2d_surfb, iemi2d_surfb, iabs2d_surfb, icont2d_surfb = read_surfb(fname)

      indx1=np.where(icont2d_surfb == 0.)
      indx2=np.where(icont2d_surfb > 0.)

      iem2d_plot = np.copy(iem2d_surfb)
      iem2d_plot[indx1] = 1.
      iem2d_plot[indx2] = iem2d_surfb[indx2]/icont2d_surfb[indx2]
      
      xp_surfb = np.zeros(shape=(nzeta_surfb,np_surfb))
      yp_surfb = np.zeros(shape=(nzeta_surfb,np_surfb))
      for i in range(0,np_surfb):
         for j in range(0,nzeta_surfb):
            xp_surfb[j][i] = p_surfb[i]*np.sin(zeta_surfb[j])
            yp_surfb[j][i] = p_surfb[i]*np.cos(zeta_surfb[j])                  
#
#------------------------------plotting----------------------------------
#
      xsize=18. #in cm
      xsize=xsize/2.54 #in inches
      aspect_ratio=1./.8
      ysize=xsize/aspect_ratio
      fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
      #windx=windx+1
      #
      #   plt.ion()
      #
      ax = fig1.subplots()
   
      titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs)
      ax.set_title(titlestr)
      ax.set_xlabel(r'$x_p$')
      ax.set_ylabel(r'$y_p$')
      ax.axis('equal') #isotropic      
      ax.set_xlim(-np.max(p_surfb),np.max(p_surfb))
      ax.set_ylim(-np.max(p_surfb),np.max(p_surfb))
      contourplot = ax.contourf(xp_surfb, yp_surfb, iem2d_plot,
                                levels=clevels,
                                extend='both',
                                cmap=cmap_lev)
      contourplot.cmap.set_over('black')
      contourplot.cmap.set_under('grey')
      contourplot.changed()
      cbar = fig1.colorbar(contourplot,ax=ax,orientation='vertical',aspect=40,ticks=ticks)
      cbar.set_label(r'$I/I_{cont}$')

      #overplot a circle
      nd=101
      xc=np.sin(np.linspace(0.,2.*np.pi))
      yc=np.cos(np.linspace(0.,2.*np.pi))
      ax.plot(xc,yc,color='black', zorder=1)
#      plt.show()


      onamea = './animation_files/spec_surfb{indx:05d}.png'.format(indx=indx)
      fig1.savefig(onamea, bbox_inches='tight')
      fig1.clf()
#     
      animation_files.append(onamea)
      animation_images.append(imageio.imread(onamea))
         
      indx=indx+1
#
   imageio.mimwrite('./animation_files/spec_surfb.gif', animation_images, fps=2, loop=1)

#
#---------------------------------------------------------------------------------
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


fnames=['../TRASH/spec_surface_00001.h5',
        '../TRASH/spec_surface_00002.h5',
        '../TRASH/spec_surface_00003.h5',
        '../TRASH/spec_surface_00004.h5',
        '../TRASH/spec_surface_00005.h5',
        '../TRASH/spec_surface_00006.h5',
        '../TRASH/spec_surface_00007.h5',
        '../TRASH/spec_surface_00008.h5',
        '../TRASH/spec_surface_00009.h5',
        '../TRASH/spec_surface_00010.h5',
        '../TRASH/spec_surface_00011.h5',
        '../TRASH/spec_surface_00012.h5',
        '../TRASH/spec_surface_00013.h5',
        '../TRASH/spec_surface_00014.h5',
        '../TRASH/spec_surface_00015.h5']


windx = 1
main(fnames=fnames,windx=windx,create_animation=False)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
