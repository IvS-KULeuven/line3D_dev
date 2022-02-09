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
#
def main(fname='../outputFILES/spec_triangles00001.dat', obs_indx=1, oname='./ps_files/spec_triangles',
         windx=1, xlim=[-3.,3.], ylim=[-3.,3.]):
#+
# NAME:
#	plot_triangles
#
# PURPOSE:
#	This procedure plots the tirangulation for the binary spec program
#
# CALLING SEQUENCE:
#
#	plot_triangles, dir
#
# INPUTS:
#	dir:	Directory, where appropriate files are stored.
#	
# KEYWORD PARAMETERS:
#       fname:  Set this keyword to the name of a file inside directory dir, 
#               from which the triangles shall be read in
#
#       obs_indx: Set this keyword to an integer, to describe from which file
#                 inside directory dir, the triangles are read in:
#                 e.g. obs_indx=1 reads from 'dir/spec_triangles_00001.dat'
#
#       oname:  Set this keyword to the output-name of a ps-file, if output
#               shall be plotted to that file
#
#       windx:  Set this keyword to an integer defining the window, to which
#               output is plotted
#
#       ylim:   Set this keyword to a 2-d array, which sets the yrange
#
#       xlim:   Set this keyword to a 2-d array, which sets the xrange
#
#       oplot:  Set this keyword (flag), to overplot in the current window
#               This keyword is not allowed, when oname is set
#
#       cindx:  Set this keyword to an integer, describing the color in which
#               the emergent profile is plotted (affects only overplotted
#               profiles)
#
#       isotropic: Set this keyword to make isotropic plot
#
#       help:   Set this keyword (flag) to show the documentation of the
#               procedure
#
# EXAMPLE:
#
#-
#
#------------------read all information from hdf5-file------------------
#
   y1,x1,y2,x2,y3,x3 = get_triangles(fname)
#
#-----------------------plotting----------------------------------------
#

   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./1.
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

   npoints=np.size(x1)

   phi=0.*np.pi/180.
   for i in range(0,npoints):
      px1 = np.cos(phi)*x1[i] - np.sin(phi)*y1[i]
      py1 = np.sin(phi)*x1[i] + np.cos(phi)*y1[i]

      px2 = np.cos(phi)*x2[i] - np.sin(phi)*y2[i]
      py2 = np.sin(phi)*x2[i] + np.cos(phi)*y2[i]

      px3 = np.cos(phi)*x3[i] - np.sin(phi)*y3[i]
      py3 = np.sin(phi)*x3[i] + np.cos(phi)*y3[i]

      ax.plot([px1,px2], [py1,py2],color='black', lw=0.3)
      ax.plot([px2,px3], [py2,py3],color='black', lw=0.3)
      ax.plot([px3,px1], [py3,py1],color='black', lw=0.3)      
      
#      ax.plot([x1[i],x2[i]], [y1[i],y2[i]],color='black', lw=0.3)
#      ax.plot([x2[i],x3[i]], [y2[i],y3[i]],color='black', lw=0.3)
#      ax.plot([x3[i],x1[i]], [y3[i],y1[i]],color='black', lw=0.3)      
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


fname='../outputFILES/spec_triangles_00001.dat'
oname='ps_files/spec_triangles01'
windx = 1
main(fname=fname,xlim=[-3.,3.], ylim=[-3.,3.],windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
