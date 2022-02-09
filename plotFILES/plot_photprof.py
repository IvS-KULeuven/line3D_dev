import sys
sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from lev_misc import readcol
#
def main(fnames=['../outputFILES/photprof_star01.dat','../outputFILES/photprof_star02.dat'], xlim=[30.,-30], ylim=[0.,1.1], oname='./ps_files/spec_photprof', windx=1):
#
   ymin=ylim[0]
   ymax=ylim[1]
   
   xmin=xlim[0]
   xmax=xlim[1]
#
#-----------------------------------------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   fig = plt.figure(windx,figsize=(xsize,ysize))
#
   plt.ion()
   plt.show()
#
   ax = fig.subplots()
#
   ax.set_xlabel(r'$x_{obs} [v_{th}]$')
   ax.set_ylabel(r'$F_{L}/F_{c}$')
   ax.set_ylim(ymin,ymax)
   ax.set_xlim(xmin,xmax)      
#
#--------------------------read everything------------------------------
#
   for fname in fnames:
      xobs, xic, xicc=readcol(fname,ncols=3,nskip=1)      
      ax.plot(xobs,xic/xicc)
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig.savefig(oname1, bbox_inches='tight')
   fig.savefig(oname2, bbox_inches='tight')   
#
#---------------------plot everything for a certain model--------------------
#
fnames=['../outputFILES/photprof_star01.dat',
        '../outputFILES/photprof_star02.dat']
#fnames=['modspec_model00.h5']

windx=1

main(fnames=fnames,xlim=[30.,-30.],ylim=[0.,1.1],windx=windx)



sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
