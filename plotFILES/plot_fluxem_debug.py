import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_vbin import *
#
def main(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[1.2,-1.2], ylim=[0.,1.5], xscalefac=1., oname='./ps_files/fluxem', windx=1):
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
   ax.set_xlabel(r'$x_{obs} [v_\infty]$')
   ax.set_ylabel(r'$F_{tot}/F_{c}$')
   ax.set_ylim(ymin,ymax)
   ax.set_xlim(xmin,xmax)      
#
#--------------------------read everything------------------------------
#
   for fname in fnames:
      nxobs, alpha, gamma, xobs, ftot, fcont, femi, fabs = get_fluxem_debug(fname)
      xobs=xobs*xscalefac
      color=next(ax._get_lines.prop_cycler)['color']
      ax.plot(xobs,ftot/fcont,label=r'$(\alpha,\gamma)=({alpha:.2f},{gamma:.2f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color=color)
      ax.plot(xobs,femi/fcont,linestyle='dotted', color=color)
      ax.plot(xobs,fabs/fcont,linestyle='--', color=color)            

   ax.legend()
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig.savefig(oname1, bbox_inches='tight')
   fig.savefig(oname2, bbox_inches='tight')   
#
#---------------------plot everything for a certain model--------------------
#
fnames=['/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/smooth/kline1d12/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/smooth/kline5d12/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/smooth/kline1d13/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/smooth/kline5d13/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/clumped/kline1d12/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/clumped/kline5d12/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/clumped/kline1d13/vmicro100/FLUXEM_DEBUG_00001.dat',
        '/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/clumped/kline5d13/vmicro100/FLUXEM_DEBUG_00001.dat',
        '../outputFILES/nico_wr3d/FLUXEM_DEBUG_00001.dat']

windx = 1

fnames=['../outputFILES/FLUXEM_DEBUG_00001.dat']

fnames=['../outputFILES/FLUXEM_DEBUG_00001.dat',
        '/STER/levin/models_florian/model_ken_0G_highres/kline1d4/vmicro1d2/snap0080/FLUXEM_DEBUG_00001.dat']

vthfid=100.
vinf=1.
main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[8000.,-8000.],ylim=[0.,2.5],windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
