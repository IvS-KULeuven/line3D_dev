import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_modelspec import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from read_vbin import *
#
def main(fnames=['../outputFILES/FLUXEM_00001.dat'], phase=None, xlim=[1.2,-1.2], ylim=[0.,360.], xlim_wave=[6500.,6600.], clim=[0.,1.5], xscalefac=1., lambda0=6562.4, vth_fid=100., oname='./ps_files/dynspec', windx=1):
#
   ymin=ylim[0]
   ymax=ylim[1]
   
   xmin=xlim[0]
   xmax=xlim[1]

   xmin_wave=xlim_wave[0]
   xmax_wave=xlim_wave[1]
   
#prepare contour plots
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1.4
   ysize=xsize/aspect_ratio

   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1

   ax1 = fig1.subplots()   
   ax2 = fig2.subplots()

   plt.ion()
   plt.show()
   
   ax1.set_xlabel(r'$x_{obs} [km/s]$')
   ax1.set_ylabel(r'$\gamma$')   
   ax1.set_ylim(ymin,ymax)
   ax1.set_xlim(xmin,xmax)

   ax2.set_xlabel(r'$\lambda [A]$')
   ax2.set_ylabel(r'$\gamma$')
   ax2.set_ylim(ymin,ymax)   
   ax2.set_xlim(xmin_wave,xmax_wave)
#
#-----------------------------------------------------------------------
#
   nxobs = 301
   ngamma = np.size(fnames)
   xobs = np.linspace(xmin,xmax,nxobs)
   xobs_wave = lambda0/(vth_fid*1.e5/cgs_clight * xobs + 1.)
   fnorm2d = np.zeros(shape=(ngamma,nxobs))

   try:
      phase_set = phase != None
      phase_set = phase_set[0]
   except:
      phase_set = False
      
   if not phase_set:
      phase = np.zeros(ngamma)

   for j in range(0,ngamma):
      fname = fnames[j]
      nxobs_in, alpha_in, gamma_in, xobs_in, ftot_in, fcont_in = get_fluxem(fname)
      xobs_wave_in = lambda0/(vth_fid*1.e5/cgs_clight * xobs_in + 1.)      
      xobs_in=xobs_in*xscalefac

      if not phase_set:
         phase[j] = gamma_in*180./np.pi

      for i in range(0,nxobs):
         iim1, ii = find_indx(xobs[i], xobs_in, nxobs_in)
         fnorm_iim1 = ftot_in[iim1]/fcont_in[iim1]
         fnorm_ii = ftot_in[ii]/fcont_in[ii]
         fnorm2d[j][i] = interpol1d_2p_lin(fnorm_iim1, fnorm_ii, xobs_in[iim1], xobs_in[ii], xobs[i])

#      ax1.plot(xobs,ftot/fcont,label=r'$(\alpha,\gamma)=({alpha:.2f},{gamma:.2f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi))
#      ax2.plot(xobs_wave,ftot/fcont,label=r'$(\alpha,\gamma)=({alpha:.2f},{gamma:.2f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi))      

#   print(fnorm2d)
#   exit()
   ax1.set_ylabel(r'$F_{tot}/F_{c}$')
#
#-------------------temperature-----------------------------------------
#
   cmin=clim[0]
   cmax=clim[1]
   dcol=(cmax-cmin)/99.
   clevels=np.arange(cmin,cmax+dcol,dcol)

   contourplot1= ax1.contourf(xobs, phase, fnorm2d,
                  levels=clevels,
                  extend='both',
                  cmap='seismic')
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig1.colorbar(contourplot1,ax=ax1,orientation='vertical')
   cbar.set_label(r'$F_{tot}/F_{c}$')
#
#------------------------output to files-------------------------------
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')


#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
#
fnames=['../outputFILES/FLUXEM_00001.dat',
        '../outputFILES/FLUXEM_00002.dat',
        '../outputFILES/FLUXEM_00003.dat',
        '../outputFILES/FLUXEM_00004.dat',
        '../outputFILES/FLUXEM_00005.dat',
        '../outputFILES/FLUXEM_00006.dat',
        '../outputFILES/FLUXEM_00007.dat',
        '../outputFILES/FLUXEM_00008.dat',
        '../outputFILES/FLUXEM_00009.dat',
        '../outputFILES/FLUXEM_00010.dat',
        '../outputFILES/FLUXEM_00011.dat',
        '../outputFILES/FLUXEM_00012.dat',
        '../outputFILES/FLUXEM_00013.dat',
        '../outputFILES/FLUXEM_00014.dat',
        '../outputFILES/FLUXEM_00015.dat',
        '../outputFILES/FLUXEM_00016.dat',
        '../outputFILES/FLUXEM_00017.dat']

windx = 1

vthfid=100.
vinf=1.
#main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[400.,-400.],xlim_wave=[6550.,6570.],clim=[0.2,1.8],vth_fid=vthfid,windx=windx)




windx = 1

vthfid=100.
vinf=1.

base_dir = 'photprof0_sline0/'
base_dir = 'photprof0_sline0_highres2/'
base_dir = 'reference_test/'
base_dir = 'photprof5_sline1_vmicro20/'
base_dir = 'photprof5_sline1_vmicro20_kline1d-2/'
base_dir = 'photprof5_sline1_vmicro20_kline1d-5/'
base_dir = 'photprof5_sline1_vmicro20_kline1d-3/'
#base_dir = 'photprof0_sline0_vmicro20_kline1d0/'
phases = np.linspace(np.pi/2, 3./2.*np.pi, 17)
phases = np.linspace(3./2.*np.pi, 1./2.*np.pi, 17)


#Argument of periapsis [deg].
omega = 27.4 #297.4
nd=33
#phases = np.linspace(omega+180.-90., omega+180.+90., nd)
phases = np.linspace(0., 360., nd)
phases = phases*np.pi/180.

phases_time = np.zeros(nd)
eccentricity = 0.23
for i in range(0,nd):
   phases_time[i] = phase_to_time(phases[i], eccentricity)
#phases = phases*180./np.pi
phases = phases_time

#base_dir = 'photprof5_sline0/'
#phases = np.linspace(0.,2.*np.pi, 17)*180./np.pi

fnames=['../outputFILES/models_olivier/'+base_dir+'phase001/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase002/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase003/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase004/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase005/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase006/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase007/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase008/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase009/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase010/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase011/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase012/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase013/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase014/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase015/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase016/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase017/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase018/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase019/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase020/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase021/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase022/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase023/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase024/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase025/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase026/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase027/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase028/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase029/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase030/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase031/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase032/FLUXEM_00001.dat',
        '../outputFILES/models_olivier/'+base_dir+'phase033/FLUXEM_00001.dat']

main(fnames=fnames, phase=phases, xscalefac=vthfid/vinf,xlim=[300.,-300.],xlim_wave=[6550.,6570.],clim=[0.,2.],ylim=[1.,0.], vth_fid=vthfid,windx=windx)





sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
