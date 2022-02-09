import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_vbin import *
#
def main(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[1.2,-1.2], xlim_wave=[6500.,6600.], ylim=[0.,1.5], xscalefac=1., lambda0=6562.4, vth_fid=100., oname='./ps_files/fluxem', windx=1, snr=None, llegend=True):
#
   ymin=ylim[0]
   ymax=ylim[1]
   
   xmin=xlim[0]
   xmax=xlim[1]

   xmin_wave=xlim_wave[0]
   xmax_wave=xlim_wave[1]   
#
#-----------------------------------------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   
   fig1 = plt.figure(windx,figsize=(xsize,ysize))
   ax1 = fig1.subplots()

   fig2 = plt.figure(windx+1,figsize=(xsize,ysize))
   ax2 = fig2.subplots()      
#
   plt.ion()
   plt.show()
#
#
#   ax1.set_xlabel(r'$x_{obs} [v_\infty]$')
   ax1.set_xlabel(r'$x_{obs} [km/s]$')
   ax1.set_ylabel(r'$F_{tot}/F_{c}$')
   ax1.set_ylim(ymin,ymax)
   ax1.set_xlim(xmin,xmax)

   ax2.set_xlabel(r'$\lambda [A]$')
   ax2.set_ylabel(r'$F_{tot}/F_{c}$')
   ax2.set_ylim(ymin,ymax)   
   ax2.set_xlim(xmin_wave,xmax_wave)
#
#--------------------------read everything------------------------------
#
   for fname in fnames:
      nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname, snr=snr)
      xobs_wave = lambda0/(vth_fid*1.e5/cgs_clight * xobs + 1.)      
      xobs=xobs*xscalefac
      ax1.plot(xobs,ftot/fcont,label=r'$(\alpha,\gamma)=({alpha:.2f},{gamma:.2f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi))
      ax2.plot(xobs_wave,ftot/fcont,label=r'$(\alpha,\gamma)=({alpha:.2f},{gamma:.2f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi))

   if llegend:
      ax1.legend(ncol=2)
      ax2.legend(ncol=2)   
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

   oname1 = oname+'_lambda.png'
   oname2 = oname+'_lambda.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')      
#
#---------------------plot everything for a certain model--------------------
#
#fnames=['../outputFILES/FLUXEM_00001.dat',
#        '../outputFILES/FLUXEM_00005.dat',
#        '../outputFILES/FLUXEM_00009.dat',
#        '../outputFILES/FLUXEM_00013.dat']
#
#windx = 1
#
#vthfid=100.
#vinf=1.
#main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[350.,-350.],xlim_wave=[6550.,6570.],ylim=[0.,6.],vth_fid=vthfid,windx=windx)


#fnames=['../outputFILES/nico_wr3d/FLUXEM_00001.dat',
#       '/lhome/levin/Postdoc/organisation/plots/models/wr3d/spherical3d/clumped/kline5d12/vmicro100/FLUXEM_00001.dat']

#fnames=['/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE/snap17_kl1d13_vmicro1d2_sresol/FLUXEM_00001.dat',
#        '/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_0.50/snap40_kl1d13_vmicro1d1_sresol/FLUXEM_00001.dat',
#        '../outputFILES/nico_wr3d/FLUXEM_00001.dat']

#fnames=['/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p02s12xa3_halpha_vrot2/FLUXEM_00001.dat',
#        '../outputFILES/FLUXEM_00001.dat',
#        '../outputFILES/FLUXEM_00005.dat',
#        '../outputFILES/FLUXEM_00009.dat']#
#
#windx = 1
#
#vthfid=100.
#vinf=1.
#main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[4000.,-4000.],lambda0=5696.,xlim_wave=[5600.,5800.],ylim=[0.,4.],vth_fid=vthfid,windx=windx)



vthfid=100.
vinf=1.
windx=1
fnames=['../outputFILES/nico_wr3d/opt_vlat0_kl1d13/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/opt_vlat0_kl5d13/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/opt_vlat0_kl1d14/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/opt_vlat1_kl1d13/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/FLUXEM_00001.dat']


fnames=['../outputFILES/nico_wr3d/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/opt_vlat1_lte/FLUXEM_00001.dat',
        '../outputFILES/nico_wr3d/opt_vlat0_kl5d13/FLUXEM_00001.dat',        
        '../outputFILES/nico_wr3d/opt_vlat0_kl1d14/FLUXEM_00001.dat']

fnames=['../outputFILES/nico_wr3d/FLUXEM_00001.dat',
        '/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap40_kl1d13_vmicro1d1_sresol/FLUXEM_00001.dat']

#main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[4000.,-4000.],lambda0=5696.,xlim_wave=[5600.,5800.],ylim=[0.,4.],vth_fid=vthfid,windx=windx)


#fnames=['../outputFILES/FLUXEM_00001.dat',
#        '/STER/levin/models_florian/model_ken_0G_highres/kline1d4/vmicro1d2/snap0080/FLUXEM_00001.dat']
#main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[8000.,-8000.],lambda0=1548.,xlim_wave=[1000.,2000.],ylim=[0.,4.],vth_fid=vthfid,windx=windx, snr=50.)

base_dir = 'photprof5_sline0/'
base_dir = 'photprof0_sline0/'
#base_dir = 'photprof0_sline0_highres/'
base_dir = ''
base_dir = 'photprof5_sline1_vmicro20_kline1d-2/'

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
        '../outputFILES/models_olivier/'+base_dir+'phase017/FLUXEM_00001.dat']

#fnames=['../outputFILES/models_olivier/reference_test/phase001/FLUXEM_00001.dat',
#        '../outputFILES/models_olivier/photprof5_sline1_vmicro20_kline1d-2/phase001/FLUXEM_00001.dat',]

        

main(fnames=fnames, xscalefac=vthfid/vinf, xlim=[400.,-400.], lambda0=6562.4,xlim_wave=[6550.,6570.],ylim=[0.,1.5],vth_fid=100.,windx=windx, llegend=False)



sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
