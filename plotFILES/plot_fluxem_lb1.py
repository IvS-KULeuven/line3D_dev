import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from const_cgs import *
from read_vbin import *
from lev_interpol1d import *
from lev_misc import readcol
from lev_fluxem import get_linecenter
#
def fluxem(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[1.2,-1.2], ylim=[0.,1.5], xscalefac=1., oname='../ps_files/fluxem', windx=1):
#
   ymin=ylim[0]
   ymax=ylim[1]
   
   xmin=xlim[0]
   xmax=xlim[1]
#
#-----------------------------------------------------------------------
#
#read observations (already postprocessed)
   xobs_shenar, fnorm1_shenar, fnorm2_shenar, fnorm3_shenar, \
      fnorm4_shenar, fnorm5_shenar, fnorm6_shenar, fnorm7_shenar, fnorm8_shenar \
      =readcol('/lhome/levin/Postdoc/papers/paperIII/pro_files/lb1_observations_halpha.dat',ncols=9,nskip=1)

   xobs_shenar=xobs_shenar*xscalefac
#
#-----------------------------------------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.7
   ysize=xsize/aspect_ratio
   fig = plt.figure(windx,figsize=(xsize,ysize))

   ax0 = fig.subplots()
#
   plt.ion()
   plt.show()
#
#   ax = fig.subplots(2,1)
#
   ax0.set_xlabel(r'$x_{\rm obs} [{\rm km/s}]$')
   ax0.set_ylabel(r'$F_{\rm tot}/F_{\rm c}$')
   ax0.set_ylim(ymin,ymax)
   ax0.set_xlim(xmin,xmax)
#
#--------------------------read everything------------------------------
#
   vth_fiducial=100.
#calculating equivalent width in range
   xobs_min=-3.0
   xobs_max=3.0   

   #first profile binary code
   fname=fnames[0]
   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
   ew0 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)
   xobs=xobs*xscalefac
   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='black',linewidth=2.)
   
   fname=fnames[1]
   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
   ew1 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
   xobs=xobs*xscalefac      
#   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='red')   

   fname=fnames[2]
   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
   ew2 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
   xobs=xobs*xscalefac      
#   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='blue')   

   fname=fnames[3]
   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
   ew3 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
   xobs=xobs*xscalefac      
   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='green')   

   fname=fnames[4]
   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
   ew4 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max) 
   xobs=xobs*xscalefac      
   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='red')

#   fname=fnames[5]
#   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
#   ew5 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
#   xobs=xobs*xscalefac      
##   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='cyan')

#   fname=fnames[6]
#   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
#   ew6 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
#   xobs=xobs*xscalefac      
#   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='green')
   
#   fname=fnames[7]
#   nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
#   ew7 = equivalent_width(xobs,ftot,fcont,xobs_min,xobs_max)   
#   xobs=xobs*xscalefac      
##   ax0.plot(xobs,ftot/fcont,label=r'$(i,\gamma)=({alpha:.1f},{gamma:.1f})$'.format(alpha=alpha*180./np.pi,gamma=gamma*180./np.pi), color='yellow')            
   ax0.legend()

#   titlestr=r'$\bar{{W}}_x={ew:.4f}$'.format(ew=(ew0+ew1+ew2+ew3+ew4)/5.)
   titlestr=r'$\bar{{W}}_x={ew:.4f}$'.format(ew=ew0)
   ax0.set_title(titlestr)

#   ax0.plot(xobs_shenar,fnorm1_shenar,color='black',linestyle='--')
   ax0.plot(xobs_shenar,fnorm2_shenar,color='black',linestyle='--', linewidth=1.5)
#   ax0.plot(xobs_shenar,fnorm3_shenar,color='black',linestyle='dotted')
#  ax0.plot(xobs_shenar,fnorm4_shenar,color='black',linestyle='dotted')
#   ax0.plot(xobs_shenar,fnorm5_shenar,color='black',linestyle='--')
#  ax0.plot(xobs_shenar,fnorm6_shenar,color='black',linestyle='dotted')
   ax0.plot(xobs_shenar,fnorm7_shenar,color='black',linestyle='dotted', linewidth=1.5)
#   ax0.plot(xobs_shenar,fnorm8_shenar,color='black',linestyle='dotted')         
   
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig.savefig(oname1, bbox_inches='tight')
   fig.savefig(oname2, bbox_inches='tight')
#
#---------------------------------------------------------------------------
#
def rv_curve(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[0.,360.], ylim=[-60.,60.], vorb1=100., oname='../ps_files/fluxem', windx=1):
#
   ymin=ylim[0]
   ymax=ylim[1]
   
   xmin=xlim[0]
   xmax=xlim[1]
#
#-----------------------------------------------------------------------
#
   vth_fiducial=1.e7
#
   nd=np.size(fnames)

   xcenter1 = np.zeros(nd)
   vcenter1 = np.zeros(nd)
   gamma_arr = np.zeros(nd)
#
   i=0
   for fname in fnames:
      nxobs, alpha, gamma, xobs, ftot, fcont = get_fluxem(fname)
      fnorm=ftot/fcont
      xctr,xobsl,xobsr,fnorml,fnormr=get_linecenter(xobs, fnorm, flim=0.2, xlim=2.)
      
      xcenter1[i] = xctr
      vcenter1[i] = -xctr*vth_fiducial/1.e5
      gamma_arr[i] = gamma*180./np.pi
#      print(i, fname, gamma_arr[i],vcenter1[i],xctr)
      i=i+1

#   exit()


#
#----------------------------theoretical curve--------------------
#
#observed h-alpha wing radial velocities
   valpha_obs = -11.2*np.sin(gamma_arr*np.pi/180.0)
   valpha_obs2 = -6.4*np.sin(gamma_arr*np.pi/180.0)
#observed b-star radial velocity
   vb_obs = 52.8*np.sin(gamma_arr*np.pi/180.0)

#for model
   vb_model = -vorb1*np.sin(gamma_arr*np.pi/180.0)*np.sin(alpha)   
#
#-----------------------------------------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.7
   ysize=xsize/aspect_ratio
   fig = plt.figure(windx,figsize=(xsize,ysize))

   ax0 = fig.subplots()
#
   plt.ion()
   plt.show()
#
#   ax = fig.subplots(2,1)
#
   ax0.set_xlabel(r'$\gamma$')
   ax0.set_ylabel(r'${\rm RV} [\rm{km/s}]$')
   ax0.set_ylim(ymin,ymax)
   ax0.set_xlim(xmin,xmax)

   titlestr=r'$i={alpha:.2f}$'.format(alpha=alpha*180./np.pi)
   ax0.set_title(titlestr)   
#
#--------------------------read everything------------------------------
#
   ax0.plot(gamma_arr, valpha_obs, color='black', linestyle='solid')
   ax0.plot(gamma_arr, valpha_obs2, color='black', linestyle='solid')
   ax0.plot(gamma_arr, vb_obs, color='black', linestyle='solid')
   
   ax0.plot(gamma_arr, vb_model, color='blue', linestyle=(0, (5, 10)), marker='x')
   ax0.plot(gamma_arr, vcenter1, color='red', linestyle=(0, (5, 10)), marker='x')
   
#   ax0.legend()
   
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig.savefig(oname1, bbox_inches='tight')
   fig.savefig(oname2, bbox_inches='tight')   
#
#---------------------------------------------------------------------------
#
def equivalent_width(xobs,ftot,fcont,xmin,xmax):
   ew = 0.
   nd=np.size(xobs)
   for i in range(1,nd):
      if xobs[i] >= xmin and xobs[i] <= xmax:
         ew = ew + 0.5*(ftot[i]/fcont[i] + ftot[i-1]/fcont[i-1] - 2.)*(xobs[i]-xobs[i-1])
   return ew
   
#
#---------------------plot everything for a certain model--------------------
#
#fnames=['../outputFILES/FLUXEM_00001.dat',
#        '../outputFILES/FLUXEM_00002.dat',
#       '../outputFILES/FLUXEM_00003.dat'
#        '/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p02s12xa3_halpha_vrot2/FLUXEM_00001.dat']
#        '/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p02s11xa3_halpha_new/FLUXEM_00053.dat']
fnames=['../outputFILES/FLUXEM_00001.dat',
         '../outputFILES/FLUXEM_00002.dat',
         '../outputFILES/FLUXEM_00003.dat',
         '/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p02s12xa3_halpha_vrot2/FLUXEM_00001.dat',
         '/lhome/levin/Postdoc/papers/paperIII/models/lb1/parameter_study/model_p02s12_mdisc1v1.06e-14_mdisc2v1.78e-22_tdisc1v15000_tdisc2v10000_vmicro1v100.0_vmicro2v50.0_slope1v1.5_slope2v2.5//FLUXEM_00001.dat']
#/lhome/levin/Postdoc/papers/paperIII/models/lb1/parameter_study/model_p02s12_mdisc1v1.06e-14_mdisc2v1.78e-22_tdisc1v8000_tdisc2v10000_vmicro1v100.0_vmicro2v50.0_slope1v1.5_slope2v2.5/

windx = 1
fluxem(fnames=fnames,xscalefac=100, xlim=[300.,-300.],ylim=[0.,6.],windx=windx,oname='./ps_files/fluxem_lb1_halpha')
#
#
#rv_curve(fnames=fnames, vorb1=-141.3,windx=windx,oname='../ps_files/rv_lb1_p02s11xa3_halpha')




sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
