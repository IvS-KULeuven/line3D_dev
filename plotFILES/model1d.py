import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_models import *
from lev_fit import*
#
def main(dir='../inputFILES', oname='./ps_files/model1d', windx=1):
#
#--------------------------read model-----------------------------------
#   
   fname=dir+'/model1d.h5'
#
   nr, r, rho1d, velr1d, t1d, vth1d = read_model1d(fname)
#
   mdot1d = 4.*np.pi*r**2*velr1d*rho1d*cgs_yr/cgs_msu
#
#get fit for mdot1d
   nr_fit = 800
   r_fit = r[nr-nr_fit:nr]
   mdot1d_fit = mdot1d[nr-nr_fit:nr]
   sigma1d=1.+np.zeros(nr_fit)   
#   acoeff, bcoeff = fit_regrl(r_fit,mdot1d_fit,sigma1d)
#   mdot1d_fit = acoeff+bcoeff*r_fit
   acoeff = fit_regrc(r_fit,mdot1d_fit,sigma1d)   
   mdot1d_fit = acoeff + np.zeros(nr_fit)
   
#beta velocity law
   vmin=5.e5
   vinf=4100.e5
   beta=0.5
   bconst=1.-(vmin/vinf)**(1./beta)
   velr1d_bvel = vinf*(1.-bconst*r[0]/r)**beta
#
#-----------------------------------------------------------------------
#

   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./1.6
   ysize=xsize/aspect_ratio
   fig1 = plt.figure(windx,figsize=(xsize,ysize))
#
   plt.ion()
   plt.show()
#
   ax = fig1.subplots(5,1)
#
#density
   ax[0].set_xlabel(r'$r [cm]$')
   ax[0].set_ylabel(r'$\rho [g/cm]$')
   ax[0].set_yscale('log')
   ax[0].plot(r,rho1d)
#
#radial velocity
   ax[1].set_xlabel(r'$r [cm]$')
   ax[1].set_ylabel(r'$v_r [km/s]$')
   ax[1].set_yscale('linear')
   ax[1].plot(r,velr1d/1.e5)
   ax[1].plot(r,velr1d_bvel/1.e5)   
#
#temperature
   ax[2].set_xlabel(r'$r [cm]$')
   ax[2].set_ylabel(r'$T [kK]$')
   ax[2].set_yscale('linear')
   ax[2].plot(r,t1d/1.e3)
#
#vthermal
   ax[3].set_xlabel(r'$r [cm]$')
   ax[3].set_ylabel(r'$v_{th} [km/s]$')
   ax[3].set_yscale('linear')
   ax[3].plot(r,vth1d/1.e5)
#
#mass loss
   ax[4].set_xlabel(r'$r [cm]$')
   ax[4].set_ylabel(r'$\dot{M} [10^{-6} M_{sun}/yr]$')
   ax[4].set_yscale('linear')
   ax[4].plot(r,mdot1d*1.e6)
   ax[4].plot(r_fit,mdot1d_fit*1.e6)   
#
#
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

#
#---------------------plot everything for a certain model--------------------
#
dir='../inputFILES'
windx = 1
main(dir=dir,windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
