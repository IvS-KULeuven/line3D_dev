import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from astropy.io import fits
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
   lam0_halpha=6562.8
   lam0_hbeta=4861.3
   lam0_hgamma=4340.5

   lam1_halpha=lam0_halpha-90.
   lam2_halpha=lam0_halpha+90.

   lam1_hbeta=lam0_hbeta-90.
   lam2_hbeta=lam0_hbeta+90.
   
   lam1_hgamma=lam0_hgamma-90.
   lam2_hgamma=lam0_hgamma+90.   


#   fname = 't03800_g+5.5_m05p04_sed.fits'
   fname = 't04250_g+1.5_m10p00_sed.fits'
   hdul = fits.open(fname)
#   hdul.info()
   header = hdul[0].header
   crval1 = np.double(header['CRVAL1'])
   cdelt1 = np.double(header['CDELT1'])
   print(header)
   print(crval1, cdelt1)
#
   ftot=hdul[0].data
   nd = np.size(ftot)
   indx = np.arange(0,nd)
   lambd = crval1 + cdelt1*indx
   lambd = 10.**lambd
#   print(lambd)

#   log10(wavelength) = CRVAL1 + CDELT1 * pixel    
#   exit()  

   indx = np.where((lambd >= lam1_halpha) & (lambd <= lam2_halpha))   
   lambd_halpha = lambd[indx]
   ftot_halpha = ftot[indx]
      
#   planck = np.zeros(nd,dtype='float64')
#   planck = 2.*cgs_planck*cgs_clight**2 / (lambd*1.e-7)**5 / (np.exp((cgs_planck*cgs_clight)/(lambd*1.e-7*cgs_kb*4250.))-1.)
   
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
   ax.set_xlabel(r'$\lambda [A]$')
   ax.set_ylabel(r'$F_{tot}$')
#   ax.set_xscale('log')
#   ax.set_yscale('log')   
#   ax.set_ylim(ymin,ymax)
#   ax.set_xlim(xmin,xmax)      
#
   ax.plot(lambd,ftot)
#   ax.plot(lambd_halpha*10.,ftot_halpha)
#   ax.plot(lambd*10.,planck)                                                              
#
#---------------------plot everything for a certain model--------------------
#
fnames=['../outputFILES/FLUXEM_00001.dat',
        '../outputFILES/FLUXEM_00002.dat',
        '../outputFILES/FLUXEM_00003.dat',
        '../outputFILES/FLUXEM_00004.dat',
        '../outputFILES/FLUXEM_00005.dat']
windx = 1

vthfid=100.
vinf=1.
main(fnames=fnames,xscalefac=vthfid/vinf,xlim=[300.,-300.],ylim=[0.,6.],windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
