import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from astropy.io import fits
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline
#
def calc_continuum(lambd,ftot):
#   
#smooth everything to get the normalization
   wgauss=1  #width of gaussian
   ngauss=5*wgauss #number of pixels for each gaussian
   xgauss=np.linspace(-(ngauss-1)/2,(ngauss-1)/2,ngauss,dtype='int')  #indices used
   kgauss=1./np.sqrt(np.pi)/np.float(wgauss) * np.exp(-(np.double(xgauss)/np.float(wgauss))**2)

   nd=np.size(ftot)
   ftot_smooth=np.copy(ftot)
   for i in np.arange(np.int((ngauss-1)/2),np.int(nd-(ngauss-1)/2)):
      indx = i + xgauss
      ftot_smooth[i] = np.sum(ftot[indx]*kgauss)
#
#calculate first derivative and second derivative
   dftot1_smooth=np.zeros(nd)
   dftot2_smooth=np.zeros(nd)
   dftot3_smooth=np.zeros(nd)   
   h=lambd[1]-lambd[0]
   for i in np.arange(np.int(ngauss),np.int(nd-ngauss)):
      dftot1_smooth[i] = np.abs((0.5*ftot_smooth[i+1]-0.5*ftot_smooth[i-1])/h)
      dftot2_smooth[i] = np.abs((ftot_smooth[i-1]-2.*ftot_smooth[i] + ftot_smooth[i+1])/h**2)
      dftot3_smooth[i] = np.abs((-0.5*ftot_smooth[i-2]+ftot_smooth[i-1]-ftot_smooth[i+1]+0.5*ftot_smooth[i+2])/h**3)

   dftot1_smooth=dftot1_smooth/np.sum(dftot1_smooth)
   dftot2_smooth=dftot2_smooth/np.sum(dftot2_smooth)
   dftot3_smooth=dftot3_smooth/np.sum(dftot3_smooth)
   dftot12_smooth = 0.1*dftot1_smooth + 0.5*dftot2_smooth + 0.4*dftot3_smooth

   indx=np.where(dftot12_smooth > 0.)
   df12_max=np.max(dftot12_smooth[indx])
   df12_min=np.min(dftot12_smooth[indx])
   indxl=np.arange(0,np.int((ngauss-1)))
   indxr=np.arange(np.int(nd-(ngauss-1)),nd)
   dftot12_smooth[indxl]=df12_max
   dftot12_smooth[indxr]=df12_max

   deltaf12 = df12_max-df12_min
   fac=np.double(1.e-3)
   for i in np.arange(0,20):
      threshold12 = df12_min + np.double(fac)*deltaf12
      indx = np.where(dftot12_smooth <= threshold12)
      if np.size(indx) > 100:
         break
      else:
         fac=fac*1.9
   if(np.size(indx) < 10):
      print('error: indx smaller 10')
      exit()
#
   lambd_res = lambd[indx]
   fcont_res = ftot[indx]

#
#calculate global continuum by linear regression
   res = stats.linregress(lambd_res, fcont_res)
   fcont = res.intercept + res.slope*lambd
   fnorm = ftot/fcont
#
#or calculate a spline
#   spl = InterpolatedUnivariateSpline(lambd, fcont)
#   fcont = spl(lambd)
#   fnorm = ftot/fcont

#   fig = plt.figure(1)
#   ax = fig.subplots(6)
#   ax[0].plot(lambd,ftot)   
#   ax[0].plot(lambd,ftot_smooth)
#   ax[0].plot(lambd,fcont)
#   ax[0].scatter(lambd_res,fcont_res,color='red',marker='.')   
#   ax[1].plot(lambd,dftot1_smooth)
#   ax[2].plot(lambd,dftot2_smooth)
#   ax[3].plot(lambd,dftot3_smooth)                                            
#   ax[4].plot(lambd,dftot12_smooth)
#   ax[4].plot([np.min(lambd),np.max(lambd)],[threshold12,threshold12])
#   ax[5].plot(lambd,fnorm)   
#   plt.ion()
#   plt.show()   
#   sdum = input("Press [q] to exit.")
#   if sdum == 'q':
#       exit()


   return fnorm
#  


def main(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[1.2,-1.2], ylim=[0.,1.5], show_plot=False, oname='./ps_files/fluxem', windx=1):
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

   for fname in fnames:
      print('working on ', fname)
      hdul = fits.open(fname)
#      hdul.info()
      header = hdul[0].header
      crval1 = np.double(header['CRVAL1'])
      cdelt1 = np.double(header['CDELT1'])
#
      ftot=hdul[0].data
      nd = np.size(ftot)
      indx = np.arange(0,nd)
      lambd = crval1 + cdelt1*indx

      indx = np.where((lambd >= lam1_halpha) & (lambd <= lam2_halpha))   
      lambd_halpha = lambd[indx]
      ftot_halpha = ftot[indx]

      indx = np.where((lambd >= lam1_hbeta) & (lambd <= lam2_hbeta))   
      lambd_hbeta = lambd[indx]
      ftot_hbeta = ftot[indx]
   
      indx = np.where((lambd >= lam1_hgamma) & (lambd <= lam2_hgamma))   
      lambd_hgamma = lambd[indx]
      ftot_hgamma = ftot[indx]
#
      fnorm_halpha = calc_continuum(lambd_halpha,ftot_halpha)
      fnorm_hbeta = calc_continuum(lambd_hbeta,ftot_hbeta)
      fnorm_hgamma = calc_continuum(lambd_hgamma,ftot_hgamma)


#output to file
      fout_hbeta=fname[0:19]+'_hbeta.dat'      
      fout_hgamma=fname[0:19]+'_hgamma.dat'
      fout_halpha=fname[0:19]+'_halpha.dat'


      if np.size(lambd_halpha) != np.size(lambd_hbeta):
         print('error: nd(halpha) ne nd(hbeta)')
         exit()
      if np.size(lambd_halpha) != np.size(lambd_hgamma):
         print('error: nd(halpha) ne nd(hgamma)')
         exit()

      file_halpha = open(fout_halpha, "w")
      file_hbeta = open(fout_hbeta, "w")
      file_hgamma = open(fout_hgamma, "w")
   
      file_halpha.write('{nd:10.0f}\n'.format(nd=np.size(lambd_halpha)))
      file_hbeta.write('{nd:10.0f}\n'.format(nd=np.size(lambd_hbeta)))
      file_hgamma.write('{nd:10.0f}\n'.format(nd=np.size(lambd_hgamma)))   
      for i in range(0,np.size(lambd_halpha)):
         file_halpha.write('{lambd:20.4f}{flux:20.6f}\n'.format(lambd=lambd_halpha[i],flux=fnorm_halpha[i]))
         file_hbeta.write('{lambd:20.4f}{flux:20.6f}\n'.format(lambd=lambd_hbeta[i],flux=fnorm_hbeta[i]))
         file_hgamma.write('{lambd:20.4f}{flux:20.6f}\n'.format(lambd=lambd_hgamma[i],flux=fnorm_hgamma[i]))            
      file_halpha.close()
      file_hbeta.close()
      file_hgamma.close()   

#------------------------------------------------------------------------------
   if show_plot:   
      xsize=18. #in cm
      xsize=xsize/2.54 #in inches
      aspect_ratio=1./.8
      ysize=xsize/aspect_ratio
      fig = plt.figure(windx,figsize=(xsize,ysize))
#
      plt.ion()
      plt.show()
#
      ax = fig.subplots(2)
#
      ax[0].set_xlabel(r'$\lambda [A]$')
      ax[0].set_ylabel(r'$F_{tot}$')
#
      ax[0].plot(lambd_halpha,ftot_halpha)
      ax[0].plot(lambd_hbeta,ftot_hbeta)
      ax[0].plot(lambd_hgamma,ftot_hgamma)      
#   
      ax[1].plot(lambd_halpha,fnorm_halpha)
      ax[1].plot(lambd_hbeta,fnorm_hbeta)
      ax[1].plot(lambd_hgamma,fnorm_hgamma)         
      ax.plot(lambd*10.,planck)


#
#---------------------plot everything for a certain model--------------------
#
#fnames=['t04250_g+1.5_m10p00_hrplc.fits']
#main(fnames=fnames,windx=windx)


fnames=[]
f = open('filelist.dat', 'r')
for line in f.readlines():
   line = line.strip()
   columns = line.split()
   for entry in columns:
      fnames.append(entry)
f.close()

main(fnames=fnames,xlim=[4000.,7000.],ylim=[0.,6.],show_plot=False)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
