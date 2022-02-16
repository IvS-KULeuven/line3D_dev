import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from astropy.io import fits
#
def main(fnames=['../outputFILES/FLUXEM_00001.dat'], xlim=[1.2,-1.2], ylim=[0.,1.5], oname='./ps_files/photprof', windx=1,show_plot=True):
#
   if show_plot:
      xsize=18. #in cm
      xsize=xsize/2.54 #in inches
      aspect_ratio=1./0.3
      ysize=xsize/aspect_ratio
      fig = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
#
      plt.ion()
      plt.show()
#
      ax = fig.subplots(1,4)

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
      hdul = fits.open(fname)
      header = hdul[0].header
#      hdul.info()
#      print(header)
      fall=hdul[0].data
      ftot = fall[1][:]
      fnorm = fall[0][:]   
      nd = np.size(ftot)
      lambd = np.linspace(3000.-8.39,18000.,nd,dtype='float64',endpoint=True)

      indx = np.where((lambd >= lam1_halpha) & (lambd <= lam2_halpha))   
      lambd_halpha = lambd[indx]
      fnorm_halpha = fnorm[indx]

      indx = np.where((lambd >= lam1_hbeta) & (lambd <= lam2_hbeta))   
      lambd_hbeta = lambd[indx]
      fnorm_hbeta = fnorm[indx]
   
      indx = np.where((lambd >= lam1_hgamma) & (lambd <= lam2_hgamma))   
      lambd_hgamma = lambd[indx]
      fnorm_hgamma = fnorm[indx]

#output to file
      fout_halpha=fname[:-8]+'_halpha.ms.dat'
      fout_hbeta=fname[:-8]+'_hbeta.ms.dat'      
      fout_hgamma=fname[:-8]+'_hgamma.ms.dat'

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

      if show_plot:
#
         ax[0].set_xlabel(r'$\lambda [A]$')
         ax[0].set_ylabel(r'$F_{tot}$')
         ax[0].set_xscale('linear')
         ax[0].set_yscale('linear')   
         ax[0].set_xlim(xmin,xmax)

         ax[1].set_xlabel(r'$\lambda [A]$')
         ax[1].set_ylabel(r'$F_{tot}$')
         ax[1].set_xscale('linear')
         ax[1].set_yscale('linear')
   
         ax[2].set_xlabel(r'$\lambda [A]$')
         ax[2].set_ylabel(r'$F_{tot}$')
         ax[2].set_xscale('linear')
         ax[2].set_yscale('linear')

         ax[3].set_xlabel(r'$\lambda [A]$')
         ax[3].set_ylabel(r'$F_{tot}$')
         ax[3].set_xscale('linear')
         ax[3].set_yscale('linear')   
#
         ax[0].plot(lambd,fnorm)
         ax[1].plot(lambd_halpha,fnorm_halpha)
         ax[2].plot(lambd_hbeta,fnorm_hbeta)
         ax[3].plot(lambd_hbeta,fnorm_hbeta)
         
   if show_plot:
      oname1 = oname+'.png'
      oname2 = oname+'.ps'
      fig.savefig(oname1, bbox_inches='tight')
      fig.savefig(oname2, bbox_inches='tight')   
#
#---------------------plot everything for a certain model--------------------
#
fnames=['6000_30_m25p00.ms.fits']
main(fnames=fnames,xlim=[4000.,7000.],ylim=[0.,6.],show_plot=True)


#fnames=[]
#f = open('filelist.dat', 'r')
#for line in f.readlines():
#   line = line.strip()
#   columns = line.split()
#   for entry in columns:
#      print(entry)
#      fnames.append(entry)
#f.close()
#
#main(fnames=fnames,xlim=[4000.,7000.],ylim=[0.,6.],show_plot=False)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
