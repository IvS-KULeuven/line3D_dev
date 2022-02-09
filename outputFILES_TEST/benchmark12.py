import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_models import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
#
def main(dir='.', windx=1, oname='./ps_files/benchmark12', xlim=[0.,10.]):
#
#plots contours of searchlight along a direction




#   fig0 = plt.figure(windx)
#   ax0 = fig0.subplots()
#   ax0.scatter([0.,1.,2.],[0.,1.,2.], color='blue', marker='.', s=1, label='3D SC')
#   plt.show()   
#   sdum = input("Press [q] to exit.")
#   if sdum == 'q':
#       exit()   
#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/benchmark12.h5'

   xic1, kcont, eps_cont, nx, ny, nz, nr, x, y, z, r, \
      itmaxc, epsmaxc_sc, epsmaxc_fvm, t1d_jo, opac1d_jo, \
      mint1d_joray, mint1d_jomom, \
      fcont1d_joray, fcont1d_jomom, \
      mask3d, t3d, opac3d, mint3d_sc, mint3d_fvm, \
      fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, \
      fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm, \
      kcontrr3d_sc, kcontthth3d_sc, kcontphiphi3d_sc, \
      kcontrth3d_sc, kcontrphi3d_sc, kcontthphi3d_sc, \
      kcontrr3d_fvm, kcontthth3d_fvm, kcontphiphi3d_fvm, \
      kcontrth3d_fvm, kcontrphi3d_fvm, kcontthphi3d_fvm = read_benchmark12(fname)

   rstar=19.*cgs_rsu
#
#normalize everything to xic1
   mint3d_sc=mint3d_sc/xic1
   mint3d_fvm=mint3d_fvm/xic1
   mint1d_joray=mint1d_joray/xic1
   mint1d_jomom=mint1d_jomom/xic1
   fcont1d_joray=fcont1d_joray/xic1
   fcont1d_jomom=fcont1d_jomom/xic1
   fcontr3d_sc=fcontr3d_sc/xic1
   fcontth3d_sc=fcontth3d_sc/xic1
   fcontphi3d_sc=fcontphi3d_sc/xic1
   fcontr3d_fvm=fcontr3d_fvm/xic1
   fcontth3d_fvm=fcontth3d_fvm/xic1
   fcontphi3d_fvm=fcontphi3d_fvm/xic1
   kcontrr3d_sc=kcontrr3d_sc/xic1
   kcontthth3d_sc=kcontthth3d_sc/xic1
   kcontphiphi3d_sc=kcontphiphi3d_sc/xic1
   kcontrth3d_sc=kcontrth3d_sc/xic1
   kcontrphi3d_sc=kcontrphi3d_sc/xic1
   kcontthphi3d_sc=kcontthphi3d_sc/xic1   
   kcontrr3d_fvm=kcontrr3d_fvm/xic1
   kcontthth3d_fvm=kcontthth3d_fvm/xic1
   kcontphiphi3d_fvm=kcontphiphi3d_fvm/xic1
   kcontrth3d_fvm=kcontrth3d_fvm/xic1
   kcontrphi3d_fvm=kcontrphi3d_fvm/xic1
   kcontthphi3d_fvm=kcontthphi3d_fvm/xic1   

#calculate optical depth along z-axis
   nradial=1
   for k in np.arange(nz-1,-1,-1):
      if z[k] == 1.:
         break
      nradial=nradial+1
   opac1d_radial=np.zeros(nradial)
   r1d_radial=np.zeros(nradial)
   tau1d_radial=np.zeros(nradial)
   for k in np.arange(0,nradial):
      kk = nz-nradial+k
      r1d_radial[k] = z[kk]
      opac1d_radial[k] = opac3d[kk][int((ny-1)/2)][int((nx-1)/2)]

   for k in np.arange(nradial-2,-1,-1):
      dtau = 0.5*(opac1d_radial[k+1]+opac1d_radial[k])*(r1d_radial[k+1]-r1d_radial[k])*rstar
      tau1d_radial[k] = tau1d_radial[k+1]+dtau
   
#define radius and calculate mint_joray on 3d grid to obtain errors
   r3d=np.zeros(shape=(nz,ny,nx))
   mint3d_joray=np.zeros(shape=(nz,ny,nx))
   for i in  np.arange(0,nx):
      for j in np.arange(0,ny):
         for k in np.arange(0,nz):
            if mask3d[k][j][i] == 1:
               r3d[k][j][i] = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)
               iim1, ii = find_indx(r3d[k][j][i], r, nr)
               mint3d_joray[k][j][i] = interpol_yp(r[iim1], r[ii], mint1d_joray[iim1], mint1d_joray[ii], r3d[k][j][i])
            elif mask3d[k][j][i] == 2:
               r3d[k][j][i] = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)
               iim1, ii = find_indx(r3d[k][j][i], r, nr)
               mint3d_joray[k][j][i] = interpol_yp(r[iim1], r[ii], mint1d_joray[iim1], mint1d_joray[ii], r3d[k][j][i])
            elif mask3d[k][j][i] == 3:
               r3d[k][j][i] = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)
               iim1, ii = find_indx(r3d[k][j][i], r, nr)
               mint3d_joray[k][j][i] = interpol_yp(r[iim1], r[ii], mint1d_joray[iim1], mint1d_joray[ii], r3d[k][j][i])
            else:
               r3d[k][j][i] = 1.0 #dummy value
#
#***********************************************************************
#
#                       radial stratification
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1.4
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize), constrained_layout=True)
   windx=windx+1   
   plt.ion()
   plt.show()
#   ax0 = fig0.subplots(2,4)
   ax0 = fig0.subplots(4,1)  
#
#-----------------------------------------------------------------------
#
#   ax0[0,0].set_xlabel(r'$R/R_\ast$')
#   ax0[0,0].set_ylabel(r'$J\cdot r^2)$')
#   ax0[0,0].plot(r,mint1d_joray*r**2, color='black', linestyle='solid', label='1d (ray by ray)')
#   ax0[0,0].plot(r,mint1d_jomom*r**2, color='black', linestyle='solid', label='1d (moments)')
#   ax0[0,0].scatter(r3d,mint3d_fvm*r3d**2, color='red', marker='.', s=1, label='3D FVM')   
#   ax0[0,0].scatter(r3d,mint3d_sc*r3d**2, color='blue', marker='.', s=1, label='3D SC')
#   ax0[0,0].legend()
#
#-----------------------------------------------------------------------
#
#   ax0[0,1].set_xlabel(r'$R/R_\ast$')
#   ax0[0,1].set_ylabel(r'$F_r\cdot r^2)$')
#   ax0[0,1].plot(r,fcont1d_joray*r**2, color='black', linestyle='solid', label='1d (ray by ray)')
#   ax0[0,1].plot(r,fcont1d_jomom*r**2, color='black', linestyle='solid', label='1d (moments)')
#   ax0[0,1].scatter(r3d,fcontr3d_fvm*r3d**2, color='red', marker='.', s=1, label='3D FVM')   
#   ax0[0,1].scatter(r3d,fcontr3d_sc*r3d**2, color='blue', marker='.', s=1, label='3D SC')
#   ax0[0,1].legend()
#
#-----------------------------------------------------------------------
#
#   indx=np.where(r3d >= 1.)
#   ax0[0,2].set_xlabel(r'$R/R_\ast$')
#   ax0[0,2].set_ylabel(r'$K_{rr} \cdot r^2)$')
#   ax0[0,2].scatter(r3d[indx],kcontrr3d_fvm[indx]*r3d[indx]**2, color='red', marker='.', s=1, label='3D FVM')   
#   ax0[0,2].scatter(r3d[indx],kcontrr3d_sc[indx]*r3d[indx]**2, color='blue', marker='.', s=1, label='3D SC')
#   ax0[0,2].legend()
#
#-----------------------------------------------------------------------
#
   indx=np.where(r3d >= 1.)
   ax0[0].set_ylim([0.,1.])
   ax0[0].set_xlabel(r'$R/R_\ast$')
#   ax0[0].set_xlabel(r'$1-1/r$')
   ax0[0].set_ylabel(r'$f_{rr}$')
#   ax0[0].scatter(r3d,kcontrr3d_fvm/mint3d_fvm, color='red', marker='.', s=1, label='3D FVM')   
   ax0[0].scatter(r3d[indx],kcontrr3d_sc[indx]/mint3d_sc[indx], color='blue', marker='.', s=1, label='3D SC')
#   ax0[0].legend()

   ax0[1].set_ylim([0.,1.])
   ax0[1].set_xlabel(r'$R/R_\ast$')
#   ax0[1].set_xlabel(r'$1-1/r$')
   ax0[1].set_ylabel(r'$f_{\theta \theta}$')
#   ax0[1].scatter(r3d,kcontrr3d_fvm/mint3d_fvm, color='red', marker='.', s=1, label='3D FVM')   
   ax0[1].scatter(r3d[indx],kcontthth3d_sc[indx]/mint3d_sc[indx], color='blue', marker='.', s=1, label='3D SC')
   ax0[1].legend()

   ax0[2].set_ylim([0.,1.])
   ax0[2].set_xlabel(r'$R/R_\ast$')
#   ax0[2].set_xlabel(r'$1-1/r$')
   ax0[2].set_ylabel(r'$f_{\phi \phi}$')
#   ax0[2].scatter(r3d,kcontrr3d_fvm/mint3d_fvm, color='red', marker='.', s=1, label='3D FVM')   
   ax0[2].scatter(r3d[indx],kcontphiphi3d_sc[indx]/mint3d_sc[indx], color='blue', marker='.', s=1, label='3D SC')
   ax0[2].legend()      

   ax0[3].set_xlabel(r'$R/R_\ast$')
   ax0[3].set_ylabel(r'$\tau_r$')
   ax0[3].plot(r1d_radial,tau1d_radial, color='black', linestyle='solid')
   ax0[3].plot([np.min(r1d_radial),np.max(r1d_radial)],[1.,1.], color='black', linestyle='dashed')   
 
#
#------------------------output to files-------------------------------
#
   oname1 = oname+'.png'
   oname2 = oname+'.ps'
   fig0.savefig(oname1, bbox_inches='tight')
   fig0.savefig(oname2, bbox_inches='tight')
#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
   

windx = 1
main(windx=windx)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
