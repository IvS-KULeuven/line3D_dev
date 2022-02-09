import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_models import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from lev_contour import *


def map3d(dir='../inputFILES', oname3d='./ps_files/model3d', rp=4., windx=1,
         clim_rho=[-12.,-8.], clim_temp=[10.e3,100.e3], clim_velr=[1.e2,4.e3], xlim=[0.,8.], ylim=[-4.,4.]):  



#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/model3d.h5'

   nr, ntheta, nphi, radius, theta, phi, \
      rho3d, velr3d, velth3d, velphi3d, \
      t3d, vth3d = read_model3d(fname)

#r in r_star
   radius=radius/radius[0]
#velocity in km/s
   velr3d=velr3d/1.e5
   velth3d=velth3d/1.e5
   velphi3d=velphi3d/1.e5
#
#-----------------------------------------------------------------------

   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1.6
   ysize=xsize/aspect_ratio

   fig1 = plt.figure(windx,figsize=(xsize,ysize))
   windx=windx+1
   plt.ion()
   plt.show()

   ax1=fig1.subplots()
#   ax1 = fig1.add_subplot(111, projection='3d')
#
#-----------------------------------------------------------------------
#
#   ax1.set_xlim(xlim)
#   ax1.set_ylim(ylim)
#   ax1.set_xlabel(r'$x$')
#   ax1.set_ylabel(r'$z$')   

#
   rho2d=np.zeros(shape=(ntheta,nphi))
   velr2d=np.zeros(shape=(ntheta,nphi))
   t2d=np.zeros(shape=(ntheta,nphi))
#
   x_coord=np.zeros(shape=(ntheta,nphi))
   y_coord=np.zeros(shape=(ntheta,nphi))
   z_coord=np.zeros(shape=(ntheta,nphi))
   theta2d = np.zeros(shape=(ntheta,nphi))
   phi2d = np.zeros(shape=(ntheta,nphi))

   iim1, ii = find_indx(rp, radius, nr)

   for k in np.arange(nphi):
      for j in np.arange(ntheta):
         theta2d[j][k] = theta[j]
         phi2d[j][k] = phi[k]
         x_coord[j][k] = rp*np.sin(theta[j])*np.cos(phi[k])
         y_coord[j][k] = rp*np.sin(theta[j])*np.sin(phi[k])
         z_coord[j][k] = rp*np.cos(theta[j])
         rho2d[j][k] = interpol1d_2p_lin(rho3d[k][j][iim1],rho3d[k][j][ii],radius[iim1],radius[ii],rp)
         t2d[j][k] = interpol1d_2p_lin(t3d[k][j][iim1],t3d[k][j][ii],radius[iim1],radius[ii],rp)
         velr2d[j][k] = interpol1d_2p_lin(velr3d[k][j][iim1],velr3d[k][j][ii],radius[iim1],radius[ii],rp)

   clevels, ticks = get_clevels(clim=clim_rho)
   cmap_lev = get_cmap('jet')
   ctitle=r'$\rho$'
   lrho2d=np.log10(rho2d)

   contourplot1 = ax1.contourf(theta2d, phi2d, lrho2d,
                  levels=clevels,
                  extend='both',
                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig1.colorbar(contourplot1,ax=ax1,orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)

########################################################################

def compare_model3d(fname1, fname2):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
   small_number=1.e-4
#
   nr_01, ntheta_01, nphi_01, radius_01, theta_01, phi_01, \
      rho3d_01, velr3d_01, velth3d_01, velphi3d_01, \
      t3d_01, vth3d_01 = read_model3d(fname1)

   nr_02, ntheta_02, nphi_02, radius_02, theta_02, phi_02, \
      rho3d_02, velr3d_02, velth3d_02, velphi3d_02, \
      t3d_02, vth3d_02 = read_model3d(fname2)

   print('{str1:10}{str2:25}{str3:25}'.format(str1='##',str2='model1',str3='model2'))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='nr',x1=nr_01,x2=nr_02))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='ntheta',x1=ntheta_01,x2=ntheta_02))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='nphi',x1=nphi_01,x2=nphi_02))

   for i in np.arange(nr_01):
      print('{str1:10}{x1:25.14f}{x2:25.14f}{x3:25.15f}'.format(str1='r(i)',x1=radius_01[i],x2=radius_02[i],x3=radius_01[i]-radius_02[i]))
      if np.abs(radius_01[i]-radius_02[i]) > small_number:
         print('error: radii different')
         exit()

   for i in np.arange(ntheta_01):
      print('{str1:10}{x1:25.14f}{x2:25.14f}{x3:25.15f}'.format(str1='theta(i)',x1=theta_01[i],x2=theta_02[i],x3=theta_01[i]-theta_02[i]))
      if np.abs(theta_01[i]-theta_02[i]) > small_number:
         print('error: thetas different')
         exit()

   for i in np.arange(nphi_01):
      print('{str1:10}{x1:25.14f}{x2:25.14f}{x3:25.15f}'.format(str1='phi(i)',x1=phi_01[i],x2=phi_02[i],x3=phi_01[i]-phi_02[i]))
      if np.abs(phi_01[i]-phi_02[i]) > small_number:
         print('error: phis different')
         exit()         

   for i in np.arange(nr_01):
      for j in np.arange(ntheta_01):
         for k in np.arange(nphi_01):
#            print(i,j,k)            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='rho(i,j,k)',x1=rho3d_01[k][j][i],x2=rho3d_02[k][j][i],x3=rho3d_01[k][j][i]-rho3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='velr(i,j,k)',x1=velr3d_01[k][j][i],x2=velr3d_02[k][j][i],x3=velr3d_01[k][j][i]-velr3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='velth(i,j,k)',x1=velth3d_01[k][j][i],x2=velth3d_02[k][j][i],x3=velth3d_01[k][j][i]-velth3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='velphi(i,j,k)',x1=velphi3d_01[k][j][i],x2=velphi3d_02[k][j][i],x3=velphi3d_01[k][j][i]-velphi3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='t(i,j,k)',x1=t3d_01[k][j][i],x2=t3d_02[k][j][i],x3=t3d_01[k][j][i]-t3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='vth(i,j,k)',x1=vth3d_01[k][j][i],x2=vth3d_02[k][j][i],x3=vth3d_01[k][j][i]-vth3d_02[k][j][i]))
#            print('')
                        
            if np.abs(rho3d_01[k][j][i]-rho3d_02[k][j][i]) > small_number:
               print('error: rhos different')
               exit()         
            if np.abs(velr3d_01[k][j][i]-velr3d_02[k][j][i]) > small_number:
               print('error: velrs different')
               exit()   
            if np.abs(velth3d_01[k][j][i]-velth3d_02[k][j][i]) > small_number:
               print('error: velths different')
#               exit()
            if np.abs(velphi3d_01[k][j][i]-velphi3d_02[k][j][i]) > small_number:
               print('error: velphis different')
#               exit()
            if np.abs(t3d_01[k][j][i]-t3d_02[k][j][i]) > small_number:
               print('error: ts different')
               exit()         
            if np.abs(vth3d_01[k][j][i]-vth3d_02[k][j][i]) > small_number:
               print('error: vths different')
               exit()         

   exit()   


########################################################################
   
def main(dir='../inputFILES', oname3d='./ps_files/model3d', thetac=0., phic=0., windx=1,
         clim_rho=[-12.,-8.], clim_temp=[10.e3,100.e3], clim_velr=[1.e2,4.e3], xlim=[0.,8.], ylim=[-4.,4.], rstar=1.):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/model3d.h5'

   nr, ntheta, nphi, radius, theta, phi, \
      rho3d, velr3d, velth3d, velphi3d, \
      t3d, vth3d = read_model3d(fname)

#r in r_star
   radius=radius/radius[0]
#velocity in km/s
   velr3d=velr3d/1.e5
   velth3d=velth3d/1.e5
   velphi3d=velphi3d/1.e5
#
#***********************************************************************
#
#                       radial stratification
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=.4
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize), constrained_layout=True)
   windx=windx+1   
   plt.ion()
   plt.show()
   ax0 = fig0.subplots(4,1)

   rho1d=np.zeros(nr)
   velr1d=np.zeros(nr)
   velth1d=np.zeros(nr)
   velphi1d=np.zeros(nr)
   t1d=np.zeros(nr)
#
   jjm1, jj = find_indx(thetac, theta, ntheta)
   kkm1, kk = find_indx(phic, phi, nphi)
#   thetac = theta[jjm1]
#   print(theta[jjm1]*180./np.pi, theta[jj]*180./np.pi, thetac*180./np.pi)
#   exit()
#
   for i in range(0,nr):
      rho1d[i] = interpol2d_4p_lin(rho3d[kkm1,jjm1,i], rho3d[kkm1][jj][i], rho3d[kk][jjm1][i], rho3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      t1d[i] = interpol2d_4p_lin(t3d[kkm1,jjm1,i], t3d[kkm1][jj][i], t3d[kk][jjm1][i], t3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velr1d[i] = interpol2d_4p_lin(velr3d[kkm1,jjm1,i], velr3d[kkm1][jj][i], velr3d[kk][jjm1][i], velr3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velth1d[i] = interpol2d_4p_lin(velth3d[kkm1,jjm1,i], velth3d[kkm1][jj][i], velth3d[kk][jjm1][i], velth3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velphi1d[i] = interpol2d_4p_lin(velphi3d[kkm1,jjm1,i], velphi3d[kkm1][jj][i], velphi3d[kk][jjm1][i], velphi3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
#
#-----------------------------------------------------------------------
#
   titlestr=r'at $\theta, \phi=$({sth:.2f},{sphi:.2f})'.format(sth=thetac*180./np.pi,sphi=phic*180./np.pi)
#
   ax0[0].set_xlabel(r'$R/R_\ast$')
   ax0[0].set_ylabel(r'$\log(\rho)$')
   ax0[0].set_title(titlestr)
   ax0[0].plot(radius,np.log10(rho1d))

   ax0[1].set_xlabel(r'$R/R_\ast$')
   ax0[1].set_ylabel(r'$\log(T)$')
   ax0[1].plot(radius,np.log10(t1d))

   ax0[2].set_xlabel(r'$R/R_\ast$')
   ax0[2].set_ylabel(r'$v_i [km/s]$')
   ax0[2].plot(radius,velr1d,color='black', label=r'$v_r$')
   ax0[2].plot(radius,velth1d,color='red', label=r'$v_\theta$')
   ax0[2].plot(radius,velphi1d,color='blue', label=r'$v_\phi$')            
   ax0[2].legend(fontsize=12)


   ax0[3].set_xlabel(r'$R/R_\ast$')
   ax0[3].set_ylabel(r'$\dot{M} [Msun/yr]$')
   ax0[3].plot(radius,4.*np.pi*rho1d*velr1d*1.e5*(radius*rstar*cgs_rsu)**2 *cgs_yr/cgs_msu,color='black')
   ax0[3].set_yscale('log')   


#   sdum = input("Press [q] to exit.")
#   if sdum == 'q':
#       exit()
   
#
#***********************************************************************
#
#                           contour plots
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1.6
   ysize=xsize/aspect_ratio

   fig1 = plt.figure(windx,figsize=(xsize,ysize))
   windx=windx+1
   fig2 = plt.figure(windx,figsize=(xsize,ysize))
   windx=windx+1
   fig3 = plt.figure(windx,figsize=(xsize,ysize))   
   windx=windx+1
   plt.ion()
   plt.show()

   ax1 = fig1.subplots()
   ax2 = fig2.subplots()
   ax3 = fig3.subplots()      
#
#-----------------------------------------------------------------------
#
   ax1.set_xlim(xlim)
   ax1.set_ylim(ylim)
   ax1.set_xlabel(r'$x$')
   ax1.set_ylabel(r'$z$')


   ax2.set_xlim(xlim)
   ax2.set_ylim(ylim)
   ax2.set_xlabel(r'$x$')
   ax2.set_ylabel(r'$z$')


   ax3.set_xlim(xlim)
   ax3.set_ylim(ylim)
   ax3.set_xlabel(r'$x$')
   ax3.set_ylabel(r'$z$')      
#
#-----------------------------------------------------------------------
#
   rho2d=np.zeros(shape=(ntheta,nr))
   velr2d=np.zeros(shape=(ntheta,nr))
   t2d=np.zeros(shape=(ntheta,nr))
#
   x_coord=np.zeros(shape=(ntheta,nr))
   y_coord=np.zeros(shape=(ntheta,nr))

   kkm1, kk = find_indx(phic, phi, nphi)

   for i in range(0,nr):
      for j in range(0,ntheta):
         x_coord[j][i] = radius[i]*np.sin(theta[j])
         y_coord[j][i] = radius[i]*np.cos(theta[j])         
         rho2d[j][i] = interpol1d_2p_lin(rho3d[kkm1][j][i],rho3d[kk][j][i],phi[kkm1],phi[kk],phic)
         t2d[j][i] = interpol1d_2p_lin(t3d[kkm1][j][i],t3d[kk][j][i],phi[kkm1],phi[kk],phic)
         velr2d[j][i] = interpol1d_2p_lin(velr3d[kkm1][j][i],velr3d[kk][j][i],phi[kkm1],phi[kk],phic)                  
#
#-------------------density---------------------------------------------
#
   lrho2d = get_logvals(rho2d)
   clevels, ticks = get_clevels(clim=clim_rho)
   cmap_lev = get_cmap('jet')   

   contourplot = ax1.contourf(x_coord, y_coord, lrho2d,
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1,orientation='vertical', aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\rho)$')                 
#
#-------------------temperature-----------------------------------------
#
#  indx = np.where(rho2d > 0.)
   cmin=clim_temp[0]
   cmax=clim_temp[1]
   dcol=(cmax-cmin)/99.
   clevels=np.arange(cmin,cmax+dcol,dcol)
   clevels, ticks = get_clevels(clim=clim_temp)
   cmap_lev = get_cmap('jet')   

   contourplot= ax2.contourf(x_coord, y_coord, t2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig2.colorbar(contourplot,ax=ax2,orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$T$')                   
#
#-------------------radial velocity-------------------------------------
#
   cmin=clim_velr[0]
   cmax=clim_velr[1]
   dcol=(cmax-cmin)/99.
   clevels=np.arange(cmin,cmax+dcol,dcol)
   clevels, ticks = get_clevels(clim=clim_velr)
   cmap_lev = get_cmap('jet')   

   contourplot= ax3.contourf(x_coord, y_coord, velr2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3,orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$v_r$')   

#
#------------------------output to files-------------------------------
#
   oname1 = oname3d+'_radial.png'
   oname2 = oname3d+'_radial.ps'
   fig0.savefig(oname1, bbox_inches='tight')
   fig0.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_rho.png'
   oname2 = oname3d+'_rho.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_temp.png'
   oname2 = oname3d+'_temp.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')
   
   oname1 = oname3d+'_velr.png'
   oname2 = oname3d+'_velr.ps'
   fig3.savefig(oname1, bbox_inches='tight')
   fig3.savefig(oname2, bbox_inches='tight')
#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################

fname1='../TRASH/model3d.h5'
fname2='../inputFILES/model3d.h5'

#compare_model3d(fname1,fname2)




dir='../inputFILES'
windx = 1

map3d(dir=dir,windx=windx,rp=4.)
windx=windx+1
main(dir=dir,windx=windx,thetac=90.*np.pi/180.,phic=315.*np.pi/180., rstar=20., clim_rho=[-16.,-12.], clim_temp=[10.e3,200.e3], clim_velr=[1.e2,2.e3], xlim=[0.,7.], ylim=[-3.5,3.5])

         
sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
