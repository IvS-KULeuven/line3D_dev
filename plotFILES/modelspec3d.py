import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_modelspec import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from lev_contour import *




########################################################################

def compare_modelspec3d(fname1, fname2):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
   small_number=1.e-4
#


   teff_01, trad_01, xnue0_01, xic1_01, rstar_01, vth_fiducial_01, vmicro_01, vmax_01, logg_01, na_01, \
      nr_01, ntheta_01, nphi_01, radius_01, theta_01, phi_01, t3d_01, opac3d_01, opalbar3d_01, velx3d_01, vely3d_01, velz3d_01, \
      sline3d_01, scont3d_01 = read_modelspec3d_spc(fname1)

   teff_02, trad_02, xnue0_02, xic1_02, rstar_02, vth_fiducial_02, vmicro_02, vmax_02, logg_02, na_02, \
      nr_02, ntheta_02, nphi_02, radius_02, theta_02, phi_02, t3d_02, opac3d_02, opalbar3d_02, velx3d_02, vely3d_02, velz3d_02, \
      sline3d_02, scont3d_02 = read_modelspec3d_spc(fname2)   


   print('{str1:10}{str2:25}{str3:25}'.format(str1='##',str2='model1',str3='model2'))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='nr',x1=nr_01,x2=nr_02))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='ntheta',x1=ntheta_01,x2=ntheta_02))
   print('{str1:10}{x1:25d}{x2:25d}'.format(str1='nphi',x1=nphi_01,x2=nphi_02))
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='teff',x1=teff_01,x2=teff_02,x3=teff_01-teff_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='trad',x1=trad_01,x2=trad_02,x3=trad_01-trad_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='xnue0',x1=xnue0_01,x2=xnue0_02,x3=xnue0_01-xnue0_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='xic1',x1=xic1_01,x2=xic1_02,x3=xic1_01-xic1_02))
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='rstar',x1=rstar_01,x2=rstar_02,x3=rstar_01-rstar_02))      
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='vth_fiducial',x1=vth_fiducial_01,x2=vth_fiducial_02,x3=vth_fiducial_01-vth_fiducial_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='vmicro',x1=vmicro_01,x2=vmicro_02,x3=vmicro_01-vmicro_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='vmax',x1=vmax_01,x2=vmax_02,x3=vmax_01-vmax_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='logg',x1=logg_01,x2=logg_02,x3=logg_01-logg_02))   
   print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='na',x1=na_01,x2=na_02,x3=na_01-na_02))   

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
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='opac(i,j,k)',x1=opac3d_01[k][j][i],x2=opac3d_02[k][j][i],x3=opac3d_01[k][j][i]-opac3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='opalbar(i,j,k)',x1=opalbar3d_01[k][j][i],x2=opalbar3d_02[k][j][i],x3=opalbar3d_01[k][j][i]-opalbar3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='scont(i,j,k)',x1=scont3d_01[k][j][i],x2=scont3d_02[k][j][i],x3=scont3d_01[k][j][i]-scont3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='sline(i,j,k)',x1=sline3d_01[k][j][i],x2=sline3d_02[k][j][i],x3=sline3d_01[k][j][i]-sline3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='velx(i,j,k)',x1=velx3d_01[k][j][i],x2=velx3d_02[k][j][i],x3=velx3d_01[k][j][i]-velx3d_02[k][j][i]))            
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='vely(i,j,k)',x1=vely3d_01[k][j][i],x2=vely3d_02[k][j][i],x3=vely3d_01[k][j][i]-vely3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='velz(i,j,k)',x1=velz3d_01[k][j][i],x2=velz3d_02[k][j][i],x3=velz3d_01[k][j][i]-velz3d_02[k][j][i]))
#            print('{str1:10}{x1:25.14e}{x2:25.14e}{x3:25.15e}'.format(str1='t(i,j,k)',x1=t3d_01[k][j][i],x2=t3d_02[k][j][i],x3=t3d_01[k][j][i]-t3d_02[k][j][i]))
#            print('')
            
            if np.abs(opac3d_01[k][j][i]-opac3d_02[k][j][i]) > small_number:
               print('error: opacs different')
               exit()
            if np.abs(opalbar3d_01[k][j][i]-opalbar3d_02[k][j][i]) > small_number:
               print('error: opalbars different')
               exit()               
            if np.abs(scont3d_01[k][j][i]-scont3d_02[k][j][i]) > small_number:
               print('error: sconts different')
               exit()
            if np.abs(sline3d_01[k][j][i]-sline3d_02[k][j][i]) > small_number:
               print('error: slines different')
               exit()
            if np.abs(velx3d_01[k][j][i]-velx3d_02[k][j][i]) > small_number:
               print('error: velxs different')
               exit()
            if np.abs(vely3d_01[k][j][i]-vely3d_02[k][j][i]) > small_number:
               print('error: velys different')
               exit()   
            if np.abs(velz3d_01[k][j][i]-velz3d_02[k][j][i]) > small_number:
               print('error: velzs different')
               exit()
            if np.abs(t3d_01[k][j][i]-t3d_02[k][j][i]) > small_number:
               print('error: ts different')
               exit()         

   exit()   
#
########################################################################
#
def main(fname='../outputFILES/modspec_model00.h5', oname3d='./ps_files/modspec', thetac=0., phic=0., windx=1,
         version='spc', clim_opalbar=[-20.,-8.], clim_sline=[-5.,-2.], clim_temp=[10.e3,100.e3], clim_velr=[1.e2,4.e3], xlim=[0.,8.], ylim=[-4.,4.]):
#
#
#+
# NAME:
#	modelspec3d
#
# PURPOSE:
#	This procedure plots contours on a slice through the 3d grid from
#	output of modelspec.eo. The slice is defined by the angle phic
#       Contours are shown for:
#          3d model: rho, v_r, t
#
#	This procedure plots the radial stratification of the 3d grid from
#	output of modelspec.eo, at a given thetac, phic
#       Plots are shown for:
#          3d model: rho, v_r, v_th, v_phi, t
#
# INPUTS:
#	
# KEYWORD PARAMETERS:
#       thetac: theta coordinate for radial stratification  
#       phic:   phi coordinate for for slice and radial stratification
#       windx:  Set this keyword to an integer, defining the window-index
#       oname:  Set this keyword to a string, defining the output-name
#               (ps-file)
#       xlim:   Set this keyword to an array of length 2, defining the xrange
#       ylim:   Set this keyword to an array of length 2, defining the yrange
#       xlim1d: Set this keyword to an array of length 2, defining the xrange of the radial stratification
#       version:  Set this keyword to 'spc' or 'cac'
#
# OUTPUTS:
#	
#------------------read all information from hdf5-file------------------
#
   if version == 'cac':
      print('version cac in modelspec3d.py to be implemented')
      exit()
   elif version == 'spc':
      teff, trad, xnue0, xic1, rstar, vth_fiducial, vmicro, vmax, logg, na, \
         nr, ntheta, nphi, radius, theta, phi, t3d, opac3d, opalbar3d, velx3d, vely3d, velz3d, \
         sline3d, scont3d = read_modelspec3d_spc(fname)
   else:
      print('error in modelspec3d.py: version needs to be set to cac or spc')
      exit()

   sr=rstar*cgs_rsu

   velr3d=np.zeros(shape=(nphi,ntheta,nr))
   velth3d=np.zeros(shape=(nphi,ntheta,nr))
   velphi3d=np.zeros(shape=(nphi,ntheta,nr))
#   for i in range(0,nr):
#      print(i)
#      for j in range(0,ntheta):
#         for k in range(0,nphi):
   for k in range(0,nphi):
#      print(k)
      sinp=np.sin(phi[k])
      cosp=np.cos(phi[k])      
      for j in range(0,ntheta):
         sint=np.sin(theta[j])
         cost=np.cos(theta[j])
         velx=velx3d[k][j][:]
         vely=vely3d[k][j][:]
         velz=velz3d[k][j][:]         
         velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
         velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
         velphi3d[k][j][:] = -velx*sinp + vely*cosp         
   
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
   xsize=xsize/1.54 #in inches
   aspect_ratio=.4
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize))
   windx=windx+1   
   plt.ion()
   plt.show()
   ax0 = fig0.subplots(3,2)

   opalbar1d=np.zeros(nr)
   opac1d=np.zeros(nr)   
   sline1d=np.zeros(nr)   
   scont1d=np.zeros(nr)
   velr1d=np.zeros(nr)
   velth1d=np.zeros(nr)
   velphi1d=np.zeros(nr)
   t1d=np.zeros(nr)
#
   jjm1, jj = find_indx(thetac, theta, ntheta)
   kkm1, kk = find_indx(phic, phi, nphi)
#
   for i in range(0,nr):
      opalbar1d[i] = interpol2d_4p_lin(opalbar3d[kkm1,jjm1,i], opalbar3d[kkm1][jj][i], opalbar3d[kk][jjm1][i], opalbar3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      opac1d[i] = interpol2d_4p_lin(opac3d[kkm1,jjm1,i], opac3d[kkm1][jj][i], opac3d[kk][jjm1][i], opac3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      sline1d[i] = interpol2d_4p_lin(sline3d[kkm1,jjm1,i], sline3d[kkm1][jj][i], sline3d[kk][jjm1][i], sline3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      scont1d[i] = interpol2d_4p_lin(scont3d[kkm1,jjm1,i], scont3d[kkm1][jj][i], scont3d[kk][jjm1][i], scont3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velr1d[i] = interpol2d_4p_lin(velr3d[kkm1,jjm1,i], velr3d[kkm1][jj][i], velr3d[kk][jjm1][i], velr3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velth1d[i] = interpol2d_4p_lin(velth3d[kkm1,jjm1,i], velth3d[kkm1][jj][i], velth3d[kk][jjm1][i], velth3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      velphi1d[i] = interpol2d_4p_lin(velphi3d[kkm1,jjm1,i], velphi3d[kkm1][jj][i], velphi3d[kk][jjm1][i], velphi3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
      t1d[i] = interpol2d_4p_lin(t3d[kkm1,jjm1,i], t3d[kkm1][jj][i], t3d[kk][jjm1][i], t3d[kk][jj][i],
                                   theta[jjm1],theta[jj],phi[kkm1],phi[kk],thetac,phic)
#
#-----------------------------------------------------------------------
#
   titlestr=r'at $\theta, \phi=$({sth:.2f},{sphi:.2f})'.format(sth=thetac*180./np.pi,sphi=phic*180./np.pi)
#
   ax0[0,0].set_xlabel(r'$R/R_\ast$')
   ax0[0,0].set_ylabel(r'$\chi_c$')
   ax0[0,1].set_yscale('linear')   
   ax0[0,0].plot(radius,opac1d)
   #ax0[0,0].set_title(titlestr)

   ax0[0,1].set_xlabel(r'$R/R_\ast$')
   ax0[0,1].set_ylabel(r'$\bar{\chi}_L$')
   ax0[0,1].set_yscale('log')   
   ax0[0,1].plot(radius,opalbar1d)

   ax0[1,0].set_xlabel(r'$R/R_\ast$')
   ax0[1,0].set_ylabel(r'$S_c$')
   ax0[1,0].set_yscale('linear')      
   ax0[1,0].plot(radius,scont1d)

   ax0[1,1].set_xlabel(r'$R/R_\ast$')
   ax0[1,1].set_ylabel(r'$S_L$')
   ax0[1,1].set_yscale('log')      
   ax0[1,1].plot(radius,sline1d)

   ax0[2,0].set_xlabel(r'$R/R_\ast$')
   ax0[2,0].set_ylabel(r'$T$')
   ax0[2,0].set_yscale('linear')         
   ax0[2,0].plot(radius,t1d)


   ax0[2,1].set_xlabel(r'$R/R_\ast$')
   ax0[2,1].set_ylabel(r'$v_i [km/s]$')
   ax0[2,1].plot(radius,velr1d,color='black', label=r'$v_r$')
   ax0[2,1].plot(radius,velth1d,color='red', label=r'$v_\theta$')
   ax0[2,1].plot(radius,velphi1d,color='blue', label=r'$v_\phi$')            
   ax0[2,1].legend(fontsize=12)
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
   fig4 = plt.figure(windx,figsize=(xsize,ysize))   
   windx=windx+1   
   plt.ion()
   plt.show()

   ax1 = fig1.subplots()
   ax2 = fig2.subplots()
   ax3 = fig3.subplots()
   ax4 = fig4.subplots()
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
   
   ax4.set_xlim(xlim)
   ax4.set_ylim(ylim)
   ax4.set_xlabel(r'$x$')
   ax4.set_ylabel(r'$z$')   
#
#-----------------------------------------------------------------------
#
#
   opalbar2d=np.zeros(shape=(ntheta,nr))
   opac2d=np.zeros(shape=(ntheta,nr))   
   sline2d=np.zeros(shape=(ntheta,nr))   
   scont2d=np.zeros(shape=(ntheta,nr))
   velr2d=np.zeros(shape=(ntheta,nr))
   velth2d=np.zeros(shape=(ntheta,nr))
   velphi2d=np.zeros(shape=(ntheta,nr))
   t2d=np.zeros(shape=(ntheta,nr))
   
#
   x_coord=np.zeros(shape=(ntheta,nr))
   y_coord=np.zeros(shape=(ntheta,nr))

   kkm1, kk = find_indx(phic, phi, nphi)

   for i in range(0,nr):
      for j in range(0,ntheta):
         x_coord[j][i] = radius[i]*np.sin(theta[j])
         y_coord[j][i] = radius[i]*np.cos(theta[j])         
         opalbar2d[j][i] = interpol1d_2p_lin(opalbar3d[kkm1][j][i],opalbar3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         opac2d[j][i] = interpol1d_2p_lin(opac3d[kkm1][j][i],opac3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         sline2d[j][i] = interpol1d_2p_lin(sline3d[kkm1][j][i],sline3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         scont2d[j][i] = interpol1d_2p_lin(scont3d[kkm1][j][i],scont3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         velr2d[j][i] = interpol1d_2p_lin(velr3d[kkm1][j][i],velr3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         velth2d[j][i] = interpol1d_2p_lin(velth3d[kkm1][j][i],velth3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         velphi2d[j][i] = interpol1d_2p_lin(velphi3d[kkm1][j][i],velphi3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
         t2d[j][i] = interpol1d_2p_lin(t3d[kkm1][j][i],t3d[kkm1][j][i],phi[kkm1],phi[kk],phic)
#
#-------------------line opacity----------------------------------------
#
   lopalbar2d = get_logvals(opalbar2d)
   clevels, ticks = get_clevels(clim=clim_opalbar)
   cmap_lev = get_cmap('jet')   

   contourplot = ax1.contourf(x_coord, y_coord, lopalbar2d,
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1,orientation='vertical', aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\bar{\chi}_L)$')
#
#-------------------line source function--------------------------------
#

   lsline2d = get_logvals(sline2d)
   clevels, ticks = get_clevels(clim=clim_sline)
   cmap_lev = get_cmap('jet')   

   contourplot = ax2.contourf(x_coord, y_coord, lsline2d,
                              levels=clevels,
                              extend='both',
                              cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig2.colorbar(contourplot,ax=ax2,orientation='vertical', aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(S_L)$')
#
#-------------------temperature-----------------------------------------
#

   clevels, ticks = get_clevels(clim=clim_temp)
   cmap_lev = get_cmap('jet')   

   contourplot= ax3.contourf(x_coord, y_coord, t2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3,orientation='vertical', aspect=40,ticks=ticks)
   cbar.set_label(r'$T$')                   
#
#-------------------radial velocity-------------------------------------
#
   clevels, ticks = get_clevels(clim=clim_velr)
   cmap_lev = get_cmap('jet')   

   contourplot= ax4.contourf(x_coord, y_coord, velr2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig4.colorbar(contourplot,ax=ax4,orientation='vertical', aspect=40,ticks=ticks)
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

   oname1 = oname3d+'_sline.png'
   oname2 = oname3d+'_sline.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_temp.png'
   oname2 = oname3d+'_temp.ps'
   fig3.savefig(oname1, bbox_inches='tight')
   fig3.savefig(oname2, bbox_inches='tight')
   
   oname1 = oname3d+'_velr.png'
   oname2 = oname3d+'_velr.ps'
   fig4.savefig(oname1, bbox_inches='tight')
   fig4.savefig(oname2, bbox_inches='tight')
#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################

fname1='../outputFILES/nico_wr3d/modspec_model00.h5'
fname2='../TRASH/modspec_model00.h5'
#compare_modelspec3d(fname1,fname2)


fname='../outputFILES/spc/modspec_model00.h5'
fname='../outputFILES/nico_wr3d/modspec_model00.h5'
#fname='/lhome/levin/Postdoc/for_jon/plots_wr3d/spherical3d/clumped/kline5d12/modspec_model00.h5'
#windx = 1
#main(fname=fname,windx=windx,thetac=np.pi/4.,phic=np.pi/4., clim_velr=[1.e2,2.e3], clim_temp=[1.e4,2.e5], clim_sline=[-2.5,-1.], clim_opalbar=[-10,-3], xlim=[0.,7.], ylim=[-3.5,3.5])


fname='../outputFILES/modspec_model00.h5'
windx=1
main(fname=fname,windx=windx,thetac=np.pi/4.,phic=np.pi/4., clim_velr=[1.e2,10.e3], clim_temp=[1.e4,5.e4], clim_sline=[-4.,-1.5], clim_opalbar=[-13,-8], xlim=[0.,7.], ylim=[-3.5,3.5])

sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
