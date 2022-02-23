import sys
sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from const_cgs import *
from read_models import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from lev_contour import *
#
def main(dir='../inputFILES', oname1d='./ps_files/model2d_rad',
         oname2d='./ps_files/model2d', windx=1, clim_rho=[-17.,-9.], clim_velr=[0.,2000.], xlim=[0.,8.],
         ylim=[-4.,4.]):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/model2d.h5'

   nr, ntheta, r, theta, rho2d, velr2d, velth2d, \
   velphi2d, t2d, vth2d = read_model2d(fname)

#velocity in km/s
   velr2d=velr2d/1.e5
   velth2d=velth2d/1.e5
   velphi2d=velphi2d/1.e5
   vth2d = vth2d/1.e5

   mdot2d = 4.*np.pi*rho2d*velr2d*1.e5*cgs_yr/cgs_msu
   for i in np.arange(0,nr):
      mdot2d[:,i] = mdot2d[:,i] * r[i]**2 

#r in r_star
   r=r/r[0]
#
#-----------------------define range------------------------------------
#
   rlog = np.log10(r-1.0)
   rlog[0] = -6.0
#
#----------------------------title-strings------------------------------
#
   params = {'legend.fontsize': 12,
              'figure.figsize': (10, 8),
              'figure.titlesize': 14,             
              'axes.labelsize': 12,
              'axes.titlesize': 12,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10}
   pylab.rcParams.update(params)

   fig1 = plt.figure(windx, constrained_layout=True)
   windx=windx+1

   fig2 = plt.figure(windx, constrained_layout=True)
   windx=windx+1   
#
   plt.ion()
   plt.show()
#
   ax1 = fig1.subplots(2,4)
   ax2 = fig2.subplots(2,4)   
#
#-----------------------------------------------------------------------
#
#plots along radial direction for several theta angles
#
   theta0=0.0
   theta1=20.0
   theta2=40.0
   theta3=60.0
   theta4=80.0
   theta5=90.0
#
   indx0, ii = find_indx(theta0*np.pi/180.0, theta, ntheta)
   iim1, indx1 = find_indx(theta1*np.pi/180.0, theta, ntheta)
   iim1, indx2 = find_indx(theta2*np.pi/180.0, theta, ntheta)
   iim1, indx3 = find_indx(theta3*np.pi/180.0, theta, ntheta)
   iim1, indx4 = find_indx(theta4*np.pi/180.0, theta, ntheta)
   indx5, ii = find_indx(theta5*np.pi/180.0, theta, ntheta)
#
#--------------v_r(r) for different theta-------------------------------
#
   ax1[0,0].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[0,0].set_ylabel(r'$v_r$ [km/s]')
   ax1[0,0].plot(rlog,velr2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[0,0].plot(rlog,velr2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[0,0].plot(rlog,velr2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[0,0].plot(rlog,velr2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[0,0].plot(rlog,velr2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[0,0].plot(rlog,velr2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[0,0].legend()

   ax2[0,0].set_xlabel(r'$R/R_\ast$')
   ax2[0,0].set_ylabel(r'$v_r$ [km/s]')
#   ax2[0,0].plot(r,velr2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
#   ax2[0,0].plot(r,velr2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
#   ax2[0,0].plot(r,velr2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
#   ax2[0,0].plot(r,velr2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
#   ax2[0,0].plot(r,velr2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[0,0].plot(r,velr2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[0,0].legend()
#
#------------v_theta(r) for different theta----------------------------
#
   ax1[0,1].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[0,1].set_ylabel(r'$v_\theta$ [km/s]')
   ax1[0,1].plot(rlog,velth2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[0,1].plot(rlog,velth2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[0,1].plot(rlog,velth2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[0,1].plot(rlog,velth2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[0,1].plot(rlog,velth2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[0,1].plot(rlog,velth2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[0,1].legend()

   ax2[0,1].set_xlabel(r'$R/R_\ast$')
   ax2[0,1].set_ylabel(r'$v_\theta$ [km/s]')
   ax2[0,1].plot(r,velth2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax2[0,1].plot(r,velth2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax2[0,1].plot(r,velth2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax2[0,1].plot(r,velth2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax2[0,1].plot(r,velth2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[0,1].plot(r,velth2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[0,1].legend()   
#
#------------v_phi(r) for different theta-------------------------------
#
   ax1[0,2].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[0,2].set_ylabel(r'$v_\phi$ [km/s]')
   ax1[0,2].plot(rlog,velphi2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[0,2].plot(rlog,velphi2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[0,2].plot(rlog,velphi2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[0,2].plot(rlog,velphi2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[0,2].plot(rlog,velphi2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[0,2].plot(rlog,velphi2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[0,2].legend()

   ax2[0,2].set_xlabel(r'$R/R_\ast$')
   ax2[0,2].set_ylabel(r'$v_\phi$ [km/s]')
   ax2[0,2].plot(r,velphi2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax2[0,2].plot(r,velphi2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax2[0,2].plot(r,velphi2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax2[0,2].plot(r,velphi2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax2[0,2].plot(r,velphi2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[0,2].plot(r,velphi2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[0,2].legend()   
#
#-----------------------density for different theta---------------------
#
   ax1[1,0].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[1,0].set_ylabel(r'$\rho \rm{[g/cm^3]}$')
   ax1[1,0].plot(rlog,rho2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[1,0].plot(rlog,rho2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[1,0].plot(rlog,rho2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[1,0].plot(rlog,rho2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[1,0].plot(rlog,rho2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[1,0].plot(rlog,rho2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[1,0].legend()
   ax1[1,0].set_yscale('log')

   ax2[1,0].set_xlabel(r'$R/R_\ast$')
   ax2[1,0].set_ylabel(r'$\rho \rm{[g/cm^3]}$')
#   ax2[1,0].plot(r,rho2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
#   ax2[1,0].plot(r,rho2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
#   ax2[1,0].plot(r,rho2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
#   ax2[1,0].plot(r,rho2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
#   ax2[1,0].plot(r,rho2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[1,0].plot(r,rho2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[1,0].legend()
   ax2[1,0].set_yscale('log')   
#
#-----------------------temperature for different theta-----------------
#
   ax1[1,1].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[1,1].set_ylabel(r'$T$ [K]')
   ax1[1,1].plot(rlog,t2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[1,1].plot(rlog,t2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[1,1].plot(rlog,t2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[1,1].plot(rlog,t2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[1,1].plot(rlog,t2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[1,1].plot(rlog,t2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[1,1].legend()

   ax2[1,1].set_xlabel(r'$R/R_\ast$')
   ax2[1,1].set_ylabel(r'$T$ [K]')
   ax2[1,1].plot(r,t2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax2[1,1].plot(r,t2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax2[1,1].plot(r,t2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax2[1,1].plot(r,t2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax2[1,1].plot(r,t2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[1,1].plot(r,t2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[1,1].legend()   
#
#-------------------thermal velocity for different theta-----------------
#
   ax1[1,2].set_xlabel(r'$\log (R/R_\ast)-1$')
   ax1[1,2].set_ylabel(r'$v_{th}$ [km/s]')
   ax1[1,2].plot(rlog,vth2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax1[1,2].plot(rlog,vth2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax1[1,2].plot(rlog,vth2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax1[1,2].plot(rlog,vth2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax1[1,2].plot(rlog,vth2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax1[1,2].plot(rlog,vth2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax1[1,2].legend()

   ax2[1,2].set_xlabel(r'$R/R_\ast$')
   ax2[1,2].set_ylabel(r'$v_{th}$ [km/s]')
   ax2[1,2].plot(r,vth2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
   ax2[1,2].plot(r,vth2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
   ax2[1,2].plot(r,vth2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
   ax2[1,2].plot(r,vth2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
   ax2[1,2].plot(r,vth2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[1,2].plot(r,vth2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[1,2].legend()

#
#-------------------mass - loss rate-------------------------------------
#
#   ax1[1,2].set_xlabel(r'$\log (R/R_\ast)-1$')
#   ax1[1,2].set_ylabel(r'$v_{th}$ [km/s]')
#   ax1[1,2].plot(rlog,vth2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
#   ax1[1,2].plot(rlog,vth2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
#   ax1[1,2].plot(rlog,vth2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
#   ax1[1,2].plot(rlog,vth2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
#   ax1[1,2].plot(rlog,vth2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
#   ax1[1,2].plot(rlog,vth2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
#   ax1[1,2].legend()

   ax2[0,3].set_xlabel(r'$R/R_\ast$')
   ax2[0,3].set_ylabel(r'$\dot{M} [Msun/yr]$')   
#  ax2[0,3].plot(r,mdot2d[:][indx0],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx0]*180./np.pi))
#  ax2[0,3].plot(r,mdot2d[:][indx1],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx1]*180./np.pi))
#  ax2[0,3].plot(r,mdot2d[:][indx2],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx2]*180./np.pi))
#  ax2[0,3].plot(r,mdot2d[:][indx3],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx3]*180./np.pi))
#  ax2[0,3].plot(r,mdot2d[:][indx4],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx4]*180./np.pi))
   ax2[0,3].plot(r,mdot2d[:][indx5],label=r'$\theta=${sth:.0f}'.format(sth=theta[indx5]*180./np.pi))
   ax2[0,3].legend()
   ax2[0,3].set_yscale('log')   

   oname1 = oname1d+'_a.png'
   oname2 = oname1d+'_a.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')

   oname1 = oname1d+'_b.png'
   oname2 = oname1d+'_b.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')   
#
#***********************************************************************
#
#                           contour plots
#
#***********************************************************************
#

   params = {'legend.fontsize': 12,
              'figure.figsize': (10, 8),
              'figure.titlesize': 14,             
              'axes.labelsize': 12,
              'axes.titlesize': 12,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10}
   pylab.rcParams.update(params)

   fig3 = plt.figure(windx,constrained_layout=True)
   ax3 = fig3.subplots()
   windx=windx+1   

   fig4 = plt.figure(windx,constrained_layout=True)
   ax4 = fig4.subplots()
   windx=windx+1   
#
#-----------------------------------------------------------------------
#
   ax3.set_xlim(xlim)
   ax3.set_ylim(ylim)
   ax3.set_xlabel(r'$x$')
   ax3.set_ylabel(r'$y$')

   ax4.set_xlim(xlim)
   ax4.set_ylim(ylim)
   ax4.set_xlabel(r'$x$')
   ax4.set_ylabel(r'$y$')      
#
#-----------------------------------------------------------------------
#
   x_coord=np.zeros(shape=(ntheta,nr))
   y_coord=np.zeros(shape=(ntheta,nr))
   for i in range(0,nr):
       for j in range(0,ntheta):
          x_coord[j][i] = r[i]*np.sin(theta[j])
          y_coord[j][i] = r[i]*np.cos(theta[j])
#
#-----------------------------------------------------------------------
#
   clevels, ticks = get_clevels(clim=clim_rho)

   lrho2d=np.log10(rho2d)   
#   indx = np.where(rho2d > 0.)
#   rho_min = np.min(rho2d[indx])
#   lrho2d = rho2d*1.
#   indx = np.where(lrho2d <= 0.)
#   lrho2d[indx] = rho_min
#   lrho2d = np.log10(lrho2d)

   cmap_lev = get_cmap('jet')
   contourplot = ax3.contourf(x_coord, y_coord, lrho2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3,orientation='vertical')
   cbar.set_label(r'$\log(\rho)$')
#
#--------------------------prepare velocity vectors---------------------
#
   dim_vxslice=15
   dim_vyslice=15
#
   xmin=xlim[0]
   xmax=0.9*xlim[1]
   ymin=0.9*ylim[0]
   ymax=0.9*ylim[1]
#
   vx_slice=grid_equi(xmin,xmax,dim_vxslice)
   vy_slice=grid_equi(ymin,ymax,dim_vyslice)
#
   vx2d=np.zeros(shape=(dim_vyslice,dim_vxslice))
   vy2d=np.zeros(shape=(dim_vyslice,dim_vxslice))
#
   for i in range(0,dim_vxslice):
      for j in range(0,dim_vyslice):
         rad=np.sqrt(vx_slice[i]**2 + vy_slice[j]**2)
         th=np.arccos(vy_slice[j]/rad)
#
         if rad > 1.:
            vtheta= interpol_bilinear(velth2d, r, theta, [rad,th])
            vrad= interpol_bilinear(velr2d, r, theta, [rad,th])
            vx2d[j][i]=vrad*np.sin(th)+vtheta*np.cos(th)
            vy2d[j][i]=vrad*np.cos(th)-vtheta*np.sin(th)
#
#set maximum allowed v to 2 (to have a nice plot)
            vlim=0.5
            vel=np.sqrt(vx2d[j,i]**2 + vy2d[j,i]**2)
            vnorm=vlim/vel
            if vel > vlim:
               vx2d[j][i] = vx2d[j][i]*vnorm
               vy2d[j][i] = vy2d[j][i]*vnorm

   ax3.quiver(vx_slice,vy_slice,vx2d,vy2d)
#
#--------------------------plot radial velocities-----------------------
#
   clevels, ticks = get_clevels(clim=clim_velr)

   cmap_lev = get_cmap('jet')
   contourplot = ax4.contourf(x_coord, y_coord, velr2d,
                             levels=clevels,
                             extend='both',
                             cmap=cmap_lev)
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax4,orientation='vertical')
   cbar.set_label(r'$v_r$ [km/s]')
#
#
#-----------------------------------------------------------------------
#                                
#
#
   oname1 = oname2d+'.png'
   oname2 = oname2d+'.ps'
   fig3.savefig(oname1, bbox_inches='tight')
   fig3.savefig(oname2, bbox_inches='tight')

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
   
dir='../inputFILES'
windx = 1
main(dir=dir,windx=windx, xlim=[0.,6.], ylim=[-3.,3.], clim_rho=[-16.,-10.], clim_velr=[0.,2500.])


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
