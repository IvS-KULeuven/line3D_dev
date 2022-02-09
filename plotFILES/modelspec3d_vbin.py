import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from const_cgs import *
from read_modelspec import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_interpol3d import *
from lev_misc import *
from lev_contour import *
#
def main(fnames=['../outputFILES/modspec_model00.h5'], versions=['v00'], oname3d='./ps_files/modspec_vbin', thetac=0., phic=0., offset=0., windx=1,
         version='spc', clim_rho = [-4.,4.], clim_opalbar=[-20.,-8.], clim_sline=[-5.,-2.], clim_temp=[10.e3,100.e3], clim_velr=[1.e2,4.e3], xlim=[0.,8.], ylim=[-4.,4.]):
#
#+
# NAME:
#	modelspec3d_vbin
#
# PURPOSE:
#	This procedure plots contours on a slice through the 3d grid from
#	output of modelspec_vbin.eo. The slice is defined by the angle phic in each coordinate system
#       Contours are shown for:
#          3d model: rho, v_r, t
#
#	This procedure plots the radial stratification of the 3d grid from
#	output of modelspec.eo, at a given thetac, phic defined in each coordinate system
#       Plots are shown for:
#          3d model: rho, v_r, v_th, v_phi, t
#
# INPUTS:
#	
# KEYWORD PARAMETERS:
#       thetac: theta coordinate for radial stratification  
#       phic:   phi coordinate for for slice and radial stratification
#       thetac, phic define normal vector for slice that will be plotted
#       offset: offset in global coordinate system [unit_length0]
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
#prepare radial plots
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=.4
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   plt.ion()
   plt.show()
   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1

   ax0 = fig0.subplots(3,3)   
   ax1 = fig1.subplots(3,3)
   
   indx_model=0
   for fname in fnames:
#
#------------------read all information from hdf5-file------------------
#
      version=versions[indx_model]
      indx_model=indx_model+1
      iline, xnue0, vth_fiducial, vmax, unit_length, \
         x01, y01, z01, vx01, vy01, vz01, rstar1, teff1, logg1, lstar1, yhe1, vrot1, vmicro1, \
         x02, y02, z02, vx02, vy02, vz02, rstar2, teff2, logg2, lstar2, yhe2, vrot2, vmicro2, \
         cs1_nr,  cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta, cs1_phi, \
         cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d, cs1_rho3d, cs1_opac3d, cs1_opalbar3d, \
         cs1_scont3d, cs1_sline3d, \
         cs2_nr,  cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta, cs2_phi, \
         cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d, cs2_rho3d, cs2_opac3d, cs2_opalbar3d, \
      cs2_scont3d, cs2_sline3d = read_modelspec3d_vbin(fname,version=version) 

      print('model ', fname)
      print('iline,xnue0,vth_fiducial,vmax,unit_length', iline,xnue0,vth_fiducial,vmax,unit_length)
      print('x01,y01,z01',x01,y01,z01)
      print('vx01,vy01,vz01',vx01,vy01,vz01)
      print('rstar1,teff1,logg1,lstar1,yhe1',rstar1,teff1,logg1,lstar1,yhe1)
      print('vrot1,vmicro1',vrot1,vmicro1)
      print('cs1_nr,cs1_ntheta,cs1_nphi',cs1_nr,cs1_ntheta,cs1_nphi)
      print('')
      print('x02,y02,z02',x02,y02,z02)
      print('vx02,vy02,vz02',vx02,vy02,vz02)
      print('rstar2,teff2,logg2,lstar2,yhe2',rstar2,teff2,logg2,lstar2,yhe2)
      print('vrot2,vmicro2',vrot2,vmicro2)
      print('cs2_nr,cs2_ntheta,cs2_nphi',cs2_nr,cs2_ntheta,cs2_nphi)
      print('')
      print('----------------------')
      print('')

      sr1=rstar1*cgs_rsu
   
      cs1_velr3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
      cs1_velth3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
      cs1_velphi3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
      for k in range(0,cs1_nphi):
         sinp=np.sin(cs1_phi[k])
         cosp=np.cos(cs1_phi[k])      
         for j in range(0,cs1_ntheta):
            sint=np.sin(cs1_theta[j])
            cost=np.cos(cs1_theta[j])
            velx=cs1_velx3d[k][j][:]
            vely=cs1_vely3d[k][j][:]
            velz=cs1_velz3d[k][j][:]         
            cs1_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
            cs1_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
            cs1_velphi3d[k][j][:] = -velx*sinp + vely*cosp         
      
   #velocity in km/s
      cs1_velr3d=cs1_velr3d/1.e5
      cs1_velth3d=cs1_velth3d/1.e5
      cs1_velphi3d=cs1_velphi3d/1.e5
   
   
   
   
      sr2=rstar2*cgs_rsu
   
      cs2_velr3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
      cs2_velth3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
      cs2_velphi3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
      for k in range(0,cs2_nphi):
         sinp=np.sin(cs2_phi[k])
         cosp=np.cos(cs2_phi[k])      
         for j in range(0,cs2_ntheta):
            sint=np.sin(cs2_theta[j])
            cost=np.cos(cs2_theta[j])
            velx=cs2_velx3d[k][j][:]
            vely=cs2_vely3d[k][j][:]
            velz=cs2_velz3d[k][j][:]         
            cs2_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
            cs2_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
            cs2_velphi3d[k][j][:] = -velx*sinp + vely*cosp         
      
#velocity in km/s
      cs2_velr3d=cs2_velr3d/1.e5
      cs2_velth3d=cs2_velth3d/1.e5
      cs2_velphi3d=cs2_velphi3d/1.e5   
#
#***********************************************************************
#
#             radial stratification star 1
#
#***********************************************************************
#
      cs1_rho1d=np.zeros(cs1_nr)   
      cs1_opalbar1d=np.zeros(cs1_nr)
      cs1_opac1d=np.zeros(cs1_nr)   
      cs1_sline1d=np.zeros(cs1_nr)   
      cs1_scont1d=np.zeros(cs1_nr)
      cs1_velr1d=np.zeros(cs1_nr)
      cs1_velth1d=np.zeros(cs1_nr)
      cs1_velphi1d=np.zeros(cs1_nr)
      cs1_t1d=np.zeros(cs1_nr)
#
      jjm1, jj = find_indx(thetac, cs1_theta, cs1_ntheta)
      kkm1, kk = find_indx(phic, cs1_phi, cs1_nphi)
   #
      for i in range(0,cs1_nr):
         cs1_rho1d[i] = interpol2d_4p_lin(cs1_rho3d[kkm1,jjm1,i], cs1_rho3d[kkm1][jj][i], cs1_rho3d[kk][jjm1][i], cs1_rho3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)      
         cs1_opalbar1d[i] = interpol2d_4p_lin(cs1_opalbar3d[kkm1,jjm1,i], cs1_opalbar3d[kkm1][jj][i], cs1_opalbar3d[kk][jjm1][i], cs1_opalbar3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_opac1d[i] = interpol2d_4p_lin(cs1_opac3d[kkm1,jjm1,i], cs1_opac3d[kkm1][jj][i], cs1_opac3d[kk][jjm1][i], cs1_opac3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_sline1d[i] = interpol2d_4p_lin(cs1_sline3d[kkm1,jjm1,i], cs1_sline3d[kkm1][jj][i], cs1_sline3d[kk][jjm1][i], cs1_sline3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_scont1d[i] = interpol2d_4p_lin(cs1_scont3d[kkm1,jjm1,i], cs1_scont3d[kkm1][jj][i], cs1_scont3d[kk][jjm1][i], cs1_scont3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_velr1d[i] = interpol2d_4p_lin(cs1_velr3d[kkm1,jjm1,i], cs1_velr3d[kkm1][jj][i], cs1_velr3d[kk][jjm1][i], cs1_velr3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_velth1d[i] = interpol2d_4p_lin(cs1_velth3d[kkm1,jjm1,i], cs1_velth3d[kkm1][jj][i], cs1_velth3d[kk][jjm1][i], cs1_velth3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_velphi1d[i] = interpol2d_4p_lin(cs1_velphi3d[kkm1,jjm1,i], cs1_velphi3d[kkm1][jj][i], cs1_velphi3d[kk][jjm1][i], cs1_velphi3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)
         cs1_t1d[i] = interpol2d_4p_lin(cs1_t3d[kkm1,jjm1,i], cs1_t3d[kkm1][jj][i], cs1_t3d[kk][jjm1][i], cs1_t3d[kk][jj][i],
                                      cs1_theta[jjm1],cs1_theta[jj],cs1_phi[kkm1],cs1_phi[kk],thetac,phic)

      cs1_mu = np.cos(cs1_theta)
      cs1_mdot1d=np.zeros(cs1_nr)
      cs1_mflux3d = cs1_rho3d*cs1_velr3d*1.e5
      for i in range(0,cs1_nr):
         rad = cs1_radius[i]*sr1
         sum = 0.
         for j in range(1,cs1_ntheta):
            dmu = cs1_mu[j-1]-cs1_mu[j]
            for k in range(1,cs1_nphi):
               dphi = cs1_phi[k]-cs1_phi[k-1]
               sum = sum + (cs1_mflux3d[k][j][i] + cs1_mflux3d[k-1][j][i] + cs1_mflux3d[k][j-1][i] + cs1_mflux3d[k-1][j-1][i]) * rad**2 * dmu*dphi/4.
         cs1_mdot1d[i] = sum * cgs_yr / cgs_msu
#
#-----------------------------------------------------------------------
#
      titlestr=r'at $\theta, \phi=$({sth:.2f},{sphi:.2f})'.format(sth=thetac*180./np.pi,sphi=phic*180./np.pi)
#
      ax0[0,0].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[0,0].set_ylabel(r'$\chi_c$')
      ax0[0,1].set_yscale('linear')   
      ax0[0,0].plot(cs1_radius,cs1_opac1d)
      #ax0[0,0].set_title(titlestr)
   
      ax0[0,1].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[0,1].set_ylabel(r'$\bar{\chi}_L$')
      ax0[0,1].set_yscale('log')   
      ax0[0,1].plot(cs1_radius,cs1_opalbar1d)
   
      ax0[0,2].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[0,2].set_ylabel(r'$\rho$')
      ax0[0,2].set_yscale('log')   
      ax0[0,2].plot(cs1_radius,cs1_rho1d)
   
      ax0[1,0].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[1,0].set_ylabel(r'$S_c$')
      ax0[1,0].set_yscale('linear')      
      ax0[1,0].plot(cs1_radius,cs1_scont1d)
   
      ax0[1,1].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[1,1].set_ylabel(r'$S_L$')
      ax0[1,1].set_yscale('log')      
      ax0[1,1].plot(cs1_radius,cs1_sline1d)
   
      ax0[1,2].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[1,2].set_ylabel(r'$\dot{M} [Msun/yr]$')
      ax0[1,2].plot(cs1_radius,cs1_mdot1d)
   
      ax0[2,0].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[2,0].set_ylabel(r'$T$')
      ax0[2,0].set_yscale('linear')         
      ax0[2,0].plot(cs1_radius,cs1_t1d)
   
   
      ax0[2,1].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax0[2,1].set_ylabel(r'$v_i [km/s]$')
      ax0[2,1].plot(cs1_radius,cs1_velr1d,color='black', label=r'$v_r$')
      ax0[2,1].plot(cs1_radius,cs1_velth1d,color='red', label=r'$v_\theta$')
      ax0[2,1].plot(cs1_radius,cs1_velphi1d,color='blue', label=r'$v_\phi$')            
      ax0[2,1].legend(fontsize=12)
#
#***********************************************************************
#
#             radial stratification star 2
#
#***********************************************************************
#
      cs2_rho1d=np.zeros(cs2_nr)   
      cs2_opalbar1d=np.zeros(cs2_nr)
      cs2_opac1d=np.zeros(cs2_nr)   
      cs2_sline1d=np.zeros(cs2_nr)   
      cs2_scont1d=np.zeros(cs2_nr)
      cs2_velr1d=np.zeros(cs2_nr)
      cs2_velth1d=np.zeros(cs2_nr)
      cs2_velphi1d=np.zeros(cs2_nr)
      cs2_t1d=np.zeros(cs2_nr)
#
      jjm1, jj = find_indx(thetac, cs2_theta, cs2_ntheta)
      kkm1, kk = find_indx(phic, cs2_phi, cs2_nphi)
#
      for i in range(0,cs2_nr):
         cs2_rho1d[i] = interpol2d_4p_lin(cs2_rho3d[kkm1,jjm1,i], cs2_rho3d[kkm1][jj][i], cs2_rho3d[kk][jjm1][i], cs2_rho3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)      
         cs2_opalbar1d[i] = interpol2d_4p_lin(cs2_opalbar3d[kkm1,jjm1,i], cs2_opalbar3d[kkm1][jj][i], cs2_opalbar3d[kk][jjm1][i], cs2_opalbar3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_opac1d[i] = interpol2d_4p_lin(cs2_opac3d[kkm1,jjm1,i], cs2_opac3d[kkm1][jj][i], cs2_opac3d[kk][jjm1][i], cs2_opac3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_sline1d[i] = interpol2d_4p_lin(cs2_sline3d[kkm1,jjm1,i], cs2_sline3d[kkm1][jj][i], cs2_sline3d[kk][jjm1][i], cs2_sline3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_scont1d[i] = interpol2d_4p_lin(cs2_scont3d[kkm1,jjm1,i], cs2_scont3d[kkm1][jj][i], cs2_scont3d[kk][jjm1][i], cs2_scont3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_velr1d[i] = interpol2d_4p_lin(cs2_velr3d[kkm1,jjm1,i], cs2_velr3d[kkm1][jj][i], cs2_velr3d[kk][jjm1][i], cs2_velr3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_velth1d[i] = interpol2d_4p_lin(cs2_velth3d[kkm1,jjm1,i], cs2_velth3d[kkm1][jj][i], cs2_velth3d[kk][jjm1][i], cs2_velth3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_velphi1d[i] = interpol2d_4p_lin(cs2_velphi3d[kkm1,jjm1,i], cs2_velphi3d[kkm1][jj][i], cs2_velphi3d[kk][jjm1][i], cs2_velphi3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)
         cs2_t1d[i] = interpol2d_4p_lin(cs2_t3d[kkm1,jjm1,i], cs2_t3d[kkm1][jj][i], cs2_t3d[kk][jjm1][i], cs2_t3d[kk][jj][i],
                                      cs2_theta[jjm1],cs2_theta[jj],cs2_phi[kkm1],cs2_phi[kk],thetac,phic)

      cs2_mdot1d=np.zeros(cs2_nr)
      cs2_mflux3d = cs2_rho3d*cs2_velr3d*1.e5
      cs2_mu = np.cos(cs2_theta)
      for i in range(0,cs2_nr):
         rad = cs2_radius[i]*sr2
         sum = 0.
         for j in range(1,cs2_ntheta):
            dmu = cs2_mu[j-1]-cs2_mu[j]
            for k in range(1,cs2_nphi):
               dphi = cs2_phi[k]-cs2_phi[k-1]
               sum = sum + (cs2_mflux3d[k][j][i] + cs2_mflux3d[k-1][j][i] + cs2_mflux3d[k][j-1][i] + cs2_mflux3d[k-1][j-1][i]) * rad**2 * dmu*dphi/4.
         cs2_mdot1d[i] = sum * cgs_yr / cgs_msu
#
#-----------------------------------------------------------------------
#
      titlestr=r'at $\theta, \phi=$({sth:.2f},{sphi:.2f})'.format(sth=thetac*180./np.pi,sphi=phic*180./np.pi)
#
      ax1[0,0].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[0,0].set_ylabel(r'$\chi_c$')
      ax1[0,0].set_yscale('linear')   
      ax1[0,0].plot(cs2_radius,cs2_opac1d)
      #ax1[0,0].set_title(titlestr)
   
      ax1[0,1].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[0,1].set_ylabel(r'$\bar{\chi}_L$')
      ax1[0,1].set_yscale('log')   
      ax1[0,1].plot(cs2_radius,cs2_opalbar1d)

      ax1[0,2].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax1[0,2].set_ylabel(r'$\rho$')
      ax1[0,2].set_yscale('log')   
      ax1[0,2].plot(cs2_radius,cs2_rho1d)      
   
      ax1[1,0].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[1,0].set_ylabel(r'$S_c$')
      ax1[1,0].set_yscale('linear')      
      ax1[1,0].plot(cs2_radius,cs2_scont1d)
   
      ax1[1,1].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[1,1].set_ylabel(r'$S_L$')
      ax1[1,1].set_yscale('log')      
      ax1[1,1].plot(cs2_radius,cs2_sline1d)

      ax1[1,2].set_xlabel(r'$R/R_\ast^{(1)}$')
      ax1[1,2].set_ylabel(r'$\dot{M} [Msun/yr]$')
      ax1[1,2].plot(cs2_radius,cs2_mdot1d)      
   
      ax1[2,0].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[2,0].set_ylabel(r'$T$')
      ax1[2,0].set_yscale('linear')         
      ax1[2,0].plot(cs2_radius,cs2_t1d)
   
   
      ax1[2,1].set_xlabel(r'$R/R_\ast^{(2)}$')
      ax1[2,1].set_ylabel(r'$v_i [km/s]$')
      ax1[2,1].plot(cs2_radius,cs2_velr1d,color='black', label=r'$v_r$')
      ax1[2,1].plot(cs2_radius,cs2_velth1d,color='red', label=r'$v_\theta$')
      ax1[2,1].plot(cs2_radius,cs2_velphi1d,color='blue', label=r'$v_\phi$')            
      ax1[2,1].legend(fontsize=12)   
#
#***********************************************************************
#
#                contour plots
#
#***********************************************************************
#
#dimension of slices
#make the grid of the slice in polar coordinates
   dim_p1=101 #301
   dim_p2=101 #301   
#note that dim_zeta = 4*uneven - 3 if 45 degree shall be plotted
   dim_zeta1=4*5-3 #4*27-3
   dim_zeta2=27*5-3 #4*27-3   
   
#prepare contour plots
   xsize=18. #in cm
   xsize=xsize/1.5 #/2.54 #in inches
   aspect_ratio=3.
   ysize=xsize/aspect_ratio
 
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig3 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig4 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)   
   windx=windx+1
   fig5 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig6 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   

   ax2 = fig2.subplots(1,3)   
   ax3 = fig3.subplots(1,3)
   ax4 = fig4.subplots(1,3)   
   ax5 = fig5.subplots(1,3)
   ax6 = fig6.subplots(1,2)   

#   ax1.set_xlim(xlim)
#   ax1.set_ylim(ylim)
#   ax1.set_xlabel(r'$x$')
#   ax1.set_ylabel(r'$z$')   
   #
#------------------read all information from hdf5-file------------------
#
   indx_model=0
   fname=fnames[indx_model]
   version=versions[indx_model]
#
   iline, xnue0, vth_fiducial, vmax, unit_length, \
      x01, y01, z01, vx01, vy01, vz01, rstar1, teff1, logg1, lstar1, yhe1, vrot1, vmicro1, \
      x02, y02, z02, vx02, vy02, vz02, rstar2, teff2, logg2, lstar2, yhe2, vrot2, vmicro2, \
      cs1_nr,  cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta, cs1_phi, \
      cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d, cs1_rho3d, cs1_opac3d, cs1_opalbar3d, \
      cs1_scont3d, cs1_sline3d, \
      cs2_nr,  cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta, cs2_phi, \
      cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d, cs2_rho3d, cs2_opac3d, cs2_opalbar3d, \
      cs2_scont3d, cs2_sline3d = read_modelspec3d_vbin(fname,version=version)
#
   sr1=rstar1*cgs_rsu
   
   cs1_velr3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   cs1_velth3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   cs1_velphi3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   for k in range(0,cs1_nphi):
      sinp=np.sin(cs1_phi[k])
      cosp=np.cos(cs1_phi[k])      
      for j in range(0,cs1_ntheta):
         sint=np.sin(cs1_theta[j])
         cost=np.cos(cs1_theta[j])
         velx=cs1_velx3d[k][j][:]
         vely=cs1_vely3d[k][j][:]
         velz=cs1_velz3d[k][j][:]         
         cs1_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
         cs1_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
         cs1_velphi3d[k][j][:] = -velx*sinp + vely*cosp
   
#velocity in km/s
   cs1_velr3d=cs1_velr3d/1.e5
   cs1_velth3d=cs1_velth3d/1.e5
   cs1_velphi3d=cs1_velphi3d/1.e5




   sr2=rstar2*cgs_rsu

   cs2_velr3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   cs2_velth3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   cs2_velphi3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   for k in range(0,cs2_nphi):
      sinp=np.sin(cs2_phi[k])
      cosp=np.cos(cs2_phi[k])      
      for j in range(0,cs2_ntheta):
         sint=np.sin(cs2_theta[j])
         cost=np.cos(cs2_theta[j])
         velx=cs2_velx3d[k][j][:]
         vely=cs2_vely3d[k][j][:]
         velz=cs2_velz3d[k][j][:]         
         cs2_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
         cs2_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
         cs2_velphi3d[k][j][:] = -velx*sinp + vely*cosp         
   
#velocity in  km/s
   cs2_velr3d=cs2_velr3d/1.e5
   cs2_velth3d=cs2_velth3d/1.e5
   cs2_velphi3d=cs2_velphi3d/1.e5   
#
#------------------------create and interpolate onto slice-----------------------
#
   nhat=np.array([np.sin(thetac)*np.cos(phic),
                  np.sin(thetac)*np.sin(phic),
                  np.cos(thetac)])

   rmin1=cs1_radius[0]
   rmin2=cs2_radius[0]
#
#--------------calculate scalar array on slice for star 1---------------
#
#offset for coordinate system 1
   offset1 = offset - nhat[0]*x01 - nhat[1]*y01-nhat[2]*z01
   offset1 = offset1*unit_length/rstar1
               
#for star 1
   cs1_rho2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_opalbar2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_opac2d=np.zeros(shape=(dim_zeta1,dim_p1))   
   cs1_sline2d=np.zeros(shape=(dim_zeta1,dim_p1))   
   cs1_scont2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_velr2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_velth2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_velphi2d=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_t2d=np.zeros(shape=(dim_zeta1,dim_p1))
   
#
   cs1_xcoord=np.zeros(shape=(dim_zeta1,dim_p1))
   cs1_ycoord=np.zeros(shape=(dim_zeta1,dim_p1))
#
#--------------------equiduistant zeta-grid-----------------------------
#
   zeta_min=0.
   zeta_max=2.*np.pi
   cs1_zeta_arr = np.linspace(zeta_min, zeta_max, dim_zeta1)
#
#---------------p-grid equidistant (future in log space)----------------
#
   p_min=rmin1
   p_max=np.max(cs1_radius)
   cs1_p_arr = np.linspace(p_min,p_max,dim_p1)


   transmat, transmat_inv = calc_transmat2(thetac,phic)

   vec_slice = np.zeros(3)
   vec_cac = np.zeros(3)
   
   for i in range(0,dim_p1):
      for j in range(0,dim_zeta1):
         cs1_xcoord[j][i] = cs1_p_arr[i]*np.sin(cs1_zeta_arr[j])
         cs1_ycoord[j][i] = cs1_p_arr[i]*np.cos(cs1_zeta_arr[j])

         vec_slice[0] = cs1_xcoord[j][i]
         vec_slice[1] = cs1_ycoord[j][i]
         vec_slice[2] = offset1
         vec_cac = np.dot(transmat,vec_slice)
         radp=np.sqrt(vec_cac[0]*vec_cac[0] + vec_cac[1]*vec_cac[1] + vec_cac[2]*vec_cac[2])                  
#need to make sure that values on axes are treated the same
         if np.abs(vec_cac[0]) < 1.e-5: vec_cac[0]=0.
         if np.abs(vec_cac[1]) < 1.e-5: vec_cac[1]=0.
         if np.abs(vec_cac[2]) < 1.e-5: vec_cac[2]=0.
         
#get angles in spherical coordinates
         thetap, phip = angles_spc(vec_cac[0], vec_cac[1], vec_cac[2])
#find indices of cube-vertices for interpolation
         ir1, ir2, it1, it2, ip1, ip2, expol = get_rtp_indx(radp, thetap, phip, cs1_radius, cs1_theta, cs1_phi, cs1_nr, cs1_ntheta, cs1_nphi)
#
#----------------------interpolate 3d arrays onto slice-----------------
#
         rad1=cs1_radius[ir1]
         rad2=cs1_radius[ir2]
         theta1=cs1_theta[it1]
         theta2=cs1_theta[it2]
         phi1=cs1_phi[ip1]
         phi2=cs1_phi[ip2]
#
#density
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_rho3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_rho2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#opalbar
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_opalbar3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_opalbar2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#opac         
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_opac3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_opac2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)

#sline
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_sline3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_sline2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#scont
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_scont3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_scont2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#temperature
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_t3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_t2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#velocity components
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_velr3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_velr2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)

         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_velth3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_velth2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)

         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs1_nr, cs1_ntheta, cs1_nphi, cs1_velphi3d, ir1, ir2, it1, it2, ip1, ip2)
         cs1_velphi2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)         

#
#--------------calculate scalar array on slice for star 2---------------
#
#for star 2

#offset for coordinate system 2
   offset2 = offset - nhat[0]*x02 - nhat[1]*y02-nhat[2]*z02
   offset2 = offset2*unit_length/rstar2

#for star 2
   cs2_rho2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_opalbar2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_opac2d=np.zeros(shape=(dim_zeta2,dim_p2))   
   cs2_sline2d=np.zeros(shape=(dim_zeta2,dim_p2))   
   cs2_scont2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_velr2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_velth2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_velphi2d=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_t2d=np.zeros(shape=(dim_zeta2,dim_p2))
   
#
   cs2_xcoord=np.zeros(shape=(dim_zeta2,dim_p2))
   cs2_ycoord=np.zeros(shape=(dim_zeta2,dim_p2))
#
#--------------------equiduistant zeta-grid-----------------------------
#
   zeta_min=0.
   zeta_max=2.*np.pi
   cs2_zeta_arr = np.linspace(zeta_min, zeta_max, dim_zeta2)
#
#---------------p-grid equidistant (future in log space)----------------
#
   p_min=rmin1
   p_max=np.max(cs2_radius)
   cs2_p_arr = np.linspace(p_min,p_max,dim_p2)


   transmat, transmat_inv = calc_transmat2(thetac,phic)
   
   vec_slice = np.zeros(3)
   vec_cac = np.zeros(3)
   
   for i in range(0,dim_p2):
      for j in range(0,dim_zeta2):
         cs2_xcoord[j][i] = cs2_p_arr[i]*np.sin(cs2_zeta_arr[j])
         cs2_ycoord[j][i] = cs2_p_arr[i]*np.cos(cs2_zeta_arr[j])

         vec_slice[0] = cs2_xcoord[j][i]
         vec_slice[1] = cs2_ycoord[j][i]
         vec_slice[2] = offset2
         vec_cac = np.dot(transmat,vec_slice)
         radp=np.sqrt(vec_cac[0]*vec_cac[0] + vec_cac[1]*vec_cac[1] + vec_cac[2]*vec_cac[2])                  
#need to make sure that values on axes are treated the same
         if np.abs(vec_cac[0]) < 1.e-5: vec_cac[0]=0.
         if np.abs(vec_cac[1]) < 1.e-5: vec_cac[1]=0.
         if np.abs(vec_cac[2]) < 1.e-5: vec_cac[2]=0.
         
#get angles in spherical coordinates
         thetap, phip = angles_spc(vec_cac[0], vec_cac[1], vec_cac[2])
         #print(cs2_p_arr[i], cs2_zeta_arr[j], thetap, phip)
#find indices of cube-vertices for interpolation
         ir1, ir2, it1, it2, ip1, ip2, expol = get_rtp_indx(radp, thetap, phip, cs2_radius, cs2_theta, cs2_phi, cs2_nr, cs2_ntheta, cs2_nphi)
#
#----------------------interpolate 3d arrays onto slice-----------------
#
         rad1=cs2_radius[ir1]
         rad2=cs2_radius[ir2]
         theta1=cs2_theta[it1]
         theta2=cs2_theta[it2]
         phi1=cs2_phi[ip1]
         phi2=cs2_phi[ip2]
#
#density
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_rho3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_rho2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)

#opalbar
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_opalbar3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_opalbar2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)


#opac
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_opac3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_opac2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)


#sline
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_sline3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_sline2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)


#scont
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_scont3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_scont2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)

#temperature
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_t3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_t2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#velr
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_velr3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_velr2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#velth
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_velth3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_velth2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)
#velphi
         vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(cs2_nr, cs2_ntheta, cs2_nphi, cs2_velphi3d, ir1, ir2, it1, it2, ip1, ip2)
         cs2_velphi2d[j][i] = trilin_complete(radp, thetap, phip, rad1, rad2, theta1, theta2, phi1, phi2, \
                                           vala, valb, valc, vald, vale, valf, valg, valh, \
                                           rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2, radp, \
                                           expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=1)         

#
#-----------------combine to global coordinate system-------------------
#
   transmat, transmat_inv = calc_transmat2(thetac,phic)
   cs0_xcoord = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_ycoord = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_zcoord = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)      
   cs0_rho2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_opalbar2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_opac2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_sline2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_scont2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)
   cs0_t2d = np.zeros(dim_p1*dim_zeta1+dim_p2*dim_zeta2)         
#
   ii=0
   for i in range(0,dim_p1):
      for j in range(0,dim_zeta1):
         vec_cyc = np.array([cs1_xcoord[j,i],cs1_ycoord[j][i],offset1])   #in cylindrical system 1
         vec_cac = np.dot(transmat,vec_cyc)                  #in carthesian system 1
         vec_cac = vec_cac*rstar1/unit_length + np.array([x01,y01,z01]) #in global system
         vec_cyc = np.dot(transmat_inv,vec_cac)                #in global cylindrical system
         cs0_xcoord[ii] = vec_cyc[0]
         cs0_ycoord[ii] = vec_cyc[1]
         cs0_rho2d[ii] = cs1_rho2d[j][i]
         cs0_opalbar2d[ii] = cs1_opalbar2d[j][i]
         cs0_opac2d[ii] = cs1_opac2d[j][i]
         cs0_sline2d[ii] = cs1_sline2d[j][i]
         cs0_scont2d[ii] = cs1_scont2d[j][i]
         cs0_t2d[ii] = cs1_t2d[j][i]
         ii=ii+1

   for i in range(0,dim_p2):
      for j in range(0,dim_zeta2):
         vec_cyc = np.array([cs2_xcoord[j][i], cs2_ycoord[j][i],offset2])   #in cylindrical system 2
         vec_cac = np.dot(transmat,vec_cyc)                  #in carthesian system 2
         vec_cac = vec_cac*rstar2/unit_length + np.array([x02,y02,z02]) #in global system
         vec_cyc = np.dot(transmat_inv,vec_cac)                #in global cylindrical system
         cs0_xcoord[ii] = vec_cyc[0]
         cs0_ycoord[ii] = vec_cyc[1]
         cs0_rho2d[ii] = cs2_rho2d[j][i]
         cs0_opalbar2d[ii] = cs2_opalbar2d[j][i]
         cs0_opac2d[ii] = cs2_opac2d[j][i]
         cs0_sline2d[ii] = cs2_sline2d[j][i]
         cs0_scont2d[ii] = cs2_scont2d[j][i]
         cs0_t2d[ii] = cs2_t2d[j][i]
         ii=ii+1
#
#create triangulation
   cs0_triang = tri.Triangulation(cs0_xcoord, cs0_ycoord)
#
#-------------------density---------------------------------------------
#
   indx = np.where(cs1_rho2d > 0.)
   if np.size(indx) == 0:
      lrho2d = cs1_rho2d
      clevels, ticks = get_clevels(clim=[0.,1.0])
      ctitle=r'$\rho$'
   if np.size(indx) > 0:
      rho_min = 1.e-10*np.min(cs1_rho2d[indx])
      lrho2d = cs1_rho2d*1.
      indx = np.where(lrho2d <= 0.)
      lrho2d[indx] = rho_min
      lrho2d = np.log10(lrho2d)
      clevels, ticks = get_clevels(clim=clim_rho)      
      ctitle=r'$\log(\rho)$'

   cmap_lev = get_cmap('jet')      
   contourplot1 = ax2[0].contourf(cs1_xcoord, cs1_ycoord, lrho2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig2.colorbar(contourplot1,ax=ax2[0],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)

   indx = np.where(cs2_rho2d > 0.)
   if np.size(indx) == 0:
      lrho2d = cs2_rho2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$\rho$'        
   if np.size(indx) > 0:
      rho_min = 1.e-10*np.min(cs2_rho2d[indx])
      lrho2d = np.copy(cs2_rho2d)
      indx = np.where(lrho2d <= 0.)
      lrho2d[indx] = 1.e-100
      lrho2d = np.log10(lrho2d)
      clevels, ticks = get_clevels(clim=clim_rho)      
      ctitle=r'$\log(\rho)$'


   cmap_lev = get_cmap('jet')      
   contourplot2 = ax2[1].contourf(cs2_xcoord, cs2_ycoord, lrho2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot2.cmap.set_over('black')
   contourplot2.cmap.set_under('grey')
   contourplot2.changed()   
#
   cbar = fig2.colorbar(contourplot2,ax=ax2[1],orientation='vertical', ticks=ticks)
   cbar.set_label(ctitle)
#
#
#
   indx = np.where(cs0_rho2d > 0.)
   if np.size(indx) == 0:
      lrho2d = cs0_rho2d
      clevels, ticks = get_clevels(clim=[0.,1.])      
      ctitle=r'$\rho$'        
   if np.size(indx) > 0:
      rho_min = 1.e-10*np.min(cs0_rho2d[indx])
      lrho2d = np.copy(cs0_rho2d)
      indx = np.where(lrho2d <= 0.)
      lrho2d[indx] = 1.e-100
      lrho2d = np.log10(lrho2d)
      clevels, ticks = get_clevels(clim=clim_rho)      
      ctitle=r'$\log(\rho)$'


   cmap_lev = get_cmap('jet')      
   contourplot3 = ax2[2].tricontourf(cs0_triang, lrho2d,
                                     levels=clevels,
                                     extend='both',
                                     cmap=cmap_lev)
   contourplot3.cmap.set_over('black')
   contourplot3.cmap.set_under('grey')
   contourplot3.changed()   
#
   cbar = fig2.colorbar(contourplot3,ax=ax2[2],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)   
#
#-------------------line opacity----------------------------------------
#
   indx = np.where(cs1_opalbar2d > 0.)
   if np.size(indx) == 0:
      lopalbar2d = cs1_opalbar2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$\bar{\chi}_L$'      
   if np.size(indx) > 0:
      opalbar_min = 1.e-10*np.min(cs1_opalbar2d[indx])
      lopalbar2d = cs1_opalbar2d*1.
      indx = np.where(lopalbar2d <= 0.)
#      print(np.min(cs1_opalbar2d))
#      print(np.min(lopalbar2d))
      lopalbar2d[indx] = opalbar_min
#      print(np.min(lopalbar2d))
      lopalbar2d = np.log10(lopalbar2d)
#      print(np.min(lopalbar2d))
      clevels, ticks = get_clevels(clim=clim_opalbar)            
      ctitle=r'$\log(\bar{\chi}_L)$'

#   print('test',np.min(lopalbar2d), opalbar_min)
#   exit()

   cmap_lev = get_cmap('jet')
   contourplot1 = ax3[0].contourf(cs1_xcoord, cs1_ycoord, lopalbar2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig3.colorbar(contourplot1,ax=ax3[0],orientation='vertical', ticks=ticks)
   cbar.set_label(ctitle)


   
   indx = np.where(cs2_opalbar2d > 0.)
   if np.size(indx) == 0:
      lopalbar2d = cs2_opalbar2d
      clevels, ticks = get_clevels(clim=[0.,1.])      
      ctitle=r'$\bar{\chi}_L$'      
   if np.size(indx) > 0:
      opalbar_min = 1.e-10*np.min(cs2_opalbar2d[indx])
      lopalbar2d = cs2_opalbar2d*1.
      indx = np.where(lopalbar2d <= 0.)
      lopalbar2d[indx] = opalbar_min
      lopalbar2d = np.log10(lopalbar2d)
      clevels, ticks = get_clevels(clim=clim_opalbar)
      ctitle=r'$\log(\bar{\chi}_L)$'


   cmap_lev = get_cmap('jet')      
   contourplot2 = ax3[1].contourf(cs2_xcoord, cs2_ycoord, lopalbar2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot2.cmap.set_over('black')
   contourplot2.cmap.set_under('grey')
   contourplot2.changed()

   cbar = fig3.colorbar(contourplot2,ax=ax3[1],orientation='vertical', ticks=ticks)
   cbar.set_label(ctitle)
#
#
   indx = np.where(cs0_opalbar2d > 0.)
   if np.size(indx) == 0:
      lopalbar2d = cs0_opalbar2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$\bar{\chi}_L$'      
   if np.size(indx) > 0:
      opalbar_min = 1.e-10*np.min(cs0_opalbar2d[indx])
      lopalbar2d = np.copy(cs0_opalbar2d)
      indx = np.where(lopalbar2d <= 0.)
      lopalbar2d[indx] = 1.e-100
      lopalbar2d = np.log10(lopalbar2d)
      clevels, ticks = get_clevels(clim=clim_opalbar)
      ctitle=r'$\log(\bar{\chi}_L)$'


   cmap_lev = get_cmap('jet')
   
   contourplot3 = ax3[2].tricontourf(cs0_triang, lopalbar2d,
                                     levels=clevels,
                                     extend='both',
                                     cmap=cmap_lev)
   contourplot3.cmap.set_over('black')
   contourplot3.cmap.set_under('grey')
   contourplot3.changed()   
#
   cbar = fig3.colorbar(contourplot3,ax=ax3[2],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)      
#
#-------------------line source function--------------------------------
#
   indx = np.where(cs1_sline2d > 0.)
   if np.size(indx) == 0:
      lsline2d = cs1_sline2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$S_L$'
   if np.size(indx) > 0:
      sline_min = 1.e-5*np.min(cs1_sline2d[indx])
      lsline2d = cs1_sline2d*1.
      indx = np.where(lsline2d <= 0.)
      lsline2d[indx] = sline_min
      lsline2d = np.log10(lsline2d)
      clevels, ticks = get_clevels(clim=clim_sline)
      ctitle=r'$\log(S_L)$'


   cmap_lev = get_cmap('jet')      
   contourplot1 = ax4[0].contourf(cs1_xcoord, cs1_ycoord, lsline2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig4.colorbar(contourplot1,ax=ax4[0],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)


   
   indx = np.where(cs2_sline2d > 0.)
   if np.size(indx) == 0:
      lsline2d = cs2_sline2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$S_L$'
   if np.size(indx) > 0:
      sline_min = 1.e-5*np.min(cs2_sline2d[indx])
      lsline2d = cs2_sline2d*1.
      indx = np.where(lsline2d <= 0.)
      lsline2d[indx] = sline_min
      lsline2d = np.log10(lsline2d)
      clevels, ticks = get_clevels(clim=clim_sline)
      ctitle=r'$\log(S_L)$'


   cmap_lev = get_cmap('jet')      
   contourplot2 = ax4[1].contourf(cs2_xcoord, cs2_ycoord, lsline2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot2.cmap.set_over('black')
   contourplot2.cmap.set_under('grey')
   contourplot2.changed()

   cbar = fig4.colorbar(contourplot2,ax=ax4[1],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)


   
   indx = np.where(cs0_sline2d > 0.)
   if np.size(indx) == 0:
      lsline2d = cs0_sline2d
      clevels, ticks = get_clevels(clim=[0.,1.])
      ctitle=r'$S_L$'
   if np.size(indx) > 0:
      sline_min = 1.e-5*np.min(cs0_sline2d[indx])
      lsline2d = np.copy(cs0_sline2d)
      indx = np.where(lsline2d <= 0.)
      lsline2d[indx] = 1.e-100
      lsline2d = np.log10(lsline2d)
      clevels, ticks = get_clevels(clim=clim_sline)
      ctitle=r'$\log(S_L)$'


   cmap_lev = get_cmap('jet')      
   contourplot3 = ax4[2].tricontourf(cs0_triang, lsline2d,
                                     levels=clevels,
                                     extend='both',
                                     cmap=cmap_lev)
   contourplot3.cmap.set_over('black')
   contourplot3.cmap.set_under('grey')
   contourplot3.changed()   
#
   cbar = fig4.colorbar(contourplot3,ax=ax4[2],orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)   
#
#-------------------temperature-----------------------------------------
#
#  indx = np.where(rho2d > 0.)
   clevels, ticks = get_clevels(clim=clim_temp)


   cmap_lev = get_cmap('jet')   
   contourplot1= ax5[0].contourf(cs1_xcoord, cs1_ycoord, cs1_t2d,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig5.colorbar(contourplot1,ax=ax5[0],orientation='vertical',ticks=ticks)
   cbar.set_label(r'$T$')


   cmap_lev = get_cmap('jet')   
   contourplot2= ax5[1].contourf(cs2_xcoord, cs2_ycoord, cs2_t2d,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)
   contourplot2.cmap.set_over('black')
   contourplot2.cmap.set_under('grey')
   contourplot2.changed()

   cbar = fig5.colorbar(contourplot2,ax=ax5[1],orientation='vertical',ticks=ticks)
   cbar.set_label(r'$T$')



   cmap_lev = get_cmap('jet')   
   contourplot3 = ax5[2].tricontourf(cs0_triang, cs0_t2d,
                                     levels=clevels,
                                     extend='both',
                                     cmap=cmap_lev)
   contourplot3.cmap.set_over('black')
   contourplot3.cmap.set_under('grey')
   contourplot3.changed()   
#
   cbar = fig5.colorbar(contourplot3,ax=ax5[2],orientation='vertical',ticks=ticks)
   cbar.set_label(r'$T$')

   sdum = input("Press [q] to exit.")
   if sdum == 'q':
       exit()
#
#-------------------radial velocity-------------------------------------
#
   clevels, ticks = get_clevels(clim=clim_velr)   

   cmap_lev = get_cmap('jet')
   contourplot1 = ax6[0].contourf(cs1_xcoord, cs1_ycoord, cs1_velr2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot2 = ax6[1].contourf(cs2_xcoord, cs2_ycoord, cs2_velr2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap='jet')
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()
   
   contourplot2.cmap.set_over('black')
   contourplot2.cmap.set_under('grey')
   contourplot2.changed()   

   cbar = fig6.colorbar(contourplot1,ax=ax6[0],orientation='vertical',ticks=ticks)
   cbar.set_label(r'$v_r$')
   cbar = fig6.colorbar(contourplot2,ax=ax6[1],orientation='vertical',ticks=ticks)
   cbar.set_label(r'$v_r$')      

#
#------------------------output to files-------------------------------
#
   oname1 = oname3d+'_star1_radial.png'
   oname2 = oname3d+'_star1_radial.ps'
   fig0.savefig(oname1, bbox_inches='tight')
   fig0.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_star2_radial.png'
   oname2 = oname3d+'_star2_radial.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')   

   oname1 = oname3d+'_rho.png'
   oname2 = oname3d+'_rho.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_opalbar.png'
   oname2 = oname3d+'_opalbar.ps'
   fig3.savefig(oname1, bbox_inches='tight')
   fig3.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_sline.png'
   oname2 = oname3d+'_sline.ps'
   fig4.savefig(oname1, bbox_inches='tight')
   fig4.savefig(oname2, bbox_inches='tight')      

   oname1 = oname3d+'_temp.png'
   oname2 = oname3d+'_temp.ps'
   fig5.savefig(oname1, bbox_inches='tight')
   fig5.savefig(oname2, bbox_inches='tight')
   
   oname1 = oname3d+'_velr.png'
   oname2 = oname3d+'_velr.ps'
   fig6.savefig(oname1, bbox_inches='tight')
   fig6.savefig(oname2, bbox_inches='tight')
#

def main_test(fname=None, version=None, clim_sline=None, clim_rho=None, clim_opalbar=None, clim_temp=None, clim_velr=None, clim_velth=None, clim_velphi=None, windx=1):
#
#------------------read all information from hdf5-file------------------
#
   iline, xnue0, vth_fiducial, vmax, unit_length, \
         x01, y01, z01, vx01, vy01, vz01, rstar1, teff1, logg1, lstar1, yhe1, vrot1, vmicro1, \
         x02, y02, z02, vx02, vy02, vz02, rstar2, teff2, logg2, lstar2, yhe2, vrot2, vmicro2, \
         cs1_nr,  cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta, cs1_phi, \
         cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d, cs1_rho3d, cs1_opac3d, cs1_opalbar3d, \
         cs1_scont3d, cs1_sline3d, \
         cs2_nr,  cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta, cs2_phi, \
         cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d, cs2_rho3d, cs2_opac3d, cs2_opalbar3d, \
      cs2_scont3d, cs2_sline3d = read_modelspec3d_vbin(fname,version=version)

   print('model ', fname)
   print('iline,xnue0,vth_fiducial,vmax,unit_length', iline,xnue0,vth_fiducial,vmax,unit_length)
   print('x01,y01,z01',x01,y01,z01)
   print('vx01,vy01,vz01',vx01,vy01,vz01)
   print('rstar1,teff1,logg1,lstar1,yhe1',rstar1,teff1,logg1,lstar1,yhe1)
   print('vrot1,vmicro1',vrot1,vmicro1)
   print('cs1_nr,cs1_ntheta,cs1_nphi',cs1_nr,cs1_ntheta,cs1_nphi)
   print('')
   print('x02,y02,z02',x02,y02,z02)
   print('vx02,vy02,vz02',vx02,vy02,vz02)
   print('rstar2,teff2,logg2,lstar2,yhe2',rstar2,teff2,logg2,lstar2,yhe2)
   print('vrot2,vmicro2',vrot2,vmicro2)
   print('cs2_nr,cs2_ntheta,cs2_nphi',cs2_nr,cs2_ntheta,cs2_nphi)
   print('')
   print('----------------------')
   print('')

   sr1=rstar1*cgs_rsu

   cs1_velr3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   cs1_velth3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   cs1_velphi3d=np.zeros(shape=(cs1_nphi,cs1_ntheta,cs1_nr))
   for k in range(0,cs1_nphi):
      sinp=np.sin(cs1_phi[k])
      cosp=np.cos(cs1_phi[k])      
      for j in range(0,cs1_ntheta):
         sint=np.sin(cs1_theta[j])
         cost=np.cos(cs1_theta[j])
         velx=cs1_velx3d[k][j][:]
         vely=cs1_vely3d[k][j][:]
         velz=cs1_velz3d[k][j][:]         
         cs1_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
         cs1_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
         cs1_velphi3d[k][j][:] = -velx*sinp + vely*cosp         
   
#velocity in km/s
   cs1_velr3d=cs1_velr3d/1.e5
   cs1_velth3d=cs1_velth3d/1.e5
   cs1_velphi3d=cs1_velphi3d/1.e5




   sr2=rstar2*cgs_rsu

   cs2_velr3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   cs2_velth3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   cs2_velphi3d=np.zeros(shape=(cs2_nphi,cs2_ntheta,cs2_nr))
   for k in range(0,cs2_nphi):
      sinp=np.sin(cs2_phi[k])
      cosp=np.cos(cs2_phi[k])      
      for j in range(0,cs2_ntheta):
         sint=np.sin(cs2_theta[j])
         cost=np.cos(cs2_theta[j])
         velx=cs2_velx3d[k][j][:]
         vely=cs2_vely3d[k][j][:]
         velz=cs2_velz3d[k][j][:]         
         cs2_velr3d[k][j][:] = velx*sint*cosp + vely*sint*sinp + velz*cost
         cs2_velth3d[k][j][:] = velx*cost*cosp + vely*cost*sinp - velz*sint            
         cs2_velphi3d[k][j][:] = -velx*sinp + vely*cosp         
   
#velocity in km/s
   cs2_velr3d=cs2_velr3d/1.e5
   cs2_velth3d=cs2_velth3d/1.e5
   cs2_velphi3d=cs2_velphi3d/1.e5   
#
#***********************************************************************
#
#prepare contour plots
   xsize=18. #in cm
   xsize=xsize/1.5 #/2.54 #in inches
   aspect_ratio=2.
   ysize=xsize/aspect_ratio
 
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig3 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig4 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   fig5 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1      

   ax2 = fig2.subplots()
   ax3 = fig3.subplots()
   ax4 = fig4.subplots()
   ax5 = fig5.subplots()      
#
#--------------calculate scalar array on slice for star 2---------------
#
#for star 2
   cs2_rho2d=np.zeros(shape=(cs2_ntheta, cs2_nr))
   cs2_velr2d=np.zeros(shape=(cs2_ntheta, cs2_nr))
   cs2_velrho2d=np.zeros(shape=(cs2_ntheta, cs2_nr))   
   cs2_velz2d=np.zeros(shape=(cs2_ntheta, cs2_nr))   
   cs2_velth2d=np.zeros(shape=(cs2_ntheta, cs2_nr))
   cs2_velphi2d=np.zeros(shape=(cs2_ntheta, cs2_nr))         
   cs2_opalbar2d=np.zeros(shape=(cs2_ntheta, cs2_nr))   
#
   cs2_xcoord=np.zeros(shape=(cs2_ntheta, cs2_nr))
   cs2_ycoord=np.zeros(shape=(cs2_ntheta, cs2_nr))

   k=11
   for i in range(0,cs2_nr):
      for j in range(0,cs2_ntheta):
         sint = np.sin(cs2_theta[j])
         cost = np.cos(cs2_theta[j])
         sinp = np.sin(cs2_phi[k])
         cosp = np.cos(cs2_phi[k])
         
         cs2_xcoord[j][i] = cs2_radius[i]*np.sin(cs2_theta[j])
         cs2_ycoord[j][i] = cs2_radius[i]*np.cos(cs2_theta[j])

         cs2_rho2d[j][i] =cs2_rho3d[k][j][i]
         cs2_opalbar2d[j][i] =cs2_opalbar3d[k][j][i]

         cs2_velr2d[j][i] = cs2_velr3d[k][j][i]         
         cs2_velth2d[j][i] = cs2_velth3d[k][j][i]
         cs2_velphi2d[j][i] = cs2_velphi3d[k][j][i]

         cs2_velrho2d[j][i] = cosp*cs2_velx3d[k][j][i]/1.e5 + sinp*cs2_vely3d[k][j][i]/1.e5
         cs2_velz2d[j][i] = cs2_velz3d[k][j][i]/1.e5

#----------------------------------------------------------------------
         
   cs2_lrho2d = get_logvals(cs2_rho2d)
   clevels, ticks = get_clevels(clim=clim_rho)
   ctitle=r'$\log(\rho)$'

   cmap_lev = get_cmap('jet')      
   contourplot1 = ax2.contourf(cs2_xcoord, cs2_ycoord, cs2_lrho2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig2.colorbar(contourplot1,ax=ax2,orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)

#----------------------------------------------------------------------
         
   clevels, ticks = get_clevels(clim=clim_velr)
   ctitle=r'$v_z$'

   cmap_lev = get_cmap('jet')      
   contourplot1 = ax3.contourf(cs2_xcoord, cs2_ycoord, cs2_velz2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig3.colorbar(contourplot1,ax=ax3,orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)   
#----------------------------------------------------------------------
         
   clevels, ticks = get_clevels(clim=clim_velth)
   ctitle=r'$v_\rho$'

   cmap_lev = get_cmap('jet')      
   contourplot1 = ax4.contourf(cs2_xcoord, cs2_ycoord, cs2_velrho2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig4.colorbar(contourplot1,ax=ax4,orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)


#----------------------------------------------------------------------
         
   clevels, ticks = get_clevels(clim=clim_velphi)
   ctitle=r'$v_\phi$'

   cmap_lev = get_cmap('jet')      
   contourplot1 = ax5.contourf(cs2_xcoord, cs2_ycoord, cs2_velphi2d,
                                  levels=clevels,
                                  extend='both',
                                  cmap=cmap_lev)
   contourplot1.cmap.set_over('black')
   contourplot1.cmap.set_under('grey')
   contourplot1.changed()

   cbar = fig5.colorbar(contourplot1,ax=ax5,orientation='vertical',ticks=ticks)
   cbar.set_label(ctitle)      

   plt.show()
   
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
   
#fnames=['../outputFILES/modspec_model00.h5']
#fnames=['/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p02s12xa3_halpha/modspec_model00.h5']
#fnames=['/lhome/levin/Postdoc/papers/paperIII/models/lb1/model_p05s10xa3_halpha/modspec_model00.h5']
#versions=['v01']
#windx = 1
#main(fnames=fnames,versions=versions,windx=windx,thetac=90.*np.pi/180.,phic=90.*np.pi/180., clim_sline=[-6.,-5.], clim_rho=[-16.,-9.], clim_opalbar=[-18.,-6.], clim_temp=[5.e3,1.e4], clim_velr=[0.,410.], offset=0.)



fnames=['../outputFILES/models_olivier/photprof0_sline0/phase001/modspec_model00.h5']
versions=['v01']
windx = 1

main_test(fname=fnames[0], version=versions[0], clim_sline=[-6.,-5.], clim_rho=[-16.,-9.], clim_opalbar=[-18.,-6.], clim_temp=[2.e3,1.e4], clim_velr=[-300.,300.], clim_velphi=[0.,100.], clim_velth=[0.,100.])

#main(fnames=fnames,versions=versions,windx=windx,thetac=90.*np.pi/180.,phic=90.*np.pi/180., clim_sline=[-6.,-5.], clim_rho=[-16.,-9.], clim_opalbar=[-18.,-6.], clim_temp=[2.e3,1.e4], clim_velr=[0.,300.], offset=0.)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
