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


########################################################################
   
def main(fname='../outputFILES/output_model00.h5', oname='./ps_files/sc3d', windx=1):
#
#plots conergence behaviour of ALI method and source functions
#
#------------------read all information from hdf5-file------------------
#
   spatial_grid3d, spatial_grid1d, input_mod_dim, opt_method, \
      opt_opal, opt_opac, opt_angint_method, opt_ltec, opt_sol2d, opt_incl_cont, \
      opt_start_cont, opt_ng_cont, opt_ait_cont, opt_incl_line, opt_start_line, \
      opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, opt_incl_gdark, \
      opt_incl_sdist, \
      kline, kcont, alpha, kappa0, eps_line, eps_cont, teff, xlogg, \
      trad, xnue0, rstar, lstar, vth_fiducial, vmicro, yhe, hei, mdot, \
      na, vmin, vmax, beta, vrot, \
      xic1, ntheta_gdark, theta_gdark, teff_gdark, xic1_gdark, \
      nx, ny, nz, x, y, z, mask3d, maskb3d, \
      opac3d, opalbar3d, velx3d, vely3d, velz3d, t3d, \
      scont3d, sline3d, ssobo3d, mint3d, mintbar3d, \
      fcontx3d, fconty3d, fcontz3d, \
      kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d, \
      itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, \
      nomega, nnue, n_x, n_y, n_z, nodes_xobs = read_sc3d(fname)

   print('----------------------------------------')
   print('{label:20s} {fname:20s}'.format(label='fname', fname=fname))
   print()
   print('--------------options-------------------')
   print('{label:20s} {value:20d}'.format(label='input_mod_dim', value=input_mod_dim))
   print('{label:20s} {value:20d}'.format(label='spatial_grid1d', value=spatial_grid1d))
   print('{label:20s} {value:20d}'.format(label='spatial_grid3d', value=spatial_grid3d))
   print('{label:20s} {value:20d}'.format(label='opt_method', value=opt_method))
   print('{label:20s} {value:20d}'.format(label='opt_opal', value=opt_opal))
   print('{label:20s} {value:20d}'.format(label='opt_angint_method', value=opt_angint_method))
   print('{label:20s} {value:20d}'.format(label='opt_sol2d', value=opt_sol2d))
   print('{label:20s} {value:20d}'.format(label='opt_incl_cont', value=opt_incl_cont))
   print('{label:20s} {value:20d}'.format(label='opt_start_cont', value=opt_start_cont))
   print('{label:20s} {value:20d}'.format(label='opt_ng_cont', value=opt_ng_cont))
   print('{label:20s} {value:20d}'.format(label='opt_ait_cont', value=opt_ait_cont))
   print('{label:20s} {value:20d}'.format(label='opt_incl_line', value=opt_incl_line))
   print('{label:20s} {value:20d}'.format(label='opt_start_line', value=opt_start_line))
   print('{label:20s} {value:20d}'.format(label='opt_ng_line', value=opt_ng_line))
   print('{label:20s} {value:20d}'.format(label='opt_ait_line', value=opt_ait_line))
   print('{label:20s} {value:20d}'.format(label='opt_alo_cont', value=opt_alo_cont))
   print('{label:20s} {value:20d}'.format(label='opt_alo_line', value=opt_alo_line))
   print()
   print('-------input/stellar parameters---------')
   print('{label:20s} {value:20.8e}'.format(label='xic1', value=xic1))
   print('{label:20s} {value:20.8e}'.format(label='kcont', value=kcont))
   print('{label:20s} {value:20.8e}'.format(label='kline', value=kline))
   print('{label:20s} {value:20.8e}'.format(label='eps_cont', value=eps_cont))
   print('{label:20s} {value:20.8e}'.format(label='eps_line', value=eps_line))
   print('{label:20s} {value:20.8e}'.format(label='teff', value=teff))
   print('{label:20s} {value:20.8e}'.format(label='trad', value=trad))
   print('{label:20s} {value:20.8e}'.format(label='xnue0', value=xnue0))
   print('{label:20s} {value:20.8e}'.format(label='rstar', value=rstar))
   print('{label:20s} {value:20.8e}'.format(label='vth_fiducial', value=vth_fiducial))
   print('{label:20s} {value:20.8e}'.format(label='vmicro', value=vmicro))
   print('{label:20s} {value:20.8e}'.format(label='yhe', value=yhe))
   print('{label:20s} {value:20.8e}'.format(label='hei', value=hei))
   print('{label:20s} {value:20.8e}'.format(label='mdot', value=mdot))
   print('{label:20s} {value:20.8e}'.format(label='na', value=na))
   print('{label:20s} {value:20.8e}'.format(label='vmin', value=vmin))
   print('{label:20s} {value:20.8e}'.format(label='vmax', value=vmax))
   print('{label:20s} {value:20.8e}'.format(label='beta', value=beta))
   print()
   print('------------dimensions------------------')
   print('{label:20s} {value:20d}'.format(label='nx', value=nx))   
   print('{label:20s} {value:20d}'.format(label='ny', value=ny))   
   print('{label:20s} {value:20d}'.format(label='nz', value=nz))   
   print('{label:20s} {value:20d}'.format(label='nomega', value=nomega))   
   print('{label:20s} {value:20d}'.format(label='nnue', value=nnue))   
   print('----------------------------------------')

#------------------prepare 3d radial array and normalization------------

#normalize everything to xic1
   sline3d=sline3d/xic1
   ssobo3d=ssobo3d/xic1
   scont3d=scont3d/xic1

   #calculate dilution factor and define radius
   sline3d_dilfac = np.zeros(shape=(nx,ny,nz))
   r3d = np.zeros(shape=(nx,ny,nz))
   for i in np.arange(nx):
      for j in np.arange(ny):
         for k in np.arange(nz):
            rad=np.sqrt(x[i]**2+y[j]**2+z[k]**2)
            if rad >= 1.:
               sline3d_dilfac[k][j][i] = dilfac(1.,rad)

            if mask3d[k][j][i] == 1:
               r3d[k][j][i] = rad
            elif mask3d[k][j][i] == 2:
               r3d[k][j][i] = rad
            elif mask3d[k][j][i] == 3:
               r3d[k][j][i] = rad
            else:
               r3d[k][j][i] = 1.
#
#----------------------prepare convergence behaviour--------------------
#
   epsmaxl_arr=np.abs(epsmaxl_arr)
   epsmaxc_arr=np.abs(epsmaxc_arr)
   iternrl_arr=np.arange(itmaxl)
   iternrc_arr=np.arange(itmaxc)
#
#-----------------------define range------------------------------------
#
   xmin=np.min(r3d)
   xmax=np.max(r3d)
   xlim=[xmin, xmax]


   indx=np.where(r3d >= 1.)   
   #indx = np.where(mask3d == 1 or mask3d == 1 or mask3d == 2 or mask3d == 3)
#
#----------------------------title-strings------------------------------
#
   if opt_opal == 0:
      titlestr=r'$k_L=${kline:10.2f}, $k_c$={kcont:10.2f}'.format(kline=kline,kcont=kcont)
   else:
      titlestr=r'$\kappa_0=${kappa0:9.5f}, $\alpha$={alpha:9.5f}, $k_c=${kcont:f9.5}'.format(kappa0=kappa0,alpha=alpha,kcont=kcont)

#
#***********************************************************************
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

   fig0 = plt.figure(windx, constrained_layout=True)
   windx=windx+1   
   plt.ion()
   plt.show()
   ax0 = fig0.subplots(3,2)

   fig0.suptitle(titlestr)
#
#continuum source function
   ax0[0,0].set_xlabel(r'$r/R_\ast$')
   ax0[0,0].set_ylabel(r'$S_c*r^2/I_c$')
   ax0[0,0].scatter(r3d[indx],scont3d[indx]*r3d[indx]**2, marker='.', s=1, color='blue')

#line source function
   ax0[0,1].set_xlabel(r'$r/R_\ast$')
   ax0[0,1].set_ylabel(r'$S_L*r^2/I_c$')
   ax0[0,1].scatter(r3d[indx],sline3d[indx]*r3d[indx]**2, marker='.', s=1, color='blue')      


#convergence behaviour continuum
   xmax=np.where(epsmaxc_arr > 0.)
   if np.size(xmax) == 0:
      ymin = 1.e-5
      ymax = 1.
      xmin=0.
      xmax=1
   else:
      xmax = np.max(xmax)      
      ymin = epsmaxc_arr[xmax]
      ymax=np.max(epsmaxc_arr)
      xmin = 0.

   ylim=[ymin,ymax]
   xlim=[xmin,xmax]
   ax0[1,0].set_xlabel(r'# iterations')
   ax0[1,0].set_ylabel(r'$((S_c^{(k-1)} - S_c^{(k)})/S_c^{(k)})_{\rm{max}}$')
   ax0[1,0].plot(iternrc_arr, epsmaxc_arr, color='blue')
   ax0[1,0].set_yscale('log')
   ax0[1,0].set_xlim(xlim)
   ax0[1,0].set_ylim(ylim)

#convergence behaviour line   
   xmax=np.where(epsmaxl_arr > 0.)
   if np.size(xmax) == 0:
      ymin = 1.e-5
      ymax = 1.
      xmin=0.
      xmax=1
   else:
      xmax = np.max(xmax)      
      ymin = np.min(epsmaxl_arr[xmax])
      ymax=np.max(epsmaxl_arr)
      xmin = 0.
   ylim=[ymin,ymax]
   xlim=[xmin,xmax]
   ax0[1,1].set_xlabel(r'# iterations')
   ax0[1,1].set_ylabel(r'$((S_L^{(k-1)} - S_L^{(k)})/S_L^{(k)})_{\rm{max}}$')
   ax0[1,1].plot(iternrl_arr, epsmaxl_arr, color='blue')
   ax0[1,1].set_yscale('log')
   ax0[1,1].set_xlim(xlim)
   ax0[1,1].set_ylim(ylim)   


#continuum opacity
   indx=np.where(opac3d > 0.)
   ylim=[np.min(opac3d[indx]),np.max(opac3d[indx])]
   ax0[2,0].set_xlabel(r'$r/R_\ast$')
   ax0[2,0].set_ylabel(r'$\chi_c$')
   ax0[2,0].scatter(r3d[indx],opac3d[indx], marker='.', s=1, color='blue')
   ax0[2,0].set_yscale('log')
   ax0[2,0].set_ylim(ylim)   

#line opacity
   indx=np.where(opalbar3d > 0.)
   ylim=[np.min(opalbar3d[indx]),np.max(opalbar3d[indx])]
   ax0[2,1].set_xlabel(r'$r/R_\ast$')
   ax0[2,1].set_ylabel(r'$\bar{\chi}_L$')
   ax0[2,1].set_yscale('log')      
   ax0[2,1].scatter(r3d,opalbar3d, marker='.', s=1, color='blue')
   ax0[2,1].set_ylim(ylim)
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

fname1='../TRASH/model3d.h5'
fname2='../inputFILES/model3d.h5'

#compare_model3d(fname1,fname2)


fname='../outputFILES/output_model00.h5'
fname='../outputFILES/test_js_lh/output_model00.h5'
fname='../outputFILES/ablation/output_model00.h5'

windx = 1
main(fname=fname,windx=windx)

         
sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
