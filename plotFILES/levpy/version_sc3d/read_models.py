import numpy as np
import h5py
#
def read_model1d(fname):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   nr = group_id.attrs['nr']
   nr = nr[0]

   group_id = file_id['coordinates']
   radius = np.array(group_id['r'])
   
   group_id = file_id['model']
   rho1d = np.array(group_id['rho'])
   velr1d = np.array(group_id['velr'])
   t1d = np.array(group_id['temperature'])
   vth1d = np.array(group_id['vth'])
   

   file_id.close()

   return nr, radius, rho1d, velr1d, t1d, vth1d

###################################################################


def read_model2d(fname):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   nr = group_id.attrs['nr']
   nr = nr[0]
   ntheta = group_id.attrs['ntheta']
   ntheta = ntheta[0]

   group_id = file_id['coordinates']
   radius = np.array(group_id['r'])
   theta = np.array(group_id['theta'])
   
   group_id = file_id['model']
   rho2d = np.array(group_id['rho'])
   velr2d = np.array(group_id['velr'])
   velth2d = np.array(group_id['velth'])
   velphi2d = np.array(group_id['velphi'])      
   t2d = np.array(group_id['temperature'])
   vth2d = np.array(group_id['vth'])
   
   file_id.close()

   return nr, ntheta, radius, theta, rho2d, velr2d, velth2d, velphi2d, t2d, vth2d


###################################################################


def read_model3d(fname):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   nr = group_id.attrs['nr']
   nr = nr[0]
   ntheta = group_id.attrs['ntheta']
   ntheta = ntheta[0]
   nphi = group_id.attrs['nphi']
   nphi = nphi[0]

   group_id = file_id['coordinates']
   radius = np.array(group_id['r'])
   theta = np.array(group_id['theta'])
   phi = np.array(group_id['phi'])
   
   group_id = file_id['model']
   rho3d = np.array(group_id['rho'])
   velr3d = np.array(group_id['velr'])
   velth3d = np.array(group_id['velth'])
   velphi3d = np.array(group_id['velphi'])      
   t3d = np.array(group_id['temperature'])
   vth3d = np.array(group_id['vth'])
   
   file_id.close()

   return nr, ntheta, nphi, radius, theta, phi, rho3d, velr3d, velth3d, velphi3d, t3d, vth3d


###################################################################


def read_model_slab3d(fname):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   nx = group_id.attrs['nx']
   nx = nx[0]
   ny = group_id.attrs['ny']
   ny = ny[0]
   nz = group_id.attrs['nz']
   nz = nz[0]

   group_id = file_id['units']
   unit_length = group_id.attrs['unit_length']
   unit_length = unit_length[0]
   unit_density = group_id.attrs['unit_density']
   unit_density = unit_density[0]
   unit_velocity = group_id.attrs['unit_velocity']
   unit_velocity = unit_velocity[0]
   unit_temperature = group_id.attrs['unit_temperature']
   unit_temperature = unit_temperature[0]   

   group_id = file_id['coordinates']
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])
   
   group_id = file_id['model']
   rho3d = np.array(group_id['rho'])
   velx3d = np.array(group_id['velx'])
   vely3d = np.array(group_id['vely'])
   velz3d = np.array(group_id['velz'])      
   tgas3d = np.array(group_id['tgas'])
   trad3d = np.array(group_id['trad'])   
   
   file_id.close()

   return nx, ny, nz, x, y, z, rho3d, velx3d, vely3d, velz3d, tgas3d, trad3d, unit_length, unit_density, unit_velocity, unit_temperature

###################################################################

def read_sc3d_v0(fname, read='all'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['options']
   spatial_grid3d = group_id.attrs['spatial_grid3d']
   spatial_grid3d = spatial_grid3d[0]
   spatial_grid1d = group_id.attrs['spatial_grid1d']
   spatial_grid1d = spatial_grid1d[0]
   input_mod_dim = group_id.attrs['input_mod_dim']
   input_mod_dim = input_mod_dim[0]
   opt_method = group_id.attrs['opt_method']
   opt_method = opt_method[0]
   opt_opal = group_id.attrs['opt_opal']
   opt_opal = opt_opal[0]
   opt_angint_method = group_id.attrs['opt_angint_method']
   opt_angint_method = opt_angint_method[0]
   opt_sol2d = group_id.attrs['opt_sol2d']
   opt_sol2d = opt_sol2d[0]
   opt_incl_cont = group_id.attrs['opt_incl_cont']
   opt_incl_cont = opt_incl_cont[0]   
   opt_start_cont = group_id.attrs['opt_start_cont']
   opt_start_cont = opt_start_cont[0]
   opt_ng_cont = group_id.attrs['opt_ng_cont']
   opt_ng_cont = opt_ng_cont[0]
   opt_ait_cont = group_id.attrs['opt_ait_cont']
   opt_ait_cont = opt_ait_cont[0]
   opt_incl_line = group_id.attrs['opt_incl_line']
   opt_incl_line = opt_incl_line[0]   
   opt_start_line = group_id.attrs['opt_start_line']
   opt_start_line = opt_start_line[0]   
   opt_ng_line = group_id.attrs['opt_ng_line']
   opt_ng_line = opt_ng_line[0]
   opt_ait_line = group_id.attrs['opt_ait_line']
   opt_ait_line = opt_ait_line[0]   
   opt_alo_cont = group_id.attrs['opt_alo_cont']
   opt_alo_cont = opt_alo_cont[0]
   opt_alo_line = group_id.attrs['opt_alo_line']
   opt_alo_line = opt_alo_line[0]
   opt_incl_gdark = group_id.attrs['opt_incl_gdark']
   opt_incl_gdark = opt_incl_gdark[0]   
   opt_incl_sdist = group_id.attrs['opt_incl_sdist']
   opt_incl_sdist = opt_incl_sdist[0]   

   group_id = file_id['dimensions']
   nx = group_id.attrs['ndxmax']
   nx = nx[0]
   ny = group_id.attrs['ndymax']
   ny = ny[0]
   nz = group_id.attrs['ndzmax']
   nz = nz[0]
   nomega = group_id.attrs['dim_omega']
   nomega = nomega[0]
   nnue = group_id.attrs['nxobs']
   nnue = nnue[0]

   group_id = file_id['input_parameters']
   kline = group_id.attrs['kline']
   kline = kline[0]
   kcont = group_id.attrs['kcont']
   kcont = kcont[0]
   alpha = group_id.attrs['alpha']
   alpha = alpha[0]
   kappa0 = group_id.attrs['kappa0']
   kappa0 = kappa0[0]
   eps_line = group_id.attrs['eps_line']
   eps_line = eps_line[0]
   eps_cont = group_id.attrs['eps_cont']
   eps_cont = eps_cont[0]
   teff = group_id.attrs['teff']
   teff = teff[0]
   xlogg = group_id.attrs['xlogg']
   xlogg = xlogg[0]
   trad = group_id.attrs['trad']
   trad = trad[0]
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]   
   rstar = group_id.attrs['rstar']
   rstar = rstar[0]
   lstar = group_id.attrs['lstar']
   lstar = lstar[0]
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]
   vmicro = group_id.attrs['vmicro']
   vmicro = vmicro[0]
   yhe = group_id.attrs['yhe']
   yhe = yhe[0]            
   hei = group_id.attrs['hei']
   hei = hei[0]
   mdot = group_id.attrs['mdot']
   mdot = mdot[0]
   na = group_id.attrs['na']
   na = na[0]
   vmin = group_id.attrs['vmin']
   vmin = vmin[0]
   vmax = group_id.attrs['vmax']
   vmax = vmax[0]
   beta = group_id.attrs['beta']
   beta = beta[0]
   vrot = group_id.attrs['vrot']
   vrot = vrot[0]


   group_id = file_id['bcondition']
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]
   ntheta_gdark = group_id.attrs['ntheta_gdark']
   ntheta_gdark = ntheta_gdark[0]      
   theta_gdark = np.array(group_id['theta_gdark'])
   teff_gdark = np.array(group_id['teff_gdark'])
   xic1_gdark = np.array(group_id['xic1_gdark'])
   
   group_id = file_id['coordinates']
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])
   
   group_id = file_id['angles']
   n_x = np.array(group_id['n_x'])
   n_y = np.array(group_id['n_y'])
   n_z = np.array(group_id['n_z'])
   
   group_id = file_id['frequencies']
   nodes_xobs = np.array(group_id['nodes_xobs'])   
    
   group_id = file_id['convergence_behaviour']
   itmaxc = group_id.attrs['itmaxc']
   itmaxc = itmaxc[0]
   devmaxc = group_id.attrs['devmaxc']
   devmaxc = devmaxc[0]
   epsmaxc_arr = np.array(group_id['epsmaxc_arr'])
   itmaxl = group_id.attrs['itmaxl']
   itmaxl = itmaxl[0]
   devmaxl = group_id.attrs['devmaxl']
   devmaxl = devmaxl[0]
   epsmaxl_arr = np.array(group_id['epsmaxl_arr'])   
   
   group_id = file_id['solution3d']
   scont3d = np.array(group_id['scont3d'])
   mint3d = np.array(group_id['mint3d'])
   fcontx3d = np.array(group_id['fcontx3d'])
   fconty3d = np.array(group_id['fconty3d'])
   fcontz3d = np.array(group_id['fcontz3d'])
   sline3d = np.array(group_id['sline3d'])
   mintbar3d = np.array(group_id['mintbar3d'])
   ssobo3d = np.array(group_id['ssobo3d'])   

   group_id = file_id['model3d']
   mask3d = np.array(group_id['mask3d'])
   maskb3d = np.array(group_id['mask_bpoint3d'])
   opac3d = np.array(group_id['opac3d'])
   opalbar3d = np.array(group_id['opalbar3d'])   
   velx3d = np.array(group_id['velx3d'])
   vely3d = np.array(group_id['vely3d'])
   velz3d = np.array(group_id['velz3d'])      
   t3d = np.array(group_id['t3d'])

   if read=='all':
      return  spatial_grid3d, spatial_grid1d, input_mod_dim, opt_method, \
         opt_opal, opt_angint_method, opt_sol2d, opt_incl_cont, \
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
         itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, \
         nomega, nnue, n_x, n_y, n_z, nodes_xobs
   if read=='options':
      return spatial_grid3d, spatial_grid1d, input_mod_dim, opt_method, \
         opt_opal, opt_angint_method, opt_sol2d, opt_incl_cont, \
         opt_start_cont, opt_ng_cont, opt_ait_cont, opt_incl_line, opt_start_line, \
         opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, opt_incl_gdark, \
         opt_incl_sdist
   if read=='input_parameters':
      return kline, kcont, alpha, kappa0, eps_line, eps_cont, teff, xlogg, \
         trad, xnue0, rstar, lstar, vth_fiducial, vmicro, yhe, hei, mdot, \
         na, vmin, vmax, beta, vrot
   if read=='bcondition':
      return xic1, ntheta_gdark, theta_gdark, teff_gdark, xic1_gdark   
   if read=='model':
      return nx, ny, nz, x, y, z, mask3d, maskb3d, \
         opac3d, opalbar3d, velx3d, vely3d, velz3d, t3d
   if read=='solution':
      return nx, ny, nz, x, y, z, scont3d, sline3d, ssobo3d, mint3d, mintbar3d, \
         fcontx3d, fconty3d, fcontz3d
   if read=='convergence':
      return itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
   if read=='directions':
      return nomega, nnue, n_x, n_y, n_z, nodes_xobs

###################################################################
   
def read_sc3d(fname, read='all'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['options']
   spatial_grid3d = group_id.attrs['spatial_grid3d']
   spatial_grid3d = spatial_grid3d[0]
   spatial_grid1d = group_id.attrs['spatial_grid1d']
   spatial_grid1d = spatial_grid1d[0]
   input_mod_dim = group_id.attrs['input_mod_dim']
   input_mod_dim = input_mod_dim[0]
   opt_method = group_id.attrs['opt_method']
   opt_method = opt_method[0]
   opt_opal = group_id.attrs['opt_opal']
   opt_opal = opt_opal[0]
   opt_opac = group_id.attrs['opt_opac']
   opt_opac = opt_opac[0]
   opt_angint_method = group_id.attrs['opt_angint_method']
   opt_angint_method = opt_angint_method[0]
   opt_ltec = group_id.attrs['opt_ltec']
   opt_ltec = opt_ltec[0]
   opt_sol2d = group_id.attrs['opt_sol2d']
   opt_sol2d = opt_sol2d[0]
   opt_incl_cont = group_id.attrs['opt_incl_cont']
   opt_incl_cont = opt_incl_cont[0]   
   opt_start_cont = group_id.attrs['opt_start_cont']
   opt_start_cont = opt_start_cont[0]
   opt_ng_cont = group_id.attrs['opt_ng_cont']
   opt_ng_cont = opt_ng_cont[0]
   opt_ait_cont = group_id.attrs['opt_ait_cont']
   opt_ait_cont = opt_ait_cont[0]
   opt_incl_line = group_id.attrs['opt_incl_line']
   opt_incl_line = opt_incl_line[0]   
   opt_start_line = group_id.attrs['opt_start_line']
   opt_start_line = opt_start_line[0]   
   opt_ng_line = group_id.attrs['opt_ng_line']
   opt_ng_line = opt_ng_line[0]
   opt_ait_line = group_id.attrs['opt_ait_line']
   opt_ait_line = opt_ait_line[0]   
   opt_alo_cont = group_id.attrs['opt_alo_cont']
   opt_alo_cont = opt_alo_cont[0]
   opt_alo_line = group_id.attrs['opt_alo_line']
   opt_alo_line = opt_alo_line[0]
   opt_incl_gdark = group_id.attrs['opt_incl_gdark']
   opt_incl_gdark = opt_incl_gdark[0]   
   opt_incl_sdist = group_id.attrs['opt_incl_sdist']
   opt_incl_sdist = opt_incl_sdist[0]   

   group_id = file_id['dimensions']
   nx = group_id.attrs['ndxmax']
   nx = nx[0]
   ny = group_id.attrs['ndymax']
   ny = ny[0]
   nz = group_id.attrs['ndzmax']
   nz = nz[0]
   nomega = group_id.attrs['dim_omega']
   nomega = nomega[0]
   nnue = group_id.attrs['nxobs']
   nnue = nnue[0]

   group_id = file_id['input_parameters']
   kline = group_id.attrs['kline']
   kline = kline[0]
   kcont = group_id.attrs['kcont']
   kcont = kcont[0]
   alpha = group_id.attrs['alpha']
   alpha = alpha[0]
   kappa0 = group_id.attrs['kappa0']
   kappa0 = kappa0[0]
   eps_line = group_id.attrs['eps_line']
   eps_line = eps_line[0]
   eps_cont = group_id.attrs['eps_cont']
   eps_cont = eps_cont[0]
   teff = group_id.attrs['teff']
   teff = teff[0]
   xlogg = group_id.attrs['xlogg']
   xlogg = xlogg[0]
   trad = group_id.attrs['trad']
   trad = trad[0]
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]   
   rstar = group_id.attrs['rstar']
   rstar = rstar[0]
   lstar = group_id.attrs['lstar']
   lstar = lstar[0]
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]
   vmicro = group_id.attrs['vmicro']
   vmicro = vmicro[0]
   yhe = group_id.attrs['yhe']
   yhe = yhe[0]            
   hei = group_id.attrs['hei']
   hei = hei[0]
   mdot = group_id.attrs['mdot']
   mdot = mdot[0]
   na = group_id.attrs['na']
   na = na[0]
   vmin = group_id.attrs['vmin']
   vmin = vmin[0]
   vmax = group_id.attrs['vmax']
   vmax = vmax[0]
   beta = group_id.attrs['beta']
   beta = beta[0]
   vrot = group_id.attrs['vrot']
   vrot = vrot[0]


   group_id = file_id['bcondition']
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]
   ntheta_gdark = group_id.attrs['ntheta_gdark']
   ntheta_gdark = ntheta_gdark[0]      
   theta_gdark = np.array(group_id['theta_gdark'])
   teff_gdark = np.array(group_id['teff_gdark'])
   xic1_gdark = np.array(group_id['xic1_gdark'])
   
   group_id = file_id['coordinates']
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])
   
   group_id = file_id['angles']
   n_x = np.array(group_id['n_x'])
   n_y = np.array(group_id['n_y'])
   n_z = np.array(group_id['n_z'])
   
   group_id = file_id['frequencies']
   nodes_xobs = np.array(group_id['nodes_xobs'])   
    
   group_id = file_id['convergence_behaviour']
   itmaxc = group_id.attrs['itmaxc']
   itmaxc = itmaxc[0]
   devmaxc = group_id.attrs['devmaxc']
   devmaxc = devmaxc[0]
   epsmaxc_arr = np.array(group_id['epsmaxc_arr'])
   itmaxl = group_id.attrs['itmaxl']
   itmaxl = itmaxl[0]
   devmaxl = group_id.attrs['devmaxl']
   devmaxl = devmaxl[0]
   epsmaxl_arr = np.array(group_id['epsmaxl_arr'])   
   
   group_id = file_id['solution3d']
   scont3d = np.array(group_id['scont3d'])
   mint3d = np.array(group_id['mint3d'])
   fcontx3d = np.array(group_id['fcontx3d'])
   fconty3d = np.array(group_id['fconty3d'])
   fcontz3d = np.array(group_id['fcontz3d'])
   kcontxx3d = np.array(group_id['kcontxx3d'])
   kcontxy3d = np.array(group_id['kcontxy3d'])
   kcontxz3d = np.array(group_id['kcontxz3d'])
   kcontyy3d = np.array(group_id['kcontyy3d'])
   kcontyz3d = np.array(group_id['kcontyz3d'])
   kcontzz3d = np.array(group_id['kcontzz3d'])
   sline3d = np.array(group_id['sline3d'])
   mintbar3d = np.array(group_id['mintbar3d'])
   ssobo3d = np.array(group_id['ssobo3d'])   

   group_id = file_id['model3d']
   mask3d = np.array(group_id['mask3d'])
   maskb3d = np.array(group_id['mask_bpoint3d'])
   opac3d = np.array(group_id['opac3d'])
   opalbar3d = np.array(group_id['opalbar3d'])   
   velx3d = np.array(group_id['velx3d'])
   vely3d = np.array(group_id['vely3d'])
   velz3d = np.array(group_id['velz3d'])      
   t3d = np.array(group_id['t3d'])
   eps_cont3d = np.array(group_id['eps_cont3d'])

   if read=='all':
      return  spatial_grid3d, spatial_grid1d, input_mod_dim, opt_method, \
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
         nomega, nnue, n_x, n_y, n_z, nodes_xobs
   if read=='options':
      return spatial_grid3d, spatial_grid1d, input_mod_dim, opt_method, \
         opt_opal, opt_opac, opt_angint_method, opt_ltec, opt_sol2d, opt_incl_cont, \
         opt_start_cont, opt_ng_cont, opt_ait_cont, opt_incl_line, opt_start_line, \
         opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, opt_incl_gdark, \
         opt_incl_sdist
   if read=='input_parameters':
      return kline, kcont, alpha, kappa0, eps_line, eps_cont, teff, xlogg, \
         trad, xnue0, rstar, lstar, vth_fiducial, vmicro, yhe, hei, mdot, \
         na, vmin, vmax, beta, vrot
   if read=='bcondition':
      return xic1, ntheta_gdark, theta_gdark, teff_gdark, xic1_gdark   
   if read=='model':
      return nx, ny, nz, x, y, z, mask3d, maskb3d, \
         opac3d, opalbar3d, velx3d, vely3d, velz3d, t3d
   if read=='solution':
      return nx, ny, nz, x, y, z, scont3d, sline3d, ssobo3d, mint3d, mintbar3d, \
         fcontx3d, fconty3d, fcontz3d, \
         kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d
   if read=='convergence':
      return itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
   if read=='directions':
      return nomega, nnue, n_x, n_y, n_z, nodes_xobs


###################################################################


def read_scslab3d(fname, read='all'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['options']
   opt_method = group_id.attrs['opt_method']
   opt_method = opt_method[0]
   opt_angint_method = group_id.attrs['opt_angint_method']
   opt_angint_method = opt_angint_method[0]
   opt_ng_cont = group_id.attrs['opt_ng_cont']
   opt_ng_cont = opt_ng_cont[0]
   opt_ait_cont = group_id.attrs['opt_ait_cont']
   opt_ait_cont = opt_ait_cont[0]
   opt_grey = group_id.attrs['opt_grey']
   opt_grey = opt_grey[0]
   opt_opac = group_id.attrs['opt_opac']
   opt_opac = opt_opac[0]

   group_id = file_id['dimensions']
   nx = group_id.attrs['nx']
   nx = nx[0]
   ny = group_id.attrs['ny']
   ny = ny[0]
   nz = group_id.attrs['nz']
   nz = nz[0]
   nomega = group_id.attrs['nomega']
   nomega = nomega[0]
   nnue = group_id.attrs['nnue']
   nnue = nnue[0]

   group_id = file_id['input_parameters']
   kcont = group_id.attrs['kcont']
   kcont = kcont[0]
   eps_cont = group_id.attrs['eps_cont']
   eps_cont = eps_cont[0]
   unit_length = group_id.attrs['unit_length']
   unit_length = unit_length[0]
   unit_density = group_id.attrs['unit_density']
   unit_density = unit_density[0]
   unit_velocity = group_id.attrs['unit_velocity']
   unit_velocity = unit_velocity[0]
   unit_temperature = group_id.attrs['unit_temperature']
   unit_temperature = unit_temperature[0]   

   group_id = file_id['coordinates']
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])

   group_id = file_id['angles']
   n_x = np.array(group_id['n_x'])
   n_y = np.array(group_id['n_y'])
   n_z = np.array(group_id['n_z'])

   group_id = file_id['frequencies']
   nodes_nue = np.array(group_id['nodes_nue'])

   group_id = file_id['convergence_behaviour']
   itmaxc = group_id.attrs['itmaxc']
   itmaxc = itmaxc[0]
   devmaxc = group_id.attrs['devmaxc']
   devmaxc = devmaxc[0]
   nconvc = group_id.attrs['nconvc']
   nconvc = nconvc[0]
   epsmaxc_arr = np.array(group_id['epsmaxc_arr'])
   
   group_id = file_id['solution3d']
   scont3d = np.array(group_id['scont3d'])
   mint3d = np.array(group_id['mint3d'])
   fcontx3d = np.array(group_id['fcontx3d'])
   fconty3d = np.array(group_id['fconty3d'])
   fcontz3d = np.array(group_id['fcontz3d'])
   kcontxx3d = np.array(group_id['kcontxx3d'])
   kcontxy3d = np.array(group_id['kcontxy3d'])
   kcontxz3d = np.array(group_id['kcontxz3d'])
   kcontyy3d = np.array(group_id['kcontyy3d'])
   kcontyz3d = np.array(group_id['kcontyz3d'])
   kcontzz3d = np.array(group_id['kcontzz3d'])
   
   group_id = file_id['model3d']
   mask3d = np.array(group_id['mask3d'])
   maskb3d = np.array(group_id['maskb3d'])
   rho3d = np.array(group_id['rho3d'])   
   opac3d = np.array(group_id['opac3d'])
   velx3d = np.array(group_id['velx3d'])
   vely3d = np.array(group_id['vely3d'])
   velz3d = np.array(group_id['velz3d'])      
   tgas3d = np.array(group_id['tgas3d'])
   trad3d = np.array(group_id['trad3d'])   
   eps_cont3d = np.array(group_id['eps_cont3d'])
   
   file_id.close()

   if read=='all':
      return opt_method, opt_angint_method, opt_ng_cont, opt_ait_cont, opt_grey, opt_opac, \
         unit_length, unit_density, unit_temperature, unit_velocity, \
         nx, ny, nz, nomega, nnue, kcont, eps_cont, x, y, z, n_x, n_y, n_z, nodes_nue, \
         itmaxc, devmaxc, nconvc, epsmaxc_arr, scont3d, mint3d, fcontx3d, fconty3d, fcontz3d, \
         kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d, mask3d, maskb3d, \
         rho3d, opac3d, velx3d, vely3d, velz3d, tgas3d, trad3d, eps_cont3d
   if read=='options':
      return opt_method, opt_angint_method, opt_ng_cont, opt_ait_cont, opt_grey, opt_opac
   if read=='model':
      return nx, ny, nz, kcont, eps_cont, x, y, z, mask3d, maskb3d, \
         rho3d, opac3d, velx3d, vely3d, velz3d, tgas3d, trad3d, eps_cont3d   
   if read=='solution':
      return nx, ny, nz, x, y, z, scont3d, mint3d, fcontx3d, fconty3d, fcontz3d, \
         kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d
   if read=='convergence':
      return itmaxc, devmaxc, nconvc, epsmaxc_arr
   if read=='directions':
      return nomega, nnue, n_x, n_y, n_z, nodes_nue
   if read=='units':
      return unit_length, unit_density, unit_temperature, unit_velocity

###################################################################

def read_surfb_slab3d(fname, read='all'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['input_parameters']
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]

   group_id = file_id['photospheric_profile']
   nxobs_fs = group_id.attrs['nxobs_fs']
   nxobs_fs = nxobs_fs[0]
   xobs_photprof = np.array(group_id['xobs'])
   xic_photprof = np.array(group_id['xic_nue'])
   xicc_photprof = np.array(group_id['xicc_nue'])      

   group_id = file_id['fluxem/fluxem00001']
   xobs_profile = np.array(group_id['xobs'])
   ftot_profile = np.array(group_id['ftot'])
   fcont_profile = np.array(group_id['fcont'])
   femi_profile = np.array(group_id['femi'])
   fabs_profile = np.array(group_id['fabs'])
   fnorm_profile = np.array(group_id['fnorm'])     

   group_id = file_id['surfb']
   nx_surfb = group_id.attrs['nx']
   nx_surfb = nx_surfb[0]
   ny_surfb = group_id.attrs['ny']
   ny_surfb = ny_surfb[0]
   nxobs_surfb = group_id.attrs['nxobs']
   nxobs_surfb = nxobs_surfb[0]   
   x_surfb = np.array(group_id['xp'])
   y_surfb = np.array(group_id['yp'])

   
   iem3d_surfb = np.zeros(shape=(nxobs_surfb,nx_surfb,ny_surfb))
   iemi3d_surfb = np.zeros(shape=(nxobs_surfb,nx_surfb,ny_surfb))
   iabs3d_surfb = np.zeros(shape=(nxobs_surfb,nx_surfb,ny_surfb))
   icont3d_surfb = np.zeros(shape=(nxobs_surfb,nx_surfb,ny_surfb))
   xobs_surfb = np.zeros(nxobs_surfb)

   for i in np.arange(0,nxobs_surfb):
      indx=i+1
      gname = 'surfb/xobs{indx:05d}'.format(indx=indx)
      group_id = file_id[gname]
      xobs_dum = group_id.attrs['xobs']
      xobs_dum = xobs_dum
      iem2d_surfb = np.array(group_id['iem_surface'])
      iemi2d_surfb = np.array(group_id['iemi_surface'])
      iabs2d_surfb = np.array(group_id['iabs_surface'])
      icont2d_surfb = np.array(group_id['icont_surface'])
      
      xobs_surfb[i]=xobs_dum
      iem3d_surfb[:,:][i]=iem2d_surfb
      iemi3d_surfb[:,:][i]=iemi2d_surfb
      iabs3d_surfb[:,:][i]=iabs2d_surfb
      icont3d_surfb[:,:][i]=icont2d_surfb
  
   file_id.close()

   if read=='all':
      return xic1, vth_fiducial, xnue0, nxobs_fs, \
         xobs_photprof, xic_photprof, xicc_photprof, \
         xobs_profile, ftot_profile, fcont_profile, femi_profile, fabs_profile, \
         nx_surfb, ny_surfb, nxobs_surfb, x_surfb, y_surfb, xobs_surfb, \
         iem3d_surfb, iemi3d_surfb, iabs3d_surfb, icont3d_surfb

###################################################################

def read_surfb(fname, read='all'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   np_surfb = group_id.attrs['np']
   np_surfb = np_surfb[0]
   nzeta_surfb = group_id.attrs['nzeta']
   nzeta_surfb = nzeta_surfb[0]

   group_id = file_id['input_parameters']
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]
   xobs = group_id.attrs['xobs']
   xobs = xobs[0]
   alpha = group_id.attrs['alpha']
   alpha = alpha[0]
   gamma = group_id.attrs['gamma']
   gamma = gamma[0]
   
   group_id = file_id['coordinates']
   p_surfb = np.array(group_id['p'])         
   zeta_surfb = np.array(group_id['zeta'])

   group_id = file_id['surfb']
   iem2d_surfb = np.array(group_id['iem_surface'])
   iemi2d_surfb = np.array(group_id['iemi_surface'])
   iabs2d_surfb = np.array(group_id['iabs_surface'])
   icont2d_surfb = np.array(group_id['icont_surface'])            
   
   file_id.close()

   if read=='all':
      return xic1, vth_fiducial, xnue0, xobs, alpha, gamma, \
         np_surfb, nzeta_surfb, p_surfb, zeta_surfb, \
         iem2d_surfb, iemi2d_surfb, iabs2d_surfb, icont2d_surfb

###################################################################

def read_benchmark12(fname, read='all', version=2):

   file_id = h5py.File(fname, 'r')

   group_id = file_id['input_parameters']
   kcont = group_id.attrs['kcont']
   kcont = kcont[0]
   eps_cont = group_id.attrs['eps_cont']
   eps_cont = eps_cont[0]
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]

   group_id = file_id['dimensions']
   nx = group_id.attrs['ndxmax']
   nx = nx[0]
   ny = group_id.attrs['ndymax']
   ny = ny[0]
   nz = group_id.attrs['ndzmax']
   nz = nz[0]
   nr = group_id.attrs['nr']
   nr = nr[0]

   group_id = file_id['coordinates']
   r = np.array(group_id['r'])
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])

   group_id = file_id['convergence_behaviour']
   itmaxc = group_id.attrs['itmaxc']
   itmaxc = itmaxc[0]
   epsmaxc_sc = np.array(group_id['epsmaxc_sc'])
   epsmaxc_fvm = np.array(group_id['epsmaxc_fvm'])

   group_id = file_id['model_cr']
   t1d_jo = np.array(group_id['t1d_jo'])
   opac1d_jo = np.array(group_id['opac1d_jo'])
   
   group_id = file_id['model3d']
   mask3d = np.array(group_id['mask3d'])
   t3d = np.array(group_id['t3d'])   
   opac3d = np.array(group_id['opac3d'])
  
   group_id = file_id['solution_cr']
   mint1d_joray = np.array(group_id['mint_joray'])
   mint1d_jomom = np.array(group_id['mint_jomom'])
   if version >= 1:
      fcont1d_joray = np.array(group_id['fcont_joray'])
      fcont1d_jomom = np.array(group_id['fcont_jomom'])
      
   
   group_id = file_id['solution3d']
   mint3d_sc = np.array(group_id['mint3d_sc'])
   mint3d_fvm = np.array(group_id['mint3d_fvm'])   
   if version >= 1:
      fcontr3d_sc = np.array(group_id['fcontr3d_sc'])
      fcontth3d_sc = np.array(group_id['fcontth3d_sc'])
      fcontphi3d_sc = np.array(group_id['fcontphi3d_sc'])      
      fcontr3d_fvm = np.array(group_id['fcontr3d_fvm'])
      fcontth3d_fvm = np.array(group_id['fcontth3d_fvm'])
      fcontphi3d_fvm = np.array(group_id['fcontphi3d_fvm'])
   if version >= 2:
      kcontrr3d_sc = np.array(group_id['kcontrr3d_sc'])
      kcontthth3d_sc = np.array(group_id['kcontthth3d_sc'])
      kcontphiphi3d_sc = np.array(group_id['kcontphiphi3d_sc'])
      kcontrth3d_sc = np.array(group_id['kcontrth3d_sc'])
      kcontrphi3d_sc = np.array(group_id['kcontrphi3d_sc'])
      kcontthphi3d_sc = np.array(group_id['kcontthphi3d_sc'])
      kcontrr3d_fvm = np.array(group_id['kcontrr3d_fvm'])
      kcontthth3d_fvm = np.array(group_id['kcontthth3d_fvm'])
      kcontphiphi3d_fvm = np.array(group_id['kcontphiphi3d_fvm'])
      kcontrth3d_fvm = np.array(group_id['kcontrth3d_fvm'])
      kcontrphi3d_fvm = np.array(group_id['kcontrphi3d_fvm'])
      kcontthphi3d_fvm = np.array(group_id['kcontthphi3d_fvm'])            

   
   if version==0:
      return xic1, kcont, eps_cont, nx, ny, nz, nr, x, y, z, r, \
         itmaxc, epsmaxc_sc, epsmaxc_fvm, t1d_jo, opac1d_jo, mint1d_joray, mint1d_jomom, \
         mask3d, t3d, opac3d, mint3d_sc, mint3d_fvm
   if version==1:
      return xic1, kcont, eps_cont, nx, ny, nz, nr, x, y, z, r, \
         itmaxc, epsmaxc_sc, epsmaxc_fvm, t1d_jo, opac1d_jo, mint1d_joray, mint1d_jomom, \
         fcont1d_joray, fcont1d_jomom, \
         mask3d, t3d, opac3d, mint3d_sc, mint3d_fvm, fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, \
         fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm
   if version==2:
      return xic1, kcont, eps_cont, nx, ny, nz, nr, x, y, z, r, \
         itmaxc, epsmaxc_sc, epsmaxc_fvm, t1d_jo, opac1d_jo, mint1d_joray, mint1d_jomom, \
         fcont1d_joray, fcont1d_jomom, \
         mask3d, t3d, opac3d, mint3d_sc, mint3d_fvm, fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, \
         fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm, \
         kcontrr3d_sc, kcontthth3d_sc, kcontphiphi3d_sc, kcontrth3d_sc, kcontrphi3d_sc, kcontthphi3d_sc, \
         kcontrr3d_fvm, kcontthth3d_fvm, kcontphiphi3d_fvm, kcontrth3d_fvm, kcontrphi3d_fvm, kcontthphi3d_fvm
