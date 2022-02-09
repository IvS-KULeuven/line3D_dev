pro get_all, fname, $
             ndxmax, ndymax, ndzmax, dim_mu, dim_phi, nxobs, $
;
             input_mod_dim, spatial_grid1d, spatial_grid3d, $
             mu_grid, phi_grid, opt_sol2d, opt_incl_cont, $
             opt_start_cont, opt_ng_cont, opt_ait_cont, $
             opt_method, opt_opal, opt_agint_method, $
;
             x, y, z, nodes_mu, nodes_phi_pair, phi_mask, nodes_xobs, $
;
             itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, $
;
             scont3d, mint3d, $
;
             mask_totreg3d, mask_innreg3d, mask_bpoint3d, t3d, opac3d, $
             help=print_help
;
;+
; NAME:
;	get_all, fname, ndxmax, ndymax, ndzmax, dim_mu, dim_phi, nxobs, $
;             input_mod_dim, spatial_grid1d, spatial_grid3d, mu_grid, phi_grid, $
;             opt_sol2d, opt_incl_cont, opt_start_cont, opt_ng_cont, opt_ait_cont, opt_method, opt_opal, opt_angint_method, $
;             x, y, z, nodes_mu, nodes_phi_pair, phi_mask, nodes_xobs, $
;             itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, $
;             scont3d, mint3d, $
;             mask_totreg3d, mask_innreg3d, mask_bpoint3d, t3d, opac3d
;
; PURPOSE:
;	This procedure reads *.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	get_all, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from sc3d.f90
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;	get_all, fname, ndxmax, ndymax, ndzmax, dim_mu, dim_phi, nxobs, $
;             input_mod_dim, spatial_grid1d, spatial_grid3d, mu_grid, phi_grid, $
;             opt_sol2d, opt_incl_cont, opt_start_cont, opt_ng_cont, opt_ait_cont, opt_method, opt_opal, opt_angint_method, $
;             x, y, z, nodes_mu, nodes_phi_pair, phi_mask, nodes_xobs, $
;             itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, $
;             scont3d, mint3d, $
;             mask_totreg3d, mask_innreg3d, mask_bpoint3d, t3d, opac3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'get_all'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in get_all: fname not specified'
   doc_library, 'get_all'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
;
;------------------------------options----------------------------------
;
   group_id = h5g_open(file_id, 'options')
      att_id=h5a_open_name(group_id, 'input_mod_dim')
         input_mod_dim=h5a_read(att_id)
         input_mod_dim=input_mod_dim(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'spatial_grid1d')
         spatial_grid1d=h5a_read(att_id)
         spatial_grid1d=spatial_grid1d(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'spatial_grid3d')
         spatial_grid3d=h5a_read(att_id)
         spatial_grid3d=spatial_grid3d(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_method')
         opt_method=h5a_read(att_id)
         opt_method=opt_method(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_opal')
         opt_opal=h5a_read(att_id)
         opt_opal=opt_opal(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_angint_method')
         opt_angint_method=h5a_read(att_id)
         opt_angint_method=opt_angint_method(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'mu_grid')
         mu_grid=h5a_read(att_id)
         mu_grid=mu_grid(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'phi_grid')
         phi_grid=h5a_read(att_id)
         phi_grid=phi_grid(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_sol2d')
         opt_sol2d=h5a_read(att_id)
         opt_sol2d=opt_sol2d(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_incl_cont')
         opt_incl_cont=h5a_read(att_id)
         opt_incl_cont=opt_incl_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_start_cont')
         opt_start_cont=h5a_read(att_id)
         opt_start_cont=opt_start_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_ng_cont')
         opt_ng_cont=h5a_read(att_id)
         opt_ng_cont=opt_ng_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_ait_cont')
         opt_ait_cont=h5a_read(att_id)
         opt_ait_cont=opt_ait_cont(0)
      h5a_close, att_id
   h5g_close, group_id
;
;-----------------------dimensions--------------------------------------
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'dim_mu')
         dim_mu=h5a_read(att_id)
         dim_mu=dim_mu(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'dim_phi')
         dim_phi=h5a_read(att_id)
         dim_phi=dim_phi(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nxobs')
         nxobs=h5a_read(att_id)
         nxobs=nxobs(0)
      h5a_close, att_id
   h5g_close, group_id
;
;-----------------------dimensions--------------------------------------
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kcont')
         kcont=h5a_read(att_id)
         kcont=kcont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kline')
         kline=h5a_read(att_id)
         kline=kline(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_cont')
         eps_cont=h5a_read(att_id)
         eps_cont=eps_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'teff')
         teff=h5a_read(att_id)
         teff=teff(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'rstar')
         rstar=h5a_read(att_id)
         rstar=rstar(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vth_fiducial')
         vth_fiducial=h5a_read(att_id)
         vth_fiducial=vth_fiducial(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmicro')
         vmicro=h5a_read(att_id)
         vmicro=vmicro(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'yhe')
         yhe=h5a_read(att_id)
         yhe=yhe(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'hei')
         hei=h5a_read(att_id)
         hei=hei(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'mdot')
         mdot=h5a_read(att_id)
         mdot=mdot(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'na')
         na=h5a_read(att_id)
         na=na(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmin')
         vmin=h5a_read(att_id)
         vmin=vmin(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'beta')
         beta=h5a_read(att_id)
         beta=beta(0)
      h5a_close, att_id
   h5g_close, group_id
;
;----------------------spatial coordinates------------------------------
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         y=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;----------------------angular grids------------------------------------
;
   group_id = h5g_open(file_id, 'angles')
      dset_id=h5d_open(group_id, 'nodes_mu')
         nodes_mu=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'nodes_phi_pair')
         nodes_phi_pair=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'phi_mask')
         phi_mask=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------frequencies-----------------------------------
;
   group_id = h5g_open(file_id, 'frequencies')
      dset_id=h5d_open(group_id, 'nodes_xobs')
         nodes_xobs=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;----------------------convergence behaviour----------------------------
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxc')
         itmaxc=h5a_read(att_id)
         itmaxc=itmaxc(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxc')
         devmaxc=h5a_read(att_id)
         devmaxc=devmaxc(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxl')
         devmaxl=h5a_read(att_id)
         devmaxl=devmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_arr')
         epsmaxc_arr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_arr')
         epsmaxl_arr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d solution-----------------------------------
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'scont3d')
         scont3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint3d')
         mint3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d model--------------------------------------
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask_totreg3d')
         mask_totreg3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mask_innreg3d')
         mask_innreg3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mask_bpoint3d')
         mask_bpoint3d=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 't3d')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
