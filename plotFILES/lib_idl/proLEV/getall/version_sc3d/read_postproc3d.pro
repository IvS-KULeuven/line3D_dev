pro read_postproc3d, fname, input_mod=input_mod, opt_opal=opt_opal, xic1=xic1, kline=kline, $
                     eps_line=eps_line, teff=teff, trad=trad, xnue0=xnue0, rstar=rstar, vth_fiducial=vth_fiducial, vmicro=vmicro, $
                     yhe=yhe, hei=hei, mdot=mdot, na=na, vmin=vmin, vmax=vmax, beta=beta, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                     dim_omega=dim_omega, nxobs=nxobs, xarr=xarr, yarr=yarr, zarr=zarr, $
                     n_x=n_x, n_y=n_y, n_z=n_z, nodes_xobs=nodes_xobs, nodesx_xnue=nodes_xnue, xic_nue=xic_nue, $
                     sline3d=sline3d, mintbar3d=mintbar3d, $
                     fluxx3d=fluxx3d, fluxy3d=fluxy3d, fluxz3d=fluxz3d, fluxr3d=fluxr3d, fluxth3d=fluxth3d, fluxphi3d=fluxphi3d, $
                     gradx3d=gradx3d, grady3d=grady3d, gradz3d=gradz3d, gradr3d=gradr3d, gradth3d=gradth3d, gradphi3d=gradphi3d, $, $
                     mask3d=mask3d, t3d=t3d, opalbar3d=opalbar3d, velx3d=velx3d, vely3d=vely3d, velz3d=velz3d, help=print_help


;+
; NAME:
;	read_postproc3d, fname
;
; PURPOSE:
;	This procedure reads *.h5 file from postproc.eo output
;
; CALLING SEQUENCE:
;
;	read_postproc3d, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from postproc.f90
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       options:
;          input_mod, opt_opal
;
;       input parameters / stellar parameters:
;          xic1, kcont, kline, eps_cont, eps_line, teff, trad, xnue0, rstar
;          vth_fiducial, vmicro, yhe, hei, mdot, na, vmin, vmax, beta
;
;       dimensions:
;          ndxmax, ndymax, ndzmax, dim_omega_omega, nxobs
;
;       spatial/angular/frequency grids:
;          xarr, yarr, zarr, n_x, n_y, n_z, nodes_xobs
;
;       3d model and solution:
;         sline3d, mintbar3d, 
;         fluxx3d, fluxy3d, fluxz3d, fluxr3d, fluxth3d, fluxphi3d
;         gradx3d, grady3d, gradz3d, gradr3d, gradth3d, gradphi3d, mask3d
;         t3d, opalbar3d, velx3d, vely3d, velz3d
;
; EXAMPLE:
;	read_postproc3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
;                  xarr=xarr, yarr=yarr, zarr=zarr, sline3d=sline3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_postproc3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_postproc3d: fname not specified'
   doc_library, 'read_postproc3d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_postproc3d: file does not exist'
   return
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
;------------------------------options----------------------------------
;
   group_id = h5g_open(file_id, 'options')
      att_id=h5a_open_name(group_id, 'input_mod')
         input_mod=h5a_read(att_id)
         input_mod=input_mod(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_opal')
         opt_opal=h5a_read(att_id)
         opt_opal=opt_opal(0)
      h5a_close, att_id
   h5g_close, group_id
;
;-------------------------boundary condition----------------------------
;
   group_id = h5g_open(file_id, 'bcondition')
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id
;
;--------------------------input_parameters-----------------------------
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kline')
         kline=h5a_read(att_id)
         kline=kline(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'teff')
         teff=h5a_read(att_id)
         teff=teff(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'trad')
         trad=h5a_read(att_id)
         trad=trad(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xnue0')
         xnue0=h5a_read(att_id)
         xnue0=xnue0(0)
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
      att_id=h5a_open_name(group_id, 'dim_omega')
         dim_omega=h5a_read(att_id)
         dim_omega=dim_omega(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nxobs')
         nxobs=h5a_read(att_id)
         nxobs=nxobs(0)
      h5a_close, att_id
   h5g_close, group_id
;
;----------------------spatial coordinates------------------------------
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         xarr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         yarr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         zarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;----------------------angular grids------------------------------------
;
   group_id = h5g_open(file_id, 'angles')
      dset_id=h5d_open(group_id, 'n_x')
         n_x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'n_y')
         n_y=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'n_z')
         n_zz=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------frequencies-----------------------------------
;
   group_id = h5g_open(file_id, 'frequencies')
      dset_id=h5d_open(group_id, 'nodes_xobs')
         nodes_xobs=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'nodes_xnue')
         nodes_xnue=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xic_nue')
         xic_nue=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d solution-----------------------------------
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'sline3d')
         sline3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mintbar3d')
         mintbar3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxx3d')
         fluxx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxy3d')
         fluxy3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxz3d')
         fluxz3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxr3d')
         fluxr3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxth3d')
         fluxth3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'fluxphi3d')
         fluxphi3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'gradx3d')
         gradx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'grady3d')
         grady3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'gradz3d')
         gradz3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'gradr3d')
         gradr3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'gradth3d')
         gradth3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'gradphi3d')
         gradphi3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d model--------------------------------------
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't3d')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar3d')
         opalbar3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velx3d')
         velx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vely3d')
         vely3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velz3d')
         velz3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
