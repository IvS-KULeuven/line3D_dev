pro read_sc3d, fname, opt_method=opt_method, opt_ng_cont=opt_ng_cont, opt_ait_cont=opt_ait_cont, $
                      opt_alo_cont=opt_alo_cont, opt_aningt_method=opt_angint_method, kcont=kcont, eps_cont=eps_cont, $
                      teff=teff, trad=trad, rstar=rstar, $
                      yhe=yhe, hei=hei, mdot=mdot, vmin=vmin, vmax=vmax, beta=beta, $
                      ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                      dim_omega=dim_omega, nnue=nnue, xarr=xarr, yarr=yarr, zarr=zarr, $
                      n_x=n_x, n_y=n_y, n_z=n_z, nodes_nue=nodes_nue, itmaxc=itmaxc, nconvc=nconvc, devmaxc=devmax, $
                      epsmaxc_arr=epsmaxc_arr, scont3d=scont3d, mint3d=mint3d, $
                      mask3d=mask3d, t3d=t3d, opac3d=opac3d, $
                      velx3d=velx3d, vely3d=vely3d, velz3d=velz3d, $
                      help=print_help


;+
; NAME:
;	read_sc3d, fname
;
; PURPOSE:
;	This procedure reads *.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_sc3d, fname
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
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       options:
;          opt_method, opt_ng_cont, opt_ait_cont, opt_alo_cont
;
;       input parameters / stellar parameters:
;          kcont, eps_cont, teff, trad, rstar
;          yhe, hei, mdot, vmin, vmax, beta
;
;       dimensions:
;          ndxmax, ndymax, ndzmax, dim_omega, nnue
;
;       spatial/angular/frequency grids:
;          xarr, yarr, zarr, nodes_nue, n_x, n_y, n_z
;
;       convergence behaviour:
;          itmaxc, devmaxc, nconvc, epsmaxc_arr
;
;       3d model and solution:
;         scont3d, mint3d
;         mask3d
;         t3d, opac3d, velx3d, vely3d, velz3d
;
; EXAMPLE:
;	read_sc3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
;                  xarr=xarr, yarr=yarr, zarr=zarr, scont3d=scont3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_sc3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_sc3d: fname not specified'
   doc_library, 'read_sc3d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_sc3d: file does not exist'
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
      att_id=h5a_open_name(group_id, 'opt_method')
         opt_method=h5a_read(att_id)
         opt_method=floor(opt_method(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_ng_cont')
         opt_ng_cont=h5a_read(att_id)
         opt_ng_cont=floor(opt_ng_cont(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_ait_cont')
         opt_ait_cont=h5a_read(att_id)
         opt_ait_cont=floor(opt_ait_cont(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_alo_cont')
         opt_alo_cont=h5a_read(att_id)
         opt_alo_cont=floor(opt_alo_cont(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'opt_angint_method')
         opt_angint_method=h5a_read(att_id)
         opt_angint_method=floor(opt_angint_method(0))
      h5a_close, att_id
   h5g_close, group_id
;
;-------------------------boundary condition----------------------------
;
;--------------------------input_parameters-----------------------------
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kcont')
         kcont=h5a_read(att_id)
         kcont=kcont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_cont')
         eps_cont=h5a_read(att_id)
         eps_cont=eps_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'teff')
         teff=h5a_read(att_id)
         teff=teff(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'trad')
         trad=h5a_read(att_id)
         trad=trad(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'rstar')
         rstar=h5a_read(att_id)
         rstar=rstar(0)
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
      att_id=h5a_open_name(group_id, 'vmin')
         vmin=h5a_read(att_id)
         vmin=vmin(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)
      h5a_close, att_id
;      att_id=h5a_open_name(group_id, 'beta')
;         beta=h5a_read(att_id)
;         beta=beta(0)
;      h5a_close, att_id
   h5g_close, group_id
;
;-----------------------dimensions--------------------------------------
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=floor(ndxmax(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=floor(ndymax(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=floor(ndzmax(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'dim_omega')
         dim_omega=h5a_read(att_id)
         dim_omega=floor(dim_omega(0))
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nnue')
         nnue=h5a_read(att_id)
         nnue=floor(nnue(0))
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
         n_z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------frequencies-----------------------------------
;
   group_id = h5g_open(file_id, 'frequencies')
      dset_id=h5d_open(group_id, 'nodes_nue')
         nodes_nue=h5d_read(dset_id)
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
      att_id=h5a_open_name(group_id, 'devmaxc')
         devmaxc=h5a_read(att_id)
         devmaxc=devmaxc(0)
         h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nconvc')
         nconvc=h5a_read(att_id)
         nconvc=nconvc(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_arr')
         epsmaxc_arr=h5d_read(dset_id)
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
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
         mask3d=floor(mask3d)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't3d')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
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
