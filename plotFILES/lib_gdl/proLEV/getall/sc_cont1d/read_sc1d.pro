pro read_sc1d, fname, nr=nr, eps_line=eps_line, kline=kline, kappa0=kappa0, alpha0=alpha0, opt_opal=opt_opal, $
               eps_cont=eps_cont, kcont=kcont, teff=teff, trad=trad, xic1=xic1, xnue0=xnue0, vmicro=vmicro, $
               vth_fiducial=vth_fiducial, mdot=mdot, $
               vmin=vmin, vmax=vmax, beta=beta, yhe=yhe, hei=hei, rstar=rstar, $ 
               itmaxl=itmaxl, devmaxl=devmaxl, epsmaxl_arr=epsmaxl_arr, $
               itmaxc=itmaxc, devmaxc=devmaxc, epsmaxc_arr=epsmaxc_arr, $
               r1d=r1d, t1d=t1d, velr1d=velr1d, opalbar1d=opalbar1d, sline1d=sline1d, $
               opac1d=opac1d, scont1d=scont1d, mint1d=mint1d, fcont1d=fcont1d, $
               help=print_help

;+
; NAME:
;	read_sc1d, fname
;
; PURPOSE:
;	This procedure reads *.h5 file from sc1d.eo output
;
; CALLING SEQUENCE:
;
;	read_sc1d, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from sc1d.f90
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
;          opt_opal
;
;       input parameters / stellar parameters:
;          xic1, kline, kappa0, alpha0, eps_line, teff, trad, xnue0
;          vth_fiducial, vmicro, yhe, hei, mdot, vmin, vmax, beta, na, rstar
;          eps_cont, kcont
;
;       dimensions:
;          nr
;
;       spatial/angular/frequency grids:
;          r1d
;
;       convergence behaviour:
;          itmaxl, devmaxl, epsmaxl_arr, 
;          itmaxc, devmaxc, epsmaxc_arr
;
;       3d model and solution:
;         sline1d, t1d, velr1d, opalbar1d, 
;         scont1d, mint1d, opac1d, fcont1d
;
;
; EXAMPLE:
;	read_sc1d, fname, nr=nr, r1d=r1d, sline1d=sline1d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_sc1d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_sc1d: fname not specified'
   doc_library, 'read_sc1d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_sc1d: file does not exist'
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
;      att_id=h5a_open_name(group_id, 'kcont')
;         kcont=h5a_read(att_id)
;         kcont=kcont(0)
;      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'alpha0')
         alpha0=h5a_read(att_id)
         alpha0=alpha0(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kappa0')
         kappa0=h5a_read(att_id)
         kappa0=kappa0(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
;      att_id=h5a_open_name(group_id, 'eps_cont')
;         eps_cont=h5a_read(att_id)
;         eps_cont=eps_cont(0)
;      h5a_close, att_id
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
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
   h5g_close, group_id
;
;----------------------spatial coordinates------------------------------
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r1d')
         r1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;----------------------convergence behaviour----------------------------
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxl')
         devmaxl=h5a_read(att_id)
         devmaxl=devmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxl_arr')
         epsmaxl_arr=h5d_read(dset_id)
      h5d_close, dset_id         
;      att_id=h5a_open_name(group_id, 'itmaxc')
;         itmaxc=h5a_read(att_id)
;         itmaxc=itmaxc(0)
;      h5a_close, att_id
;      att_id=h5a_open_name(group_id, 'devmaxc')
;         devmaxc=h5a_read(att_id)
;         devmaxc=devmaxl(0)
;      h5a_close, att_id
;      dset_id=h5d_open(group_id, 'epsmaxc_arr')
;         epsmaxc_arr=h5d_read(dset_id)
;      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d solution-----------------------------------
;
   group_id = h5g_open(file_id, 'solution1d')
      dset_id=h5d_open(group_id, 'sline1d')
         sline1d=h5d_read(dset_id)
      h5d_close, dset_id
;      dset_id=h5d_open(group_id, 'scont1d')
;         scont1d=h5d_read(dset_id)
;      h5d_close, dset_id
;      dset_id=h5d_open(group_id, 'mint1d')
;         mint1d=h5d_read(dset_id)
;      h5d_close, dset_id
;      dset_id=h5d_open(group_id, 'fcont1d')
;         fcont1d=h5d_read(dset_id)
;      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d model--------------------------------------
;
   group_id = h5g_open(file_id, 'model1d')
      dset_id=h5d_open(group_id, 't1d')
         t1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d')
         opalbar1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr1d')
         velr1d=h5d_read(dset_id)
      h5d_close, dset_id
;      dset_id=h5d_open(group_id, 'opac1d')
;         opac1d=h5d_read(dset_id)
;      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
;
end
