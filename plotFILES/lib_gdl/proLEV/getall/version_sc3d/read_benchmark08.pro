pro read_benchmark08, dir=dir, opt_opal=opt_opal, xic1=xic1, kline=kline, $
                      eps_line=eps_line, rstar=rstar, vmax=vmax, nr=nr, $
                      rarr=rarr, itmaxl=itmaxl, epsmaxl_sc, epsmaxl_fvm, $
                      t1d=t1d, opalbar1d=opalbar1d, velr1d=velr1d, sline1d_sc=sline1d_sc, sline1d_fvm=sline1d_fvm, $
                      ssobo1d_cr=ssobo1d_cr, sline1d_jo=sline1d_jo, ssobo1d_jo=ssobo1d_jo, help=print_help
;
;+
; NAME:
;	read_benchmark08
;
; PURPOSE:
;	This procedure reads benchmark08.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_benchmark08
;
; INPUTS:
;
;       all inputs correspond to the variables from sc3d.f90
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       dir:    Set this keyword to the directory, where benchmark08.h5 shall
;               be read
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       options:
;          opt_opal
;
;       input parameters / stellar parameters:
;          kline, alpha, kappa0, eps_line, xic1, vmax, rstar
;
;       dimensions:
;          nr
;
;       spatial/angular/frequency grids:
;          rarr
;
;       convergence behaviour:
;          itmaxl, epsmaxl_sc, epsmaxl_fvm
;
;       model and solution on central ray:
;         t1d, t1d_jo, opalbar1d, opalbar1d_jo, velr1d, sline1d_sc,
;         sline1d_fvm, ssobo1d_cr, ssobo1d_jo, sline1d_jo
;
; EXAMPLE:
;	read_benchmark08, dir='.', nr=nr, rarr=rarr, sline1d_sc=sline1d_sc
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_benchmark08'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark08.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_benchmark08: file does not exist'
   return
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'opt_opal')
         opt_opal=h5a_read(att_id)
         opt_opal=opt_opal(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kline')
         kline=h5a_read(att_id)
         kline=kline(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'alpha')
         alpha=h5a_read(att_id)
         alpha=alpha(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kappa0')
         kappa0=h5a_read(att_id)
         kappa0=kappa0(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)*1.d5
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'rstar')
         rstar=h5a_read(att_id)
         rstar=rstar(0)
         sr=rstar*!rsu
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         rarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxl_sc')
         epsmaxl_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_fvm')
         epsmaxl_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model_cr')
      dset_id=h5d_open(group_id, 't1d')
         t1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't1d_jo')
         t1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d')
         opalbar1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_jo')
         opalbar1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr1d')
         velr1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'sline_sc')
         sline1d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline_fvm')
         sline1d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
     dset_id=h5d_open(group_id, 'ssobo1d_cr')
        ssobo1d_cr=h5d_read(dset_id)
     h5d_close, dset_id
     dset_id=h5d_open(group_id, 'sline_jo')
         sline1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
     dset_id=h5d_open(group_id, 'ssobo1d_jo')
         ssobo1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
