pro read_model2d_cac, fname, ndxmax=ndxmax, ndymax=ndymax, xarr=xarr, yarr=yarr, rho2d=rho2d, $
                         velr2d=velr2d, velth2d=velth2d, velphi2d=velphi2d, t2d=t2d, $
                         vth2d=vth2d, help=print_help
;+
; NAME:
;	read_model2d_cac, fname
;
; PURPOSE:
;	This procedure reads model2d.h5 file in cartesian coordinates
;
; CALLING SEQUENCE:
;
;	read_model2d, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;          ndxmax, ndymax, xarr, yarr, rho2d, velr2d, velth2d, velphi2d, t2d, vth2d
;
; EXAMPLE:
;	read_model2d, fname, ndxmax=ndxmax, ndymax=ndymax, xarr=xarr, yarr=yarr, rho2d=rho2d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_model2d_cac'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_model2d_cac: fname not specified'
   doc_library, 'read_model2d_cac'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_model2d_cac: file does not exist'
   return
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
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
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         xarr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         yarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model')
      dset_id=h5d_open(group_id, 'rho')
         rho2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr')
         velr2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velth')
         velth2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velphi')
         velphi2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'temperature')
         t2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth')
         vth2d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
;
end
