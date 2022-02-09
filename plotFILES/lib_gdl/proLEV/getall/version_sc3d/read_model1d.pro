pro read_model1d, fname, nr=nr, radius=radius, rho1d=rho1d, velr1d=velr1d, t1d=t1d, vth1d=vth1d, help=print_help
;
;+
; NAME:
;	read_model1d, fname
;
; PURPOSE:
;	This procedure reads model1d.h5 file, that has been used as input for sc3d.eo
;
; CALLING SEQUENCE:
;
;	read_model1d, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from model.eo
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;          nr, radius, rho1d, velr1d, t1d, vth1d
;
; EXAMPLE:
;	read_model1d, fname, nr=nr, radius=radius, rho1d=rho1d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_model1d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_model1d: fname not specified'
   doc_library, 'read_model1d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_model1d: file does not exist'
   return
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
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
         radius=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model')
      dset_id=h5d_open(group_id, 'rho')
         rho1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr')
         velr1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'temperature')
         t1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth')
         vth1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   h5f_close, file_id
;
end
