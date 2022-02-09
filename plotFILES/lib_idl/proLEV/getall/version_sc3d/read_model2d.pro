pro read_model2d, fname, nr=nr, ntheta=ntheta, radius=radius, theta=theta, rho2d=rho2d, $
                         velr2d=velr2d, velth2d=velth2d, velphi2d=velphi2d, t2d=t2d, $
                         vth2d=vth2d, help=print_help
;+
; NAME:
;	read_model2d, fname
;
; PURPOSE:
;	This procedure reads model2d.h5 file, that has been used as input for sc3d.eo
;
; CALLING SEQUENCE:
;
;	read_model2d, fname
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
;          nr, ntheta, radius, theta, rho2d, velr2d, velth2d, velphi2d, t2d, vth2d
;
; EXAMPLE:
;	read_model2d, fname, nr=nr, radius=radius, theta=theta, rho2d=rho2d, 
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_model2d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_model2d: fname not specified'
   doc_library, 'read_model2d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_model2d: file does not exist'
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
      att_id=h5a_open_name(group_id, 'ntheta')
         ntheta=h5a_read(att_id)
         ntheta=ntheta(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         radius=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'theta')
         theta=h5d_read(dset_id)
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
