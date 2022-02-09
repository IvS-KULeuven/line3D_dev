pro read_model3d, fname, nr=nr, ntheta=ntheta, nphi=nphi, radius=radius, theta=theta, phi=phi, rho3d=rho3d, $
                         velr3d=velr3d, velth3d=velth3d, velphi3d=velphi3d, t3d=t3d, $
                         vth3d=vth3d, help=print_help
;+
; NAME:
;	read_model3d, fname
;
; PURPOSE:
;	This procedure reads model3d.h5 file, that has been used as input for sc3d.eo
;
; CALLING SEQUENCE:
;
;	read_model3d, fname
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
;          nr, ntheta, nphi, radius, theta, phi, rho3d, velr3d, velth3d, velphi3d, t3d, vth3d
;
; EXAMPLE:
;	read_model3d, fname, nr=nr, radius=radius, theta=theta, phi=phi, rho3d=rho3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_model3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_model3d: fname not specified'
   doc_library, 'read_model3d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_model3d: file does not exist'
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
      att_id=h5a_open_name(group_id, 'nphi')
         nphi=h5a_read(att_id)
         nphi=nphi(0)
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
      dset_id=h5d_open(group_id, 'phi')
         phi=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model')
      dset_id=h5d_open(group_id, 'rho')
         rho3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr')
         velr3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velth')
         velth3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velphi')
         velphi3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'temperature')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth')
         vth3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
;
end
