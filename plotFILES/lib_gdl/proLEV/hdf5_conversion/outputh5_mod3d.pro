pro outputh5_mod3d, fname=fname, radius=radius, theta=theta, phi=phi, $
                  rho3d=rho3d, temp3d=temp3d, velr3d=velr3d, $
                  velth3d=velth3d, velphi3d=velphi3d, vth3d=vth3d, help=print_help
;
;+
; NAME:
;       outputh5_mod3d
;
; PURPOSE:
;       saves 3d atmospheric structure as h5 file that can be read in
;       by fortran routine model.eo
;
; CALLING SEQUENCE:
;       output_mod3d
;
; INPUTS:
;
; KEYWORDS:
;       fname:   file-name where everything shall be stored
;       radius:  radial grid (in cm)
;       theta:   co-latitude
;       phi:     azimuth
;       rho3d:   density in 3d (in g/cm^3)
;       temp3d:  temperature in 3d (in K)
;       velr3d:    radial velocity in 3d (in cm/s)
;       velth3d:   polar velocity in 3d  (in cm/s)
;       velphi3d:  azimuthal velocity in 3d  (in cm/s)
;       vth3d:     thermal velocity in 3d  (in cm/s)
;
; EXAMPLE:
;       outputh5_mod3d, fname='test.h5', radius=radius, theta=theta, rho3d=rho3d
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'outputh5_mod3d'
   return
endif
;
if(not keyword_set(fname)) then fname='model3d.h5'
if(not keyword_set(radius)) then print, 'error in outputh5_mod3d: radius not specified'
if(not keyword_set(theta)) then print, 'error in outputh5_mod3d: theta not specified'
nr=n_elements(radius)
nt=n_elements(theta)
np=n_elements(phi)
if(nr lt 2) then print, 'error in outputh5_mod3d: radius not specified'
if(nt lt 2) then print, 'error in outputh5_mod3d: theta not specified'
if(np lt 2) then print, 'error in outputh5_mod3d: phi not specified'
;
if(not keyword_set(rho3d)) then rho3d=fltarr(nr,nt,np)*0.d0
if(not keyword_set(temp3d)) then temp3d=fltarr(nr,nt,np)*0.d0
if(not keyword_set(velr3d)) then velr3d=fltarr(nr,nt,np)*0.d0
if(not keyword_set(velth3d)) then velth3d=fltarr(nr,nt,np)*0.d0
if(not keyword_set(velphi3d)) then velphi3d=fltarr(nr,nt,np)*0.d0
if(not keyword_set(vth3d)) then begin
   vth3d=fltarr(nr,nt,np)*0.d0
   for i=0, nr-1 do begin
      for j=0, nt-1 do begin
         for k=0, np-1 do begin
            vth3d(i,j,k) = v_thermal(temp3d(i,j,k),12,vmicro=1.d7)
         endfor
      endfor
   endfor
endif

;
;-----------------------------------------------------------------------
;
file_id = h5f_create(fname)
;
   group_id = h5g_create(file_id,'dimensions')
      dtype_id = h5t_idl_create(nr)
         dspace_id = h5s_create_simple(1)
            dset_id = h5a_create(group_id,'nr', dtype_id, dspace_id)
               h5a_write, dset_id, nr
            h5a_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
      dtype_id = h5t_idl_create(nt)
         dspace_id = h5s_create_simple(1)
            dset_id = h5a_create(group_id,'ntheta', dtype_id, dspace_id)
               h5a_write, dset_id, nt
            h5a_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
      dtype_id = h5t_idl_create(nt)
         dspace_id = h5s_create_simple(1)
            dset_id = h5a_create(group_id,'nphi', dtype_id, dspace_id)
               h5a_write, dset_id, np
            h5a_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
   h5g_close, group_id
;
   group_id = h5g_create(file_id,'coordinates')
      dtype_id = h5t_idl_create(radius)
         dspace_id = h5s_create_simple([nr])
            dset_id = h5d_create(group_id,'r', dtype_id, dspace_id)
               h5d_write, dset_id, radius
            h5d_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
      dtype_id = h5t_idl_create(theta)
         dspace_id = h5s_create_simple([nt])
            dset_id = h5d_create(group_id,'theta', dtype_id, dspace_id)
               h5d_write, dset_id, theta
            h5d_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
      dtype_id = h5t_idl_create(phi)
         dspace_id = h5s_create_simple([np])
            dset_id = h5d_create(group_id,'phi', dtype_id, dspace_id)
               h5d_write, dset_id, phi
            h5d_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
   h5g_close, group_id
;
   group_id = h5g_create(file_id,'model')
      dtype_id = h5t_idl_create(rho3d)
         dspace_id = h5s_create_simple([nr,nt,np])
            dset_id = h5d_create(group_id,'rho', dtype_id, dspace_id)
               h5d_write, dset_id, rho3d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velr', dtype_id, dspace_id)
               h5d_write, dset_id, velr3d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velth', dtype_id, dspace_id)
               h5d_write, dset_id, velth3d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velphi', dtype_id, dspace_id)
               h5d_write, dset_id, velphi3d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'temperature', dtype_id, dspace_id)
               h5d_write, dset_id, temp3d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'vth', dtype_id, dspace_id)
               h5d_write, dset_id, vth3d
            h5d_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
   h5g_close, group_id
      
h5f_close, file_id
;
end
