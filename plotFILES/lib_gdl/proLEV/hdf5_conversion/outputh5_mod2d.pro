pro outputh5_mod2d, fname=fname, radius=radius, theta=theta, $
                  rho2d=rho2d, temp2d=temp2d, velr2d=velr2d, $
                  velth2d=velth2d, velphi2d=velphi2d, vth2d=vth2d, help=print_help
;
;+
; NAME:
;       outputh5_mod2d
;
; PURPOSE:
;       saves 2d atmospheric structure as h5 file that can be read in
;       by fortran routine sc3d.eo
;
; CALLING SEQUENCE:
;       output_mod2d
;
; INPUTS:
;
; KEYWORDS:
;       fname:   file-name where everything shall be stored
;       radius:  radial grid (in rstar)
;       theta:   co-latitude
;       rho2d:   density in 2d (in g/cm^3)
;       temp2d:  temperature in 2d (in K)
;       velr2d:    radial velocity in 2d (in cm/s)
;       velth2d:   polar velocity in 2d  (in cm/s)
;       velphi2d:  azimuthal velocity in 2d  (in cm/s)
;       vth2d:     thermal velocity in 2d  (in cm/s)
;
; EXAMPLE:
;       outputh5_mod2d, fname='test.h5', radius=radius, theta=theta, rho2d=rho2d
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'outputh5_mod2d'
   return
endif
;
if(not keyword_set(fname)) then fname='model2d.h5'
if(not keyword_set(radius)) then print, 'error in outputh5_mod2d: radius not specified'
if(not keyword_set(theta)) then print, 'error in outputh5_mod2d: theta not specified'
nr=n_elements(radius)
nt=n_elements(theta)
if(nr lt 2) then print, 'error in outputh5_mod2d: radius not specified'
if(nt lt 2) then print, 'error in outputh5_mod2d: theta not specified'
;
if(not keyword_set(rho2d)) then rho2d=fltarr(nr,nt)*0.d0
if(not keyword_set(temp2d)) then temp2d=fltarr(nr,nt)*0.d0
if(not keyword_set(velr2d)) then velr2d=fltarr(nr,nt)*0.d0
if(not keyword_set(velth2d)) then velth2d=fltarr(nr,nt)*0.d0
if(not keyword_set(velphi2d)) then velphi2d=fltarr(nr,nt)*0.d0
if(not keyword_set(vth2d)) then begin
   vth2d=fltarr(nr,nt)*0.d0
   for i=0, nr-1 do begin
      for j=0, nt-1 do begin
         vth2d(i,j) = v_thermal(temp2d(i,j),12,vmicro=1.d7)
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
   h5g_close, group_id
;
   group_id = h5g_create(file_id,'model')
      dtype_id = h5t_idl_create(rho2d)
         dspace_id = h5s_create_simple([nr,nt])
            dset_id = h5d_create(group_id,'rho', dtype_id, dspace_id)
               h5d_write, dset_id, rho2d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velr', dtype_id, dspace_id)
               h5d_write, dset_id, velr2d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velth', dtype_id, dspace_id)
               h5d_write, dset_id, velth2d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'velphi', dtype_id, dspace_id)
               h5d_write, dset_id, velphi2d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'temperature', dtype_id, dspace_id)
               h5d_write, dset_id, temp2d
            h5d_close, dset_id
            dset_id = h5d_create(group_id,'vth', dtype_id, dspace_id)
               h5d_write, dset_id, vth2d
            h5d_close, dset_id
         h5s_close, dspace_id
      h5t_close, dtype_id
   h5g_close, group_id
      
h5f_close, file_id
;
end
