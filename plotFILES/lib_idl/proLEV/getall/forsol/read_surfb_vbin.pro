pro read_surfb_vbin, fname=fname, npoints=npoints, points_xcoord=points_xcoord, points_ycoord=points_ycoord, $
                     xic1=xic1, xic2=xic2, xobs=xobs, alpha=alpha, gamma=gamma, int2d_tot=int2d_tot, int2d_abs=int2d_abs, int2d_emi=int2d_emi, help=print_help


;+
; NAME:
;	read_surfb_vbin
;
; PURPOSE:
;	This procedure reads surface-brightnes *.h5 file from spec_vbin.eo
;
; CALLING SEQUENCE:
;
;	read_surfb_vbin
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
;       fname:  Set this keyword to the filename
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       dimensions:
;          npoints
;
;       coordinates:
;          points_xcoord, points_ycoord
;
;       parameters:
;          xic1, xic2, xobs, alpha, gamma
;
;       surface brigthness:
;          int2d, int2d_abs, int2d_emi
;
; EXAMPLE:
;	read_surfb, fname='spec_surface.h5', np=np, nzeta=nzeta, p=p, $
;                   zeta=zeta, int2d=int2d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_surfb'
   return
endif
;
if(not keyword_set(fname)) then fname='spec_surface.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_surfb: file does not exist'
   return
endif
;
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
         h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic2')
         xic2=h5a_read(att_id)
         xic2=xic2(0)
      h5a_close, att_id         
      att_id=h5a_open_name(group_id, 'alpha')
         alpha=h5a_read(att_id)
         alpha=alpha(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'gamma')
         gamma=h5a_read(att_id)
         gamma=gamma(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xobs')
         xobs=h5a_read(att_id)
         xobs=xobs(0)
      h5a_close, att_id
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'npoints')
         npoints=h5a_read(att_id)
         npoints=npoints(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'points_xcoord')
         points_xcoord=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'points_ycoord')
         points_ycoord=h5d_read(dset_id)
      h5d_close, dset_id         
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'surfaceb')
      dset_id=h5d_open(group_id, 'iem_surface')
         int2d_tot=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'iabs_surface')
         int2d_abs=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'iemi_surface')
         int2d_emi=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
