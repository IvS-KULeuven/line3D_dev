pro read_modelspec3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, xcoord=xcoord, ycoord=ycoord, zcoord=zcoord, mask3d=mask3d, $
                             velx3d=velx3d, vely3d=vely3d, velz3d=velz3d, t3d=t3d, opac3d=opac3d, opalbar3d=opalbar3d, $
                             scont3d=scont3d, sline3d=sline3d, $
                             teff=teff, trad=trad, xnue0=xnue0, rstar=rstar, lstar=lstar, vth_fiducial=vth_fiducial, $
                             vmicro=vmicro, vmax=vmax, vrot=vrot, logg=logg, na=na, $
                             help=print_help
;+
; NAME:
;	read_modelspec3d, fname
;
; PURPOSE:
;	This procedure reads modspec.h5 file, that has been used as input for
;	spec.eo
;
; CALLING SEQUENCE:
;
;	read_modelspec3d, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from modelspec.eo
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;          ndxmax, ndymax, ndzmax, xcoord, ycoord, zcoord, velx3d, vely3d, velz3d, t3d,
;          opac3d, opalbar3d, scont3d, sline3d, mask3d
;
;          teff, trad, xnue0, rstar, lstar, vth_fiducial, $
;          vmicro, vmax, vrot, logg, na, $
;
; EXAMPLE:
;	read_modelspec3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax,
;	x=x, y=y, z=z, sline3d=sline3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_modelspec3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_modelspec3d: fname not specified'
   doc_library, 'read_modelspec3d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_modelspec3d: file does not exist'
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
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
   h5g_close, group_id
;
   print, 'group dimensions read in'
;
   group_id = h5g_open(file_id, 'input_parameters')
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
      att_id=h5a_open_name(group_id, 'lstar')
         lstar=h5a_read(att_id)
         lstar=lstar(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vth_fiducial')
         vth_fiducial=h5a_read(att_id)
         vth_fiducial=vth_fiducial(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmicro')
         vmicro=h5a_read(att_id)
         vmicro=vmicro(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vrot')
         vrot=h5a_read(att_id)
         vrot=vrot(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'logg')
         logg=h5a_read(att_id)
         logg=logg(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'na')
         na=h5a_read(att_id)
         na=na(0)
      h5a_close, att_id
   h5g_close, group_id

   print, 'group input_parameters read in'      

   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         xcoord=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         ycoord=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         zcoord=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   print, 'group coordinates read in'      
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
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
      dset_id=h5d_open(group_id, 't3d')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar3d')
         opalbar3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   print, 'group model3d read in'   
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'scont3d')
         scont3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline3d')
         sline3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   print, 'group solution3d read in'

;
h5f_close, file_id
;
;
;
end
