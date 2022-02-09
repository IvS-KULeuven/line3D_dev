pro read_modelspec3d_spc, fname, teff=teff, xnue0=xnue0, rstar=rstar, vth_fiducial=vth_fiducial, vmicro=vmicro, $
                          vmax=vmax, vmin=vmin, xlogg=xlogg, na=na, xic1=xic1, nr=nr, ntheta=ntheta, nphi=nphi, $
                          radius=radius, theta=theta, phi=phi, $
                          velx3d=velx3d, vely3d=vely3d, velz3d=velz3d, t3d=t3d, opac3d=opac3d, opalbar3d=opalbar3d, $
                          scont3d=scont3d, sline3d=sline3d, help=print_help
;+
; NAME:
;	read_modelspec3d_spc, fname
;
; PURPOSE:
;	This procedure reads modspec.h5 file, that has been used as input for
;	spec.eo (spherical coordinates)
;
; CALLING SEQUENCE:
;
;	read_modelspec3d_spc, fname
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
;          nr, ntheta, nphi, radius, theta, phi, velx3d, vely3d, velz3d, t3d,
;          opac3d, opalbar3d, scont3d, sline3d
;
; EXAMPLE:
;	read_modelspec3d_spc, fname, nr=nr, radius=radius, theta=theta, phi=phi, sline3d=sline3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_modelspec3d_spc'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_modelspec3d: fname not specified'
   doc_library, 'read_modelspec3d_spc'
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
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'teff')
         teff=h5a_read(att_id)
         teff=teff(0)
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
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmin')
         vmin=h5a_read(att_id)
         vmin=vmin(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'logg')
         xlogg=h5a_read(att_id)
         xlogg=xlogg(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'na')
         na=h5a_read(att_id)
         na=na(0)
      h5a_close, att_id
   h5g_close, group_id
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
   group_id = h5g_open(file_id, 'bcondition')
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model3d')
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
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'scont3d')
         scont3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline3d')
         sline3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
;
end
