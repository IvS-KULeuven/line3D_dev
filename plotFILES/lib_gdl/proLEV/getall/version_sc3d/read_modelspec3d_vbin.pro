pro read_modelspec3d_vbin, fname, x01=x01, y01=y01, z01=z01, x02=x02, y02=y02, z02=z02, $
                           vx01=vx01, vy01=vy01, vz01=vz01, vx02=vx02, vy02=vy02, vz02=vz02, $
                           rstar1=rstar1, teff1=teff1, logg1=logg1, lstar1=lstar1, yhe1=yhe1, vrot1=vrot1, $
                           rstar2=rstar2, teff2=teff2, logg2=logg2, lstar2=lstar2, yhe2=yhe2, vrot2=vrot2, $
                           vth_fiducial=vth_fiducial, vmicro1=vmicro1, vmicro2=vmicro2, vmax=vmax, xnue0=xnue0, na=na, $
                           unit_length=unit_length, $
                           cs1_nr=cs1_nr, cs1_ntheta=cs1_ntheta, cs1_nphi=cs1_nphi, $
                           cs2_nr=cs2_nr, cs2_ntheta=cs2_ntheta, cs2_nphi=cs2_nphi, $                           
                           cs1_radius=cs1_radius, cs1_theta=cs1_theta, cs1_phi=cs1_phi, $
                           cs2_radius=cs2_radius, cs2_theta=cs2_theta, cs2_phi=cs2_phi, $
                           cs1_scont3d=cs1_scont3d, cs1_sline3d=cs1_sline3d, cs1_t3d=cs1_t3d, cs1_rho3d=cs1_rho3d, cs1_opac3d=cs1_opac3d, $
                           cs1_opalbar3d=cs1_opalbar3d, cs1_velx3d=cs1_velx3d, cs1_vely3d=cs1_vely3d, cs1_velz3d=cs1_velz3d, $
                           cs2_scont3d=cs2_scont3d, cs2_sline3d=cs2_sline3d, cs2_t3d=cs2_t3d, cs2_rho3d=cs2_rho3d, cs2_opac3d=cs2_opac3d, $
                           cs2_opalbar3d=cs2_opalbar3d, cs2_velx3d=cs2_velx3d, cs2_vely3d=cs2_vely3d, cs2_velz3d=cs2_velz3d, $
                           version=version, $
                           help=print_help
;+
; NAME:
;	read_modelspec3d_vbin, fname
;
; PURPOSE:
;	This procedure reads modspec.h5 file, that has been used as input for
;	spec_vbin.eo (spherical coordinates), binary version
;
; CALLING SEQUENCE:
;
;	read_modelspec3d_vbin, fname
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from modelspec_vbin.eo
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;          cs1_nr, cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta,
;          cs1_phi, cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d,
;          cs1_opac3d, cs1_opalbar3d, cs1_scont3d, cs1_sline3d
;          cs2_nr, cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta,
;          cs2_phi, cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d,
;          cs2_opac3d, cs2_opalbar3d, cs2_scont3d, cs2_sline3d

;
; EXAMPLE:
;	read_modelspec3d_vbin, fname, cs1_nr=cs1_nr,
;	cs1_radius=cs1_radius, cs1_theta=cs1_theta, cs1_phi=cs1_phi,
;	cs1_sline3d=cs1_sline3d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_modelspec3d_vbin'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_modelspec3d_vbin: fname not specified'
   doc_library, 'read_modelspec3d_vbin'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_modelspec3d_vbin: file does not exist'
   return
endif

if(not keyword_set(version)) then version='v00'
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'x01')
         x01=h5a_read(att_id)
         x01=x01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'y01')
         y01=h5a_read(att_id)
         y01=y01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'z01')
         z01=h5a_read(att_id)
         z01=z01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vx01')
         vx01=h5a_read(att_id)
         vx01=vx01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vy01')
         vy01=h5a_read(att_id)
         vy01=vy01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vz01')
         vz01=h5a_read(att_id)
         vz01=vz01(0)
      h5a_close, att_id
         
      att_id=h5a_open_name(group_id, 'x02')
         x02=h5a_read(att_id)
         x02=x02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'y02')
         y02=h5a_read(att_id)
         y02=y02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'z02')
         z02=h5a_read(att_id)
         z02=z02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vx02')
         vx02=h5a_read(att_id)
         vx02=vx02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vy02')
         vy02=h5a_read(att_id)
         vy02=vy02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vz02')
         vz02=h5a_read(att_id)
         vz02=vz02(0)
      h5a_close, att_id
         
      att_id=h5a_open_name(group_id, 'rstar1')
         rstar1=h5a_read(att_id)
         rstar1=rstar1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'rstar2')
         rstar2=h5a_read(att_id)
         rstar2=rstar2(0)
      h5a_close, att_id

      att_id=h5a_open_name(group_id, 'teff1')
         teff1=h5a_read(att_id)
         teff1=teff1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'teff2')
         teff2=h5a_read(att_id)
         teff2=teff2(0)
      h5a_close, att_id

      att_id=h5a_open_name(group_id, 'logg1')
         logg1=h5a_read(att_id)
         logg1=logg1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'logg2')
         logg2=h5a_read(att_id)
         logg2=logg2(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'lstar1')
         lstar1=h5a_read(att_id)
         lstar1=lstar1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'lstar2')
         lstar2=h5a_read(att_id)
         lstar2=lstar2(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'yhe1')
         yhe1=h5a_read(att_id)
         yhe1=yhe1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'yhe2')
         yhe2=h5a_read(att_id)
         yhe2=yhe2(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'yhe1')
         yhe1=h5a_read(att_id)
         yhe1=yhe1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'yhe2')
         yhe2=h5a_read(att_id)
         yhe2=yhe2(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vrot1')
         vrot1=h5a_read(att_id)
         vrot1=vrot1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vrot2')
         vrot2=h5a_read(att_id)
         vrot2=vrot2(0)
      h5a_close, att_id
         
      att_id=h5a_open_name(group_id, 'xnue0')
         xnue0=h5a_read(att_id)
         xnue0=xnue0(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vth_fiducial')
         vth_fiducial=h5a_read(att_id)
         vth_fiducial=vth_fiducial(0)
      h5a_close, att_id
      if(version eq 'v00') then begin   
         att_id=h5a_open_name(group_id, 'vmicro')
            vmicro=h5a_read(att_id)
            vmicro=vmicro(0)
            vmicro1=vmicro
            vmicro2=vmicro
         h5a_close, att_id
      endif
      if(version eq 'v01') then begin   
         att_id=h5a_open_name(group_id, 'vmicro1')
            vmicro1=h5a_read(att_id)
            vmicro1=vmicro1(0)
            h5a_close, att_id
         att_id=h5a_open_name(group_id, 'vmicro2')
            vmicro2=h5a_read(att_id)
            vmicro2=vmicro2(0)
         h5a_close, att_id
         vmicro=0.d0
      endif      
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'na')
         na=h5a_read(att_id)
         na=na(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'unit_length')
         unit_length=h5a_read(att_id)
         unit_length=unit_length(0)
      h5a_close, att_id
   h5g_close, group_id
;
;----------------------------star1-------------------------------------
;   
   group1_id = h5g_open(file_id, 'star1')
   group2_id = h5g_open(group1_id, 'dimensions')   
      att_id=h5a_open_name(group2_id, 'nr')
         cs1_nr=h5a_read(att_id)
         cs1_nr=cs1_nr(0)
      h5a_close, att_id
      att_id=h5a_open_name(group2_id, 'ntheta')
         cs1_ntheta=h5a_read(att_id)
         cs1_ntheta=cs1_ntheta(0)
      h5a_close, att_id
      att_id=h5a_open_name(group2_id, 'nphi')
         cs1_nphi=h5a_read(att_id)
         cs1_nphi=cs1_nphi(0)
      h5a_close, att_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'coordinates')
      dset_id=h5d_open(group2_id, 'r')
         cs1_radius=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'theta')
         cs1_theta=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'phi')
         cs1_phi=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'model3d')
      dset_id=h5d_open(group2_id, 'velx3d')
         cs1_velx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'vely3d')
         cs1_vely3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'velz3d')
         cs1_velz3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 't3d')
         cs1_t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'rho3d')
         cs1_rho3d=h5d_read(dset_id)
      h5d_close, dset_id         
      dset_id=h5d_open(group2_id, 'opac3d')
         cs1_opac3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'opalbar3d')
         cs1_opalbar3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'solution3d')
      dset_id=h5d_open(group2_id, 'scont3d')
         cs1_scont3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'sline3d')
         cs1_sline3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
   h5g_close, group1_id 
;
;----------------------------star2-------------------------------------
;   
   group1_id = h5g_open(file_id, 'star2')
   group2_id = h5g_open(group1_id, 'dimensions')   
      att_id=h5a_open_name(group2_id, 'nr')
         cs2_nr=h5a_read(att_id)
         cs2_nr=cs2_nr(0)
      h5a_close, att_id
      att_id=h5a_open_name(group2_id, 'ntheta')
         cs2_ntheta=h5a_read(att_id)
         cs2_ntheta=cs2_ntheta(0)
      h5a_close, att_id
      att_id=h5a_open_name(group2_id, 'nphi')
         cs2_nphi=h5a_read(att_id)
         cs2_nphi=cs2_nphi(0)
      h5a_close, att_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'coordinates')
      dset_id=h5d_open(group2_id, 'r')
         cs2_radius=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'theta')
         cs2_theta=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'phi')
         cs2_phi=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'model3d')
      dset_id=h5d_open(group2_id, 'velx3d')
         cs2_velx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'vely3d')
         cs2_vely3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'velz3d')
         cs2_velz3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 't3d')
         cs2_t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'rho3d')
         cs2_rho3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'opac3d')
         cs2_opac3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'opalbar3d')
         cs2_opalbar3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
;
   group2_id = h5g_open(group1_id, 'solution3d')
      dset_id=h5d_open(group2_id, 'scont3d')
         cs2_scont3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group2_id, 'sline3d')
         cs2_sline3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group2_id
   h5g_close, group1_id      
;
h5f_close, file_id
;
;
;
end
