pro read_vh1, fname=fname, tstep=tstep, radius=radius, theta=theta, $
              rho2d=rho2d, temp2d=temp2d, accr2d=accr2d, $
              accth2d=accth2d, accphi2d=accphi2d, velr2d=velr2d, $
              velth2d=velth2d, velphi2d=velphi2d, help=print_help


;+
; NAME:
;	read_vh1
;
; PURPOSE:
;	This procedure reads hdf4 files from output of the vh1-code
;
; CALLING SEQUENCE:
;
;	read_vh1
;y
; INPUTS:
;
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       fname:  Set this keyword to the file-name (without extension), that
;               needs to be read in
;
;       tstep:   time step that shall be read in
;       radius:  radius
;       theta:   co-latitude
;       rho2d:   rho2d in g/cm^3
;       temp2d:  temperature in K
;       velr2d:    radial velocity in cm/s
;       velth2d:   latitudinal velocity in cm/s
;       velphi2d:  polar velocity in cm/s
;       accr2d:    radial acceleration in cm/s^2
;       accth2d:   latitudinal acceleration in cm/s^2
;       accphi2d:  polar acceleration in cm/s^2
;
; EXAMPLE:
;	read_vh1, fname='zp0_', radius=radius, theta=theta, rho2d=rho2d, velr2d=velr2d
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_vh1'
   return
endif
;
if(not keyword_set(tstep)) then begin
   print, 'error in read_vh1: timestep not specified'
   doc_library, 'read_vh1'
   stop
endif
;
if(not keyword_set(fname)) then begin
   print, 'error in read_vh1: fname keyword not specified'
   doc_library, 'read_vh1'
   stop
endif
;
fname_rho=fname+'d.'+string(tstep, format='(i4)')
fname_ar=fname+'g.'+string(tstep, format='(i4)')
fname_at=fname+'h.'+string(tstep, format='(i4)')
fname_ap=fname+'i.'+string(tstep, format='(i4)')
fname_temp=fname+'t.'+string(tstep, format='(i4)')
fname_velr=fname+'u.'+string(tstep, format='(i4)')
fname_velth=fname+'v.'+string(tstep, format='(i4)')
fname_velphi=fname+'w.'+string(tstep, format='(i4)')
;
if(file_test(fname_rho) ne 1) then begin
   print, 'error in read_vh1: file for rho does not exist'
   stop
endif
if(file_test(fname_ar) ne 1) then begin
   print, 'error in read_vh1: file for radial accel. does not exist'
   stop
endif
if(file_test(fname_at) ne 1) then begin
   print, 'error in read_vh1: file for latitudinal accel. does not exist'
   stop
endif
if(file_test(fname_ap) ne 1) then begin
   print, 'error in read_vh1: file for polar accel. does not exist'
   stop
endif
if(file_test(fname_velr) ne 1) then begin
   print, 'error in read_vh1: file for radial velocity does not exist'
   stop
endif
if(file_test(fname_velth) ne 1) then begin
   print, 'error in read_vh1: file for latitudinal velocity does not exist'
   stop
endif
if(file_test(fname_velphi) ne 1) then begin
   print, 'error in read_vh1: file for polar velocity does not exist'
   stop
endif
;
;--------------------------read density---------------------------------
;
fid= hdf_sd_start(fname_rho, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_rho
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_rho
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, rho2d
rho2d=10.d0^rho2d
;
ndr_rho=n_elements(radius_rho)
ndt_rho=n_elements(theta_rho)
;
;----------------------read radial line acceleration--------------------
;
fid= hdf_sd_start(fname_ar, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_ar
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_ar
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, accr2d
;
ndr_ar=n_elements(radius_ar)
ndt_ar=n_elements(theta_ar)
;
;----------------------read latitudinal line acceleration---------------
;
fid= hdf_sd_start(fname_at, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_at
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_at
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, accth2d
;
ndr_at=n_elements(radius_at)
ndt_at=n_elements(theta_at)
;
;----------------------read polar line acceleration---------------------
;
fid= hdf_sd_start(fname_ap, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_ap
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_ap
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, accphi2d
;
ndr_ap=n_elements(radius_ap)
ndt_ap=n_elements(theta_ap)
;
;----------------------read temperature---------------------------------
;
fid= hdf_sd_start(fname_temp, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_temp
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_temp
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, temp2d
;
ndr_temp=n_elements(radius_temp)
ndt_temp=n_elements(theta_temp)
;
;----------------------read radial velocity-----------------------------
;
fid= hdf_sd_start(fname_velr, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_vrad
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_vrad
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, velr2d
;
ndr_vrad=n_elements(radius_vrad)
ndt_vrad=n_elements(theta_vrad)
;
;----------------------read latitudinal velocity------------------------
;
fid= hdf_sd_start(fname_velth, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_vtheta
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_vtheta
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, velth2d
;
ndr_vtheta=n_elements(radius_vtheta)
ndt_vtheta=n_elements(theta_vtheta)
;
;----------------------read polar velocity------------------------------
;
fid= hdf_sd_start(fname_velphi, /read)
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim0')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, theta_vphi
;
dset_indx =hdf_sd_nametoindex(fid,'fakeDim1')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, radius_vphi
;
dset_indx =hdf_sd_nametoindex(fid,'Data-Set-2')
dset_id = hdf_sd_select(fid, dset_indx)
hdf_sd_getdata, dset_id, velphi2d
;
ndr_vphi=n_elements(radius_vphi)
ndt_vphi=n_elements(theta_vphi)
;
;-----------------------test if same grid has been used-----------------
;----------------------------and interpolate if not---------------------
;
linterp=0
if (ndr_rho ne ndr_ar) then linterp=1
if (ndr_rho ne ndr_at) then linterp=1
if (ndr_rho ne ndr_ap) then linterp=1
if (ndr_rho ne ndr_temp) then linterp=1
if (ndr_rho ne ndr_vrad) then linterp=1
if (ndr_rho ne ndr_vtheta) then linterp=1
if (ndr_rho ne ndr_vphi) then linterp=1
if (ndt_rho ne ndt_ar) then linterp=1
if (ndt_rho ne ndt_at) then linterp=1
if (ndt_rho ne ndt_ap) then linterp=1
if (ndt_rho ne ndt_temp) then linterp=1
if (ndt_rho ne ndt_vrad) then linterp=1
if (ndt_rho ne ndt_vtheta) then linterp=1
if (ndt_rho ne ndt_vphi) then linterp=1

;test if grids are the same
if(linterp eq 0) then begin
   for i=0, ndr_rho-1 do begin
      if(abs(radius_rho(i)-radius_ar(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_at(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_ap(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_temp(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_vrad(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_vtheta(i)) gt 1.d-10) then linterp=1
      if(abs(radius_rho(i)-radius_vphi(i)) gt 1.d-10) then linterp=1
   endfor
   for i=0, ndt_rho-1 do begin
      if(abs(theta_rho(i)-theta_ar(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_at(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_ap(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_temp(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_vrad(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_vtheta(i)) gt 1.d-10) then linterp=1
      if(abs(theta_rho(i)-theta_vphi(i)) gt 1.d-10) then linterp=1
   endfor
endif

if(linterp eq 1) then begin
   print, 'error in read_vh1: grids are not the same'
   print, '-> need to interpolate on common grid (to be implemented)'
   stop
endif
;
radius=radius_rho
theta=theta_rho
;
end
