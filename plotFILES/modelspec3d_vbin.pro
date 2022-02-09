pro modelspec3d_vbin, fname=fname, alpha=alpha, gamma=gamma, offset=offset, windx=windx, oname=oname, $
                      xlim=xlim, ylim=ylim, print_help=print_help
;
;+
; NAME:
;	modelspec3d_vbin
;
; PURPOSE:
;	This procedure plots contours on a slice through the 3d grid from
;	output of modelspec_vbin.eo. The slice is defined by the
;	angles alpha, gamma and an offset
;
; CALLING SEQUENCE:
;	model3dspec_vbin, fname=fname
;
; INPUTS:
;	
; KEYWORD PARAMETERS:
;       alpha:  latitude in global coordinate system
;       gamma:  azimuth in global coordinate system  
;       offset: offset in global coordinate system [unit_length0]
;       windx:  Set this keyword to an integer, defining the window-index
;       oname:  Set this keyword to a string, defining the output-name
;               (ps-file)
;       xlim:   Set this keyword to an array of length 2, defining the xrange
;       ylim:   Set this keyword to an array of length 2, defining the yrange
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	
; EXAMPLE:
;	modelspec3d_vbin, fname='modspec_model00.h5'
;       modelspec3d_vbin, fname='modspec_model00.h5', alpha=90, gamma=90
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
pdefault=!p
;
if(keyword_set(print_help)) then begin
   doc_library, 'modelspec3d_vbin'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in modelspec3d_vbin: fname not specified'
   doc_library, 'modelspec3d_vbin'
   stop
endif
;
if(not keyword_set(windx)) then windx=0
;
if(not keyword_set(alpha)) then alpha=0.d0
if(not keyword_set(gamma)) then gamma=0.d0
alpha=alpha*!pi/180
gamma=gamma*!pi/180.
;
;-----------------------------------------------------------------------
;
fname='/lhome/levin/Postdoc/line3D/outputFILES/modspec_model00.h5'
;
read_modelspec3d_vbin, fname, $
         x01=x01, y01=y01, z01=z01, x02=x02, y02=y02, z02=z02, $
         vx01=vx01, vy01=vy01, vz01=vz01, vx02=vx02, vy02=vy02, vz02=vz02, $
         rstar1=unit_length1, teff1=teff1, logg1=logg1, lstar1=lstar1, yhe1=yhe1, vrot1=vrot1, $
         rstar2=unit_length2, teff2=teff2, logg2=logg2, lstar2=lstar2, yhe2=yhe2, vrot2=vrot2, $
         vth_fiducial=vth_fiducial, vmax=vmax, xnue0=xnue0, na=na, $
         unit_length=unit_length0, $
         cs1_nr=cs1_nr, cs1_ntheta=cs1_ntheta, cs1_nphi=cs1_nphi, $
         cs2_nr=cs2_nr, cs2_ntheta=cs2_ntheta, cs2_nphi=cs2_nphi, $                           
         cs1_radius=cs1_radius, cs1_theta=cs1_theta, cs1_phi=cs1_phi, $
         cs2_radius=cs2_radius, cs2_theta=cs2_theta, cs2_phi=cs2_phi, $
         cs1_scont3d=cs1_scont3d, cs1_sline3d=cs1_sline3d, cs1_t3d=cs1_t3d, cs1_opac3d=cs1_opac3d, $
         cs1_opalbar3d=cs1_opalbar3d, cs1_velx3d=cs1_velx3d, cs1_vely3d=cs1_vely3d, cs1_velz3d=cs1_velz3d, $
         cs2_scont3d=cs2_scont3d, cs2_sline3d=cs2_sline3d, cs2_t3d=cs2_t3d, cs2_opac3d=cs2_opac3d, $
         cs2_opalbar3d=cs2_opalbar3d, cs2_velx3d=cs2_velx3d, cs2_vely3d=cs2_vely3d, cs2_velz3d=cs2_velz3d, version='v01'           
;
;------------------------------opacity----------------------------------
;
if(keyword_set(oname)) then oname='modelspec3d_vbin_opalbar.ps'
;
clim=[-14.d0,-3.d0]
logscale=1
;
titlestr1='star 1'
titlestr2='star 2'
ctitlestr=textoidl('log(<\chi>)')


if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
contour_slice3d_vbin, cs1_radius, cs1_theta, cs1_phi, cs1_opalbar3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, $
                      cs2_radius, cs2_theta, cs2_phi, cs2_opalbar3d, cs2_velx3d, cs2_vely3d, cs2_velz3d, $
                      x01, y01, z01, x02, y02, z02, vx01, vy01, vz01, vx02, vy02, vz02, $
                      unit_length0, unit_length1, unit_length2, $
                      offset=offset, alpha=alpha, gamma=gamma, $
                      logscale=logscale, xlim=xlim, ylim=ylim, clim=clim, $
                      titlestr1=titlestr1, titlestr2=titlestr2, ctitlestr=ctitlestr, isoc=isoc

;
;-----------------------------------------------------------------------
;
!p=pdefault

end
