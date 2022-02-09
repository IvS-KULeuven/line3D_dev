function conv_xobs_lamb, xobs, lam0=lam0, vth_fiducial=vth_fiducial, help=print_help
;
;+
; NAME:
;	conv_xobs_lamb
;
; PURPOSE:
;	This function converts from xobs-regime to wavelength-regime.
;
; CALLING SEQUENCE:
;	Result = conv_xobs_lamb(xobs)
;
; INPUTS:
;	xobs:	Frequency shift from line center in units of fiducial
;        	doppler-width
;	
; KEYWORD PARAMETERS:
;	lam0:	Set this keyword to the rest-wavelength of transition in
;	        Angstrom (default: Halpha-transition)
;	vth_fiducial: Set this keyword to the fiducial thermal velocitiy in
;                     km/s (default: 1000.)
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; OUTPUTS:
;	Wavelength in Angstrom.
;
; EXAMPLE:
;	lamb = conv_xobs_lamb(10., lam0=6562., vth_fiducial=100.)
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'conv_xobs_lamb'
   return, 0
endif
;
;------------------set default values-----------------------------------
;
if(not keyword_set(lam0)) then begin
   print, 'rest wavelength not specified, H-alpha is assumed'
   print, '   lam0=', !lam_halpha
   lam0=!lam_halpha
endif
;
if(not keyword_set(vth_fiducial)) then begin
   print, 'vth_fiducial not specified, use default'
   print, '   vth_fiducial=', 1.d3
   vth_fiducial=1.d3
endif
;
;-----------------------------------------------------------------------
;
lamb=lam0/(xobs*vth_fiducial*1.d5/!cgs_clight + 1.d0)
;
return, lamb
;
end
