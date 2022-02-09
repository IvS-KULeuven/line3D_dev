function conv_lamb_xobs, lamb, lam0=lam0, vth_fiducial=vth_fiducial, help=print_help
;
;+
; NAME:
;	conv_lamb_xobs
;
; PURPOSE:
;	This function converts from wavelength-regime to xobs-regime.
;
; CALLING SEQUENCE:
;	Result = conv_lamb_xobs(lamb)
;
; INPUTS:
;	lamb:	Wavelength in Angstrom
;	
; KEYWORD PARAMETERS:
;	lam0:	Set this keyword to the rest-wavelength of transition in
;	        Angstrom (default: Halpha-transition)
;	vth_fiducial: Set this keyword to the fiducial thermal velocitiy in
;                     km/s (default: 1000.)
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; OUTPUTS:
;	Frequency shift from line center in units of fiducial doppler-width.
;
; EXAMPLE:
;	xobs = conv_lamb_xobs(6000., lam0=6562., vth_fiducial=100.)
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'conv_lamb_xobs'
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
xobs = (lam0/lamb - 1.d0)*!cgs_clight/vth_fiducial/1.d5
;
return, xobs
;
end
