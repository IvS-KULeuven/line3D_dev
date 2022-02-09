function rho_spherical, rad, velr, xmloss, rstar, help=print_help
;
;+
; NAME:
;       rho_spherical
;
; PURPOSE:
;       This function calculates density at a given radius for spherically
;       symmetric models

;
; CALLING SEQUENCE:
;       Result = rho_spherical(rad, velr, xmloss, rstar)
;
; INPUTS:
;       rad:    Radius in r_star.
;       velr:   radial velocity in cgs
;       xmloss: mass loss rate in m_sun/yr
;       rstar:  stellar radius in r_sun
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = rho_spherical(2., 100.d5, 1.d-6, 19.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'rho_spherical'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
sr = rstar*!rsu
xmloss_cgs = xmloss*!msu/365.25d0/24.d0/3600.d0
;
rho_spherical = xmloss_cgs/(4.d0*!pi*rad^2 * sr^2 * velr)

return, rho_spherical

end
