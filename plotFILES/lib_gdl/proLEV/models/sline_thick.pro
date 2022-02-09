function sline_thick, rad, beta, vmin, vmax, help=print_help
;
;+
; NAME:
;       sline_thick
;
; PURPOSE:
;       This function calculates the line source function in optically
;       thick spherical symmetric atmospheres, in units of the core-intensity
;       I_c.
;       A beta-velocity law is assumed.
;       (see Owocki&Rybicki 1985, eq. 26)
;
; CALLING SEQUENCE:
;       Result = sline_thick(rad, beta, vmin, vmax)
;
; INPUTS:
;       rad:    Radius in r_star.
;       vmin:   Minimum velocity.
;       vmax:   Maximum velocity (v_inf).
;       beta:   velocity-law exponent.
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = sline_thick([1.,1.5,2.], 10., 2000., 1.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'sline_thick'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
mu_star=sqrt(1.d0-1.d0/rad/rad)
;
;for rmin=1.d0
cb=(1.d0-(vmin/vmax)^(1.d0/beta))
;
sigma=beta/(rad/cb - 1.d0) - 1.d0
;
ca=(1.d0-mu_star^3)/(1.d0-mu_star)/3.d0
;
sline=0.5d0*(1.d0-mu_star) * (1.d0+sigma*ca)/(1.d0+sigma/3.d0)

return, sline

end
