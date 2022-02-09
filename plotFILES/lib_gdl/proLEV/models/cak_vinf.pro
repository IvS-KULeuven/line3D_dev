function cak_vinf, lstar, mstar, rstar, alpha, help=print_help
;
;+
; NAME:
;       cak_vinf
;
; PURPOSE:
;       This function calculates Mdot for CAK model
;
; CALLING SEQUENCE:
;       Result = cak_vinf(lstar, mstar, rstar, alpha)
;
; INPUTS:
;       lstar:    Luminosity in cgs.
;       mstar:   Mass in cgs
;       rstar:   radius in cgs
;       alpha:   line-force parameter
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = cak_vinf(1.d5*!lsun_cgs, 4.*!msun_cgs, 5.*!rsu, 0.6)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'cak_vinf'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;calculate eddington factor
;
;thomson opacity
kappae = 0.4d0

gammae = kappae*lstar/4./!pi/!cgs_clight/!cgs_grav/mstar
;
vesc = sqrt(2.d0*!cgs_grav * mstar / rstar * (1.d0-gammae))
;
vinf = vesc*sqrt(alpha/(1.d0-alpha))
;

return, vinf

end
