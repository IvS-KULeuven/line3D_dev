function cak_mdot, lstar, mstar, alpha, qbar, help=print_help
;
;+
; NAME:
;       cak_mdot
;
; PURPOSE:
;       This function calculates Mdot for CAK model
;
; CALLING SEQUENCE:
;       Result = cak_mdot(lstar, mstar, alpha, qbar)
;
; INPUTS:
;       lstar:    Luminosity in cgs.
;       mstar:   Mass in cgs
;       alpha, qbar:   line-force parameters
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = cak_mdot(1.d5*!lsun_cgs, 4.*!msun_cgs, 0.6, 2000.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'cak_mdot'
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
fdum1 = (qbar*gammae/(1.d0-gammae))^((1.d0-alpha)/alpha)
;
;calculate mdot
mdot = lstar/!cgs_clight^2 * alpha/(1.d0-alpha) * fdum1
;
;
;in msun/yr
mdot = mdot/!msu * !yr

return, mdot

end
