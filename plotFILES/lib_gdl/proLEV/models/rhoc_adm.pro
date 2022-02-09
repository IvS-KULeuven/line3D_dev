function rhoc_adm, rad, delta, mu, help=print_help
;
;+
; NAME:
;       rhoc_adm
;
; PURPOSE:
;       This function calculates density/density_c_star at a given radius,
;       smoothing parameter delta and mu for ADM-model (cooled downlow component)
;
; CALLING SEQUENCE:
;       Result = rhoc_adm(rad, delta, mu)
;
; INPUTS:
;       rad:    Radius in r_star.
;       delta:  Smoothing length parameter delta in r_star.
;       mu:     Cosine of latitude-angle.
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = rhoc_adm(1.5d0, 0.3d0, [-1.,0.,1.])
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'rhoc_adm'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
rhoc = sqrt(rad-1.d0+mu^2)*sqrt(1.d0+3.d0*mu^2) / $
      sqrt(mu^2+delta^2/rad^2)/(4.d0*rad-3.d0+3.d0*mu^2)/rad^2

return, rhoc

end
