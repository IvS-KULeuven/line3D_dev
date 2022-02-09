function rhow_adm, rad, mu, help=print_help
;
;+
; NAME:
;       rhow_adm
;
; PURPOSE:
;       This function calculates density/density_w_star at a given radius
;       and mu for ADM-model (wind upflow component)
;
; CALLING SEQUENCE:
;       Result = rhow_adm(rad, mu)
;
; INPUTS:
;       rad:    Radius in r_star.
;       mu:     Cosine of latitude-angle.
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = rhow_adm(1.5d0, [-1.,0.,1.])
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'rhow_adm'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
rhow = sqrt(rad-1.d0+mu^2)*sqrt(1.d0+3.d0*mu^2) / $
      sqrt(rad - 1)/(4.d0*rad-3.d0+3.d0*mu^2)/rad^1.5d0

return, rhow

end
