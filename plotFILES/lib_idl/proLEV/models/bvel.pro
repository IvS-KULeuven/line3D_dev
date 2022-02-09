function bvel, vmin, vinf, beta, rad, help=print_help
;
;+
; NAME:
;       bvel
;
; PURPOSE:
;       This function calculates a beta velocity law for given radii
;
; CALLING SEQUENCE:
;       Result = bvel(vmin, vinf, beta, rad)
;
; INPUTS:
;       rad:    Radius in r_star.
;       vmin:   Minimum velocity in arbitrary units
;       vmax:   Maximum velcoity in same units as vmin
;       beta:   -
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = bvel(1., 2700., 1., [1.,2.,3.,4.])
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'bvel'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
b=1.d0-(vmin/vinf)^(1.d0/beta)
vel=vinf*(1.d0-b/rad)^beta

return, vel

end
