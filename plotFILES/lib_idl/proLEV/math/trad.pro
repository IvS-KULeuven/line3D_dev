function trad, xnue, J, help=print_help
;
;+
; NAME:
;       trad
;
; PURPOSE:
;       This function calculates the radiation temperature with Ansatz:
;           J_nue = B_nue(T_rad)
;
; CALLING SEQUENCE:
;       Result = trad(xnue, J)
;
; INPUTS:
;       xnue:   Considered frequency in Hz
;       J:      Mean intensity at given frequency in cgs
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; OUTPUTS:
;       Radiation temperature in K
;
; EXAMPLE:
;       t = trad(4.57d14, 1.d-2)
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'trad'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
term1 = !cgs_planck * xnue / !cgs_kb
term2 = 2.d0 * !cgs_planck * xnue^3.d0 / !cgs_clight / !cgs_clight / J + 1.d0
term2 = alog(term2)
;
trad = term1 / term2
;
return, trad
;
end function trad
