function trad_dilfac, xnue, r, J, help=print_help
;
;+
; NAME:
;       trad_dilfac
;
; PURPOSE:
;       This function calculates the radiation temperature with Ansatz:
;           J_nue = W*B_nue(T_rad)
;
; CALLING SEQUENCE:
;       Result = trad_dilfac(xnue, r, J)
;
; INPUTS:
;       xnue:   Considered frequency in Hz
;       r:      Radial position in the atmosphere in stellar radii
;       J:      Mean intensity at given frequency in cgs
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; OUTPUTS:
;       Radiation temperature in K
;
; EXAMPLE:
;       t = trad(4.57d14, 3., 1.d-2)
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
dilfac = 0.5d0*(1.d0-sqrt(1.d0-1.d0/r/r))
;
term1 = !cgs_planck * xnue / !cgs_kb
term2 = 2.d0 * !cgs_planck * xnue^3.d0 * dilfac / !cgs_clight / !cgs_clight / J + 1.d0
term2 = alog(term2)
;
trad = term1 / term2
;
return, trad
;
end function trad_dilfac
