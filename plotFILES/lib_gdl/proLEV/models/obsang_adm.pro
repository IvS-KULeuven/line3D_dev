function obsang_adm, phi, inc, obl, help=print_help
;
;+
; NAME:
;       obsang_adm
;
; PURPOSE:
;       This function calculates the angle between magnetic pole axis and
;       observer for a given phase, inclination and obliquity.
;
; CALLING SEQUENCE:
;       Result = obsang_adm(phi, inc, obl)
;
; INPUTS:
;       phi:    Phase-angle in rad.
;       inc:    Inclination angle in rad.
;       obl:    Obliquity in rad.
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = obsang_adm(0., 30., 67.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'obsang_adm'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
obsang = acos(sin(obl)*cos(phi)*sin(inc) + cos(obl)*cos(inc))

return, obsang

end
