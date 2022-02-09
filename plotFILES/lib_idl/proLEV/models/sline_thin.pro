function sline_thin, rad, help=print_help
;
;+
; NAME:
;       sline_thin
;
; PURPOSE:
;       This function calculates the line source function in optically
;       thin atmospheres, in units of the core-intensity I_c
;
; CALLING SEQUENCE:
;       Result = sline_thin(rad)
;
; INPUTS:
;       rad:    Radius in r_star.
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = sline_thin(1.5d0)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'sline_thin'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
sline=0.5d0*(1.d0-sqrt(1.d0-1.d0/rad/rad))

return, sline

end
