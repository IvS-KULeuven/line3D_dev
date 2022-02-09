function bnue, nue, T, help=print_help
;
;+
; NAME:
;	bnue
;
; PURPOSE:
;	This function calculates the Planck-function at a given frequency
;       and temperature in cgs.
;
;
; CALLING SEQUENCE:
;	Result = bnue(nue, T)
;
; INPUTS:
;	nue:	Frequency in Hz
;	T:	Temperature in K
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = bnue(4.57d14, 36000.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'bnue'
   return, 0
endif
;
if(n_params() lt 2) then begin
   print, 'Syntax - BNUE(nue, T)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
const1 = 2.d0 * !cgs_planck / !cgs_clight / !cgs_clight
const2 = !cgs_planck / !cgs_kb
;
bnue = 0.d0
bnue = const1 * nue^3.d0 / (exp(const2*nue/T) - 1.d0)
;
return, bnue
;
end
