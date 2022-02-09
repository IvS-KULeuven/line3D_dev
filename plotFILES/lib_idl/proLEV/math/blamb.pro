function blamb, lambda, T, help=print_help
;
;+
; NAME:
;	blamb
;
; PURPOSE:
;	This function calculates the Planck-function at a given wavelength
;       and temperature in cgs.
;
;
; CALLING SEQUENCE:
;	Result = bnue(lambda, T)
;
; INPUTS:
;	lambda:	Wavelength in Angstrom
;	T:	Temperature in K
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = blamb(1548.d0, 36000.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'blamb'
   return, 0
endif
;
if(n_params() lt 2) then begin
   print, 'Syntax - blamb(lambda, T)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
const1 = 2.d0 * !cgs_planck * !cgs_clight^2
const2 = !cgs_planck *!cgs_clight / !cgs_kb
;
;convert lambda to cm
lamb=lambda*1.d-8
;
blamb = 0.d0
blamb = const1 / lamb^5.d0 / (exp(const2/lamb/T) - 1.d0)
;
return, blamb
;
end
