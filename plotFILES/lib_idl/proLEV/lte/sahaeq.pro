function sahaeq, nue, nel, T, eion, gip1, gi, help=print_help
;
;+
; NAME:
;	sahaeq
;
; PURPOSE:
;	This function calculates the ionization balance in LTE from
;	the Saha equation (n_i+1/n_i)
;
; CALLING SEQUENCE:
;	Result = sahaeq(nue, nel, T)
;
; INPUTS:
;	nue:	Frequency in Hz
;	T:	Temperature in K
;       nel:    electron density in 1/cm^3
;       eion:   ionization energy in erg
;       gip1,gi:   degeneracy level of upper and lower stage
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = sahaeq(1.93798D15, 1.d10, 40.d3, 1, 1)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'sahaeq'
   return, 0
endif
;
if(n_params() lt 4) then begin
   print, 'Syntax - sahaeq(nue, nel, T, eion)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
const = 2.d0 / !cgs_clight^3
;
saha=const *nue^3/nel * gip1/gi * exp(-eion/!cgs_kb/T)
;
return, saha
;
end
