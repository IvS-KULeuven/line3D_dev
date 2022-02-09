function sahaeq, nel, T, eion, gip1, gi, help=print_help
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
;	T:	Temperature in K
;       nel:    electron density in 1/cm^3
;       eion:   ionization energy in erg
;       gip1,gi:   degeneracy level of upper and lower stage
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = sahaeq(1.d10, 40.d3, 8.7177d-11, 1, 1)
;
;ionization energy HeII: 54.4177650 eV ; 8.7177260e-11 erg
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'sahaeq'
   return, 0
endif
;
if(n_params() lt 5) then begin
   print, 'Syntax - sahaeq(nel, T, eion, gip1, gip)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;electron de broglie wavelength
lam=sqrt(!cgs_planck^2/2.d0/!pi/!cgs_me/!cgs_kb/T)
const = 2.d0 / lam^3
;
saha = 2.d0/lam^3 /nel * gip1/gi * exp(-eion/!cgs_kb/T)
;
return, saha
;
end
