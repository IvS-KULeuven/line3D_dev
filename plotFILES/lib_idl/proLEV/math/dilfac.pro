function dilfac, rstar, r, help=print_help
;
;+
; NAME:
;	dilfac
;
; PURPOSE:
;	This function calculates the Dilution-factor for a given radius.
;
; CALLING SEQUENCE:
;	Result = dilfac(rstar, r)
;
; INPUTS:
;	rstar:	Stellar radius in arbitrary units
;	r:	Radius, in same units as rstar
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = dilfac(1., 12.)
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'dilfac'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
if(r lt rstar) then begin
   print, 'error in dilfac: radius less than stellar radius'
   stop
endif
;
x=(rstar/r)^2
;
dilfac=0.5d0*(1.d0-sqrt(1.d0-x))
;
return, dilfac
;
end
