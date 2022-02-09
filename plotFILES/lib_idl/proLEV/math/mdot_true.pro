function mdot_true, mdot_zero, ralfven, help=print_help
;
;+
; NAME:
;	mdot_true
;
; PURPOSE:
;	This function calculates the true mass loss rate for given input
;       (spherical) mass loss rate and alfven radius
;
; CALLING SEQUENCE:
;	Result = mdot_true(mdot_zero, ralfven)
;
; INPUTS:
;	mdot_zero: (spherical) mass loss rate in arbitrary units
;	ralfven:   alfven radius in stellar radii
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = mdot_true(1., 3.5)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'mdot_true'
   return, 0
endif
;
if(n_params() lt 2) then begin
   print, 'Syntax - mdot_true(mdot_zero, ralfven)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
s3=sqrt(3.d0)
muc=sqrt(1.d0-1.d0/ralfven)
mdot=(1.d0-atan(s3)/s3-muc+atan(s3*muc)/s3)*4.d0/3.d0

;mdot_petit=1.d0-muc
;
return, mdot
;
end
