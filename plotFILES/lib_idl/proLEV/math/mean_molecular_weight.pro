function mean_molecular_weight, y, ih, ihe, mass=mass, help=print_help
;
;+
; NAME:
;	mean_molecular_weight
;
; PURPOSE:
;	This function calculates the mean molecular weight for a
;	hydrogen-helium-plasma.
;
; CALLING SEQUENCE:
;	Result = mean_molecular_weight(z,ih,ihe)
;
; INPUTS:
;	y:	helium abundance (either by mass: y=nhe*mhe/(nh*mh+nhe*mhe)
;                                   or by number: y=nhe/nh
;	ih:	ionization degree of hydrogen
;       ihe:    ioniyation degree of helium (# of free electrons for
;       each helium atom)
;	
; KEYWORD PARAMETERS:
;       mass:   Set this keyword if input abundance shall be given by mass-fraction
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = mean_molecular_weight(0.1,1.,2.)
;	Result = mean_molecular_weight(0.3,1.,2.,/mass)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'mean_molecular_weight'
   return, 0
endif
;
if(n_params() lt 3) then begin
   print, 'Syntax - mean_molecular_weight(y,ih,ihe)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(mass)) then begin
   x=1.d0-y
   mu=1.d0/(x*(1.d0+ih)+y/4.d0*(1.d0+ihe))
   return, mu
endif else begin
   mu=(1.d0+4.d0*y)/(1.d0+ih+(1.d0+ihe)*y)
   return, mu
endelse

;
end
