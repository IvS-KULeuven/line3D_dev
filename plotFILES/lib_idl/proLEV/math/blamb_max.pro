function blamb_max, T, help=print_help
;
;+
; NAME:
;	blamb_max
;
; PURPOSE:
;	This function calculates Wien's displacement law at a given and
;	temperature in Angstrom
;
;
; CALLING SEQUENCE:
;	Result = blamb_max(T)
;
; INPUTS:
;	T:	Temperature in K
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = blamb_max(36000.)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'blamb_max'
   return, 0
endif
;
if(n_params() lt 1) then begin
   print, 'Syntax - blamb_max(T)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;xmax from maximizing planck function with x=hc/lambda/k/T
xmax=4.965114232d0
;
lamb_max = !cgs_planck*!cgs_clight/!cgs_kb/xmax/T
;print, !cgs_planck*!cgs_clight/!cgs_kb/xmax
;
;convert lambda to angstrom
lamb_max=lamb_max*1.d8
;
return, lamb_max
;
end
