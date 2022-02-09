function v_thermal, t, na, vmicro=vmicro, help=print_help
;
;+
; NAME:
;	v_thermal
;
; PURPOSE:
;	This function calculates the thermal velocity
;       for a given temperature and atomic number in cm/s
;
;
; CALLING SEQUENCE:
;	Result = v_thermal(t, na)
;
; INPUTS:
;	t:      Temperature in K
;	na:	Atomic number
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;	help:	Set this keyword to turbulent velocity in cm/s.
;
; EXAMPLE:
;	Result = v_therma(43.d3, 12, vmicro=1.d5)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'v_thermal'
   return, 0
endif
;
if(n_params() lt 2) then begin
   print, 'Syntax - v_thermal(t, na)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
if(not keyword_set(vmicro)) then begin
   vmicro=0.d0
endif
;
vth = sqrt(2.d0*!cgs_kb * t / float(na) / !cgs_mp + vmicro^2 )
;
return, vth
;
end
