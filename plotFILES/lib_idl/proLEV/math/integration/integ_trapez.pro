function integ_trapez, x, f, help=print_help
;
;-----------------------------------------------------------------------
;+
; NAME:
;       integ_trapez
;
; PURPOSE:
;       This function calculates the integral of a function f given on a
;       x-grid, using the trapezoidal rule
;
; CALLING SEQUENCE:
;       Result = integ_trapez(x,f)
;
; INPUTS:
;       x:    x coordinates
;       f:    corresponding function values
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = integ_trapez([0.,1.],[0.,1.])
;-
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'integ_trapez'
   return, 0
endif
;
nd=n_elements(x)
if(nd lt 2) then begin
   print, 'error integ_trapez: at minimum 2 values required'
endif
if(n_elements(f) ne nd) then begin
   print, 'error integ_trapez: number elements of f(x) and x different'
endif
;
sum=0.d0
for i=1, nd-1 do begin
   sum=sum+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
endfor
;
return, sum
;
;
end function integ_trapez
