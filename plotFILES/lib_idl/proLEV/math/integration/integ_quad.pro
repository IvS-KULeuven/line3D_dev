function integ_quad, x, f, help=print_help
;
;-----------------------------------------------------------------------
;+
; NAME:
;       integ_quad
;
; PURPOSE:
;       This function calculates the integral of a function f given on a
;       x-grid, by approximating the function as piecewise quadratic function
;
; CALLING SEQUENCE:
;       Result = integ_quad(x,f)
;
; INPUTS:
;       x:    x coordinates
;       f:    corresponding function values
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = integ_quad([0.,1.,2.],[0.,1.,4.])
;-
;-----------------------------------------------------------------------
;
print, 'integ_quad to be debugged'
stop

if(keyword_set(print_help)) then begin
   doc_library, 'integ_quad'
   return, 0
endif
;
nd=n_elements(x)
if(nd lt 3) then begin
   print, 'error integ_quad: at minimum 3 values required'
   stop
endif
if(n_elements(f) ne nd) then begin
   print, 'error integ_quad: number elements of f(x) and x different'
   stop
endif
if(nd mod 3 ne 0) then begin
   print, 'error integ_quad: number of data points need to be a multiple of 3'
   stop
endif
;
sum=0.d0
for i=1, nd-1, 3 do begin
   dxi=x(i)-x(i-1)
   dxip1=x(i+1)-x(i)
   dx=dxi+dxip1
   a=(2.d0*dxi-dxip1)*dx/6.d0/dxi
   b=dx^3/6.d0/dxi/dxip1
   c=(2.d0*dxip1-dxi)*dx/6.d0/dxip1
   sum=sum+a*f(i-1)+b*f(i)+c*f(i+1)
endfor
;
return, sum
;
;
end function integ_quad
