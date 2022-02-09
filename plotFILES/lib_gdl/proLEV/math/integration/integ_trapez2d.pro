function integ_trapez2d, x, y, f, help=print_help
;
;-----------------------------------------------------------------------
;+
; NAME:
;       integ_trapez2d
;
; PURPOSE:
;       This function calculates the 2d integral of a function f given on a
;       regular x,y-grid, using the trapezoidal rule
;
; CALLING SEQUENCE:
;       Result = integ_trapez2d(x,f)
;
; INPUTS:
;       x:    x coordinates
;       y:    y coordinates
;       f:    corresponding function values (2d array)
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = integ_trapez2d(x, y, f)
;-
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'integ_trapez2d'
   return, 0
endif
;
nx=n_elements(x)
ny=n_elements(y)
if(nx lt 2) then begin
   print, 'error integ_trapez2d: at minimum 2 x-values required'
   stop
endif
if(ny lt 2) then begin
   print, 'error integ_trapez2d: at minimum 2 y-values required'
   stop
endif
if(n_elements(f) ne nx*ny) then begin
   print, 'error integ_trapez2d: number elements of f(x,y) and x, y different'
   stop
endif
;
;perfrom x-integral on each y-level
sumx=fltarr(ny)*0.d0
for j=0, ny-1 do begin
   sumx(j)=0.d0
   for i=1, nx-1 do begin
      sumx(j) = sumx(j) + 0.5d0*(f(i-1,j)+f(i,j))*(x(i)-x(i-1))
   endfor
endfor
;
;perform y-integral
sum=0.d0
for i=1, ny-1 do begin
   sum = sum + 0.5d0*(sumx(i-1)+sumx(i-1))*(y(i)-y(i-1))
endfor
;
return, sum
;
;
end
