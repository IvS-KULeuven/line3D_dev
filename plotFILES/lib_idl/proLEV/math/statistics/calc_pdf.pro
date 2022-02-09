pro calc_pdf, nx, x, x_mid, pdf_x, p_x
;
;+
; NAME:
;	calc_pdf
;
; PURPOSE:
;	This procedure calculates the probability density function for an
;	input grid, where:
;             pdf_x at (x(i)+x(i-1)/2. is 1.d0/(x(i)-x(i-1))/(nx-1)
;       
;
; CALLING SEQUENCE:
;	calc_pdf, nx, x, x_mid, pdf_x, p_x
;
; INPUTS:
;	nx:	number of coordinates of grid
;	x:      coordinates of input grid
;
; KEYWORDS:
;
; OUTPUTS:
;      x_mid     mid-points of coordinates (where pdf lives)
;      pdf_x     probability density function
;      p_x       probability of finding x in each interval
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	calc_pdf, nx, x, x_mid, pdf_x, p_x
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'calc_pdf'
   return
endif
;
x_mid = fltarr(nx-1)*0.d0
pdf_x = fltarr(nx-1)*0.d0
p_x = fltarr(nx-1)*0.d0
;
for i=0, nx-2 do begin
   x_mid(i) = (x(i+1)+x(i))/2.d0
   pdf_x(i) = 1.d0/(x(i+1)-x(i))/(nx-1)
   p_x(i) = pdf_x(i)*(x(i+1)-x(i))
endfor
;
;normalize pdf (numerical reasons)
pdf_x=pdf_x/total(p_x)
;
end
