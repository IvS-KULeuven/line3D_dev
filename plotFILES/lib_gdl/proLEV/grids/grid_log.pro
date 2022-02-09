pro grid_log, nx, x, xmin, xmax, help=print_help
;
;+
; NAME:
;	grid_log
;
; PURPOSE:
;	This procedure calculates a logarithmic 1d grid
;
; CALLING SEQUENCE:
;	grid_equi, nx, x, xmin, xmax
;
; INPUTS:
;	nx:	dimension
;       xmin, xmax:  boundaries of grid
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	x: x-grid
;
; EXAMPLE:
;
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'grid_log'
   return
endif
;
;-----------------------------------------------------------------------
;
if(xmax/xmin le 0.d0) then begin
   print, 'error in grid_log: xmax/xmin le 0.'
   stop
endif

del = alog10(xmax/xmin)/float(nx-1)
x=fltarr(nx)*0.d0
x(0) = xmin
for i=1, nx-1 do begin
   x(i) = x(i-1)*10.d0^del
endfor
;
;
end
