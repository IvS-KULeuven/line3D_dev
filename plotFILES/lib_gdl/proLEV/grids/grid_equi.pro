pro grid_equi, nx, x, xmin, xmax, help=print_help
;
;+
; NAME:
;	grid_equi
;
; PURPOSE:
;	This procedure calculates an equidistant 1d grid
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
;	down-scaled 2d-array, arr2d_new
;       corresponding downscaled coordinates, r_new, phi_new
;
; EXAMPLE:
;	downscale_arr2d, r, phi, arr2d, r_new, phi_new, r_new, phi_new, arr2d_new
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'grid_equi'
   return
endif
;
;-----------------------------------------------------------------------
;
x=xmin+findgen(nx)*(xmax-xmin)/(nx-1)
;
;
end
