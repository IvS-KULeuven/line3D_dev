pro grid_polar, r, phi, arr2d, nr_new, nphi_new, r_new, phi_new, arr2d_new, help=print_help
;
;+
; NAME:
;	grid_polar
;
; PURPOSE:
;	This procedure calculates a new 2d-array out of a given 2d-array
;       but with smaller dimensions, given by nr_new, nphi_new
;
; CALLING SEQUENCE:
;	grid_polar, r, phi, arr2d, nr_new, nphi_new, r_new, phi_new, arr2d_new
;
; INPUTS:
;	arr2d:	input 2d-array
;       r:      r-coordinate of the input array
;       phi:    phi-coordinate of the input array
;       nr_new:   dimension of new r-array
;       nphi_new: dimension of new phi-array
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
   doc_library, 'downscale_arr2d'
   return
endif
;
;-----------------------------------------------------------------------
;
nr=n_elements(r)
nphi=n_elements(phi)
;
r_new=fltarr(nr_new)
phi_new=fltarr(nphi_new)
;
arr2d_new=fltarr(nr_new,nphi_new)
;
;------------------------calculate new polar grids----------------------
;
nr_core=floor(0.25*nr_new)
nr_ncore=nr_new-nr_core
;
;equidistant inside core
del=1.d0/float(nr_core-1)
for i=0, nr_core-1 do begin
   r_new(i)=i*del
endfor
;
;log-log in outer part
;r_new(nr_new-nr_ncore-1)=1.d0
;r_new(nr_new-nr_ncore)=1.04d0
;del=alog10(alog10(max(r))/alog10(1.04d0))/float(nr_ncore-1)
;for i=nr_new-nr_ncore+1, nr_new-1 do begin
;   r_new(i)=r_new(i-1)^(10.d0^del)
;endfor
;
;or log in outer part
r_new(nr_new-nr_ncore-1)=1.d0
del=alog10(max(r)/1.0d0)/float(nr_ncore)
for i=nr_new-nr_ncore, nr_new-1 do begin
   r_new(i)=r_new(i-1)*10.d0^del
endfor

;
;angle grid
phi_new=min(phi)+findgen(nphi_new)*(max(phi)-min(phi))/float(nphi_new-1)
;
;-----------------------------now, perform interpolation----------------
;
for i=0, nr_new-1 do begin
   for j=0, nphi_new-1 do begin
      interpol_bilinear, arr2d, r, phi, [r_new(i),phi_new(j)], valp
      arr2d_new(i,j)=valp
   endfor
endfor
;
;
end
