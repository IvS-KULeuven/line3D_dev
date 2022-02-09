pro downsize_arr2d, r, phi, arr2d, nr_new, nphi_new, r_new, phi_new, arr2d_new, arr2db=arr2db, arr2d2_new=arr2d2_new, help=print_help
;
;+
; NAME:
;	downsize_arr2d
;
; PURPOSE:
;	This procedure calculates a new 2d-array out of a given 2d-array
;       but with smaller dimensions, given by nr_new, nphi_new
;
; CALLING SEQUENCE:
;	downsize_arr2d, r, phi, arr2d, nr_new, nphi_new, r_new, phi_new, arr2d_new
;
; INPUTS:
;	arr2d:	input 2d-array
;       r:      r-coordinate of the input array
;       phi:    phi-coordinate of the input array
;       nr_new:   dimension of new r-array
;       nphi_new: dimension of new phi-array
;	
; KEYWORD PARAMETERS:
;       arr2db: Set this keyword to an additional 2d-array that shall be downsized
;       arr2d2_new: New downsized 2d array corresponding to arr2db
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	down-scaled 2d-array, arr2d_new
;       corresponding downsized coordinates, r_new, phi_new
;
; EXAMPLE:
;	downsize_arr2d, r, phi, arr2d, r_new, phi_new, r_new, phi_new, arr2d_new
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'downsize_arr2d'
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
n_fixed=2 ;fixed inner and outer boundary
;
;------------------------create new radius-grid-------------------------
;
;number of used coordinates for averaging
nav=floor(float(nr)/float(nr_new-n_fixed))
if(nav eq 0) then begin
   print, 'error in downsize_arr2d: nav is zero => need larger input grid, or smaller output grid, nr'
endif
;
nrest_in=nr
nrest_out=nr_new-n_fixed
;
k=0
for i=0, nr_new-n_fixed-1 do begin
   xav=0.d0
   for j=0, nav-1 do begin
      xav=xav+r(k)
      k=k+1
   endfor
   xav=xav/float(nav)
   r_new(i)=xav
;
;adapt nav in order that outer points are all resolved
   nrest_in=nrest_in-nav
   nrest_out=nrest_out-1
   if(nrest_out ne 0) then begin
      nav=floor(float(nrest_in)/float(nrest_out))
   endif
;
endfor
;
r_new(nr_new-2)=r(0)
r_new(nr_new-1)=r(nr-1)
;
r_new=r_new(sort(r_new))
;
;------------------------create new phi-grid----------------------------
;
;number of used coordinates for averaging
nav=floor(float(nphi)/float(nphi_new-n_fixed))
if(nav eq 0) then begin
   print, 'error in downsize_arr2d: nav is zero => need larger input grid, or smaller output grid, nphi'
endif
;
if(nphi_new eq nphi) then begin
   phi_new=phi
endif else begin
;
   nrest_in=nphi
   nrest_out=nphi_new-n_fixed
;
   k=0
   for i=0, nphi_new-n_fixed-1 do begin
      xav=0.d0
      for j=0, nav-1 do begin
         xav=xav+phi(k)
         k=k+1
      endfor
      xav=xav/float(nav)
      phi_new(i)=xav
;
;adapt nav in order that outer points are all resolved
      nrest_in=nrest_in-nav
      nrest_out=nrest_out-1
      if(nrest_out ne 0) then begin
         nav=floor(float(nrest_in)/float(nrest_out))
      endif
;
   endfor
;
   phi_new(nphi_new-2)=phi(0)
   phi_new(nphi_new-1)=phi(nphi-1)
;
   phi_new=phi_new(sort(phi_new))
;
endelse
;
;------------------------2d interpolation-------------------------------
;
arr2d_new=fltarr(nr_new,nphi_new)*0.d0
;
if(keyword_set(arr2db)) then begin
;interpolate arr2d and arr2db
   arr2d2_new=fltarr(nr_new,nphi_new)*0.d0
;
   for i=0, nr_new-1 do begin
      for j=0, nphi_new-1 do begin
         interpol_bilinear, arr2d, r, phi, [r_new(i),phi_new(j)], valp
         arr2d_new(i,j)=valp
         interpol_bilinear, arr2db, r, phi, [r_new(i),phi_new(j)], valp
         arr2d2_new(i,j)=valp
      endfor
   endfor
;
endif else begin
;interpolate arr2d
;
   for i=0, nr_new-1 do begin
      for j=0, nphi_new-1 do begin
         interpol_bilinear, arr2d, r, phi, [r_new(i),phi_new(j)], valp
         arr2d_new(i,j)=valp
      endfor
   endfor
;
endelse
;
;
;
end
