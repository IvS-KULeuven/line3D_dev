pro errors3d, num3d, theo3d, mask3d, erri, errm, err_max, devm, help=print_help
;
;+
; NAME:
;	errors3d
;
; PURPOSE:
;	This procedure calculates the deviations of a 3d-grid from a
;	(theoretical) 3d-grid
;       
;
; CALLING SEQUENCE:
;	errors3d, num3d, theo3d, mask3d, erri, errm, devm
;
; INPUTS:
;	num3d:	numerical solution in 3d
;	theo3d: theoretical solution in 3d
;       mask3d: mask where values are actually calculated in 3d
;
; KEYWORDS:
;
; OUTPUTS:
;      erri:    The relative error at each point
;      errm:    The mean (relative) error of all points
;      err_max: The maximum relative error of all points
;      devm:    The standard deviation of errors from the mean error
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	errors3d, f_interp, f_theo, mask3d, error3d, error_mean, deviation_mean
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'errors3d'
   return
endif
;
;-----------------------------------------------------------------------
;
snum3d=size(num3d)
stheo3d=size(theo3d)
if(not array_equal(snum3d,stheo3d)) then begin
   print, 'error in errors3d: size of input arrays does not match'
endif
;
nx=snum3d(1)
ny=snum3d(2)
nz=snum3d(3)
;
erri=fltarr(nx,ny,nz)*0.d0
ratioi=fltarr(nx,ny,nz)*0.d0
;
;total number of points
ntot=0L
errm=0.d0
ratiom=0.d0
;
for i=0, nx-1 do begin
   for j=0, ny-1 do begin
      for k=0, nz-1 do begin
         if(mask3d(i,j,k) eq 1) then begin
            sn=num3d(i,j,k)
            st=theo3d(i,j,k)
            error = (sn-st)/st
            erri(i,j,k) = error
            error=abs(error)
            errm=errm+error
;
            ratio = sn/st
            ratioi(i,j,k) = ratio
            ratiom = ratiom+ratio
;
            ntot=ntot+1L
         endif
      endfor
   endfor
endfor
;
errm=errm/ntot
ratiom=ratiom/ntot
devm=0.d0
;
for i=0, nx-1 do begin
   for j=1, ny-1 do begin
      for k=0, nz-1 do begin
         if(mask3d(i,j,k) eq 1) then begin
            devm=devm+(erri(i,j,k)-errm)^2
         endif
      endfor
   endfor
endfor
;
devm=sqrt(devm/float(ntot-1))
;
err_max=max(erri)
ratio_max=max(ratioi)
ratio_min=min(ratioi)
;
;
print, '--------CALCULATING ERROR FOR 3D-GRID------'
print, 'mean relative error', errm
print, 'max  relative error', err_max
print, 'standard deviation estimator of mean error', devm
print, ' '
print, 'mean ratio (num/theo)', ratiom
print, 'max  ratio (num/theo)', ratio_max
print, 'min  ratio (num/theo)', ratio_min
;
;
;
end


