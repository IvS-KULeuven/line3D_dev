pro errors2d, num2d, theo2d, erri, errm, err_max, devm, help=print_help
;
;+
; NAME:
;	errors2d
;
; PURPOSE:
;	This procedure calculates the deviations of a 2d-grid from a
;	(theoretical) 2d-grid
;       
;
; CALLING SEQUENCE:
;	errors2d, num2d, theo2d, erri, errm, devm
;
; INPUTS:
;	num2d:	numerical solution in 2d
;	theo2d: theoretical solution in 2d
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
;	errors2d, f_interp, f_theo, error2d, error_mean, deviation_mean
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'errors2d'
   return
endif
;
;-----------------------------------------------------------------------
;
if(not array_equal(num2d, theo2d)) then begin
   print, 'error in errors2d: size of input arrays does not match'
endif
;
nx=size(num2d)
nx=nx(1)
ny=size(num2d)
ny=ny(1)
;
erri=fltarr(nx,ny)*0.d0
ratioi=fltarr(nx,ny)*0.d0
;
;total number of points
ntot=0L
errm=0.d0
ratiom=0.d0
;
for i=0, nx-1 do begin
   for j=0, ny-1 do begin
      sn=num2d(i,j)
      st=theo2d(i,j)
      error = abs(sn-st)/st
      erri(i,j) = error
      errm=errm+error
;
      ratio = sn/st
      ratioi(i,j) = ratio
      ratiom = ratiom+ratio
;
      ntot=ntot+1L
   endfor
endfor
;
errm=errm/ntot
ratiom=ratiom/ntot
devm=0.d0
;
for i=0, nx-1 do begin
   for j=1, ny-1 do begin
      devm=devm+(erri(i,j)-errm)^2
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
print, '--------CALCULATING ERROR FOR 2D-GRID------'
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


