pro errors1d, x_theo, f_theo, x_num, f_num, erri, errm, err_max, devm, help=print_help, uniform=uniform
;
;+
; NAME:
;	errors1d
;
; PURPOSE:
;	This procedure calculates the deviations of two 1d arrays
;       
;
; CALLING SEQUENCE:
;	errors1d, x_theo, f_theo, x_num, f_num, erri, errm, devm
;
; INPUTS:
;	x_theo:	coordinates of theoretical grid
;	f_theo: exact function values
;	x_num:	coordinates of numerical grid
;	f_num:  solution of numerical method
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
;       uniform:  Set this keyword to calculate the errors on a uniform grid
;                 (theoretical and numerical solution obtained by linear inteprolation)
;	help:	  Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	errors1d, x_theo, f_theo, x_num, f_num, errors1d, error_mean,
;	error_max, deviation_mean
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'errors1d'
   return
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(uniform)) then begin
;
;-----------interpolate solution on a uniformly distributed grid--------
;
   nnum = n_elements(x_num)
   ntheo = n_elements(x_theo)
;
   nd=1000
;
   xmin = x_num(0)
   xmax = x_num(nnum-1)
;
   x_theo2 = xmin + findgen(nd)*(xmax-xmin)/(nd-1)
   f_num2 = fltarr(nd)
   f_theo2 = fltarr(nd)
;
   for i=0, nd-1 do begin
;inteprolate numerical solution onto new grid
      find_indx, x_theo2(i), x_num, nnum, iim1, ii
      f_num2(i) = interpol_ypfct(x_num(iim1), x_num(ii), f_num(iim1), f_num(ii), x_theo2(i))
;interpolate theoretical solution onto new grid
      find_indx, x_theo2(i), x_theo, ntheo, iim1, ii
      f_theo2(i) = interpol_ypfct(x_theo(iim1), x_theo(ii), f_theo(iim1), f_theo(ii), x_theo2(i))
   endfor
;
endif else begin
;
;---------interpolate theoretical solution onto numerical grid----------
;
   nd = n_elements(x_num)
   ntheo = n_elements(x_theo)
;
   f_num2 = f_num
   f_theo2 = fltarr(nd)
;
   for i=0, nd-1 do begin
      find_indx, x_num(i), x_theo, ntheo, iim1, ii
      f_theo2(i) = interpol_ypfct(x_theo(iim1), x_theo(ii), f_theo(iim1), f_theo(ii), x_num(i))
   endfor
;
endelse
;
;-----------------------------------------------------------------------
;
erri=fltarr(nd)*0.d0
ratioi=fltarr(nd)*0.d0
;
;
;total number of points
errm=0.d0
ratiom=0.d0
;
;note: neglect ghost-zones
for i=1, nd-1 do begin
   error = abs(f_theo2(i)-f_num2(i))/abs(f_theo2(i))
   erri(i) = error
   errm=errm+error
;
   ratio = f_num2(i)/f_theo2(i)
   ratioi(i) = ratio
   ratiom = ratiom+ratio
endfor
;
errm=errm/nd
ratiom=ratiom/nd
devm=0.d0
;
for i=1, nd-1 do begin
   devm=devm+(erri(i)-errm)^2
endfor
;
devm=sqrt(devm/float(nd-1))
;
err_max=max(erri)
ratio_max=max(ratioi)
ratio_min=min(ratioi(where(ratioi gt 0.)))
;
;
print, '--------CALCULATING ERROR FOR 1D-GRID------'
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


