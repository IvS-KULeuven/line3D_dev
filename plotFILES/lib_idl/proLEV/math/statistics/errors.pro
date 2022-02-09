pro errors, sol3d, x, y, z, sol1d, r1d, erri, errm, err_max, devm, erri2=erri2, mask3d=mask3d, help=print_help
;
;+
; NAME:
;	errors
;
; PURPOSE:
;	This procedure calculates the deviations of a 3d-grid from a 1d-grid
;       
;
; CALLING SEQUENCE:
;	errors, sol3d, sol1d, r1d, erri, errm, devm
;
; INPUTS:
;	sol3d:	Solution of a model in 3D
;	x,y,z:  Coordinates of the solution-grid
;       sol1d:  Solution of a model in 1D
;       r1d:    Radii of the 1D-solutions (same units as x,y,z)
;
; KEYWORDS:
;
; OUTPUTS:
;      erri:    The absolute relative error at each point
;      erri2:   The relative error at each point
;      errm:    The mean (relative) error of all points
;      err_max: The maximum relative error of all points
;      devm:    The standard deviation of errors from the mean error
;      mask3d:  Mask where errors shall be calculated (0 if not)
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	errors, sline3d, x, y, z, sline1d, r1d, error3d, error_mean, deviation_mean
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'errors'
   return
endif
;
;-----------------------------------------------------------------------
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
ndzmax=n_elements(z)
nr=n_elements(r1d)
;
erri=fltarr(ndxmax,ndymax,ndzmax)*0.d0
erri2=fltarr(ndxmax,ndymax,ndzmax)*0.d0
ratioi=fltarr(ndxmax,ndymax,ndzmax)*0.d0
;
rmin=min(r1d)
rmax=x(ndxmax-2)
;
;total number of points
ntot=0L
errm=0.d0
ratiom=0.d0
;
;print, ndxmax, ndymax, ndzmax, nr
;print, sol1d
;print, r1d
;
;note: neglect ghost-zones
for i=1, ndxmax-2 do begin
   for j=1, ndymax-2 do begin
      for k=1, ndzmax-2 do begin
         rad=sqrt(x(i)^2 + y(j)^2 + z(k)^2)
;
         if(rad ge rmin and rad le rmax) then begin
            if(keyword_set(mask3d)) then begin
               if(mask3d(i,j,k) ne 0) then begin
                  xnum=sol3d(i,j,k)
                  find_indx, rad, r1d, nr, indx1, indx2
                  interpol_yp, r1d(indx1), r1d(indx2), sol1d(indx1), sol1d(indx2), rad, xtheo
                  error = abs(xtheo-xnum)/xtheo
                  erri(i,j,k) = error
                  erri2(i,j,k) = (xnum-xtheo)/xtheo
                  errm=errm+error
;
                  ratio = xnum/xtheo
                  ratioi(i,j,k) = ratio
                  ratiom = ratiom+ratio
                  ntot=ntot+1L
               endif
            endif else begin
               xnum=sol3d(i,j,k)
               find_indx, rad, r1d, nr, indx1, indx2
               interpol_yp, r1d(indx1), r1d(indx2), sol1d(indx1), sol1d(indx2), rad, xtheo
               error = abs(xtheo-xnum)/xtheo
               erri(i,j,k) = error
               erri2(i,j,k) = (xnum-xtheo)/xtheo
               errm=errm+error
;
               ratio = xnum/xtheo
               ratioi(i,j,k) = ratio
               ratiom = ratiom+ratio
               ntot=ntot+1L
            endelse
         endif
;
      endfor
   endfor
endfor
;
errm=errm/ntot
ratiom=ratiom/ntot
devm=0.d0
;
for i=1, ndxmax-2 do begin
   for j=1, ndymax-2 do begin
      for k=1, ndzmax-2 do begin
         rad=sqrt(x(i)^2 + y(j)^2 + z(k)^2)
;
         if(rad ge rmin) then begin
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
ratio_min=min(ratioi(where(ratioi gt 0.)))
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


