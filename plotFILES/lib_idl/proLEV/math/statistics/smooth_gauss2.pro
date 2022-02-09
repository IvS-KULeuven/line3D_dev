pro smooth_gauss2, x, y, ysmooth, sigma=sigma, ngauss=ngauss
;
;smooth out a distribution y using a gaussian filter of width sigma
;over ngauss points
;note: works also for non-equidistant steps
;
nd = n_elements(x)
;
if(not keyword_set(sigma)) then sigma=1.
if(not keyword_set(ngauss)) then begin
;rough estimate for number of points to be included
   ngauss=0
   sum=0.d0
   for i=1, nd-1 do begin
      sum=sum+(x(i)-x(i-1))
      ngauss=ngauss+1
      if(sum gt 5.d0*sigma) then break
   endfor
;to be on safe side, take double the points
   ngauss=2*ngauss-1
endif
;
ysmooth=y
;
for i=0, nd-1 do begin
   jstart = max([0, i-floor(ngauss/2.)])
   jend = min([nd-1, i+floor(ngauss/2.)])
;
   sum=0.d0
   norm=0.d0
   for j=jstart+1, jend do begin
      gauss_iim1 = exp(-((x(i)-x(j-1))/sigma)^2)
      gauss_ii = exp(-((x(i)-x(j))/sigma)^2)
      f_iim1 = y(j-1)
      f_ii = y(j)
      dx = x(j)-x(j-1)
      sum = sum + 0.5d0*(gauss_iim1*f_iim1 + gauss_ii*f_ii)*dx
      norm = norm + 0.5d0*(gauss_iim1 + gauss_ii)*dx
   endfor
   ysmooth(i) = sum/norm
;   print, 't2', sum, norm
endfor
;
;
end
