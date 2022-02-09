pro smooth_gauss, x, y, ysmooth, sigma=sigma, ngauss=ngauss
;
;smooth out a distribution y using a gaussian filter of width sigma
;over ngauss points
;note: works only for equidistant steps
;
nd = n_elements(x)
;
if(not keyword_set(sigma)) then sigma=1.
if(not keyword_set(ngauss)) then ngauss=101
;
ysmooth=y
;
for i=0, nd-1 do begin
   jstart = max([0, i-floor(ngauss/2.)])
   jend = min([nd-1, i+floor(ngauss/2.)])
;
   sum=0.d0
   norm=0.d0
   for j=jstart, jend do begin
      gauss = exp(-((x(i)-x(j))/sigma)^2)
      sum = sum + y(j)*gauss
      norm = norm + gauss
   endfor
   ysmooth(i) = sum/norm
;   print, 't1', sum*(x(1)-x(0)), norm*(x(1)-x(0))
endfor
;
;
end
