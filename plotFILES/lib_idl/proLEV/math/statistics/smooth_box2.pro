pro smooth_box2, x, y, ysmooth, sigma=sigma, ngauss=ngauss
;
;smooth out a distribution y using a box filter of width sigma
;over ngauss points
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
   xlim1 = x(i)-sigma/2.d0
   xlim2 = x(i)+sigma/2.d0
;
   sum=0.d0
   norm=0.d0
   for j=jstart+1, jend do begin
      gauss_iim1 = 1.d0
      if(x(j-1) lt xlim1) then gauss_iim1 = 0.d0
      if(x(j-1) gt xlim2) then gauss_iim1 = 0.d0

      gauss_ii = 1.d0
      if(x(j) lt xlim1) then gauss_ii = 0.d0
      if(x(j) gt xlim2) then gauss_ii = 0.d0
      
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
