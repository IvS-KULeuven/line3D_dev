pro convert3dto3d, xin, yin, zin, arr3din, xout, yout, zout, arr3dout
;
;interpolate values on an input 3d grid onto an output 3d grid
;
nxin=n_elements(xin)
nyin=n_elements(yin)
nzin=n_elements(zin)

nxout=n_elements(xout)
nyout=n_elements(yout)
nzout=n_elements(zout)

arr3dout=fltarr(nxout,nyout,nzout)

for i=0, nxout-1 do begin
   find_indx, xout(i), xin, nxin, ii, iim1
   for j=0, nyout-1 do begin
      find_indx, yout(j), yin, nyin, jj, jjm1
      for k=0, nzout-1 do begin
         find_indx, zout(k), zin, nzin, kk, kkm1
         coeff3d_8p_lin, xin(iim1), xin(ii), yin(jjm1), yin(jj), zin(kkm1), zin(kk), xout(i), yout(j), zout(k), $
                         acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
;
         arr3dout(i,j,k) = acoeff*arr3din(iim1,jjm1,kkm1) + bcoeff*arr3din(ii,jjm1,kkm1) + $
                           ccoeff*arr3din(iim1,jj,kkm1) + dcoeff*arr3din(ii,jj,kkm1) + $
                           ecoeff*arr3din(iim1,jjm1,kk) + fcoeff*arr3din(ii,jjm1,kk) + $
                           gcoeff*arr3din(iim1,jj,kk) + hcoeff*arr3din(ii,jj,kk)
      endfor
   endfor
endfor





end
