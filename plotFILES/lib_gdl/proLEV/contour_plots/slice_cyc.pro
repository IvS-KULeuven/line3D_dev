pro slice_cyc, alpha, gamma, offset, np, nzeta, $
               rmax, rmin, x, y, z, ndxmax, ndymax, ndzmax, $
               arr3d, x_slice, y_slice, $
               arr2d
;
;calculate slice in polar coordinates coordinates
;
;on input:  alpha, gamma:           viewing angles
;           offset:                 offset of the slice from center
;           np, nzeta:              dimensions of the slice
;           rmax, rmin:             limits of information region
;           x, y, z:                coordinates of the original 3d-coordinate
;                                   system
;           p_max:                  maximum p, until which slice is calculated
;           ndxmax, ndymax, ndumax: dimensions of the 3d-coordinate system
;           arr3d:                  values to be interpolated onto slice
;
;on output: x_slice, y_slice:       x,y-coordinates of the slice (both have
;                                   dimension (dim_p, dim_zeta)
;           arr2d:                  values on slice
;
p_min=0.01d0
p_max=max([max(x), max(y), max(z)])
;
zeta_min=0.d0
zeta_max=2.d0*!pi
;
;--------------------equiduistant zeta-grid-----------------------------
;
zeta_arr = zeta_min + findgen(nzeta)*(zeta_max - zeta_min)/(nzeta-1)
;
;---------------p-grid equidistant in logspace--------------------------
;
p_arr = p_min + findgen(np)*(p_max-p_min)/(np-1)
;
del = (alog10(p_max)-alog10(p_min))/(np-2)
;print, del
p_arr(0)=0.d0
p_arr(1)=alog10(p_min)

for i=2, np-1 do begin
   p_arr(i) = p_arr(i-1) + del
   p_arr(i-1) = 10.d0^p_arr(i-1)
;   print, p_arr(i), del
endfor
p_arr(np-1)=10^p_arr(i-1)
;
;-----------------------------------------------------------------------
;
arr2d=fltarr(np, nzeta)*0.d0
x_slice=fltarr(np, nzeta)*0.d0
y_slice=fltarr(np, nzeta)*0.d0
;
;------------------interpolate arr3d onto the slice-------------------
;
;calculate transformation matrix from coordinates on slice
;to coordinates on 3d-grid
calc_transmat2, alpha, gamma, transmat
;
vec_slice=fltarr(3)*0.d0
vec_cac=fltarr(3)*0.d0
;
for i=0, np-1 do begin
   for j=0, nzeta-1 do begin
      x_slice(i,j) = p_arr(i) * cos(zeta_arr(j))
      y_slice(i,j) = p_arr(i) * sin(zeta_arr(j))
      vec_slice(0) = x_slice(i,j)
      vec_slice(1) = y_slice(i,j)
      vec_slice(2) = offset
      vec_cac=transmat#vec_slice
;
;need to make sure that values on axes are treated the same
      if(abs(vec_cac(0)) lt 1.d-5) then vec_cac(0)=0.
      if(abs(vec_cac(1)) lt 1.d-5) then vec_cac(1)=0.
      if(abs(vec_cac(2)) lt 1.d-5) then vec_cac(2)=0.

;
      radp=sqrt(vec_cac(0)*vec_cac(0) + vec_cac(1)*vec_cac(1) + vec_cac(2)*vec_cac(2))
      if(radp lt rmin or radp gt rmax) then begin
         arr2d(i,j) = 0.d0
      endif else begin
;find indices of cube-vertices for interpolation
         get_xyz_indx, vec_cac(0), vec_cac(1), vec_cac(2), x, y, z, $
                       ndxmax, ndymax, ndzmax, indx_x1, indx_x2, $
                       indx_y1, indx_y2, indx_z1, indx_z2, expol
;store all spatial values on cube-vertices for interpolation and check if
;interpolation in logspace is allowed
         get_xyz_values1, vec_cac(0), vec_cac(1), vec_cac(2),  x, y, z, $
                          ndxmax, ndymax, ndzmax, $
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                          x1, x2, y1, y2, z1, z2, $
                          rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                          llogx, llogy, llogz
;
;----------------------interpolate arr3d--------------------------------
;
;store all physical values on cube-vertices for interpolation
         get_xyz_values2, ndxmax, ndymax, ndzmax, arr3d, $
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                          vala, valb, valc, vald, vale, valf, valg, valh, llogf

;following statements can switch to pure linear interpolation
;         llogx=0
;         llogy=0
;         llogz=0
;         llogf=0
;here, can decide if values shall be interpolated by function*r^2
         lr2=1
;actual interpolation
        trilin_complete, vec_cac(0), vec_cac(1), vec_cac(2), $
                         x1, x2, y1, y2, z1, z2, $
                         vala, valb, valc, vald, vale, valf, valg, valh, $
                         rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                         expol, llogx, llogy, llogz, llogf, lr2, valp
;
         arr2d(i,j)=valp
;
      endelse
;
   endfor
endfor


end
