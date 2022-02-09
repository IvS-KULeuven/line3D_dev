pro slice_cac, alpha, gamma, offset, dim_xslice, dim_yslice, $
               rmax, rmin, x, y, z, ndxmax, ndymax, ndzmax, $
               arr3d, x_slice, y_slice, $
               arr2d
;
;calculate slice in carthesian (rectangular) coordinates
;
;on input:  alpha, gamma:           viewing angles
;           offset:                 offset of the slice from center
;           dim_xslice, dim_yslice: dimensions of the slice
;           rmax, rmin:             limits of information region
;           x, y, z:                coordinates of the original 3d-coordinate
;                                   system
;           ndxmax, ndymax, ndumax: dimensions of the 3d-coordinate system
;           arr3d:                  scalar values to be interpolated onto slice
;
;on output: x_slice, y_slice:       coordinates of the slice
;           arr2d:                  values on slice
;
xslice_min=min([min(x), min(y), min(z)])
xslice_max=max([max(x), max(y), max(z)])
;
x_slice = xslice_min + findgen(dim_xslice)*(xslice_max-xslice_min)/(dim_xslice-1)
y_slice = xslice_min + findgen(dim_xslice)*(xslice_max-xslice_min)/(dim_xslice-1)
;
;
;make equidistant in log-space
;del = (alog10(rmax)-alog10(0.01))/(dim_xslice/2-1)
;dum = fltarr(dim_xslice/2)
;dum(0) = alog10(0.01d0)
;for i=1, dim_xslice/2-1 do begin
;   dum(i) = dum(i-1)+del
;endfor

;for i=0, dim_xslice/2-1 do begin
;   dum(i) = 10.d0^dum(i)
;endfor
;
;x_slice=fltarr(dim_xslice)*0.d0
;for i=0, dim_xslice/2-1 do begin
;   x_slice(i) = -dum(dim_xslice/2-1 - i)
;   x_slice(dim_xslice - 1 - i) = dum( dim_xslice/2 - 1 - i)
;endfor

;y_slice=x_slice

;stop
;
;-----------------------------------------------------------------------
;
arr2d=fltarr(dim_xslice, dim_yslice)*0.d0
;
;----------------interpolate everything onto the slice------------------
;
;calculate transformation matrix from coordinates on slice
;to coordinates on 3d-grid
calc_transmat2, alpha, gamma, transmat
;
vec_slice=fltarr(3)*0.d0
vec_cac=fltarr(3)*0.d0
;
for i=0, dim_xslice-1 do begin
   for j=0, dim_yslice-1 do begin
      vec_slice(0) = x_slice(i)
      vec_slice(1) = y_slice(j)
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
;------------------------interpolation of arr3d-------------------------
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

;for i=0, dim_xslice-1 do begin
;   print, arr2d(i, dim_yslice/2), x_slice(i), y_slice(dim_yslice/2)
;endfor

;stop

end
