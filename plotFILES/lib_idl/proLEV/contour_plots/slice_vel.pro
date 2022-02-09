pro slice_vel, alpha, gamma, offset, nx, ny, rmax, rmin, $
               xslice_max, xslice_min, x, y, z, ndxmax, ndymax, ndzmax, $
               velx3d, vely3d, velz3d, x_slice, y_slice, $
               velx2d, vely2d
;
;calculate velocity vectors on a slice 
;
;on input:  alpha, gamma:           viewing angles
;           offset:                 offset of the slice from center
;           nx, ny:                 dimensions of the slice
;           rmax, rmin:             limits of information region
;           x, y, z:                coordinates of the original 3d-coordinate
;                                   system
;           xlice_min, xslice_max:  range of the slice
;           ndxmax, ndymax, ndumax: dimensions of the 3d-coordinate system
;           velx3d, vely3d, velz3d: values to be interpolated onto slice
;
;on output: x_slice, y_slice:       x,y-coordinates of the slice
;           velx2d:                 x-component of velocity along x_slice
;           vely2d:                 y-component of velocity along y_slice



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
;
;---------------equidistant spatial grid--------------------------------
;
x_slice = xslice_min + findgen(nx) * (xslice_max-xslice_min)/(nx-1)
y_slice = xslice_min + findgen(ny) * (xslice_max-xslice_min)/(ny-1)
;
;-----------------------------------------------------------------------
;
velx2d = fltarr(nx,ny) * 0.d0
vely2d = fltarr(nx,ny) * 0.d0
;
;----------------interpolate everything onto the slice------------------
;
;calculate transformation matrix: slice-coordinates -> 3d-coordinates
calc_transmat2, alpha, gamma, transmat
;
;calculate transformation matrix: 3d-coordinates -> slice_coordinates
transmat_trans=transpose(transmat)
;
vec_slice=fltarr(3)*0.d0
vec_cac=fltarr(3)*0.d0
vel_slice=fltarr(3)*0.d0
vel_cac=fltarr(3)*0.d0
;
for i=0, nx-1 do begin
   for j=0, ny-1 do begin
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
         velx2d(i,j) = 0.d0
         vely2d(i,j) = 0.d0
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
;------------------interpolation of velocity components-----------------
;
;x-components
         get_xyz_values2, ndxmax, ndymax, ndzmax, velx3d, $
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                          vala, valb, valc, vald, vale, valf, valg, valh, llogf
;         llogx=0
;         llogy=0
;         llogz=0
;         llogf=0
;here, can decide if values shall be interpolated by function*r^2
         lr2=0
;
        velx=0.d0
        trilin_complete, vec_cac(0), vec_cac(1), vec_cac(2), $
                         x1, x2, y1, y2, z1, z2, $
                         vala, valb, valc, vald, vale, valf, valg, valh, $
                        rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                         expol, llogx, llogy, llogz, llogf, lr2, velx
;
;y-components
         get_xyz_values2, ndxmax, ndymax, ndzmax, vely3d, $
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                          vala, valb, valc, vald, vale, valf, valg, valh, llogf
;         llogx=0
;         llogy=0
;         llogz=0
;         llogf=0
;here, can decide if values shall be interpolated by function*r^2
         lr2=0
;
        vely=0.d0
        trilin_complete, vec_cac(0), vec_cac(1), vec_cac(2), $
                         x1, x2, y1, y2, z1, z2, $
                         vala, valb, valc, vald, vale, valf, valg, valh, $
                        rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                         expol, llogx, llogy, llogz, llogf, lr2, vely
;
;z-components
         get_xyz_values2, ndxmax, ndymax, ndzmax, velz3d, $
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                          vala, valb, valc, vald, vale, valf, valg, valh, llogf
;         llogx=0
;         llogy=0
;         llogz=0
;         llogf=0
;here, can decide if values shall be interpolated by function*r^2
         lr2=0
;
        velz=0.d0
        trilin_complete, vec_cac(0), vec_cac(1), vec_cac(2), $
                         x1, x2, y1, y2, z1, z2, $
                         vala, valb, valc, vald, vale, valf, valg, valh, $
                         rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                         expol, llogx, llogy, llogz, llogf, lr2, velz
;
         vel_cac(0)=velx
         vel_cac(1)=vely
         vel_cac(2)=velz
;
;transformation of velocity components to slice-coordinates
         vel_slice=transmat_trans#vel_cac
;
         velx2d(i,j)=vel_slice(0)
         vely2d(i,j)=vel_slice(1)

;         print, vec_slice(0), vec_slice(1), vec_slice(2), vec_cac(0), vec_cac(1), vec_cac(2), velx, vely, velz, vel_slice(0), vel_slice(1)
;
      endelse
;
   endfor
endfor

end
