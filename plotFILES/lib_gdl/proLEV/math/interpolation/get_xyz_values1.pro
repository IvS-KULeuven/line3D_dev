pro get_xyz_values1, x_coord, y_coord, z_coord, xarr, yarr, zarr, ndxmax, ndymax, ndzmax, $
                     indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                     x1, x2, y1, y2, z1, z2, $
                     rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                     llogx, llogy, llogz
;
;-----------------------------------------------------------------------
;
;   get coordinates and radii of a cube
;   with vertices given by indx_x1, ... indx_z2
;
;input: dimension of arrays x,y,z: ndxmax, ndymax, ndzmax
;       arrays x,y,z
;       coordinates of a point inside cube: x_coord, y_coord, z_coord
;       indices of cube-vertices: indx_x1, ... indx_z2
;
;output: grid value: x1, x2, y1, y2, z1, z2
;        radii of vertices and of point: rada, ... radh, radp
;        flags to decide if logarithmic interpolation is allowed:
;               llogx, llogy, llogz
;
;-----------------------------------------------------------------------
;
x1=xarr(indx_x1)
x2=xarr(indx_x2)
y1=yarr(indx_y1)
y2=yarr(indx_y2)
z1=zarr(indx_z1)
z2=zarr(indx_z2)
;
rada = sqrt(x1*x1 + y1*y1 + z1*z1)
radb = sqrt(x2*x2 + y1*y1 + z1*z1)
radc = sqrt(x1*x1 + y1*y1 + z2*z2)
radd = sqrt(x2*x2 + y1*y1 + z2*z2)
rade = sqrt(x1*x1 + y2*y2 + z1*z1)
radf = sqrt(x2*x2 + y2*y2 + z1*z1)
radg = sqrt(x1*x1 + y2*y2 + z2*z2)
radh = sqrt(x2*x2 + y2*y2 + z2*z2)
;
radp = sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord)
;
;-------------check if logarithmic interpolation is allowed-------------
;
if(x1*x2 le 0.d0) then begin
   llogx=0.
endif else begin
   if(x1*x_coord le 0.d0) then begin
      llogx=0
   endif else begin
      llogx=1
   endelse
endelse
;
if(y1*y2 le 0.d0) then begin
   llogy=0
endif else begin
   if(y1*y_coord le 0.d0) then begin
      llogy=0
   endif else begin
      llogy=1
   endelse
endelse
;
if(z1*z2 le 0.d0) then begin
   llogz=0
endif else begin
   if(z1*z_coord le 0.d0) then begin
      llogz=0
   endif else begin
      llogz=1
   endelse
endelse
;
;
;
end
