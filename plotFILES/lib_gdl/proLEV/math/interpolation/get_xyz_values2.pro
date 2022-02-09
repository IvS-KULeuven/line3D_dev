pro get_xyz_values2, ndxmax, ndymax, ndzmax, yvalue3d, $
                     indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                     vala, valb, valc, vald, vale, valf, valg, valh, llogf
;
;-----------------------------------------------------------------------
;
;   get physical quantities of a cube 
;   with vertices given by indx_x1, ... indx_z2
;
;input: dimension of 3-d array: ndxmax, ndymax, ndzmax
;       physical value on grid: yvalue3d
;       indices of cube-vertices: indx_x1, ... indx_z2
;
;output: physical value on cube vertices: vala ... valh
;        flag to decide if logarithmic interpolation is allowed:
;               llogf
;
;-----------------------------------------------------------------------
;
vala = yvalue3d(indx_x1, indx_y1, indx_z1)
valb = yvalue3d(indx_x2, indx_y1, indx_z1)
valc = yvalue3d(indx_x1, indx_y1, indx_z2)
vald = yvalue3d(indx_x2, indx_y1, indx_z2)
vale = yvalue3d(indx_x1, indx_y2, indx_z1)
valf = yvalue3d(indx_x2, indx_y2, indx_z1)
valg = yvalue3d(indx_x1, indx_y2, indx_z2)
valh = yvalue3d(indx_x2, indx_y2, indx_z2)
;
if(vala le 0.d0 or valb le 0.d0 or valc le 0.d0 or vald le 0.d0 or vale le 0.d0 or valf le 0.d0 or valg le 0.d0 or valh le 0.d0) then begin
   llogf=0
endif else begin
   llogf=1
endelse
;
end
