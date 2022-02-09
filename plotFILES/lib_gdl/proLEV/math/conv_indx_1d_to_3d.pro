pro conv_indx_1d_to_3d, indx_1d, ndx, ndy, indx_x, indx_y, indx_z
;
;--------calculation of 3-d indices corresponding to 1-d index----------
;------------with: 1.: 1-d index increases by 1 along x-----------------
;------------------2.: 1-d index increases by nd_x along y--------------
;------------------3.: 1-d index increases by nd_x*nd_y along z---------
;
indx_x = indx_1d mod ndx
indx_y = ((indx_1d - indx_x)/ndx) mod ndy
indx_z = ((indx_1d - indx_x)/ndx - indx_y) / ndy
;
;
end
  
