pro conv_indx_1d_to_2d, indx_1d, ndx, indx_x, indx_y
;
;--------calculation of 2-d indices corresponding to 1-d index----------
;------------with: 1.: 1-d index increases by 1 along x-----------------
;------------------2.: 1-d index increases by nd_x along y--------------
;
indx_x = indx_1d mod ndx
indx_y = ((indx_1d - indx_x)/ndx)
;
;
end
  
