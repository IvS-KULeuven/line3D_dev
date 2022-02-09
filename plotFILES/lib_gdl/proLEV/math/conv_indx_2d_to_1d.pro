pro conv_indx_2d_to_1d, indx_x, indx_y, ndx, indx_1d
;
;--------calculation of 1-d index corresponding to 2-d indices----------
;------------with: 1.: 1-d index increases by 1 along x-----------------
;------------------2.: 1-d index increases by nd_x along y--------------
;
;
indx_1d = indx_x + ndx * indx_y
;
end
