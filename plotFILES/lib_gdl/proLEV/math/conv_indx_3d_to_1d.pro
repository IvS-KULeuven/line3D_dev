pro conv_indx_3d_to_1d, indx_x, indx_y, indx_z, ndx, ndy, indx_1d
;
;--------calculation of 1-d index corresponding to 3-d indices----------
;------------with: 1.: 1-d index increases by 1 along x-----------------
;------------------2.: 1-d index increases by nd_x along y--------------
;------------------3.: 1-d index increases by nd_x*nd_y along z---------
;
;
indx_1d = indx_x + ndx * indx_y + ndx * ndy * indx_z
;
end
