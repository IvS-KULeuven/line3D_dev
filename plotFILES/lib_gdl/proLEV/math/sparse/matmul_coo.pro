pro matmul_coo, data, col_indx, row_indx, x_vec, y_vec, nd, nnz
;
;-----------------------format: coordinate list-------------------------
;
;input: x_vec: input - vector which shall be multiplied
;       data: storage of matrix row by row
;       col_indx: column indices for each data-elements
;       row_indx: row indices fo each data-elements
;       nd: dimension x and y vectors
;       nnz: number of non-zero-matrix elements
;
;output: y_vec = matrix * x_vec
;
;-----------------------------------------------------------------------
;
y_vec=fltarr(nd)*0.d0
;
for i=0L, nnz-1 do begin
   y_vec(row_indx(i)) = y_vec(row_indx(i)) + data(i)*x_vec(col_indx(i))
endfor
;
end
