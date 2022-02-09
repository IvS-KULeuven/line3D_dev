pro newit_jor_coo, data, col_indx, row_indx, diag, b_vec, x_old, x_new, wor, nd, nnz
;
x_new=fltarr(nd)*0.d0
;
for i=0L, nnz-1 do begin
   x_new(row_indx(i)) = x_new(row_indx(i)) + data(i)*x_old(col_indx(i))
endfor
;
x_new=x_old-wor*x_new/diag - wor*b_vec/diag
;
end
