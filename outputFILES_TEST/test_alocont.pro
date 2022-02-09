pro test_alocont, dir=dir, windx=windx, oname=oname
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/alo_cont.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_alo')
         nd_alo=h5a_read(att_id)
         nd_alo=nd_alo(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_diag')
         nd_diag=h5a_read(att_id)
         nd_diag=nd_diag(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'alo')
      dset_id=h5d_open(group_id, 'alocont_data')
         alocont_data=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'alocont_data_diag')
         alocont_data_diag=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'alocont_rowindx')
         alocont_rowindx=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'alocont_colindx')
         alocont_colindx=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'b_vec')
         b_vec=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
alocont_colindx=alocont_colindx-1
alocont_rowindx=alocont_rowindx-1
;
;-----------------------define range------------------------------------
;

xlim=[1,ndxmax*ndymax*ndzmax]
ylim=[ndxmax*ndymax*ndzmax,1]

;x1=7000
;x2=13000
;xlim=[x1,x2]
;ylim=[x2,x1]
xlim=[72000.,72200.]
ylim=[74400.,74100.]
;
;----------------------------title-strings------------------------------
;
titlestr='ALO'
xtitlestr=textoidl('m (col)')
ytitlestr=textoidl('n (row)')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/alo_cont.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx;, xsize=950, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
plot, [0,0],[0,0], $
      xrange=xlim, yrange=ylim, psym=3, /xs, /ys, /isotropic
for i=0, nd_alo-1 do begin
   if(alocont_data(i) lt 0.d0) then begin
      oplot, [alocont_colindx(i),alocont_colindx(i)], [alocont_rowindx(i),alocont_rowindx(i)], $
             psym=3, color=ci_red
   endif
   if(alocont_data(i) gt 0.d0) then begin
      oplot, [alocont_colindx(i),alocont_colindx(i)], [alocont_rowindx(i),alocont_rowindx(i)], $
             psym=3, color=ci_green
   endif

endfor
;
;if(keyword_set(oname)) then begin
;   legend, [lstr1, lstr2, lstr3, lstr4], $
;           psym=[0,0,3,3], $
;           linestyle=[0,2,0,0], $
;           color=[ci_black, ci_black,ci_blue,ci_red], $
;           textcolor=ci_black, $
;           /right_legend
;endif else begin
;   legend, [lstr1, lstr2, lstr3, lstr4], $
;           psym=[0,0,3,3], $
;           linestyle=[0,2,0,0], $
;           color=[ci_white,ci_white,ci_blue,ci_red], $
;           textcolor=ci_white, $
;           /right_legend
;endelse
;
;
;--------------------------test inversion-------------------------------
;
sol_vec=fltarr(nd_diag)*0.d0
jsor_coo, alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, b_vec, nd_diag, nd_alo, sol_vec


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro jsor_coo, data, col_indx, row_indx, diag, b_vec, nd, nnz, sol_vec
;
;-----------solves a linear system a_mat * x + b_vec = 0----------------
;--------------------with jacobi-iteration------------------------------
;
;   input: data, col_indx, row_indx: matrix of linear system in coo-sparse-matrix-storage
;          diag: diagonal of matrix
;          b_vec: rhs of linear system
;          nd: dimension of linear system
;          nnz: nummber of zero elements
;
;   output: sol_vec: solution of the linear system: sol_vec=a_mat^-1 * (-b_vec)
;
;-----------------------------------------------------------------------
;
itmax=10
;
for i=0, itmax-1 do begin
;
;-----------------------------------------------------------------------
;
   newit_jor_coo, data, col_indx, row_indx, diag, b_vec, sol_vec, sol_vec_new, nd, nnz


   window, i
   plot, findgen(nd), sol_vec_new
;
;   call calc_dev(sol_vec, sol_vec_new, nd, eps_max)
;
;   sol_vec=sol_vec_new
;
;   if(abs(eps_max).lt.dev_max) then
;      write(*,*) "convergence after iteration no. ", i
;      write(*,*) "max (dev): ", eps_max
;      write(*,*)
;      exit
;   else if(i.eq.itmax) then
;      write(*,*) 'no convergence after iteration no. ', i
;      stop
;   end if
;
endfor
;
matmul_coo, data, col_indx, row_indx, sol_vec, y_vec, nd, nnz
print, y_vec
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro newit_jor_coo, data, col_indx, row_indx, diag, b_vec, x_old, x_new, nd, nnz
;
;-----------------------------------------------------------------------
;
x_new=fltarr(nd)*0.d0
;
for i=0, nnz-1 do begin
   x_new(row_indx(i)) = x_new(row_indx(i)) + data(i)*x_old(col_indx(i))
endfor
;
x_new=x_old-x_new/diag - b_vec/diag
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro matmul_coo, data, col_indx, row_indx, x_vec, y_vec, nd, nnz
;
;-----------------------format: coordinate list-------------------------
;
;input: x_vec: input - vector which shall be multiplied
;       data: storage of matrix row by row
;       col_indx: column indices for each data-elements
;       row_indx: row indices fo each data-elements
;       nd: dimension of nxn matrix
;       nnz: number of non-zero-matrix elements
;
;output: y_vec = matrix * x_vec
;
;-----------------------------------------------------------------------
;
y_vec=fltarr(nd)*0.d0
;
for i=0, nnz-1 do begin
   y_vec(row_indx(i)) = y_vec(row_indx(i)) + data(i)*x_vec(col_indx(i))
endfor
;
end
