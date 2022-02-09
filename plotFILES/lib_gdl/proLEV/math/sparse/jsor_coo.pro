pro jsor_coo, data, col_indx, row_indx, diag, b_vec, nd, nnz, sol_vec, tol=tol, itmax=itmax, wor=wor
;
;-----------solves a linear system a_mat * x + b_vec = 0----------------
;------------with jacobi-iteration-using over-relaxation----------------
;----------or gauss-seidel-iteration using over-relaxation--------------
;
;   input: data, col_indx, row_indx: matrix of linear system in coo-sparse-matrix-storage
;          diag: diagonal of matrix
;          b_vec: rhs of linear system
;          nd: dimension of linear system
;          nnz: nummber of zero elements
;          wor: over-relaxation-factor
;          opt_sor: true if sor shall be used
;
;   output: sol_vec: solution of the linear system: sol_vec=a_mat^-1 * (-b_vec)
;
;-----------------------------------------------------------------------
;
if(not keyword_set(itmax)) then itmax=200000L
if(not keyword_set(tol)) then tol=1.d-10
if(not keyword_set(wor)) then wor=1.d0
;
solvec_ng=fltarr(4,nd)*0.d0
solvec_new=fltarr(nd)*0.d0
;
;-----------------------------------------------------------------------
;
s1=3L
s2=4L
s3=5L
s4=6L
s5=7L
ng_const=8L
;
indx_max=0L
eps_max=0.d0
;
;-----------------------------------------------------------------------
;
for i=1L, itmax do begin
;calculate new iterate
   newit_jor_coo, data, col_indx, row_indx, diag, b_vec, sol_vec, solvec_new, wor, nd, nnz

;calculate deviation
   calc_dev, sol_vec, solvec_new, nd, eps_max
;
;store new iterate
   sol_vec=solvec_new
;
;-------------------------ng-extrapolation------------------------------
;
   if(i eq s1) then begin
      solvec_ng(0,*)=sol_vec
      s1=s1+ng_const
   endif
   if(i eq s2) then begin
      solvec_ng(1,*)=sol_vec
      s2=s2+ng_const
   endif
   if(i eq s3) then begin
      solvec_ng(2,*)=sol_vec
      s3=s3+ng_const
   endif
   if(i eq s4) then begin
      solvec_ng(3,*)=sol_vec   
      s4=s4+ng_const
   endif
;
   if(i eq s5) then begin
      ng_expol1d, solvec_ng, nd
;
      sol_vec=solvec_ng(0,*)
      s5=s5+ng_const
   endif
;
;-----------------------------------------------------------------------
;
   if(abs(eps_max) lt tol) then begin
      print, "convergence after iteration no. ", i
      print, "max (dev): ", eps_max
      print, ""
      return
   endif
   if(i eq itmax) then begin
      print, 'no convergence after iteration no. ', i
   endif
;
endfor
;
;
;
end
