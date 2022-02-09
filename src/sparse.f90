!
!***********************************************************************
!        routine for matrix - vector multiplication in coo format
!      routine for matrix inversion (jacobi iteration) in coo format
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine matmul_coo(data, col_indx, row_indx, x_vec, y_vec, nd, nnz)
!
!-----------------------format: coordinate list-------------------------
!
!input: x_vec: input - vector which shall be multiplied
!       data: storage of matrix row by row
!       col_indx: column indices for each data-elements
!       row_indx: row indices fo each data-elements
!       nd: dimension of nxn matrix
!       nnz: number of non-zero-matrix elements
!
!output: y_vec = matrix * x_vec
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd, nnz
real(dp), dimension(nnz), intent(in) :: data
integer(i4b), dimension(nnz), intent(in) :: col_indx, row_indx
real(dp), dimension(nd), intent(in) :: x_vec
real(dp), dimension(nd) :: y_vec
!
! ... local scalars
integer(i4b) :: i, j
!
!-----------------------------------------------------------------------
!
y_vec=0.d0
!
do i=1, nnz
   y_vec(row_indx(i)) = y_vec(row_indx(i)) + data(i)*x_vec(col_indx(i))
enddo
!
end subroutine matmul_coo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine jsor_coo(data, col_indx, row_indx, diag, b_vec, w, nd, nnz, opt_sor, sol_vec)
!
!-----------solves a linear system a_mat * x + b_vec = 0----------------
!------------with jacobi-iteration-using over-relaxation----------------
!----------or gauss-seidel-iteration using over-relaxation--------------
!
!   input: data, col_indx, row_indx: matrix of linear system in coo-sparse-matrix-storage
!          diag: diagonal of matrix
!          b_vec: rhs of linear system
!          nd: dimension of linear system
!          nnz: nummber of zero elements
!          w: over-relaxation-factor
!          opt_sor: true if sor shall be used
!
!   output: sol_vec: solution of the linear system: sol_vec=a_mat^-1 * (-b_vec)
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd, nnz
real(dp), intent(in) :: w
logical, intent(in) :: opt_sor
real(dp), dimension(nd) :: sol_vec
real(dp), dimension(nd), intent(in) :: b_vec, diag
real(dp), dimension(nnz), intent(in) :: data
integer(i4b), dimension(nnz), intent(in) ::  col_indx, row_indx
!
! ... local scalars
integer(i4b) :: i, err
integer(i4b) :: indx_max, s1, s2, s3, s4, s5
integer(i4b), parameter :: itmax=200000
integer(i4b), parameter :: ng_const=8               !ng-extrapolation is performed at each 8th iterate
real(dp) :: eps_max
real(dp), parameter :: dev_max=1.d-10
!
! ... local arrays
real(dp), dimension(:), allocatable :: sol_vec_new
real(dp), dimension(:,:), allocatable :: solvec_ng
real(dp), dimension(:), allocatable :: data_crs
integer(i4b), dimension(:), allocatable :: row_crs, col_crs
!
! ... local logicals
!
!
!------------transform to crs storage format if sor-iteration is used---
!
if(opt_sor) then
   allocate(data_crs(nnz), stat=err)
      if(err.ne.0) stop 'allocation error jsor_coo: data_crs'
   allocate(col_crs(nnz), stat=err)
      if(err.ne.0) stop 'allocation error jsor_coo: col_indx_crs'
   allocate(row_crs(nd+1), stat=err)
      if(err.ne.0) stop 'allocation error jsor_coo: row_ptr_crs'
!
   call coocsr(nd,nnz,data,row_indx,col_indx,data_crs,col_crs,row_crs)
endif
!
allocate(solvec_ng(4,nd), stat=err)
   if(err.ne.0) stop 'allocation error jsor_coo: solvec_ng'
!
allocate(sol_vec_new(nd), stat=err)
   if(err.ne.0) stop 'allocation error jsor_coo: sol_vec_new'
!
!-----------------------------------------------------------------------
!
solvec_ng=0.d0
s1=3
s2=4
s3=5
s4=6
s5=7
!
indx_max=0
eps_max=0.d0
!
!-----------------------------------------------------------------------
!
do i=1, itmax
!
!-----------------------------------------------------------------------
!
   if(opt_sor) then
      call newit_sor_crs(data_crs, col_crs, row_crs, diag, b_vec, sol_vec, sol_vec_new, w, nd, nnz)
   else
      call newit_jor_coo(data, col_indx, row_indx, diag, b_vec, sol_vec, sol_vec_new, w, nd, nnz)
   endif
!
   call calc_dev(sol_vec, sol_vec_new, nd, eps_max)
!
   sol_vec=sol_vec_new
!
!-------------------------ng-extrapolation------------------------------
!
   if(i.eq.s1) then
      solvec_ng(1,:)=sol_vec
      s1=s1+ng_const
   elseif(i.eq.s2) then
      solvec_ng(2,:)=sol_vec
      s2=s2+ng_const
   elseif(i.eq.s3) then
      solvec_ng(3,:)=sol_vec
      s3=s3+ng_const
   elseif(i.eq.s4) then
      solvec_ng(4,:)=sol_vec   
      s4=s4+ng_const
   endif
!
   if(i.eq.s5) then
      call ng_expol1d(solvec_ng, nd)
!      call ait_expol1d(solvec_ng, nd)
      sol_vec=solvec_ng(1,:)
      s5=s5+ng_const
   endif
!
!-----------------------------------------------------------------------
!
   if(abs(eps_max).lt.dev_max) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      write(*,*)
      exit
   else if(i.eq.itmax) then
      write(*,*) 'no convergence after iteration no. ', i
      stop
   end if
!
enddo
!
!
!
end subroutine jsor_coo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine newit_jor_coo(data, col_indx, row_indx, diag, b_vec, x_old, x_new, w, nd, nnz)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd, nnz
real(dp), intent(in) :: w
!
real(dp), dimension(nd) :: x_new
real(dp), dimension(nd), intent(in) :: x_old, b_vec, diag
real(dp), dimension(nnz), intent(in) :: data
integer(i4b), dimension(nnz), intent(in) :: col_indx, row_indx
!
! ... local scalars
integer(i4b) :: i
real(dp) :: sum
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
x_new=0.d0
!
do i=1, nnz
   x_new(row_indx(i)) = x_new(row_indx(i)) + data(i)*x_old(col_indx(i))
enddo
!
x_new=x_old-w*x_new/diag - w*b_vec/diag
!
!
end subroutine newit_jor_coo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine newit_sor_crs(data, col_indx, row_ptr, diag, b_vec, x_old, x_new, w, nd, nnz)
!
use prog_type
!
implicit none
!
! ... arguments
!
! ... arguments
integer(i4b), intent(in) :: nd, nnz
real(dp), intent(in) :: w
!
real(dp), dimension(nd) :: x_new
real(dp), dimension(nd), intent(in) :: x_old, b_vec, diag
real(dp), dimension(nnz), intent(in) :: data
integer(i4b), dimension(nnz), intent(in) :: col_indx
integer(i4b), dimension(nd+1), intent(in) :: row_ptr
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: sum0, sum1
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
x_new=0.d0
!
do i=1, nd
   sum0=0.d0
   sum1=0.d0
!
   do j=row_ptr(i), row_ptr(i+1)-1
      if(col_indx(j).lt.i) then
         sum0 = sum0 + data(j) * x_new(col_indx(j))
      else
         sum1 = sum1 + data(j) * x_old(col_indx(j))
      endif
   enddo
!   
   x_new(i) = x_old(i)  - (sum0 + sum1 + b_vec(i)) * w / diag(i)
!
enddo
!
!
!
end subroutine newit_sor_crs
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!
!from sparskit package: http://www-users.cs.umn.edu/~saad/software/sparskit/
use prog_type
use fund_const

!
!----------------------------------------------------------------------- 
real(dp) :: a(*), ao(*), x
integer(i4b) :: ir(*),jc(*),jao(*),iao(*)
!-----------------------------------------------------------------------
!  coordinate     to   compressed sparse row 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!--------- 
! nrow	= dimension of the matrix 
! nnz	= number of nonzero elements in matrix
! a,
! ir, 
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column 
!	  number. the order of the elements is arbitrary. 
!
! on return:
!----------- 
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao 
! 	continung the real values, jao containing the column indices, 
!	and iao being the pointer to the beginning of the row, 
!	in arrays ao, jao.
!
! notes:
!------ this routine is not in place.  see coicsr
!
!------------------------------------------------------------------------

      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
! determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of coocsr ------------------------------------------- 
!----------------------------------------------------------------------- 
      end subroutine coocsr

