!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark10_solution
!
!----------------------3d-searchlight beam test-------------------------
!--------------setting opacity and source-functions to zero-------------
!-----------calculating intensities for given direction n_y, n_z--------
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z
use dime3d, only: int3d, scont3d, opac3d, opalbar3d, sline3d, imask_bpoint3d, alocont_o_nn3d, &
                  alocont_nn3d_tmp, normalization3d_tmp, mint3d_tmp, fcontx3d_tmp, &
                  fconty3d_tmp, fcontz3d_tmp
use angles, only: n_x, n_y, n_z
use bcondition, only: xic1, xic1_gdark
use mod_benchmark, only: int3d_sc, int3d_fvm, nn_y, nn_z
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: err
real(dp) :: nn_x
!
! ... local arrays
!
nn_x=nn_y**2+nn_z**2
if(nn_x.ge.1.d0) stop 'error in benchmark10_solution: nn_y^2 + nn_z^2 ge 1'
nn_x = sqrt(1.d0-nn_x)
!
allocate(int3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark10_solution'
allocate(int3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark10_solution'
allocate(alocont_nn3d_tmp(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark10_solution: alocont_nn3d_tmp'
allocate(normalization3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark10_solution: normalization3d_tmp'
allocate(mint3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark10_solution: mint3d_tmp'
allocate(fcontx3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark10_solution: fcontx3d_tmp'
allocate(fconty3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark10_solution: fconty3d_tmp'
allocate(fcontz3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark10_solution: fcontz3d_tmp'
!
!
!
xic1=1.d0
xic1_gdark=xic1_gdark/xic1_gdark(1)  !normalized to polar value
!
scont3d=1.d-10
opac3d=1.d-10
!
n_x = nn_x
n_y = nn_y
n_z = nn_z
!

!short characteristics solution
!call fsc_cont3d(1,1)
call fsc_cont3d_lin(1)
int3d_sc=int3d
!
call ffvm_cont3d(1)
int3d_fvm=int3d
!
!
end subroutine benchmark10_solution
