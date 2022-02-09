!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark02_solution
!
!----------------------2d-searchlight beam test-------------------------
!--------------setting opacity and source-functions to zero-------------
!--------------calculating intensities for given angle n_z--------------
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z
use dime3d, only: int3d, scont3d, opac3d, opalbar3d, sline3d, imask_bpoint3d, imask_innreg3d
use angles, only: dim_mu, nodes_mu
use bcondition, only: xic1
use mod_benchmark, only: int2d_sc, int2d_fvm, n_z
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: err
!
!for output to hdf5
!
! ... local arrays
!
!
allocate(int2d_sc(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark02_solution'
allocate(int2d_fvm(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark02_solution'
!
where(imask_innreg3d.ne.0)
   scont3d=1.d-10
endwhere
opac3d=1.d-10
!
nodes_mu=n_z
!
!short characteristics solution
call fsc_cont2d(1)
int2d_sc=int3d(:,ndymax/2+1,:)
!
call ffvm_cont2d(1)
int2d_fvm=int3d(:,ndymax/2+1,:)
!
!
end subroutine benchmark02_solution
