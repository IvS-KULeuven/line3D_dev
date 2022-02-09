subroutine benchmark03_solution
!
!----------------------2d-searchlight beam test-------------------------
!--------------setting opacity and source-functions to zero-------------
!--------------calculating intensities for different angles-------------
!
use prog_type
use dime3d, only: int3d, scont3d, opac3d, ndxmax, ndymax, ndzmax
use angles, only: dim_mu, nodes_mu
use bcondition, only: xic1
use mod_benchmark, only: int2dsc_angdep, int2dfvm_angdep
!
implicit none
!
! ... local scalars
integer(i4b) :: indx_mu, err
!
!
allocate(int2dsc_angdep(ndxmax,ndzmax,dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark03_solution'
allocate(int2dfvm_angdep(ndxmax,ndzmax,dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark03_solution'
!
xic1=1.d0
!
scont3d=1.d-10
opac3d=1.d-10
!
do indx_mu=1, dim_mu
   write(*,*) 'calculating (mu)', nodes_mu(indx_mu)
!
   call fsc_cont2d(indx_mu)
   int2dsc_angdep(:,:,indx_mu)=int3d(:,ndymax/2+1,:)
!
   call ffvm_cont2d(indx_mu)
   int2dfvm_angdep(:,:,indx_mu)=int3d(:,ndymax/2+1,:)
!
enddo
!
!
end subroutine benchmark03_solution
