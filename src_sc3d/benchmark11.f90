subroutine benchmark11_solution
!
!--------------setting opacity and source-functions to zero-------------
!--------------calculating intensities for different angles-------------
!-----------integrating intensities over complete solid angle-----------
!-----------------------for 2d radiative transfer-----------------------
!------------------------(see subroutine mint_2d)-----------------------
!
use prog_type
use dimecr, only: n1d_cr, r1d_cr, norm1d_cr
use dime3d, only: mint3d, int3d, opac3d, scont3d, normalization3d, imask3d, ndxmax, ndymax, ndzmax, x, y, z, &
                  fcontx3d, fconty3d, fcontz3d, opalbar3d, alocont_nn3d, alocont_o_nn3d, &
                  fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, mint3d_tmp, alocont_nn3d_tmp, normalization3d_tmp
use angles, only: dim_omega
use bcondition, only: xic1, xic2, xic1_gdark
use mod_benchmark, only: mint3d_theo, mint3d_sc, mint3d_fvm, r1d_angdep, n1d_angdep, &
                         intsc_angdep2, intfvm_angdep2, xcoord_angdep2, ycoord_angdep2, zcoord_angdep2, &
                         fcontr3d_theo, fcontr3d_sc, fcontr3d_fvm, fcontth3d_theo, fcontth3d_sc, &
                         fcontth3d_fvm, fcontphi3d_theo, fcontphi3d_sc, fcontphi3d_fvm
use fund_const, only: pi
use params_input, only: rmax
use mod_interp1d, only: find_index
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, oindx, err
integer(i4b) :: nr_dum, iim2, iim1, ii, iip1, n1d_angdep2
real(dp) :: dilfac, fdum, rad, theta, phi
!
! ... local arrays
integer(i4b), dimension(:), allocatable :: indxx_angdep, indxy_angdep, indxz_angdep
!
!-----------------------------------------------------------------------
!
allocate(mint3d_theo(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(mint3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(mint3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
!
allocate(fcontr3d_theo(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontr3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontr3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontth3d_theo(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontth3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontth3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontphi3d_theo(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontphi3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(fcontphi3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
!
!for angular dependence of intensity at 6 different radii
n1d_angdep2=6
n1d_angdep=n1d_angdep2*7
!
allocate(xcoord_angdep2(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(ycoord_angdep2(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(zcoord_angdep2(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(indxx_angdep(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(indxy_angdep(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(indxz_angdep(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
!
allocate(intsc_angdep2(n1d_angdep, dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
allocate(intfvm_angdep2(n1d_angdep, dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
!
allocate(r1d_angdep(n1d_angdep2), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark11'
!
!
!
!define radii at which intensity shall be stored
do i=1, n1d_angdep2
   r1d_angdep(n1d_angdep2+1-i) = 1.d0+(rmax-1.d0)/(2.d0**i)
enddo
r1d_angdep=(/3.09d0,3.5d0,4.01d0,4.57d0,5.29d0,9.05d0/)
!
!inensities as function of angle at points on x-axis
do i=1, n1d_angdep2
   xcoord_angdep2(i) = r1d_angdep(i)
   ycoord_angdep2(i) = -0.01d0
   zcoord_angdep2(i) = -0.01d0
enddo
!inensities as function of angle at points on y-axis
do i=1, n1d_angdep2
   xcoord_angdep2(n1d_angdep2+i) = -0.01d0
   ycoord_angdep2(n1d_angdep2+i) = r1d_angdep(i)
   zcoord_angdep2(n1d_angdep2+i) = -0.01d0
enddo
!inensities as function of angle at points on z-axis
do i=1, n1d_angdep2
   xcoord_angdep2(2*n1d_angdep2+i) = -0.01d0
   ycoord_angdep2(2*n1d_angdep2+i) = -0.01d0
   zcoord_angdep2(2*n1d_angdep2+i) = r1d_angdep(i)
enddo
!inensities as function of angle along 45 degree in x-y plane
do i=1, n1d_angdep2
   xcoord_angdep2(3*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
   ycoord_angdep2(3*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
   zcoord_angdep2(3*n1d_angdep2+i) = -0.01d0
enddo
!inensities as function of angle along 45 degree in x-z plane
do i=1, n1d_angdep2
   xcoord_angdep2(4*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
   ycoord_angdep2(4*n1d_angdep2+i) = -0.01d0
   zcoord_angdep2(4*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
enddo
!inensities as function of angle along 45 degree in y-z plane
do i=1, n1d_angdep2
   xcoord_angdep2(5*n1d_angdep2+i) = -0.01d0
   ycoord_angdep2(5*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
   zcoord_angdep2(5*n1d_angdep2+i) = r1d_angdep(i)/sqrt(2.d0)
enddo
!inensities as function of angle along n=(1/3,1/3,1/3)
do i=1, n1d_angdep2
   xcoord_angdep2(6*n1d_angdep2+i) = r1d_angdep(i)/sqrt(3.d0)
   ycoord_angdep2(6*n1d_angdep2+i) = r1d_angdep(i)/sqrt(3.d0)
   zcoord_angdep2(6*n1d_angdep2+i) = r1d_angdep(i)/sqrt(3.d0)
enddo
!
!corresponding indices
do i=1, n1d_angdep
   call find_index(xcoord_angdep2(n1d_angdep+1-i), x, ndxmax, iim2, iim1, ii, iip1)
   xcoord_angdep2(n1d_angdep+1-i) = x(ii)
   indxx_angdep(n1d_angdep+1-i) = ii
   call find_index(ycoord_angdep2(n1d_angdep+1-i), y, ndymax, iim2, iim1, ii, iip1)
   ycoord_angdep2(n1d_angdep+1-i) = y(ii)
   indxy_angdep(n1d_angdep+1-i) = ii
   call find_index(zcoord_angdep2(n1d_angdep+1-i), z, ndzmax, iim2, iim1, ii, iip1)
   zcoord_angdep2(n1d_angdep+1-i) = z(ii)
   indxz_angdep(n1d_angdep+1-i) = ii
enddo
!
!-----------------------------------------------------------------------
!
xic1=1.d0
xic2=0.d0
xic1_gdark=1.d0
!
opac3d=1.d-10
opalbar3d=1.d-10
where(imask3d.ne.4) scont3d=1.d-10
!
!------------calculating intensities for given nodes--------------------
!
mint3d=0.d0
mint3d_sc=0.d0
mint3d_fvm=0.d0
normalization3d=0.d0
fcontx3d=0.d0
fconty3d=0.d0
fcontz3d=0.d0
fcontr3d_sc=0.d0
fcontth3d_sc=0.d0
fcontphi3d_sc=0.d0
fcontr3d_fvm=0.d0
fcontth3d_fvm=0.d0
fcontphi3d_fvm=0.d0
!
write(*,*) '------------------calculating intensities for all mu, phi----------------------'
write(*,*)
!
if(allocated(int3d)) deallocate(int3d)
if(allocated(alocont_o_nn3d)) deallocate(alocont_o_nn3d)
!
!-------------------------short characteristics solution----------------
!
write(*,*) 'short characteristics solution'
write(*,*)
!
!$omp parallel &
!$omp private(err, oindx,i,j,k,l)
!
!-------------allocation of global (threadprivate) arrays---------------
!
allocate(mint3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'

allocate(alocont_nn3d_tmp(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(normalization3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fcontx3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fconty3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fcontz3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(alocont_o_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(int3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
!
!$omp do schedule(static)
do oindx=1, dim_omega
   write(*,*) 'calculating omega', oindx, dim_omega
!   call fsc_cont3d(oindx)
   call fsc_cont3d_lin(oindx)
!   call ffvm_cont3d(oindx)
!
   do l=1, n1d_angdep
      do i=1, ndxmax
         if(indxx_angdep(l).eq.i) then
            do j=1, ndymax
               if(indxy_angdep(l).eq.j) then
                  do k=1, ndzmax
                     if(indxz_angdep(l).eq.k) then
                        intsc_angdep2(l,oindx) = int3d(i,j,k)
                     endif 
                  enddo
               endif
            enddo
         endif
      enddo
   enddo
!
enddo
!$omp enddo
!
!------------------------add up temporary arrays------------------------
!
!$omp critical
   alocont_nn3d = alocont_nn3d+alocont_nn3d_tmp
   normalization3d = normalization3d + normalization3d_tmp
   mint3d = mint3d + mint3d_tmp
   fcontx3d = fcontx3d + fcontx3d_tmp
   fconty3d = fconty3d + fconty3d_tmp
   fcontz3d = fcontz3d + fcontz3d_tmp
!$omp end critical
!
!------------deallocation of global (threadprivate) arrays--------------
!
   deallocate(alocont_nn3d_tmp)
   deallocate(normalization3d_tmp)
   deallocate(mint3d_tmp)
   deallocate(fcontx3d_tmp)
   deallocate(fconty3d_tmp)
   deallocate(fcontz3d_tmp)
   deallocate(alocont_o_nn3d)
   deallocate(int3d)
!
!$omp end parallel
!
!
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               rad = x(i)**2+y(j)**2+z(k)**2
               call get_angles_spc(x(i), y(j), z(k), theta, phi)
               fcontr3d_sc(i,j,k) = fcontx3d(i,j,k)*sin(theta)*cos(phi) + fconty3d(i,j,k)*sin(theta)*sin(phi) + fcontz3d(i,j,k)*cos(theta)
               fcontth3d_sc(i,j,k) = fcontx3d(i,j,k)*cos(theta)*cos(phi) + fconty3d(i,j,k)*cos(theta)*sin(phi) - fcontz3d(i,j,k)*sin(theta)
               fcontphi3d_sc(i,j,k) = -fcontx3d(i,j,k)*sin(phi) + fconty3d(i,j,k)*cos(phi)
            case default
         end select
      enddo
   enddo
enddo
!
mint3d_sc=mint3d


!j=ndymax/2+1
!k=ndzmax/2+1
!do i=1, ndzmax
!   fdum=abs(x(i))**2 / xic1
!   write(*,'(7es20.8)') x(i), mint3d_sc(i,j,k)*fdum, fcontr3d_sc(i,j,k)*fdum, fcontr3d_sc(i,j,k)
!enddo
!write(*,*)
!
!-------------------------finite volume method solution-----------------
!
mint3d=0.d0
fcontx3d=0.d0
fconty3d=0.d0
fcontz3d=0.d0
normalization3d=0.d0
!
write(*,*) 'finite volume solution'
write(*,*)
!
!
!$omp parallel &
!$omp private(err, oindx,i,j,k,l)
!
!-------------allocation of global (threadprivate) arrays---------------
!
allocate(mint3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(alocont_nn3d_tmp(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(normalization3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fcontx3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fconty3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(fcontz3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(alocont_o_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
allocate(int3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error benchmark11'
!
!$omp do schedule(static)
do oindx=1, dim_omega
   write(*,*) 'calculating omega', oindx, dim_omega
!
   call ffvm_cont3d(oindx)
!
   do l=1, n1d_angdep
      do i=1, ndxmax
         if(indxx_angdep(l).eq.i) then
            do j=1, ndymax
               if(indxy_angdep(l).eq.j) then
                  do k=1, ndzmax
                     if(indxz_angdep(l).eq.k) then
                        intfvm_angdep2(l,oindx) = int3d(i,j,k)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
   enddo
!
enddo
!$omp enddo
!
!------------------------add up temporary arrays------------------------
!
!$omp critical
   alocont_nn3d = alocont_nn3d+alocont_nn3d_tmp
   normalization3d = normalization3d + normalization3d_tmp
   mint3d = mint3d + mint3d_tmp
   fcontx3d = fcontx3d + fcontx3d_tmp
   fconty3d = fconty3d + fconty3d_tmp
   fcontz3d = fcontz3d + fcontz3d_tmp
!$omp end critical
!
!------------deallocation of global (threadprivate) arrays--------------
!
   deallocate(alocont_nn3d_tmp)
   deallocate(normalization3d_tmp)
   deallocate(mint3d_tmp)
   deallocate(fcontx3d_tmp)
   deallocate(fconty3d_tmp)
   deallocate(fcontz3d_tmp)
   deallocate(alocont_o_nn3d)
   deallocate(int3d)
!
!$omp end parallel
!
!write(*,'(a5, 10a20)') '#', 'j(x+)', 'j(x-)', 'j(y+)', 'j(y-)', 'j(z+)', 'j(z-)', 'norm(x+)', 'norm(y+)', 'norm(z+)'
!do i=ndxmax/2+1, ndxmax
!   write(*,fmt='(i5, 10(e20.8))') i-ndxmax/2, z(i), &
!                                   mint3d(i,ndymax/2+1,ndzmax/2+1), &
!                                   mint3d(i,ndymax/2+1,ndzmax/2+1), &
!                                   mint3d(ndxmax/2+1,i,ndzmax/2+1), &
!                                   mint3d(ndxmax/2+1,i,ndzmax/2+1), &
!                                   mint3d(ndxmax/2+1,ndymax/2+1,i), &
!                                   mint3d(ndxmax/2+1,ndymax/2+1,i), &
!                           normalization3d(i,ndymax/2+1,ndzmax/2+1), &
!                           normalization3d(ndxmax/2+1,i,ndzmax/2+1), &
!                           normalization3d(ndxmax/2+1,ndymax/2+1,i)
!end do
!
mint3d_fvm=mint3d
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               rad = x(i)**2+y(j)**2+z(k)**2
               call get_angles_spc(x(i), y(j), z(k), theta, phi)
               fcontr3d_fvm(i,j,k) = fcontx3d(i,j,k)*sin(theta)*cos(phi) + fconty3d(i,j,k)*sin(theta)*sin(phi) + fcontz3d(i,j,k)*cos(theta)
               fcontth3d_fvm(i,j,k) = fcontx3d(i,j,k)*cos(theta)*cos(phi) + fconty3d(i,j,k)*cos(theta)*sin(phi) - fcontz3d(i,j,k)*sin(theta)
               fcontphi3d_fvm(i,j,k) = -fcontx3d(i,j,k)*sin(phi) + fconty3d(i,j,k)*cos(phi)
            case default
         end select
      enddo
   enddo
enddo
!
!---------------------theoretical solution------------------------------
!
mint3d_theo=0.d0
fcontr3d_theo=0.d0
fcontth3d_theo=0.d0
fcontphi3d_theo=0.d0
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               rad = sqrt(x(i)**2+y(j)**2+z(k)**2)
               dilfac=1.d0-sqrt(1.d0-1.d0/rad**2)

               mint3d_theo(i,j,k) = 0.5d0*xic1*dilfac
               fcontr3d_theo(i,j,k) = xic1/4.d0/rad**2
            case default
         end select
      enddo
   enddo
enddo
!
write(*,'(7a20)') 'x', 'J(theo)*r^2/I_c', 'J(3d-sc)*r^2/I_c', 'J(3d-fvm)*r^2/I_c', 'F(theo)*r^2/I_c', 'F(3d-sc)*r^2/I_c', 'F(3d-fvm)*r^2)/I_c'
j=ndymax/2+1
k=ndzmax/2+1
do i=1, ndzmax
   fdum=abs(x(i))**2 / xic1
   write(*,'(7es20.8)') x(i), mint3d_theo(i,j,k)*fdum, mint3d_sc(i,j,k)*fdum, mint3d_fvm(i,j,k)*fdum, &
                        fcontr3d_theo(i,j,k)*fdum, fcontr3d_sc(i,j,k)*fdum, fcontr3d_fvm(i,j,k)*fdum
enddo
write(*,*)
!
write(*,'(7a20)') 'y', 'J(theo)*r^2/I_c', 'J(3d-sc)*r^2/I_c', 'J(3d-fvm)*r^2/I_c', 'F(theo)*r^2/I_c', 'F(3d-sc)*r^2/I_c', 'F(3d-fvm)*r^2)/I_c'
i=ndxmax/2+1
k=ndzmax/2+1
do j=1, ndzmax
   fdum=abs(y(j))**2 / xic1
   write(*,'(7es20.8)') y(j), mint3d_theo(i,j,k)*fdum, mint3d_sc(i,j,k)*fdum, mint3d_fvm(i,j,k)*fdum, &
                        fcontr3d_theo(i,j,k)*fdum, fcontr3d_sc(i,j,k)*fdum, fcontr3d_fvm(i,j,k)*fdum
enddo
write(*,*)
!
write(*,'(7a20)') 'z', 'J(theo)*r^2/I_c', 'J(3d-sc)*r^2/I_c', 'J(3d-fvm)*r^2/I_c', 'F(theo)*r^2/I_c', 'F(3d-sc)*r^2/I_c', 'F(3d-fvm)*r^2)/I_c'
i=ndxmax/2+1
j=ndymax/2+1
do k=1, ndzmax
   fdum=abs(z(k))**2 / xic1
   write(*,'(7es20.8)') z(k), mint3d_theo(i,j,k)*fdum, mint3d_sc(i,j,k)*fdum, mint3d_fvm(i,j,k)*fdum, &
                        fcontr3d_theo(i,j,k)*fdum, fcontr3d_sc(i,j,k)*fdum, fcontr3d_fvm(i,j,k)*fdum
enddo
write(*,*)
!
end subroutine benchmark11_solution
