subroutine test_pray
!
   use prog_type
   use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, sline_ray, profile_ray, opac_ray, scont_ray
   use mod_spectests, only: iin
!
   implicit none
!
! ... local scalars
   real(dp) :: iem, iem_c, iem_theo, ierr, iemi, iabs
!
!
!

   write(*,*) '--------------------testing p-ray (formal solution)----------------------------'
   write(*,*)
!
   nz_ray=50
!
   call allocate_fs1d
!
   call setup_pray_test
!
   call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin, iem, iem_c, iemi, iabs)
!
   call formal_ray_theo(iem_theo)
!
   call print_err(iem, iem_theo, ierr)
!
end subroutine test_pray
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_pray_test
!
!----------------set up one arbitrary ray (for test purposes)-----------
!
   use prog_type
   use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, &
      profile_ray, velz_ray
   use mod_spectests, only: opt_test_fs, zmin, zmax, opac_const, opac_max, &
      scont_const
!
   implicit none
!
! ... arguments
!
! ... local scalars
   integer(i4b) :: i
!
!---------------------------make z_ray grid-----------------------------
!
   do i=1, nz_ray
      z_ray(i) = zmin + (i-1.d0) * (zmax-zmin) / (nz_ray-1)
   enddo
!
!-----------------------------------------------------------------------
!
   select case(opt_test_fs)
    case(0)
      opalbar_ray=0.d0
      opac_ray=opac_const
      profile_ray=0.d0
      sline_ray=0.d0
      scont_ray=0.d0
      velz_ray=0.d0
    case(1)
      do i=1, nz_ray
         opac_ray(i)=opac_max/(z_ray(i)+1.d0)/(z_ray(i)+1.d0)
      enddo
      opalbar_ray=0.d0
      profile_ray=0.d0
      sline_ray=0.d0
      scont_ray=0.d0
      velz_ray=0.d0
    case(2)
      do i=1, nz_ray
         opac_ray(i)=opac_max/(z_ray(i)+1.d0)/(z_ray(i)+1.d0)/(z_ray(i)+1.d0)
      enddo
      opalbar_ray=0.d0
      profile_ray=0.d0
      sline_ray=0.d0
      scont_ray=0.d0
      velz_ray=0.d0
    case(3)
      opac_ray=opac_const
      opalbar_ray=0.d0
      profile_ray=0.d0
      scont_ray=scont_const
      sline_ray=0.d0
      velz_ray=0.d0
    case default
      stop 'opt_test_fs not defined'

   end select
!
!-----------------------------------------------------------------------
!
   write(*,'(a4, 7(a20))') '#', 'z [r_star]', 'opac [1/r_star]', 'opal', 'profile-fct', 's_cont', 's_line', 'velz [cm/s]'
   do i=1, nz_ray
      write(*,'(i4, 7(e20.8))') i, z_ray(i), opac_ray(i), opalbar_ray(i), profile_ray(i), scont_ray(i), sline_ray(i), velz_ray(i)
   enddo
   write(*,*)
!
!
end subroutine setup_pray_test
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine formal_ray_theo(iem)
!
!----------------calculate analytycal solution fo test problems---------
!
   use prog_type
   use mod_spectests, only: opt_test_fs, zmin, zmax, opac_const, opac_max, &
      scont_const, iin
!
   implicit none
!
! ... arguments
   real(dp) :: iem
!
! ... local scalars
   integer(i4b) :: i
!
! ... local functions
   real(dp) :: fct_the1, fct_the2, fct_the3
!
!-----------------------------------------------------------------------
!
   select case(opt_test_fs)
    case(0)
      iem=fct_the1(iin, zmax, zmin, opac_const, 0.d0)
    case(1)
      iem=fct_the2(iin, zmax, zmin, opac_max, 0.d0)
    case(2)
      iem=fct_the3(iin, zmax, zmin, opac_max, 0.d0)
    case(3)
      iem=fct_the1(iin, zmax, zmin, opac_const, scont_const)
    case default
      stop 'opt_test_fs not defined'
   end select
!
end subroutine formal_ray_theo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function fct_the1(iin, zmax, zmin, opac_const, scont_const)
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: iin, zmax, zmin, opac_const, scont_const
   real(dp) :: fct_the1
!
   fct_the1 = iin * exp(-opac_const*(zmax-zmin)) + scont_const * (1.d0-exp(-opac_const*(zmax-zmin)))
!
end function fct_the1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function fct_the2(iin, zmax, zmin, opac_max, scont_const)
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: iin, zmax, zmin, opac_max, scont_const
   real(dp) :: fct_the2
!
! ... local scalars
   real(dp) :: dtau
!
   dtau= opac_max * (1.d0/(zmax+1.d0) - 1.d0/(zmin+1.d0))
!
   fct_the2 = iin * exp(dtau) + scont_const * (1.d0-exp(dtau))
!
end function fct_the2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function fct_the3(iin, zmax, zmin, opac_max, scont_const)
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: iin, zmax, zmin, opac_max, scont_const
   real(dp) :: fct_the3
!
! ... local scalars
   real(dp) :: dtau
!
   dtau= 0.5d0 * opac_max * (1.d0/(zmax+1.d0)/(zmax+1.d0) - 1.d0/(zmin+1.d0)/(zmin+1.d0))
!
   fct_the3 = iin * exp(dtau) + scont_const * (1.d0-exp(dtau))
!
end function fct_the3
