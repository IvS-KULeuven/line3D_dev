module mod_opal
!
   use prog_type
   use fund_const
   use mod_iline, only: yhe_lteopal, dir_opal
   
   implicit none
!
   integer(i4b), parameter :: nrho_opal=19
   integer(i4b), parameter :: ntemp_opal=70
!
   real(dp), dimension(nrho_opal) :: rho_opal
   real(dp), dimension(ntemp_opal) :: temp_opal
   real(dp), dimension(ntemp_opal,nrho_opal) :: kappa_opal
!
!
contains

   subroutine get_opal_table()
      !
      !
      !... arguments
      ! real(dp), intent(in) :: yhe_mass
      !
      ! ... local scalars
      integer(i4b) :: i, indx_opal
      real(dp) :: weight_yhe, dist, fdum
      !
      ! ... local characters
      ! character(len=6) :: fname
      character :: chdum
      !
      ! ... local logicals
      logical :: lcheck
      !

      !
      !-----------------------------------------------------------------------
      !
      !
      inquire(file=trim(dir_opal), exist=lcheck)
      !
      if(.not.lcheck) then
         write(*,*) 'error in get_opal_table: file "', trim(dir_opal), '" does not exist'
         stop
      endif
      !
      write(*,*) '--------------------------reading OPAL tables----------------------------------'
      write(*,*) 'reading from file: ', trim(dir_opal)
      write(*,*)
      !
      !-----------------------------------------------------------------------
      !
      kappa_opal=zero
      !
      !read in corresponding kappa
      open(1, file=trim(dir_opal))
      !
      !skip first 4 lines
      do i = 1, 4
         read(1,*)
      enddo

      !read rho
      read(1,*) chdum, rho_opal
      !   write(*,*) chdum !fdum, rho_opal
      read(1,*)
      !
      !read termperature and kappas
      do i=1, ntemp_opal
         read(1,'(f4.2,19f7.3)') temp_opal(i), kappa_opal(i,:)
      enddo

      close(1)
      !
      !write(*,*) yhe, weight_yhe, indx_opal, yhe_opal(indx_opal), fname
      !
      !write(*,*) rho_opal
      !write(*,*) temp_opal
      !
      !do i=1, ntemp_opal
      !   write(*,'(19f10.3)') 10.d0**kappa_opal(i,:)
      !enddo
      !stop 'go on in get_opal_table'

   end subroutine get_opal_table

!
end module mod_opal
