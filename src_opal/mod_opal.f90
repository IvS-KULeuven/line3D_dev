module mod_opal
!
   use prog_type
   use fund_const
!
   implicit none
!
   integer(i4b), parameter :: nyhe_opal=8
   integer(i4b), parameter :: nrho_opal=19
   integer(i4b), parameter :: ntemp_opal=70
!
   real(dp), dimension(nyhe_opal), parameter :: yhe_opal=(/0.00999d0, 0.1d0, 0.28d0, 0.59d0, 0.61d0, 0.98d0, 0.99d0, 1.d0 /)
   real(dp), dimension(nrho_opal) :: rho_opal
   real(dp), dimension(ntemp_opal) :: temp_opal
   real(dp), dimension(ntemp_opal,nrho_opal) :: kappa_opal
!
   character(len=25), parameter :: dir_opal='opal_tables'
!
contains

   subroutine get_opal_table(yhe)
      !
      !
      !... arguments
      real(dp), intent(in) :: yhe
      !
      ! ... local scalars
      integer(i4b) :: i, indx_opal
      real(dp) :: weight_yhe, dist, fdum
      !
      ! ... local characters
      character(len=6) :: fname
      character :: chdum
      !
      ! ... local logicals
      logical :: lcheck
      !
      !-----------------------------------------------------------------------
      !
      !find nearest yhe in opal table
      weight_yhe = 1.d10
      indx_opal = 1
      do i=1, nyhe_opal
         dist = (yhe - yhe_opal(i))**2
         if(dist.lt.weight_yhe) then
            indx_opal=i
            weight_yhe=dist
         endif
      enddo
      !
      write(fname,'(a1,i5.5)') 'Y', int(10000*yhe_opal(indx_opal))
      !
      inquire(file=trim(dir_opal)//'/'//fname, exist=lcheck)
      !
      if(.not.lcheck) then
         write(*,*) 'error in get_opal_table: file "', trim(dir_opal)//'/'//fname, '" does not exist'
         stop
      endif
      !
      write(*,*) '--------------------------reading OPAL tables----------------------------------'
      write(*,*) 'reading from file: ', trim(dir_opal)//'/'//fname
      write(*,*)
      !
      !-----------------------------------------------------------------------
      !
      kappa_opal=zero
      !
      !read in corresponding kappa
      open(1, file=trim(dir_opal)//'/'//fname)
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
