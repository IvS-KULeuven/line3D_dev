module mod_lte
!
!module to read in LTE occupation numbers from Luka's code
!
   use prog_type
   use fund_const
   use mod_iline, only: element_lu, element_ll, dir_lte
!
   implicit none
!
   ! integer(i4b), parameter :: nyhe_lte=8
   integer(i4b) :: nrho_lte, ntemp_lte, nlines_lte
   real(dp) :: glower
!
   ! real(dp), dimension(nyhe_lte), parameter :: yhe_lte=(/0.00999d0, 0.1d0, 0.28d0, 0.59d0, 0.61d0, 0.98d0, 0.99d0, 1.d0 /)
   real(dp), dimension(:), allocatable :: rho_lte
   real(dp), dimension(:), allocatable :: temp_lte
   real(dp), dimension(:,:), allocatable :: nlower_lte
!
contains

   subroutine get_lte_table()
      !
      !
      !yhe_mass is mass fraction  m_he/m_tot
      !
      !
      !... arguments
      ! real(dp), intent(in) :: yhe_mass
      !
      ! ... local scalars
      integer(i4b) :: i, idum, indx, err, lu, lu_dum, lu_dum2
      real(dp) :: fdum, gf, lambda, gf_dum, lambda_dum
      ! real(dp) :: weight_yhe, dist, 
      !
      ! ... local characters
      ! character(len=6) :: fname
      character :: chdum
      !
      ! ... local logicals
      logical :: lcheck
      !
      !-----------------------------------------------------------------------
      !
      inquire(file=TRIM(dir_lte), exist=lcheck)
      !
      if(.not.lcheck) then
         write(*,*) 'error in get_lte_table: file "', TRIM(dir_lte), '" does not exist'
         stop
      endif
      !
      !-----------------------------------------------------------------------
      !
      write(*,*) '---------------------------reading LTE tables----------------------------------'
      write(*,*) 'reading from file: ', TRIM(dir_lte)
      write(*,*)
      !
      !read in corresponding nlower
      open(1, file=TRIM(dir_lte))

      !skip first lines
      do i=1, 2
         read(1,*)
      enddo

      !read nrho and ntemp and nlines
      read(1,*) nlines_lte, chdum
      read(1,*) ntemp_lte, chdum
      read(1,*) nrho_lte, chdum

      allocate(nlower_lte(ntemp_lte,nrho_lte), stat=err)
      allocate(rho_lte(nrho_lte), stat=err)
      allocate(temp_lte(ntemp_lte), stat=err)

      nlower_lte=zero
      rho_lte=zero
      temp_lte=zero

      read(1,*) glower, chdum

      !skip 2 lines
      read(1,*)
      read(1,*)

      !read line data
      do i=1, nlines_lte
         read(1,*) idum, lu_dum, gf_dum, lambda_dum
         if(element_lu.eq.lu_dum) then
            gf=gf_dum
            lambda=lambda_dum
            lu_dum2=lu_dum
         endif
      enddo

      if(element_lu.ne.lu_dum2) stop 'error in get_lte_tables: lu not found'

      !skip three lines
      do i=1, 3
         read(1,*)
      enddo

      !read temperature
      read(1,*) fdum, temp_lte
      !
      !read rhos and nlower
      do i=1, nrho_lte
         !       read(1,'(f4.2,19f7.3)') rho_lte(i), nlower_lte(:,i)
         read(1,*) rho_lte(i), nlower_lte(:,i)
      enddo

      close(1)
      !
      !write(*,*) yhe, weight_yhe, indx_lte, yhe_lte(indx_lte), fname
      !
      !write(*,*) rho_lte
      !write(*,*) temp_lte
      !
      !do i=1, ntemp_lte
      !   write(*,'(19f10.3)') 10.d0**nlower_lte(i,:)
      !enddo
      !stop 'go on in get_lte_table'

   end subroutine get_lte_table

!
end module mod_lte
