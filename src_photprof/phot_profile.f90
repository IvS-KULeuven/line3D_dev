!***********************************************************************
!***********************************************************************
!
!    subroutines to read, and interpolate the photospheric profile
!
!***********************************************************************
!***********************************************************************
!
!***********************************************************************
!***********************************************************************
!
!     old routines (to be debugged for new formal solver versions)
!
!!***********************************************************************
!***********************************************************************
!
subroutine get_photprof(iline, teff, xlogg, yhe, nxobs, nodes_xnue, xic_nue, xicc_nue, xnue0, xic1)
!
!----------------read photospheric profile from external----------------
!--------------interpolate external profile onto own grid---------------
!
!input: teff, xlogg, yhe: stellar parameter to decide which profile is read in
!       nxbos:        dimension of xnue-array
!       nodes_xnue:   xnue-array (own frequency grid)
!       xnue0:        central transition frequency
!       xic1:         b(t_rad) (since external profile is normalized to that value)
!       iline:        number of line to be considered
!
!output: xic_nue:  photospheric profile
!       xicc_nue:  photospheric continuum
!
   use prog_type
   use fund_const, only: cgs_clight
   use photprof_ext, only: fname1, fname2, dirphotprofp_v0
   use mod_interp1d, only: interpol_lin, interpol_yp
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, yhe
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xnue0, xic1
   real(dp), dimension(nxobs), intent(in) :: nodes_xnue
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
   integer(i4b), parameter :: nxobs_ext = 46
   integer(i4b) :: i, err
   real(dp) :: dum_xnue0
   real(dp) :: xnue_ext_max1, xnue_ext_min1, xnue_ext_max2, xnue_ext_min2, xlogg1, xlogg2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: xic_nue1, xic_nue2
   real(dp), dimension(nxobs_ext) :: xnue_ext1, xnue_ext2, xic_ext1, xic_ext2
!
! ... local logicals
   logical :: check1, check2
! ... local functions
!
   write(*,*) '------------getting photospheric profile---------------'
!
!---------------------read photospheric profile-------------------------
!
   call get_filename_photprof(teff, xlogg, yhe, fname1, fname2, xlogg1, xlogg2)
!
!
!check if files exist (remember: maybe other syntax
!   when compiled with gfortran)
   inquire(file=trim(dirphotprofp_v0)//'/'//fname1, exist=check1)
   inquire(file=trim(dirphotprofp_v0)//'/'//fname2, exist=check2)
!
   if(.not.check1) then
      write(*,*) 'error in get_photprof: file "', trim(dirphotprofp_v0)//'/'//fname1, '" does not exist'
      stop
   endif
   if(.not.check2) then
      write(*,*) 'error in get_photprof: file "', trim(dirphotprofp_v0)//'/'//fname2, '" does not exist'
      stop
   endif
!
   open(1, file=dirphotprofp_v0//'/'//fname1, form='formatted')
   open(2, file=dirphotprofp_v0//'/'//fname2, form='formatted')
!
!write(*,*) nxobs_ext
!read in external profile
   do i=1, nxobs_ext
!      write(*,*) i
      read(1, *) xnue_ext1(nxobs_ext+1-i), xic_ext1(nxobs_ext+1-i)
      read(2, *) xnue_ext2(nxobs_ext+1-i), xic_ext2(nxobs_ext+1-i)
   enddo
!
   close(1)
   close(2)
!
   write(*,*) 'used files:'
   write(*,*) '   ', trim(dirphotprofp_v0)//'/'//fname1
   write(*,*) '   ', trim(dirphotprofp_v0)//'/'//fname2
   write(*,*)
!
!------------transform to own units (frequency and profile in cgs)------
!
   xnue_ext1=cgs_clight/xnue_ext1/1.d-8
   xnue_ext2=cgs_clight/xnue_ext2/1.d-8
!
   xic_ext1=xic_ext1*xic1
   xic_ext2=xic_ext2*xic1
!
!----------interpolate external profiles onto own frequency grid--------
!-------since external profiles may have different frequency grids------
!
   allocate(xic_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xic_nue1'
   allocate(xic_nue2(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xic_nue2'
!
!interpolation in logspace
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnue_ext1), xic_nue1, log10(xic_ext1))
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnue_ext2), xic_nue2, log10(xic_ext2))
   xic_nue1=10.d0**xic_nue1
   xic_nue2=10.d0**xic_nue2
!
!-------------------interpolate to input xlogg--------------------------
!
   xnue_ext_max1=maxval(xnue_ext1)
   xnue_ext_max2=maxval(xnue_ext2)
   xnue_ext_min1=minval(xnue_ext1)
   xnue_ext_min2=minval(xnue_ext2)
!
   do i=1, nxobs
      if(xnue_ext_max1.lt.nodes_xnue(i).or.xnue_ext_max2.lt.nodes_xnue(i)) then
         xic_nue(i) = xic1
      else if(xnue_ext_min1.gt.nodes_xnue(i).or.xnue_ext_min2.gt.nodes_xnue(i)) then
         xic_nue(i) = xic1
      else
         xic_nue(i) = interpol_yp(xlogg1, xlogg2, xic_nue1(i), xic_nue2(i), xlogg)
      endif
   enddo
!
!-----------------------------------------------------------------------
!
!constant continuum
   xicc_nue=xic1
!
!--------------------output for debug reasons---------------------------
!
!open(1, file='trash/prof_ext1.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnue_ext1(i), xic_ext1(i)
!   enddo
!close(1)
!!
!open(1, file='trash/prof_ext2.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnue_ext2(i), xic_ext2(i)
!   enddo
!close(1)
!!
!open(1, file='trash/prof_own.dat', form='formatted')
!   do i=1, nxobs
!      write(1,'(4es20.8)') nodes_xnue(i), xic_nue1(i), xic_nue2(i), xic_nue(i)
!   enddo
!close(1)
!stop
!
!-----------------------------------------------------------------------
!
!
!
end subroutine get_photprof
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_filename_photprof(teff, xlogg, yhe, fname1, fname2, g1, g2)
!
   use prog_type
   use mod_sort, only: bubblesort
   use mod_interp1d, only: find_index, interpol_lin
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, yhe
   real(dp) :: g1, g2
   character(len=*) :: fname1, fname2
!
! ... local scalars
   integer(i4b), parameter :: nmod_ext=111, nxobs_ext=46
   integer(i4b) :: i, iim2, iim1, iip1, indx, err
   integer(i4b) :: ntemp_eq, nyhe_eq
   real(dp) :: tphot
   integer(i4b) :: tphot_int, xlogg_int
!
! ... local arrays


!catalogue stores all awailable flux-profiles (#, teff, yhe, logg)
   real(dp), dimension(nmod_ext, 3) :: catalogue
   real(dp), dimension(:), allocatable :: logg_catalogue
!
!--------------set the catalogue of given external profiles-------------
!
   catalogue= reshape((/ 2.00d4, 2.50d4, 2.50d4, 2.50d4, 2.50d4, 2.50d4, 2.50d4, 2.50d4, &
      2.75d4, 2.75d4, 2.75d4, 2.75d4, 2.75d4, 2.75d4, 3.00d4, 3.00d4, &
      3.00d4, 3.00d4, 3.00d4, 3.00d4, 3.00d4, 3.00d4, 3.00d4, 3.00d4, &
      3.00d4, 3.25d4, 3.25d4, 3.25d4, 3.25d4, 3.25d4, 3.25d4, 3.25d4, &
      3.25d4, 3.25d4, 3.25d4, 3.50d4, 3.50d4, 3.50d4, 3.50d4, 3.50d4, &
      3.50d4, 3.50d4, 3.50d4, 3.50d4, 3.50d4, 3.75d4, 3.75d4, 3.75d4, &
      3.75d4, 3.75d4, 3.75d4, 3.75d4, 3.75d4, 3.75d4, 4.00d4, 4.00d4, &
      4.00d4, 4.00d4, 4.00d4, 4.00d4, 4.00d4, 4.00d4, 4.00d4, 4.25d4, &
      4.25d4, 4.25d4, 4.25d4, 4.25d4, 4.25d4, 4.25d4, 4.25d4, 4.50d4, &
      4.50d4, 4.50d4, 4.50d4, 4.50d4, 4.50d4, 4.50d4, 4.50d4, 4.70d4, &
      4.75d4, 4.75d4, 4.75d4, 4.75d4, 4.75d4, 4.75d4, 5.00d4, 5.00d4, &
      5.00d4, 5.00d4, 5.00d4, 5.00d4, 5.25d4, 5.25d4, 5.25d4, 5.25d4, &
      5.25d4, 5.25d4, 5.50d4, 5.50d4, 5.50d4, 5.50d4, 5.75d4, 5.75d4, &
      5.75d4, 5.75d4, 5.75d4, 6.00d4, 6.00d4, 6.00d4, 6.00d4, &
      0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, 0.2d0, 0.2d0, 0.2d0, &
      0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, 0.2d0, &
      0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, &
      0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, &
      0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, &
      0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.1d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, &
      0.1d0, 0.1d0, 0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.1d0, 0.1d0, 0.2d0, 0.2d0, 0.2d0, 0.1d0, 0.1d0, 0.2d0, &
      0.2d0, 0.1d0, 0.1d0, &
      2.50d0, 2.75d0, 3.00d0, 3.50d0, 2.75d0, 3.00d0, 3.50d0, 4.00d0, &
      2.75d0, 3.00d0, 3.50d0, 2.75d0, 3.00d0, 3.50d0, 2.80d0, 3.00d0, &
      3.20d0, 3.50d0, 3.70d0, 4.00d0, 3.00d0, 3.20d0, 3.50d0, 3.70d0, &
      4.00d0, 3.00d0, 3.20d0, 3.50d0, 3.70d0, 4.00d0, 3.00d0, 3.20d0, &
      3.50d0, 3.70d0, 4.00d0, 3.20d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, &
      3.20d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, 3.20d0, 3.50d0, 3.70d0, &
      4.00d0, 4.50d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, 3.30d0, 3.50d0, &
      3.70d0, 4.00d0, 4.50d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, 3.50d0, &
      3.70d0, 4.00d0, 4.50d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, 3.50d0, &
      3.70d0, 4.00d0, 4.50d0, 3.50d0, 3.70d0, 4.00d0, 4.50d0, 3.60d0, &
      3.70d0, 4.00d0, 4.50d0, 3.70d0, 4.00d0, 4.50d0, 3.70d0, 4.00d0, &
      4.50d0, 3.75d0, 4.00d0, 4.50d0, 3.80d0, 4.00d0, 4.20d0, 4.40d0, &
      4.00d0, 4.20d0, 4.00d0, 4.20d0, 4.00d0, 4.20d0, 4.00d0, 4.20d0, &
      4.40d0, 4.20d0, 4.40d0, 4.20d0, 4.40d0, 4.20d0, 4.40d0 /), shape(catalogue))
!
!---------find correct temperatures and correct yhe in catalogue--------
!
   write(*,*) 'looking for photospheric profiles with'
   write(*,'(a10,es20.8)') 'teff', teff
   write(*,'(a10,es20.8)') 'xlogg', xlogg
   write(*,'(a10,es20.8)') 'yhe', yhe
   write(*,*)
!same as in jos routine: take profile at slightly lower value than teff,
!                        because profiles are in steps of 2.5d3 kelvin
   ntemp_eq=0
   nyhe_eq=0
   tphot = teff - mod(teff, 2500.d0)
   if (mod(teff,2500.d0).gt.1250.) tphot = tphot + 2500.d0
!
   do i=1, nmod_ext
      if(tphot.eq.catalogue(i,1)) then
         ntemp_eq=ntemp_eq+1
         if(yhe.eq.catalogue(i,2)) then
            nyhe_eq=nyhe_eq+1
         endif
      endif
   enddo
!
   if(ntemp_eq.eq.0) stop 'error get_filename_photprof: no photospheric profile given for used teff'
   if(nyhe_eq.lt.2) stop 'error get_filename_photprof: less than two photospheric profiles given for yhe, teff'
!needs to be at least two, since interpolation in logg will be performed
!
!------------------------get correct log-g:-----------------------------
!
!two files (which are lower and upper bounds of given logg)
!
   allocate(logg_catalogue(nyhe_eq), stat=err)
   if(err.ne.0) stop 'error get_filename_photprof: logg_catalogue'
!
   indx=1
!
   do i=1, nmod_ext
      if(tphot.eq.catalogue(i,1).and.yhe.eq.catalogue(i,2)) then
         logg_catalogue(indx)=catalogue(i,3)
         indx=indx+1
      endif
   enddo
!
!sort array to find lower and larger logg
   call bubblesort(logg_catalogue, nyhe_eq)
!
   call find_index(xlogg, logg_catalogue, nyhe_eq, iim2, iim1, indx, iip1)
!
   g1=logg_catalogue(indx-1)
   g2=logg_catalogue(indx)
!
!--------------------------store filenames------------------------------
!
   if(yhe.eq.0.1) then
      write(fname1, '(a1, i5, a1, i3)') 't', int(tphot), 'n', int(100.*g1)
      write(fname2, '(a1, i5, a1, i3)') 't', int(tphot), 'n', int(100.*g2)
   else
      write(fname1, '(a1, i5, a1, i3)') 't', int(tphot), 'd', int(100.*g1)
      write(fname2, '(a1, i5, a1, i3)') 't', int(tphot), 'd', int(100.*g2)
   endif
!
!
!
end subroutine get_filename_photprof
!
!***********************************************************************
!***********************************************************************
!
!                 to read in profiles of a. herrero
!
!***********************************************************************
!***********************************************************************
!
subroutine get_photprof_herrero(iline, teff, xlogg, yhe, nxobs, nodes_xnue, xic_nue, xicc_nue, xnue0, trad)
!
!----------------read photospheric profile from external----------------
!--------------interpolate external profile onto own grid---------------
!
!input: teff, xlogg, yhe: stellar parameter to decide which profile is read in
!       nxbos:        dimension of xnue-array
!       nodes_xnue:   xnue-array (own frequency grid)
!       xnue0:        central transition frequency
!       trad:         radiation temperature
!       iline:        index of line transition to be considered
!
!output: xic_nue:  photospheric profile
!        xicc_nue: photospheric continuum
!
   use prog_type
   use fund_const, only: cgs_clight
   use photprof_ext, only: fnamep_iijjkk, fnamec_iijjkk, &
      dirphotprofp_herrero, dirphotprofc_herrero, &
      dirphotprofp_herrero1, dirphotprofc_herrero1, &
      dirphotprofp_herrero2, dirphotprofc_herrero2
   use mod_interp1d, only: interpol_lin
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, yhe
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xnue0, trad
   real(dp), dimension(nxobs), intent(in) :: nodes_xnue
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
   integer(i4b), parameter :: nxobs_ext=118 !herrero profiles have 118 grid points
   integer(i4b) :: i, err
   real(dp) :: dum_xnue0, xic1
   real(dp) :: xnue_ext_max1, xnue_ext_min1, xnue_ext_max2, xnue_ext_min2, xlogg1, xlogg2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: xicc_nue1, xicp_nue1
   real(dp), dimension(nxobs_ext) :: xnuep_ext1, xnuec_ext1, xicp_ext1, xicc_ext1
!
! ... local characters
   character(len=100) :: header
!
! ... local logicals
   logical :: checkc1, checkp1
!
! ... local functions
   real(dp) :: bnue
!
   write(*,*) '------------getting photospheric profile---------------'
!
!choose line from linelist
   select case(iline)
    case(1)
      dirphotprofp_herrero = dirphotprofp_herrero1
      dirphotprofc_herrero = dirphotprofc_herrero1
    case(2)
      dirphotprofp_herrero = dirphotprofp_herrero2
      dirphotprofc_herrero = dirphotprofc_herrero2
    case default
      stop 'error in get_photprof_herrero: iline not properly specified'
   end select
!
!---------------------read photospheric profile-------------------------
!
   call get_filenames_photprof_herrero(fnamep_iijjkk, fnamec_iijjkk, teff, xlogg, yhe)
!
!check if files exist (remember: maybe other syntax
!   when compiled with gfortran)
   inquire(file=trim(dirphotprofp_herrero)//'/'//fnamep_iijjkk, exist=checkp1)
   inquire(file=trim(dirphotprofc_herrero)//'/'//fnamec_iijjkk, exist=checkc1)
!
   if(.not.checkp1) then
      write(*,*) 'error in get_photprof_herrero: file "', trim(dirphotprofp_herrero)//'/'//fnamep_iijjkk, '" does not exist'
      stop
   endif
   if(.not.checkc1) then
      write(*,*) 'error in get_photprof_herrero: file "', trim(dirphotprofc_herrero)//'/'//fnamec_iijjkk, '" does not exist'
      stop
   endif
!
   write(*,*) 'used files:'
   write(*,*) '   ', trim(dirphotprofp_herrero)//'/'//fnamep_iijjkk
   write(*,*) '   ', trim(dirphotprofc_herrero)//'/'//fnamec_iijjkk
   write(*,*)
!
!
   open(1, file=trim(dirphotprofp_herrero)//'/'//fnamep_iijjkk, form='formatted')
   open(2, file=trim(dirphotprofc_herrero)//'/'//fnamec_iijjkk, form='formatted')
!
!read in external profile
   read(1,*) header
   read(1,*) header
   read(2,*) header
   read(2,*) header
   !
!   write(*,*) nxobs_ext
   do i=1, nxobs_ext
!      write(*,*) i
      read(1,*) xnuep_ext1(nxobs_ext+1-i), xicp_ext1(nxobs_ext+1-i)
      read(2,*) xnuec_ext1(nxobs_ext+1-i), xicc_ext1(nxobs_ext+1-i)
   enddo
!
   close(1)
   close(2)
!
!------------transform to own units (frequency and profile in cgs)------
!
   xnuep_ext1 = cgs_clight/xnuep_ext1/1.d-8
   xnuec_ext1 = cgs_clight/xnuec_ext1/1.d-8
!
!normalize continuum flux to xic1 (note: herrero continuum are eddington fluxes => factor of 4)
   xic1 = bnue(xnue0,trad)
   xicc_ext1=xicc_ext1*4.d0/xic1
!
!----------interpolate external profiles onto own frequency grid--------
!-------since external profiles may have different frequency grids------
!
   allocate(xicp_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicp_nue1'
   allocate(xicc_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicc_nue1'
!
!interpolation in logspace
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuec_ext1), xicc_nue1, log10(xicc_ext1))
   xicc_nue1=10.d0**xicc_nue1

   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuep_ext1), xicp_nue1, log10(xicp_ext1))
   xicp_nue1=10.d0**xicp_nue1
!
!-----------------------------todo--------------------------------------
!-----------include different profiles that are interpolated------------
!
!including a wavelength dependent continuum needs to be debugged!!!
   xicc_nue=xic1!1.d0!xicc_nue1
   xic_nue=xicp_nue1*xicc_nue!!*xicc_nue1
!
!
!--------------------output for debug reasons---------------------------
!
!open(1, file='TRASH/prof_ext1.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnuep_ext1(i), xicp_ext1(i)
!   enddo
!close(1)
!!
!!
!open(1, file='TRASH/prof_own.dat', form='formatted')
!   do i=1, nxobs
!      write(1,'(4es20.8)') nodes_xnue(i), xic_nue(i)
!   enddo
!close(1)!
!
!stop
!
!
!
end subroutine get_photprof_herrero
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_filenames_photprof_herrero(fnamep, fnamec, teff, xlogg, yhe)
!
   use prog_type
   use photprof_ext, only: dirphotprofp_herrero, dirphotprofc_herrero
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, yhe
   character(len=6) :: fnamep, fnamec
!
! ... local scalars
   integer(i4b), parameter :: nmod_ext=606, nyhe_ext=6, nteff_ext=17, nxlogg_ext=11
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1, indx, i
   real(dp) :: teff_iim1, teff_ii, xlogg_jjm1, xlogg_jj, yhe_kkm1, yhe_kk, teff2, xlogg2, yhe2
   real(dp) :: w, wi
!
! ... local logicals
!
! ... local arrays
   character(len=4), dimension(nmod_ext) :: catalogue_fname
   real(dp), dimension(nteff_ext) :: catalogue_teff
   real(dp), dimension(nyhe_ext) :: catalogue_yhe
   real(dp), dimension(nxlogg_ext) :: catalogue_xlogg
!
! ... local characters
   character(len=100) :: header

!catalogue stores all awailable flux-profiles (#, teff, yhe, logg)
   real(dp), dimension(nmod_ext, 3) :: catalogue
   real(dp), dimension(:), allocatable :: logg_catalogue
!
!--------------set the catalogue of given external profiles-------------

   catalogue_fname = (/ 't100', 't110', 't120', 't130', 't140', 't150', 't160', 't170', &
      't180', 't190', 't200', 't210', 't220', 't230', 't240', 't250', &
      't260', 't270', 't280', 't290', 't500', 't510', 'taav', 'taaw', &
      'taax', 'tabd', 'tabe', 'tabf', 'tabg', 'tabh', 'tabi', 'tabj', &
      'tabk', 'tabl', 'tabm', 'tabn', 'tabo', 'tabp', 'tabq', 'tabr', &
      'tabs', 'tabv', 'tabw', 'taby', 'taca', 'tacb', 'tacc', 'tacd', &
      'tace', 'tacf', 'tacg', 'taci', 'tacj', 'tack', 'tacl', 'tacm', &
      'tacn', 'taco', 'tacp', 'tacq', 'tacs', 'tact', 'tacu', 'tacv', &
      'tacw', 'tacx', 'tacy', 'tacz', 'tada', 'tadb', 'tadc', 'tade', &
      'tadf', 'tadg', 'tadh', 'tadi', 'tadj', 'tadk', 'tadl', 'tadm', &
      'tadn', 'tado', 'tadp', 'tadq', 'tadr', 'tads', 'tadt', 'tadu', &
      'tadv', 'tadw', 'tadx', 'tady', 'tadz', 'taea', 'taeb', 'taec', &
      'taed', 'taee', 'taeg', 'taei', 'taem', 'taen', 'taeo', 'taep', &
      'taeq', 'taes', 'taet', 'taeu', 'taew', 'taex', 'taey', 'taez', &
      'tafa', 'tafc', 'tafd', 'taff', 'tafg', 'tafi', 'tafj', 'tafk', &
      'tafl', 'tafm', 'tafn', 'tafo', 'tafp', 'tafq', 'tafr', 'tafs', &
      'taft', 'tafu', 'tafv', 'tafw', 'tafy', 'tafz', 'tage', 'tagf', &
      'tagh', 'tagi', 'tagj', 'tagk', 'tagl', 'tagm', 'tagn', 'tagp', &
      'tagq', 'tagr', 'tags', 'tagt', 'tagu', 'tagv', 'tagw', 'tagx', &
      'tagz', 'taha', 'tahb', 'tahc', 'tahd', 'tahe', 'tahf', 'tahg', &
      'tahh', 'tahi', 'tahj', 'tahk', 'tahl', 'tahm', 'tahn', 'taho', &
      'tahp', 'tahq', 'tahr', 'tahs', 'taht', 'tahu', 'tahv', 'tahw', &
      'tahy', 'tahz', 'taia', 'taib', 'taic', 'taid', 'taie', 'taif', &
      'taig', 'taih', 'taij', 'taik', 'tail', 'taim', 'tain', 'taip', &
      'tair', 'tais', 'taiv', 'taiy', 'taiz', 'taja', 'tajb', 'tajc', &
      'taje', 'tajf', 'taji', 'tajj', 'tajl', 'tajm', 'tajo', 'tajp', &
      'tajr', 'tajs', 'tajt', 'taju', 'tajv', 'tajw', 'tajx', 'tajy', &
      'tajz', 'taka', 'takb', 'takc', 'take', 'takf', 'takg', 'taki', &
      'takk', 'takl', 'takm', 'takn', 'takp', 'takq', 'takr', 'taks', &
      'takt', 'taku', 'takv', 'takw', 'takx', 'taky', 'takz', 'tala', &
      'talb', 'talc', 'tald', 'tale', 'talf', 'talh', 'tali', 'talj', &
      'talk', 'talm', 'taln', 'talo', 'talp', 'talq', 'talr', 'tals', &
      'talt', 'talu', 'talv', 'talw', 'talx', 'taly', 'talz', 'tama', &
      'tamb', 'tamc', 'tamd', 'tame', 'tamf', 'tamg', 'tamh', 'tami', &
      'tamj', 'tamk', 'taml', 'tamm', 'tamn', 'tamo', 'tamp', 'tamq', &
      'tamr', 'tams', 'tamt', 'tamu', 'tamw', 'tamx', 'tamy', 'tamz', &
      'tana', 'tanc', 'tand', 'tanf', 'tang', 'tanh', 'tani', 'tanj', &
      'tank', 'tanl', 'tanm', 'tann', 'tano', 'tanp', 'tanq', 'tanr', &
      'tans', 'tant', 'tanu', 'tanv', 'tanw', 'tanx', 'tany', 'tanz', &
      'taoa', 'taob', 'taoc', 'taod', 'taoe', 'taof', 'taog', 'taoh', &
      'taoi', 'taoj', 'taok', 'taol', 'taom', 'taon', 'taoo', 'taoq', &
      'taor', 'taos', 'taov', 'taow', 'taox', 'tapa', 'tapb', 'tapc', &
      'tapd', 'tape', 'tapf', 'tapg', 'taph', 'tapi', 'tapq', 'tapr', &
      'taps', 'tapt', 'tapu', 'tapv', 'tapw', 'tapy', 'tapz', 'taqa', &
      'taqb', 'taqc', 'taqd', 'taqf', 'taqg', 'taqh', 'taqi', 'taqj', &
      'taqk', 'taql', 'taqm', 'taqn', 'taqq', 'taqr', 'taqs', 'taqt', &
      'taqu', 'taqv', 'taqw', 'taqx', 'taqy', 'taqz', 'tara', 'tarb', &
      'tarc', 'tard', 'tare', 'tarf', 'targ', 'tarh', 'tari', 'tarj', &
      'tarl', 'tarn', 'tarp', 'tarr', 'tars', 'tart', 'taru', 'tarv', &
      'tarw', 'tarx', 'tary', 'tarz', 'tasa', 'tasb', 'tasc', 'tasd', &
      'tase', 'tasf', 'tasg', 'tash', 'tasi', 'tasj', 'task', 'tasl', &
      'tasm', 'tasn', 'taso', 'tb11', 'tb22', 'tb33', 'tb44', 'tb55', &
      'tb66', 'tbbe', 'tbbf', 'tbbg', 'tbbh', 'tbbi', 'tbbj', 'tbbl', &
      'tbbm', 'tbbn', 'tbbo', 'tbbp', 'tbbq', 'tbbr', 'tbbs', 'tbbw', &
      'tbbx', 'tbca', 'tbcb', 'tbcc', 'tbcd', 'tbce', 'tbcf', 'tbcg', &
      'tbci', 'tbcj', 'tbck', 'tbcl', 'tbcm', 'tbcn', 'tbcp', &
      'tbcq', 'tbct', 'tbcu', 'tbcv', 'tbcw', 'tbcx', 'tbcy', 'tbcz', &
      'tbda', 'tbdb', 'tbdc', 'tbdd', 'tbde', 'tbdf', 'tbdg', 'tbdh', &
      'tbdi', 'tbdj', 'tbdk', 'tbdl', 'tbdm', 'tbdn', 'tbdo', 'tbdp', &
      'tbdq', 'tbdr', 'tbds', 'tbdt', 'tbdu', 'tbdv', 'tbdw', 'tbdx', &
      'tbdy', 'tbdz', 'tbea', 'tbeb', 'tbec', 'tbed', 'tbee', 'tbeg', &
      'tbei', 'tbej', 'tbem', 'tben', 'tbep', 'tbeq', 'tber', 'tbet', &
      'tbeu', 'tbex', 'tbez', 'tbfa', 'tbfd', 'tbfg', 'tbfj', 'tbga', &
      'tbgb', 'tbgc', 'tbgd', 'tbge', 'tbgf', 'tbgg', 'tbgh', 'tbgi', &
      'tbgj', 'tbgk', 'tbgl', 'tbgm', 'tbgn', 'tbgo', 'tbgq', 'tbgr', &
      'tbgs', 'tbgt', 'tbgw', 'tbgz', 'tbha', 'tbhb', 'tbhe', 'tbhf', &
      'tbhg', 'tbhh', 'tbhi', 'tbhj', 'tbhl', 'tbhm', 'tbho', 'tbhp', &
      'tbhq', 'tbhr', 'tbhs', 'tbht', 'tbhu', 'tbhv', 'tbhw', 'tbhx', &
      'tbhy', 'tbia', 'tbib', 'tbic', 'tbid', 'tbie', 'tbif', 'tbig', &
      'tbih', 'tbii', 'tbij', 'tbik', 'tbil', 'tbim', 'tbin', 'tbio', &
      'tbip', 'tbiq', 'tbir', 'tbis', 'tbit', 'tbiu', 'tbiv', 'tbiw', &
      'tbix', 'tbiy', 'tbiz', 'tbja', 'tbjc', 'tbje', 'tbjf', 'tbji', &
      'tbjj', 'tbjl', 'tbjm', 'tbjn', 'tbjo', 'tbjp', 'tbjq', 'tbjr', &
      'tbjs', 'tbjt', 'tbju', 'tbjv', 'tbjw', 'tbjx', 'tbjy', 'tbjz', &
      'tbka', 'tbkb', 'tbkc', 'tbkd', 'tbke', 'tbkf', 'tcca', 'tccb', &
      'tccc', 'tccf', 'tccg', 'tcch', 'tccj', 'tneu', 'tnne' /)
!
!
!find that profile with lowest difference to input parameters
   w=1.d10
   do i=1, nmod_ext
!
      open(1, file=trim(dirphotprofp_herrero)//trim('/p_'//catalogue_fname(i)))
      read(1,*) header
      read(1,*) teff2, xlogg2, yhe2
      close(1)

      wi = (teff-teff2)**2 + (xlogg-xlogg2)**2 + (yhe-yhe2)**2
      if(wi.lt.w) then
         w=wi
         indx=i
      endif
   enddo
!
   fnamep = 'p_'//catalogue_fname(indx)
   fnamec = 'c_'//catalogue_fname(indx)
!
!
!
end subroutine get_filenames_photprof_herrero
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_photprof_fastwind(iline, teff, xlogg, yhe, nxobs, nodes_xnue, xic_nue, xicc_nue, xnue0, trad)
!
!----------------read photospheric profile from external----------------
!--------------interpolate external profile onto own grid---------------
!
!input: teff, xlogg, yhe: stellar parameter to decide which profile is read in
!       nxbos:        dimension of xnue-array
!       nodes_xnue:   xnue-array (own frequency grid)
!       xnue0:        central transition frequency
!       trad:         radiation temperature
!       iline:        index of line transition to be considered
!
!note: at the moment, only two models present
!
!output: xic_nue:  photospheric profile
!        xicc_nue: photospheric continuum
!
   use prog_type
   use fund_const, only: cgs_clight
   use photprof_ext, only: dirphotprofp_fastwind
   use mod_interp1d, only: interpol_lin
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, yhe
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xnue0, trad
   real(dp), dimension(nxobs), intent(in) :: nodes_xnue
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
   integer(i4b), parameter :: nxobs_ext=161 !fastwind profiles have 161 grid points
   integer(i4b) :: i, err, idum1
   real(dp) :: dum_xnue0, xic1, fdum1, fdum2, fdum3, fdum4, fdum5
   real(dp) :: xnue_ext_max1, xnue_ext_min1, xnue_ext_max2, xnue_ext_min2, xlogg1, xlogg2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: xicc_nue1, xicp_nue1
   real(dp), dimension(nxobs_ext) :: xnuep_ext1, xnuec_ext1, xicp_ext1, xicc_ext1
!
! ... local characters
   character(len=100) :: header
   character(len=100) :: fnamep
!
! ... local logicals
   logical :: checkc1, checkp1
!
! ... local functions
   real(dp) :: bnue
!
   write(*,*) '------------getting photospheric profile---------------'
!
!---------------------read photospheric profile-------------------------
!
   call get_filenames_photprof_fastwind(iline, fnamep, teff, xlogg, yhe)
!
!check if files exist (remember: maybe other syntax
!   when compiled with gfortran)
   inquire(file=trim(dirphotprofp_fastwind)//'/'//trim(fnamep), exist=checkp1)
!
   if(.not.checkp1) then
      write(*,*) 'error in get_photprof_fastwind: file "', trim(dirphotprofp_fastwind)//'/'//trim(fnamep), '" does not exist'
      stop
   endif
!
   write(*,*) 'used files:'
   write(*,*) '   ', trim(dirphotprofp_fastwind)//'/'//trim(fnamep)
   write(*,*)
!
!
   open(1, file=trim(dirphotprofp_fastwind)//'/'//trim(fnamep), form='formatted')
!
!write(*,'(i4,F12.6,F16.2,e20.6,F12.6,F16.6)')  22, 0.761257, 6527.06, 0.569482E-07, 0.995975, 0.995975
!read in external profile
!   write(*,*) nxobs_ext
   do i=1, nxobs_ext
!      write(*,*) i
      read(1,*) idum1, fdum1, fdum2, fdum3, fdum4, fdum5
      xnuep_ext1(nxobs_ext+1-i) = fdum2
      xnuec_ext1(nxobs_ext+1-i) = fdum2
      xicp_ext1(nxobs_ext+1-i) = fdum4
!      write(*,'(2i4,F12.6,F16.2,e20.6,F12.6,F16.6)') i, idum1, fdum1, fdum2, fdum3, fdum4, fdum5 !xnuep_ext1(nxobs_ext+1-i), xicp_ext1(nxobs_ext+1-i)
   enddo
!
   close(1)
!
!------------transform to own units (frequency and profile in cgs)------
!
   xnuep_ext1 = cgs_clight/xnuep_ext1/1.d-8
   xnuec_ext1 = cgs_clight/xnuec_ext1/1.d-8
!
!normalize (or set) continuum flux to xic1
   xic1 = bnue(xnue0,trad)
   xicc_ext1 = xic1/xic1
!
!----------interpolate external profiles onto own frequency grid--------
!-------since external profiles may have different frequency grids------
!
   allocate(xicp_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicp_nue1'
   allocate(xicc_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicc_nue1'
!
!interpolation in logspace
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuec_ext1), xicc_nue1, log10(xicc_ext1))
   xicc_nue1=10.d0**xicc_nue1

   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuep_ext1), xicp_nue1, log10(xicp_ext1))
   xicp_nue1=10.d0**xicp_nue1
!
!-----------------------------todo--------------------------------------
!-----------include different profiles that are interpolated------------
!
!including a wavelength dependent continuum needs to be debugged!!!
   xicc_nue=xic1!1.d0!xicc_nue1
   xic_nue=xicp_nue1*xicc_nue!!*xicc_nue1
!
!
!--------------------output for debug reasons---------------------------
!
!open(1, file='TRASH/prof_ext1.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnuep_ext1(i), xicp_ext1(i)
!   enddo
!close(1)
!!
!!
!open(1, file='TRASH/prof_own.dat', form='formatted')
!   do i=1, nxobs
!      write(1,'(4es20.8)') nodes_xnue(i), xic_nue(i)/xic1
!   enddo
!close(1)
!
!stop
!
!
!
end subroutine get_photprof_fastwind
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_filenames_photprof_fastwind(iline, fnamep, teff, xlogg, yhe)
!
   use prog_type
   use photprof_ext, only: dirphotprofp_fastwind
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: iline
   real(dp), intent(in) :: teff, xlogg, yhe
   character(len=100) :: fnamep
!
! ... local scalars
   integer(i4b), parameter :: nmod_ext=2
   integer(i4b) :: i, indx
   real(dp) :: teff2, xlogg2, yhe2
   real(dp) :: w, wi
!
! ... local logicals
!
! ... local arrays
   character(len=100), dimension(nmod_ext) :: catalogue_fname
!
! ... local characters
   character(len=100) :: header

!catalogue stores all awailable flux-profiles (#)
   real(dp), dimension(nmod_ext) :: catalogue
   real(dp), dimension(nmod_ext):: teff_catalogue, yhe_catalogue, logg_catalogue
!
!--------------set the catalogue of given external profiles-------------
!
   teff_catalogue=(/18.d3, 12.7d3 /)
   yhe_catalogue=(/0.1d0, 0.2d0 /)
   logg_catalogue=(/4.d0, 3.1d0 /)

   select case(iline)
    case(1)
      catalogue_fname = (/ 'lb1_be/OUT.HALPHA_VT005', 'lb1_st/OUT.HALPHA_VT005' /)
    case(2)
      catalogue_fname = (/ 'lb1_be/OUT.HBETA_VT005', 'lb1_st/OUT.HBETA_VT005' /)
    case default
      stop 'error in get_filenames_photprof_fastwind: iline not properly specified'
   end select
!
!
!find that profile with lowest difference to input parameters
   w=1.d10
   do i=1, nmod_ext
!
      teff2 = teff_catalogue(i)
      yhe2 = yhe_catalogue(i)
      xlogg2 = logg_catalogue(i)
      wi = (teff-teff2)**2 + (xlogg-xlogg2)**2 + (yhe-yhe2)**2
      if(wi.lt.w) then
         w=wi
         indx=i
      endif
   enddo
!
   fnamep = catalogue_fname(indx)
!
!
!
end subroutine get_filenames_photprof_fastwind
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_photprof_coelho05(iline, teff, xlogg, fehe, aenh, nxobs, nodes_xnue, xic_nue, xicc_nue, xnue0, trad)
!
!----------------read photospheric profile from external----------------
!--------------interpolate external profile onto own grid---------------
!
!input: teff, xlogg, fehe, aenh: stellar parameter to decide which profile is read in
!       fehe:         metallicity Fe/He
!       aenh:         alpha-element enhancement
!       nxbos:        dimension of xnue-array
!       nodes_xnue:   xnue-array (own frequency grid)
!       xnue0:        central transition frequency
!       trad:         radiation temperature
!       iline:        index of line transition to be considered
!
!note: at the moment, only two models present
!
!output: xic_nue:  photospheric profile
!        xicc_nue: photospheric continuum
!
   use prog_type
   use fund_const, only: cgs_clight
   use photprof_ext, only: dirphotprofp_coelho05
   use mod_interp1d, only: interpol_lin
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, fehe, aenh
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xnue0, trad
   real(dp), dimension(nxobs), intent(in) :: nodes_xnue
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
   integer(i4b) :: i, err, idum1
   real(dp) :: dum_xnue0, xic1, fdum1, fdum2, fdum3, fdum4, fdum5
   real(dp) :: xnue_ext_max1, xnue_ext_min1, xnue_ext_max2, xnue_ext_min2, xlogg1, xlogg2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: xicc_nue1, xicp_nue1
!
!external profile
   integer(i4b) :: nxobs_ext
   real(dp), dimension(:), allocatable :: xnuep_ext1, xnuec_ext1, xicp_ext1, xicc_ext1
!
! ... local characters
   character(len=100) :: header
   character(len=100) :: fnamep
!
! ... local logicals
   logical :: checkc1, checkp1
!
! ... local functions
   real(dp) :: bnue
!
   write(*,*) '------------getting photospheric profile---------------'
!
!---------------------read photospheric profile-------------------------
!
   call get_filenames_photprof_coelho05(iline, fnamep, teff, xlogg, fehe, aenh)
!
!check if files exist (remember: maybe other syntax when compiled with gfortran)
   inquire(file=trim(dirphotprofp_coelho05)//'/'//trim(fnamep), exist=checkp1)
!
   if(.not.checkp1) then
      write(*,*) 'error in get_photprof_coelho05: file "', trim(dirphotprofp_coelho05)//'/'//trim(fnamep), '" does not exist'
      stop
   endif
!
   write(*,*) 'used files:'
   write(*,*) '   ', trim(dirphotprofp_coelho05)//'/'//trim(fnamep)
   write(*,*)
!
!
   open(1, file=trim(dirphotprofp_coelho05)//'/'//trim(fnamep), form='formatted')
!
!write(*,'(i4,F12.6,F16.2,e20.6,F12.6,F16.6)')  22, 0.761257, 6527.06, 0.569482E-07, 0.995975, 0.995975
!read in external profile
!   write(*,*) nxobs_ext
   read(1,*) nxobs_ext
   allocate(xnuep_ext1(nxobs_ext), stat=err)
   allocate(xnuec_ext1(nxobs_ext), stat=err)
   allocate(xicp_ext1(nxobs_ext), stat=err)
   allocate(xicc_ext1(nxobs_ext), stat=err)
   do i=1, nxobs_ext
      read(1,*) fdum1, fdum2
      xnuep_ext1(nxobs_ext+1-i) = fdum1
      xnuec_ext1(nxobs_ext+1-i) = fdum1
      xicp_ext1(nxobs_ext+1-i) = fdum2
!      write(*,'(2i4,F12.6,F16.2,e20.6,F12.6,F16.6)') i, idum1, fdum1, fdum2, fdum3, fdum4, fdum5 !xnuep_ext1(nxobs_ext+1-i), xicp_ext1(nxobs_ext+1-i)
   enddo
!
   close(1)
!
!------------transform to own units (frequency and profile in cgs)------
!
   xnuep_ext1 = cgs_clight/xnuep_ext1/1.d-8
   xnuec_ext1 = cgs_clight/xnuec_ext1/1.d-8
!
!!normalize (or set) continuum flux to xic1
   xic1 = bnue(xnue0,trad)
   xicc_ext1 = xic1/xic1
!
!----------interpolate external profiles onto own frequency grid--------
!-------since external profiles may have different frequency grids------
!
   allocate(xicp_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicp_nue1'
   allocate(xicc_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicc_nue1'
!
!   write(*,*) nodes_xnue
!   write(*,*) xnuec_ext1
!   stop
!interpolation in logspace
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuec_ext1), xicc_nue1, log10(xicc_ext1))
   xicc_nue1=10.d0**xicc_nue1
!
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuep_ext1), xicp_nue1, log10(xicp_ext1))
   xicp_nue1=10.d0**xicp_nue1

!
!-----------------------------todo--------------------------------------
!-----------include different profiles that are interpolated------------
!
!including a wavelength dependent continuum needs to be debugged!!!
   xicc_nue=xic1!1.d0!xicc_nue1
   xic_nue=xicp_nue1*xicc_nue!!*xicc_nue1
!
!--------------------output for debug reasons---------------------------
!
!open(1, file='TRASH/prof_ext1.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnuep_ext1(i), xicp_ext1(i)
!   enddo
!close(1)
!!
!!
!open(1, file='TRASH/prof_own.dat', form='formatted')
!   do i=1, nxobs
!      write(1,'(4es20.8)') nodes_xnue(i), xic_nue(i)/xic1
!   enddo
!close(1)
!
!stop
!
!
!
end subroutine get_photprof_coelho05
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_filenames_photprof_coelho05(iline, fnamep, teff, logg, fehe, aenh)
!
!catalogue for photospheric profiles from coelho05
!on inpu: teff, logg, [Fe/He], alpha-enhancemend
!
   use prog_type
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: iline
   real(dp), intent(in) :: teff, logg, fehe, aenh
   character(len=100) :: fnamep
!
! ... local scalars
   integer(i4b), parameter :: naenh_ext=2, nfehe_ext=8, nteff_ext=15, nlogg_ext=11
   integer(i4b) :: i, j, k, l, indxi, indxj, indxk, indxl

   real(dp) :: w, wi, w_aenh, w_fehe, w_logg, w_teff
   real(dp) :: aenh2, fehe2, logg2, teff2

!integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1, indx, i
!real(dp) :: w, wi
!
! ... local logicals
!
! ... local arrays
   real(dp), dimension(nteff_ext) :: catalogue_teff
   real(dp), dimension(nfehe_ext) :: catalogue_fehe
   real(dp), dimension(naenh_ext) :: catalogue_aenh
   real(dp), dimension(nlogg_ext) :: catalogue_logg
!
! ... local characters
   character(len=100) :: header
   character(len=2) :: string_logg
   character(len=3) :: string_fehe, string_aenh
   character(len=4) :: string_teff
!
!--------------set the catalogue of given external profiles-------------
!
   catalogue_teff = (/ 3500.d0, 3750.d0, 4000.d0, 4250.d0, 4500.d0, 4750.d0, 5000.d0, &
      5250.d0, 5500.d0, 5750.d0, 6000.d0, 6250.d0, 6500.d0, 6750.d0, 7000.d0 /)
   catalogue_logg = (/ 0.d0, 0.5d0, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 3.5d0, 4.d0, 4.5d0, 5.d0 /)
   catalogue_fehe = (/ -2.5d0, -2.d0, -1.5d0, -1.d0, -0.5d0, 0.d0, 0.2d0, 0.5d0 /)
   catalogue_aenh = (/ 0.d0, 0.4d0 /)
!
!
!find files with Teff, logg, Fe/He and alpha enhancemend yielding smallest difference to input parameters
   w=1.d10
   do i=1, nfehe_ext
      w_fehe = (fehe-catalogue_fehe(i))**2
      do j=1, naenh_ext
         w_aenh = (aenh-catalogue_aenh(j))**2
         do k=1, nteff_ext
            w_teff = (teff-catalogue_teff(k))**2
            do l=1, nlogg_ext
               w_logg = (logg-catalogue_logg(l))**2

               wi = w_fehe + w_aenh + w_teff + w_logg
!
               if(wi.lt.w) then
                  w=wi
                  indxi=i
                  indxj=j
                  indxk=k
                  indxl=l
               endif
            enddo
         enddo
      enddo
   enddo
!
   fehe2 = catalogue_fehe(indxi)
   aenh2 = catalogue_aenh(indxj)
   teff2 = catalogue_teff(indxk)
   logg2 = catalogue_logg(indxl)
!
   if(fehe2.lt.0.d0) then
      write(string_fehe,'(a,i2.2)') 'm', int(-fehe2*10.)
   else
      write(string_fehe,'(a,i2.2)') 'p', int(fehe2*10.)
   endif
!
   write(string_aenh,'(a,i2.2)') 'p', int(aenh2*10.)
   write(string_teff, '(i4.4)') int(teff2)
   write(string_logg, '(i2.2)') int(logg2*10.)



   select case(iline)
    case(1)
      fnamep = string_teff//'_'//string_logg//'_'//string_fehe//string_aenh//'_halpha.ms.dat'
    case(2)
      fnamep = string_teff//'_'//string_logg//'_'//string_fehe//string_aenh//'_hbeta.ms.dat'
    case default
      write(*,*) 'error in get_filenames_photprof_coelho05:'
      write(*,*) 'iline', iline, 'not implemented yet'
      stop
   end select
!
!fnamep = 'p_'//catalogue_fname(indx)
!fnamec = 'c_'//catalogue_fname(indx)
!
!
!
end subroutine get_filenames_photprof_coelho05
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_photprof_coelho14(iline, teff, xlogg, fehe, aenh, nxobs, nodes_xnue, xic_nue, xicc_nue, xnue0, trad)
!
!----------------read photospheric profile from external----------------
!--------------interpolate external profile onto own grid---------------
!
!input: teff, xlogg, fehe, aenh: stellar parameter to decide which profile is read in
!       fehe:         metallicity Fe/He
!       aenh:         alpha-element enhancement
!       nxbos:        dimension of xnue-array
!       nodes_xnue:   xnue-array (own frequency grid)
!       xnue0:        central transition frequency
!       trad:         radiation temperature
!       iline:        index of line transition to be considered
!
!note: at the moment, only two models present
!
!output: xic_nue:  photospheric profile
!        xicc_nue: photospheric continuum
!
   use prog_type
   use fund_const, only: cgs_clight
   use photprof_ext, only: dirphotprofp_coelho14
   use mod_interp1d, only: interpol_lin
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: teff, xlogg, fehe, aenh
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xnue0, trad
   real(dp), dimension(nxobs), intent(in) :: nodes_xnue
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
   integer(i4b) :: i, err, idum1
   real(dp) :: dum_xnue0, xic1, fdum1, fdum2, fdum3, fdum4, fdum5
   real(dp) :: xnue_ext_max1, xnue_ext_min1, xnue_ext_max2, xnue_ext_min2, xlogg1, xlogg2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: xicc_nue1, xicp_nue1
!
!external profile
   integer(i4b) :: nxobs_ext
   real(dp), dimension(:), allocatable :: xnuep_ext1, xnuec_ext1, xicp_ext1, xicc_ext1
!
! ... local characters
   character(len=100) :: header
   character(len=100) :: fnamep
!
! ... local logicals
   logical :: checkc1, checkp1
!
! ... local functions
   real(dp) :: bnue
!
   write(*,*) '------------getting photospheric profile---------------'
!
!---------------------read photospheric profile-------------------------
!
   call get_filenames_photprof_coelho14(iline, fnamep, teff, xlogg, fehe, aenh)
!
!check if files exist (remember: maybe other syntax when compiled with gfortran)
   inquire(file=trim(dirphotprofp_coelho14)//'/'//trim(fnamep), exist=checkp1)
!
   if(.not.checkp1) then
      write(*,*) 'error in get_photprof_coelho14: file "', trim(dirphotprofp_coelho14)//'/'//trim(fnamep), '" does not exist'
      stop
   endif
!
   write(*,*) 'used files:'
   write(*,*) '   ', trim(dirphotprofp_coelho14)//'/'//trim(fnamep)
   write(*,*)
!
!
   open(1, file=trim(dirphotprofp_coelho14)//'/'//trim(fnamep), form='formatted')
!
!write(*,'(i4,F12.6,F16.2,e20.6,F12.6,F16.6)')  22, 0.761257, 6527.06, 0.569482E-07, 0.995975, 0.995975
!read in external profile
!   write(*,*) nxobs_ext
   read(1,*) nxobs_ext
   allocate(xnuep_ext1(nxobs_ext), stat=err)
   allocate(xnuec_ext1(nxobs_ext), stat=err)
   allocate(xicp_ext1(nxobs_ext), stat=err)
   allocate(xicc_ext1(nxobs_ext), stat=err)
   do i=1, nxobs_ext
      read(1,*) fdum1, fdum2
      xnuep_ext1(nxobs_ext+1-i) = fdum1
      xnuec_ext1(nxobs_ext+1-i) = fdum1
      xicp_ext1(nxobs_ext+1-i) = fdum2
!      write(*,'(2i4,F12.6,F16.2,e20.6,F12.6,F16.6)') i, idum1, fdum1, fdum2, fdum3, fdum4, fdum5 !xnuep_ext1(nxobs_ext+1-i), xicp_ext1(nxobs_ext+1-i)
   enddo
!
   close(1)
!
!------------transform to own units (frequency and profile in cgs)------
!
   xnuep_ext1 = cgs_clight/xnuep_ext1/1.d-8
   xnuec_ext1 = cgs_clight/xnuec_ext1/1.d-8
!
!!normalize (or set) continuum flux to xic1
   xic1 = bnue(xnue0,trad)
   xicc_ext1 = xic1/xic1
!
!----------interpolate external profiles onto own frequency grid--------
!-------since external profiles may have different frequency grids------
!
   allocate(xicp_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicp_nue1'
   allocate(xicc_nue1(nxobs), stat=err)
   if(err.ne.0) stop 'allocation error get_photprof: xicc_nue1'
!
!   write(*,*) nodes_xnue
!   write(*,*) xnuec_ext1
!   stop
!interpolation in logspace
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuec_ext1), xicc_nue1, log10(xicc_ext1))
   xicc_nue1=10.d0**xicc_nue1
!
   call interpol_lin(nxobs, nxobs_ext, log10(nodes_xnue), log10(xnuep_ext1), xicp_nue1, log10(xicp_ext1))
   xicp_nue1=10.d0**xicp_nue1

!
!-----------------------------todo--------------------------------------
!-----------include different profiles that are interpolated------------
!
!including a wavelength dependent continuum needs to be debugged!!!
   xicc_nue=xic1!1.d0!xicc_nue1
   xic_nue=xicp_nue1*xicc_nue!!*xicc_nue1
!
!--------------------output for debug reasons---------------------------
!
!open(1, file='TRASH/prof_ext1.dat', form='formatted')
!   do i=1, nxobs_ext
!      write(1,*) xnuep_ext1(i), xicp_ext1(i)
!   enddo
!close(1)
!!
!!
!open(1, file='TRASH/prof_own.dat', form='formatted')
!   do i=1, nxobs
!      write(1,'(4es20.8)') nodes_xnue(i), xic_nue(i)/xic1
!   enddo
!close(1)
!
!
!
!
end subroutine get_photprof_coelho14
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_filenames_photprof_coelho14(iline, fnamep, teff, logg, fehe, aenh)
!
!catalogue for photospheric profiles from coelho14
!on inpu: teff, logg, [Fe/He], alpha-enhancemend
!
   use prog_type
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: iline
   real(dp), intent(in) :: teff, logg, fehe, aenh
   character(len=100) :: fnamep
!
! ... local scalars
   integer(i4b), parameter :: naenh_ext=2, nfehe_ext=4, nteff_ext=55, nlogg_ext=12
   integer(i4b) :: i, j, k, l, indxi, indxj, indxk, indxl

   real(dp) :: w, wi, w_aenh, w_fehe, w_logg, w_teff
   real(dp) :: aenh2, fehe2, logg2, teff2

!integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1, indx, i
!real(dp) :: w, wi
!
! ... local logicals
!
! ... local arrays
   real(dp), dimension(nteff_ext) :: catalogue_teff
   real(dp), dimension(nfehe_ext) :: catalogue_fehe
   real(dp), dimension(naenh_ext) :: catalogue_aenh
   real(dp), dimension(nlogg_ext) :: catalogue_logg
!
! ... local characters
   character(len=100) :: header
   character(len=5) :: string_logg
   character(len=3) :: string_fehe, string_aenh
   character(len=6) :: string_teff
!
!--------------set the catalogue of given external profiles-------------
!
   catalogue_teff = (/ 3000.d0, 3200.d0, 3400.d0, 3600.d0, 3800.d0, 4000.d0, 4250.d0, 4500.d0, 4750.d0, &
      5000.d0, 5250.d0, 5500.d0, 5750.d0, 6000.d0, 6250.d0, 6500.d0, 6750.d0, &
      8000.d0, 8250.d0, 8500.d0, 8750.d0, 9000.d0, 9250.d0, 9500.d0, 9750.d0, &
      1000.d0, 8250.d0, 8500.d0, 8750.d0, 9000.d0, 9250.d0, 9500.d0, 9750.d0, &
      10000.d0, 10250.d0, 10500.d0, 10750.d0, 11000.d0, 11250.d0, 11500.d0, 11750.d0, &
      12000.d0, 13000.d0, 14000.d0, 15000.d0, 16000.d0, 17000.d0, 18000.d0, 19000.d0, &
      20000.d0, 21000.d0, 22000.d0, 23000.d0, 24000.d0, 25000.d0 /)
   catalogue_logg = (/ -0.5d0, 0.d0, 0.5d0, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 3.5d0, 4.d0, 4.5d0, 5.d0 /)
   catalogue_fehe = (/ -1.d0, -0.5d0, 0.d0, 0.2d0 /)
   catalogue_aenh = (/ 0.d0, 0.4d0 /)
!
!
!find files with Teff, logg, Fe/He and alpha enhancemend yielding smallest difference to input parameters
   w=1.d10
   do i=1, nfehe_ext
      w_fehe = (fehe-catalogue_fehe(i))**2
      do j=1, naenh_ext
         w_aenh = (aenh-catalogue_aenh(j))**2
         do k=1, nteff_ext
            w_teff = (teff-catalogue_teff(k))**2
            do l=1, nlogg_ext
               w_logg = (logg-catalogue_logg(l))**2

               wi = w_fehe + w_aenh + w_teff + w_logg
!
               if(wi.lt.w) then
                  w=wi
                  indxi=i
                  indxj=j
                  indxk=k
                  indxl=l
               endif
            enddo
         enddo
      enddo
   enddo
!
   fehe2 = catalogue_fehe(indxi)
   aenh2 = catalogue_aenh(indxj)
   teff2 = catalogue_teff(indxk)
   logg2 = catalogue_logg(indxl)
!
   if(fehe2.lt.0.d0) then
      write(string_fehe,'(a,i2.2)') 'm', int(-fehe2*10.)
   else
      write(string_fehe,'(a,i2.2)') 'p', int(fehe2*10.)
   endif
!
   write(string_aenh,'(a,i2.2)') 'p', int(aenh2*10.)
   write(string_teff, '(a,i5.5)') 't', int(teff2)
!
   if(logg2.lt.0.d0) then
      write(string_logg,'(a,f3.1)') 'g-', -logg2
   else
      write(string_logg,'(a,f3.1)') 'g+', logg2
   endif



   select case(iline)
    case(1)
      fnamep = string_teff//'_'//string_logg//'_'//string_fehe//string_aenh//'_halpha.dat'
    case(2)
      fnamep = string_teff//'_'//string_logg//'_'//string_fehe//string_aenh//'_hbeta.dat'
    case default
      write(*,*) 'error in get_filenames_photprof_coelho14:'
      write(*,*) 'iline', iline, 'not implemented yet'
      stop
   end select
!
!fnamep = 'p_'//catalogue_fname(indx)
!fnamec = 'c_'//catalogue_fname(indx)
!
!
end subroutine get_filenames_photprof_coelho14

!***********************************************************************
!***********************************************************************
!
!                       general routines
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_photprof(iline, nxobs, xic_nue, xicc_nue, xic1)
!
!----------set photospheric profile to constant value: xic1-------------
!----------set photospheric continuum to constant value: xic1-----------
!
   use prog_type
   use fund_const, only: cgs_clight
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: nxobs, iline
   real(dp), intent(in) :: xic1
   real(dp), dimension(nxobs) :: xic_nue, xicc_nue
!
! ... local scalars
!
! ... local functions
!

   write(*,*) '----------------------no photospheric profile is included----------------------'
   write(*,*)
!
   xic_nue=xic1
   xicc_nue=xic1
!
end subroutine calc_photprof
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_photprof(outputdir, fname, nxobs, xobs, photprof, xicc_nue)
!
   use prog_type
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: nxobs
   real(dp), dimension(nxobs), intent(in) :: xobs, photprof, xicc_nue
   character(len=*), intent(in) :: outputdir, fname
!
! ... local scalars
   integer(i4b) :: i
!
!photospheric profile
   open(1, file=trim(outputdir)//'/'//trim(fname), form='formatted')
   write(1, '(i5)') nxobs
   do i=1, nxobs
      write(1,'(3(es20.8))') xobs(i), photprof(i), xicc_nue(i)
   enddo
   close(1)
!
end subroutine output_photprof
