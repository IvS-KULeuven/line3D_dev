module mod_iline

  use prog_type
  use fund_const, only: cgs_clight
  
  implicit none
  
  !index for the line
  integer(i4b) :: iline

  !transition frequency
  real(dp) :: xnue0
  !
  !atomic mass in proton mass
  integer(i4b) :: na
  !
  !upper and lower level statistical weight
  real(dp) :: gl, gu
  !
  !oscillated strength
  real(dp) :: flu
  !
  !gf value (not in log-space)
  real(dp) :: gf
  !
  !atomic number z, ionizationg stage i (in roman numbering), lower level ll, and upper level lu
  integer(i4b) :: element_z, element_i, element_ll, element_lu
  !
  !file name of the lte table if needed to be read in
  character(len=12) :: fname_lte
  character(len=255) :: lte_file_dir
  !
  !line-strength parameters (if required)
  real(dp) :: kline, alpha, kappa0  
  !
contains

  subroutine get_iline(iline)
    !
    !get all required data for a given line described by integer iline
    !note: whenever you change something in here,
    !      you also need to change corresponding entries in
    !         function get_opalbar
    !         subroutines get_photprof, get_photprof_herrero, calc_photprof
    !
    !input: iline - line identifier
    !
    ! ... arguments
    integer(i4b), intent(in) :: iline
    !
    ! ... local scalars
    integer(i4b) :: i, nlines, idum, lu_dum
    real(dp) :: lambda, fdum
    !
    ! ... local characters
    character(len=2) :: str_z, str_i, str_ll
    character :: chdum
    !
    ! ... local logicals
    logical :: lcheck
    !
    !
    !
    write(*,*) '---get_iline: getting xnue0, na, z, i, ll, lu, gl, gu, flu----'
    write(*,'(a10,i20)') 'iline', iline

    select case(iline)
       case(0)
          write(*,*) 'Line from LTE tables'

          !read in line list (only second row)
          write(*,*) "reading input line properties from : in_linelist.dat"
          open(1, file='in_linelist.dat')
             read(1,*)
             read(1,*) element_z, element_i, element_ll, element_lu
             read(1,*) lte_file_dir
          close(1)

          write(str_z,'(i2.2)') element_z
          write(str_i,'(i2.2)') element_i
          write(str_ll,'(i2.2)') element_ll
          fname_lte = str_z//'_'//str_i//'_'//str_ll//'.txt'

          write(*,*) 'reading line data from', TRIM(lte_file_dir)//'/'//fname_lte

          inquire(file=TRIM(lte_file_dir)//'/'//fname_lte, exist=lcheck)

          if(.not.lcheck) then
             write(*,*) 'error in get_iline: file "', TRIM(lte_file_dir)//'/'//fname_lte, '" does not exist'
             stop
          endif


          open(1, file=TRIM(lte_file_dir)//'/'//fname_lte)

             !skip one line
             read(1,*)
          
             !read atomic mass
             read(1,*) fdum, chdum
             na = nint(fdum)
             !
             !number of lines
             read(1,*) nlines, chdum
             !
             !skip two lines
             read(1,*)
             read(1,*)
             !
             !read statistical weight of lower level
             read(1,*) gl, chdum

             !skip two lines
             read(1,*)
             read(1,*)

             !read line data
             do i=1, nlines
                read(1,*) idum, lu_dum, gf, lambda
                if(element_lu.eq.lu_dum) then
                   element_lu = lu_dum
                   exit
                endif
             enddo

             if(element_lu.ne.lu_dum) stop 'error in get_iline: lu not found'

          close(1)

          !since lambda in angstrom
          xnue0 = cgs_clight/lambda*1.d8
          flu = gf/gl

       case(1)
          !hydrogen 2->3
          write(*,*) 'Halpha model'
          na=1
          element_z = 1
          element_i = 1
          element_ll = 2
          element_lu = 3
          gl = 8.d0
          gu = 18.d0
          flu = 6.4108d-1
          gf = gl*flu
          xnue0= 4.5680294d14          
          fname_lte='01_01_02.txt'
       case(2)
          !hydrogen 2->4
          write(*,*) 'Hbeta'
          na=1
          element_z = 1
          element_i = 1
          element_ll = 2
          element_lu = 4
          gl = 8.d0
          gu = 32.d0
          flu = 1.1938d-1
          gf = gl*flu
          xnue0= 6.1668776d14
          fname_lte='01_01_02.txt'
       case(10)
          !CIV (1s2,2s) -> (1s2,2p,spin 3/2)
          write(*,*) 'CIV resonance line'
          na=12
          element_z = 6
          element_i = 4
          element_ll = 1
          element_lu = 2          
          gl = 2.d0
          gu = 4.d0
          flu = 1.9d-1
          gf = gl*flu
          xnue0 = 1.93798d15
         !  fname_lte='06_04_01.txt'
       case(11)
          !CIII (1s2,2s,3p) -> (1s2,2s,3d)
          write(*,*) 'CIII 5696 line'
          na=12
          element_z = 6
          element_i = 3
          element_ll = 9
          element_lu = 12
          gl = 3.d0
          gu = 5.d0
          flu = 3.46d-1
          gf = gl*flu
          xnue0 = 5.2632103d14
          fname_lte='06_03_09.txt'
       case default
          stop 'error in get_iline: iline not properly specified'
    end select
    !
    write(*,'(a20,i20)') 'atomic number', element_z
    write(*,'(a20,i20)') 'ionization stage', element_i
    write(*,'(a20,i20)') 'lower level', element_ll
    write(*,'(a20,i20)')'upper level', element_lu
    write(*,'(a20,es20.8)') 'gl', gl
    write(*,'(a20,es20.8)') 'gu', gu
    write(*,'(a20,es20.8)') 'flu', flu
    write(*,'(a20,es20.8)') 'gf-value', gf
    write(*,'(a20,i20)') 'na', na
    write(*,'(a20,es20.8)') 'xnue0', xnue0
    write(*,'(a20,a20)') 'fname_lte', fname_lte
    write(*,*)
    !
    !
    !
  end subroutine get_iline

end module mod_iline
