module mod_iline
   use prog_type
   use fund_const, only: cgs_clight
   use mod_mforce

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

   integer(i4b), parameter :: nyhe_lteopal=1
   real(dp), dimension(nyhe_lteopal), parameter :: yhe_lteopal=(/0.00999d0/) !, 0.1d0, 0.28d0, 0.59d0, 0.61d0, 0.98d0, 0.99d0, 1.d0 /)
   !
   !file name of the lte table if needed to be read in
   character(len=255) :: dir_lte
   character(len=255) :: dir_opal
   !
   !line-strength parameters (if required)
   real(dp) :: kline, alpha, kappa0  
   !
   contains

      subroutine get_iline(iline, yhe_mass, X_mass, Y_mass)
         !
         !get all required data for a given line described by integer iline
         !note: whenever you change something in here,
         !      you also need to change corresponding entries in
         !         function get_opalbar
         !         subroutines get_photprof, get_photprof_herrero, calc_photprof
         !
         !input: iline - line identifier
         !       yhe_mass - helium mass fraction 
         !
         ! ... arguments
         integer(i4b), intent(in) :: iline
         real(dp), intent(in) :: yhe_mass
         real(dp), optional, intent(in) :: X_mass, Y_mass
         !
         ! ... local scalars
         integer(i4b) :: i, nlines, idum, lu_dum
         real(dp) :: lambda, fdum
         real(dp) :: weight_yhe, dist
         integer(i4b) :: indx_lteopal
         !
         ! ... local characters
         character(len=2) :: str_z, str_i, str_ll
         character :: chdum
         character(len=20) :: yname = ' '
         !
         ! ... local logicals
         logical :: lcheck
         !
         !
         !
         write(*,*) '---get_iline: getting xnue0, na, z, i, ll, lu, gl, gu, flu----'
         write(*,'(a10,i20)') 'iline', iline

         !find nearest yhe in LTE OPAL table
         weight_yhe = 1.d10
         indx_lteopal = 1
         do i=1, nyhe_lteopal
            dist = (yhe_mass - yhe_lteopal(i))**2
            if(dist.lt.weight_yhe) then
               indx_lteopal=i
               weight_yhe=dist
            endif
         enddo
         !
         write(yname,'(a1,i5.5)') 'Y', int(100000*yhe_lteopal(indx_lteopal))
         !
         !Get the destination ot the data location (check if the enviroment was set)
         CALL get_environment_variable("DIR_OPAL", dir_opal)
         IF (LEN_TRIM(dir_opal) .EQ. 0) STOP '---- Error: DIR_OPAL not set -----'
         dir_opal = TRIM(dir_opal)//'/'//TRIM(yname)

         select case(iline)
            !!!!!!!!!!!!!!!!!!!!!!!!!
            case(0)
               write(*,*) 'Line from LTE tables'

               ! Get the destination ot the data location (check if the enviroment was set)
               CALL get_environment_variable("DIR_LTE", dir_lte)
               IF (LEN_TRIM(dir_lte) .EQ. 0) STOP '---- Error: DIR_LTE not set -----'
               dir_lte = TRIM(dir_lte)//'/'//TRIM(yname)

               inquire(file="./in_line_lte.dat", exist=lcheck)
               !
               if(.not.lcheck) then
                  write(*,*) 'error in line: file ./in_line.dat does not exist'
                  stop
               endif

               !read in line (only second row)
               write(*,*) "reading input line properties from : in_line.dat"
               open(1, file='./in_line.dat')
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*) element_z, element_i, element_ll, element_lu
               close(1)
               !
               write(str_z,'(i2.2)') element_z
               write(str_i,'(i2.2)') element_i
               write(str_ll,'(i2.2)') element_ll
               dir_lte = TRIM(dir_lte)//'/'//TRIM(str_z//'_'//str_i//'_'//str_ll//'.txt')

               write(*,*) 'reading line data from', TRIM(dir_lte)

               inquire(file=TRIM(dir_lte), exist=lcheck)

               if(.not.lcheck) then
                  write(*,*) 'error in get_iline: file "', TRIM(dir_lte), '" does not exist'
                  stop
               endif


               open(1, file=TRIM(dir_lte))

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
            !!!!!!!!!!!!!!!!!!!!!!!!!
            case(1)
               inquire(file="./in_line.dat", exist=lcheck)
               !
               if(.not.lcheck) then
                  write(*,*) 'error in line: file ./in_line.dat does not exist'
                  stop
               endif

               write(*,*) "reading input line properties from : in_line.dat"
               open(1, file='./in_line.dat')
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*) element_z, element_i, element_ll, element_lu
                  read(1,*) gl, gu, flu, xnue0 , na
               close(1)       
               gf = gl*flu
            !!!!!!!!!!!!!!!!!!!!!!!!!
            case(111)
               inquire(file="./in_line.dat", exist=lcheck)
               !
               if(.not.lcheck) then
                  write(*,*) 'error in line: file ./in_line.dat does not exist'
                  stop
               endif

               !read in line (only second row)
               write(*,*) "reading input line properties from : in_line.dat"
               open(1, file='./in_line.dat')
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*)
                  read(1,*) element_z, element_i, element_ll, element_lu
               close(1)
            !
               write(str_z,'(i2.2)') element_z
               write(str_i,'(i2.2)') element_i
               write(str_ll,'(i2.2)') element_ll   
               write(*,*) 'using MFORCE for Z ',str_z, '  I ', str_i

               call Initialise_NLTE(X_mass, Y_mass)
               call Get_NLTE_line(element_z, element_i, element_ll, element_lu, gl, gu, gf, lambda)
               ! write(*,*) element_z, element_i, element_ll, element_lu, gl, gu, gf, lambda

               flu =  gf/gl
               xnue0 = cgs_clight/lambda*1.d8 
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
         write(*,*)
         !
         !
         !
      end subroutine get_iline

end module mod_iline
