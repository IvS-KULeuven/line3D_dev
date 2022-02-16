!
!-----------------------------------------------------------------------
!
!   calculates or reads 1d model atmospheres that can be read in
!               spec.eo to calculate line profiles
!
!-----------------------------------------------------------------------
!
program modspec
!
use prog_type
use options_modspec, only: input_mod
!
implicit none
!
! ... local scalars
!
call read_input
!
select case(input_mod)
   case(0)
      call read_model1da
      call output1d
   case(1,2)
      call calc_model1d_sobolev
      call output1d
   case(3,4,5)
      call read_model1db
      call output1d
   case(6)
      call calc_model1d_petrenz
      call output1d
   case(7,8,9)
      call read_model3db
      call output3d
   case(10)
      call calc_model2d_petrenz
      call output3d_spc
   case(11)
      call read_model3d
      call output3d
   case(12)
      call calc_model3d_spc3d
      call output3d_spc
   case(13)
      call calc_model3d_spc3dc
      call output3d_spc
   case(14)
      call calc_model3d_spc3db
      call output3d_spc
   case(15)
      call calc_model1d_luka
      call output1d
   case(16)
      call calc_model3d_be
      call output3d_spc
   case(17)
      call calc_model3d_ileyk
      call output3d_spc
   case(18)
      call calc_model3d_windlte
      call output3d_spc
   case(19)
      call calc_model3d_nicowr3d
      call output3d_spc
   case(20)
      call calc_test_model_js_lh
      call output3d_spc
   case default
      stop 'error in modspec: input_mod not specified'
end select
!

end program modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!

use prog_type
use fund_const, only: rsu
use options_modspec, only: indat_file, input_file, input_file2, output_file, input_mod
use params_modspec, only: teff, trad, tmin, xlogg, rstar, rmax, tmin, xmloss, vmin, &
                          vmax, beta, yhe, yhe_mass, hei, vrot, vmicro, vth_fiducial, &
                          sr, rmin, eps_line, unit_length, xic2
use fund_const, only: rsu, pi, cgs_grav, xmsu
use mod_opal
use mod_lte
use mod_iline, only: iline, get_iline, na, kline, alpha, kappa0, fname_lte
!
!
implicit none
!
! ... local scalars
real(dp) :: fdum
!
! ... local characters
!
! ... local functions
!
! ... namelist
namelist / input_options / input_file, input_file2, output_file, input_mod
namelist / input_model / teff, trad, xlogg, rstar, rmax, tmin, xmloss, vmin, vmax, vrot, &
                         vmicro, vth_fiducial, beta, yhe, hei
namelist / input_line / iline, eps_line, kline, alpha, kappa0
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) indat_file
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')

!read options, i.e. which model shall be used
   rewind 1
   read(1, nml=input_options)

!read 1d model parameters
   rewind 1
   read(1, nml=input_model)
!
!read line strength parameters etc
   rewind 1
   read(1, nml=input_line)
!
close(1)
!

call get_iline(iline)
!
sr=rstar*rsu
unit_length=sr
tmin=tmin*teff
rmin=1.d0
!
!velocities in cm
vth_fiducial=vth_fiducial*1.d5
vmicro=vmicro*1.d5
vmin=vmin*1.d5
vmax=vmax*1.d5
vrot=vrot*1.d5 !in cm/s
!
!calculate mass fraction where yhe = n_he/n_h = yhe/(1-yhe), where yhe=n_he/(n_h+n_he)
yhe_mass = one / (one + one/four/yhe)
!
!opal tables
call get_opal_table(yhe_mass)
!
!lte tables
if(iline.eq.0) then
   call get_lte_table(yhe_mass)
endif
!
!
!default xic2
xic2=zero
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model1da
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu
use options_modspec, only: input_file
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d, sline1d, scont1d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax
use mod_iline, only: xnue0
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: lambda, fdum
!
! ... local characters
character(len=500) :: cdum
!
! ... local logcials
logical :: ldum
!
!
write(*,*) '--------------------reading input from JOs program-----------------------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
open(1, file=trim(input_file), form='formatted')
   read(1,*) cdum
   read(1,'(i5, 6es20.8, i5, 12es20.8)') n1d, fdum, fdum, fdum, fdum, fdum, fdum, &
                                         ldum, fdum, teff, xic1, lambda, vth_fiducial, fdum, fdum, &
                                         vmax, fdum, fdum, fdum, sr
!
   xnue0 = cgs_clight/lambda/1.d-8
   trad = log(2.d0*cgs_planck*xnue0**3/cgs_clight**2/xic1 + 1.d0)
   trad = cgs_planck*xnue0/trad/cgs_kb
   rstar = sr/rsu
!
   allocate(r1d(n1d))
   allocate(velr1d(n1d))
   allocate(opac1d(n1d))
   allocate(opalbar1d(n1d))
   allocate(t1d(n1d))
   allocate(vth1d(n1d))
   allocate(sline1d(n1d))
   allocate(scont1d(n1d))
!
   read(1,*) cdum
!
   do i=1, n1d
      read(1,'(11es20.8)') r1d(i), t1d(i), fdum, fdum, opac1d(i), opalbar1d(i), velr1d(i), &
                            sline1d(i), fdum, fdum, scont1d(i)
      vth1d(i) = vth_fiducial
   enddo
close(1)
!
end subroutine read_model1da
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model1db
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu
use options_modspec, only: input_file, input_mod
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d, sline1d, scont1d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro
use mod_iline, only: xnue0, na
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: lambda, fdum
!
! ... local characters
character(len=500) :: cdum
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r, dims_iter
!
!
write(*,*) '----------------reading input from benchmark08: sc or fvm solution-------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)

         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         
         vmax=vmax*1.d5

         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, n1d, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      dims_r = (/ n1d /)
      allocate(r1d(n1d))
      allocate(velr1d(n1d))
      allocate(opac1d(n1d))
      allocate(opalbar1d(n1d))
      allocate(t1d(n1d))
      allocate(vth1d(n1d))
      allocate(sline1d(n1d))
      allocate(scont1d(n1d))

      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r1d, dims_r, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model_cr', group_id, err)
         call h5dopen_f(group_id, 't1d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t1d, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar1d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar1d, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velr1d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr1d, dims_r, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'solution_cr', group_id, err)
         select case(input_mod)
            case(3)
               write(*,*) 'reading sc solution'
               write(*,*)
               call h5dopen_f(group_id, 'sline_sc', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline1d, dims_r, err)
               call h5dclose_f(dset_id, err)
            case(4)
               write(*,*) 'reading fvm solution'
               write(*,*)
               call h5dopen_f(group_id, 'sline_fvm', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline1d, dims_r, err)
               call h5dclose_f(dset_id, err)
            case(5)
               write(*,*) 'reading JOs solution (on own grid)'
               write(*,*)
               call h5dopen_f(group_id, 'sline_jo', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline1d, dims_r, err)
               call h5dclose_f(dset_id, err)
            case default
               stop 'error in read_model1db: input_mod not properly specified'
         end select
!
end subroutine read_model1db         
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3db
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu
use options_modspec, only: input_file, input_mod
use dime_modspec, only: ndxmax, ndymax, ndzmax, n1d, r1d, sline1d, x, y, z, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line
use hdf5
use mod_interp1d, only: find_index, interpol_yp
use mod_iline, only: alpha, kappa0, kline, na, xnue0
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2, iim1, ii, iip1
real(dp) :: lambda, fdum, rad
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_r
integer(hsize_t), dimension(3) :: dims_3d
!
!
write(*,*) '----------------reading input from benchmark13: sc or fvm solution-------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
 
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'opt_opal', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kline, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'alpha', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kappa0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'eps_line', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)         
         vmax=vmax*1.d5
         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, n1d, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_r = (/ n1d /)
         call h5aopen_f(group_id, 'ndxmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_x = (/ ndxmax /)
         call h5aopen_f(group_id, 'ndymax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_y = (/ ndymax /)
         call h5aopen_f(group_id, 'ndzmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_z = (/ ndzmax /)
      call h5gclose_f(group_id, err)

      dims_3d = (/ ndxmax, ndymax, ndzmax /)

      allocate(x(ndxmax), stat=err)
      allocate(y(ndymax), stat=err)
      allocate(z(ndzmax), stat=err)
      allocate(imask3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(t3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opac3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opalbar3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velx3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(vely3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velz3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(sline3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(scont3d(ndxmax,ndymax,ndzmax), stat=err)

      allocate(r1d(n1d), stat=err)
      allocate(sline1d(n1d), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r1d, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 'mask3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_integer, imask3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opac3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velx3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)

      scont3d=0.d0
      select case(input_mod)
         case(7)
            write(*,*) 'reading sc solution'
            write(*,*)
            call h5gopen_f(file_id, 'solution3d', group_id, err)
               call h5dopen_f(group_id, 'sline3d_sc', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5gclose_f(group_id, err)
         case(8)
            write(*,*) 'reading fvm solution'
            write(*,*)
            call h5gopen_f(file_id, 'solution3d', group_id, err)
               call h5dopen_f(group_id, 'sline3d_fvm', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5gclose_f(group_id, err)
         case(9)
            write(*,*) 'reading jos 1d solution and interpolation'
            write(*,*)
            call h5gopen_f(file_id, 'solution_cr', group_id, err)
               call h5dopen_f(group_id, 'sline1d_jo', dset_id, err)
                  call h5dread_f(dset_id, h5t_native_double, sline1d, dims_r, err)
               call h5dclose_f(dset_id, err)
            call h5gclose_f(group_id, err)
            do i=1, ndxmax
               do j=1, ndymax
                  do k=1, ndzmax
                     select case(imask3d(i,j,k))
                        case(1,2,3)
                           rad=sqrt(x(i)**2+y(j)**2+z(k)**2)
                           call find_index(rad, r1d, n1d, iim2, iim1, ii, iip1)
                           sline3d(i,j,k) = interpol_yp(r1d(iim1), r1d(ii), sline1d(iim1), sline1d(ii), rad)
                        case default
                           sline3d(i,j,k) = 0.d0
                     end select
                  enddo
               enddo
            enddo
         case default
            stop 'error in read_model3db: input_mod not properly specified'
         end select
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!
end subroutine read_model3db
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3d
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu
use options_modspec, only: input_file, input_mod
use dime_modspec, only: ndxmax, ndymax, ndzmax, x, y, z, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot
use hdf5
use mod_iline, only: alpha, kappa0, kline, xnue0, na
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2, iim1, ii, iip1
real(dp) :: lambda, fdum, rad
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
real(dp) :: interpol_yp
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
integer(hsize_t), dimension(3) :: dims_3d
!
!
write(*,*) '----------------reading standard 3d input: cartesian coordinates---------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
 
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'options', group_id, err)
         call h5aopen_f(group_id, 'opt_opal', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kline, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'eps_line', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)         
         vmax=vmax*1.d5
         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         sr=rstar*rsu
         call h5aopen_f(group_id, 'lstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xlogg', attr_id, err)
            if(err.eq.-1) then
               write(*,*) 'error reading h5 file: use previously read xlogg'
             else
               call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            endif
         call h5aopen_f(group_id, 'vrot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         vrot=vrot*1.d5 !in cm/s
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'bcondition', group_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'ndxmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_x = (/ ndxmax /)
         call h5aopen_f(group_id, 'ndymax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_y = (/ ndymax /)
         call h5aopen_f(group_id, 'ndzmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_z = (/ ndzmax /)
      call h5gclose_f(group_id, err)

      dims_3d = (/ ndxmax, ndymax, ndzmax /)

      allocate(x(ndxmax), stat=err)
      allocate(y(ndymax), stat=err)
      allocate(z(ndzmax), stat=err)
      allocate(imask3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(t3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opac3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opalbar3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velx3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(vely3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velz3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(sline3d(ndxmax,ndymax,ndzmax), stat=err)
      allocate(scont3d(ndxmax,ndymax,ndzmax), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 'mask3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_integer, imask3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opac3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
!convert opalbar and opac to 1/cm
         opac3d=opac3d/sr
         opalbar3d=opalbar3d/sr
         call h5dopen_f(group_id, 'velx3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)

      scont3d=0.d0
      call h5gopen_f(file_id, 'solution3d', group_id, err)
         call h5dopen_f(group_id, 'sline3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!calculate new vmax (adapted to maximum absolute velocities occuring in the atmosphere
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         fdum = sqrt(velx3d(i,j,k)**2+vely3d(i,j,k)+velz3d(i,j,k)**2)
         if(fdum.gt.vmax) vmax=fdum
      enddo
   enddo
enddo
!
!
end subroutine read_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_spc3d
!
!read in standard 3d input from sc3d.eo (cartesian coordinates)
!read in any other input file in 3d (spherical coordinates)
!calculate the opacity of 3d input file (spherical coordinates)
!interpolate the source functions from cartesian grid onto the spherical grid
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use hdf5
use mod_interp1d, only: find_index, interpol_yp
use mod_interp3d, only: interpol3d_8p_lin
use mod_opacities, only: opalbar_model_kline, opalbar_model_hamann
use mod_iline, only: alpha, kappa0, kline, xnue0, na
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
integer(i4b) :: ndxmax, ndymax, ndzmax, nr_modext, ntheta_modext, nphi_modext
real(dp) :: xcoord, ycoord, zcoord, fdum, rad, scont, sline, rho, opalbar, mdot
real(dp) :: sint, cost, sinp, cosp
real(dp) :: kline_model, alpha_model, kappa0_model
!
! ... local arrays
real(dp), dimension(:), allocatable :: x, y, z
real(dp), dimension(:,:,:), allocatable :: velx3d_cac, vely3d_cac, velz3d_cac, &
                                           opac3d_cac, opalbar3d_cac, sline3d_cac, &
                                           scont3d_cac, t3d_cac
integer, dimension(:,:,:), allocatable :: imask3d_cac
real(dp), dimension(:), allocatable :: r_modext, theta_modext, phi_modext
real(dp), dimension(:,:,:), allocatable :: velr3d_modext, velth3d_modext, velphi3d_modext, &
                                           rho3d_modext, t3d_modext

real(dp), dimension(:,:,:), allocatable :: velr3d_flip, velth3d_flip, velphi3d_flip, rho3d_flip, t3d_flip, sline3d_flip, scont3d_flip
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_rad, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims3d_cac, dims3d_spc
!
!
write(*,*) '----------------reading standard 3d input: cartesian coordinates---------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'options', group_id, err)
         call h5aopen_f(group_id, 'opt_opal', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kline_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(kline.ne.kline_model.and.opt_opal.eq.0) then
            write(*,*) 'warning in calc_model3d_spc3db: kline ne kline_model'
            write(*,*) 'kline used for upcoming calculations (from indat file)', kline
            write(*,*) 'kline used within sc3d:', kline_model
            write(*,*)
         endif
         call h5aopen_f(group_id, 'alpha', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, alpha_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(alpha.ne.alpha_model.and.opt_opal.eq.1) then
            write(*,*) 'warning in calc_model3d_spc3db: alpha ne alpha_model'
            write(*,*) 'alpha used for upcoming calculations (from indat file)', alpha
            write(*,*) 'alpha used within sc3d:', alpha_model
            write(*,*)
         endif            
         call h5aopen_f(group_id, 'kappa0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kappa0_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(kappa0.ne.kappa0_model.and.opt_opal.eq.1) then
            write(*,*) 'warning in calc_model3d_spc3db: kappa0 ne kappa0_model'
            write(*,*) 'kappa0 used for upcoming calculations (from indat file)', kappa0
            write(*,*) 'kappa0 used within sc3d:', kappa0_model
            write(*,*)
         endif            
         call h5aopen_f(group_id, 'eps_line', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)         
         vmax=vmax*1.d5
         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'yhe', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'hei', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, hei, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'mdot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xmloss, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         mdot = xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         sr=rstar*rsu
         call h5aopen_f(group_id, 'lstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xlogg', attr_id, err)
            if(err.eq.-1) then
               write(*,*) 'error reading h5 file: use previously read xlogg'
             else
               call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            endif
         call h5aopen_f(group_id, 'vrot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         vrot=vrot*1.d5 !in cm/s
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'bcondition', group_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'ndxmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_x = (/ ndxmax /)
         call h5aopen_f(group_id, 'ndymax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_y = (/ ndymax /)
         call h5aopen_f(group_id, 'ndzmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_z = (/ ndzmax /)
      call h5gclose_f(group_id, err)

      dims3d_cac = (/ ndxmax, ndymax, ndzmax /)

      allocate(x(ndxmax), stat=err)
      allocate(y(ndymax), stat=err)
      allocate(z(ndzmax), stat=err)
      allocate(imask3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(t3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opac3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opalbar3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velx3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(vely3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velz3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(sline3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(scont3d_cac(ndxmax,ndymax,ndzmax), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 'mask3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_integer, imask3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opac3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opac3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
!convert opalbar and opac to 1/cm
         opac3d_cac=opac3d_cac/sr
         opalbar3d_cac=opalbar3d_cac/sr
         call h5dopen_f(group_id, 'velx3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)

      scont3d_cac=0.d0
      call h5gopen_f(file_id, 'solution3d', group_id, err)
         call h5dopen_f(group_id, 'sline3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, sline3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!calculate new vmax (adapted to maximum absolute velocities occuring in the atmosphere
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         fdum = sqrt(velx3d_cac(i,j,k)**2+vely3d_cac(i,j,k)+velz3d_cac(i,j,k)**2)
         if(fdum.gt.vmax) vmax=fdum
      enddo
   enddo
enddo
!
!-----------------------------------------------------------------------
!
write(*,*) '------------reading 3d atmospheric structure in spherical coordinates----------'
write(*,*) 'input_file: ', trim(input_file2)
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'nphi', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nphi_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)

call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims_phi=(/nphi_modext/)
dims3d_spc=(/ nr_modext, ntheta_modext, nphi_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(phi_modext(nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(rho3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(t3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velr3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velth3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velphi3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
!normalize radius to radius 1
      r_modext=r_modext/r_modext(1) 
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'phi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, phi_modext, dims_phi, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!------------------flip evereything (for debug reasons)-----------------
!
!allocate(velr3d_flip(nr_modext,ntheta_modext,nphi_modext))
!allocate(velth3d_flip(nr_modext,ntheta_modext,nphi_modext))
!allocate(velphi3d_flip(nr_modext,ntheta_modext,nphi_modext))
!allocate(rho3d_flip(nr_modext,ntheta_modext,nphi_modext))
!allocate(t3d_flip(nr_modext,ntheta_modext,nphi_modext))
!allocate(scont3d_flip(ndxmax,ndymax,ndzmax))
!allocate(sline3d_flip(ndxmax,ndymax,ndzmax))
!
!do i=1, ntheta_modext
!   velr3d_flip(:,i,:)=velr3d_modext(:,ntheta_modext+1-i,:)
!   velth3d_flip(:,i,:)= -velth3d_modext(:,ntheta_modext+1-i,:)
!   velphi3d_flip(:,i,:)=velphi3d_modext(:,ntheta_modext+1-i,:)
!   rho3d_flip(:,i,:)=rho3d_modext(:,ntheta_modext+1-i,:)
!   t3d_flip(:,i,:)=t3d_modext(:,ntheta_modext+1-i,:)
!enddo
!!
!do i=1, ndzmax
!   scont3d_flip(:,:,i)=scont3d(:,:,ndzmax-i+1)
!   sline3d_flip(:,:,i)=sline3d(:,:,ndzmax-i+1)
!enddo
!!
!velr3d_modext=velr3d_flip
!velth3d_modext=velth3d_flip
!velphi3d_modext=velphi3d_flip
!rho3d_modext=rho3d_flip
!t3d_modext=t3d_flip
!
!scont3d=scont3d_flip
!sline3d=sline3d_flip
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------interpolating source function and setting up opacities------------'
write(*,*)
!
nr=nr_modext
ntheta=ntheta_modext
nphi=nphi_modext
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'

r=r_modext
theta=theta_modext
phi=phi_modext
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!        
         sint=sin(theta_modext(j))
         cost=cos(theta_modext(j))
         sinp=sin(phi_modext(k))
         cosp=cos(phi_modext(k))
!        
         xcoord=r_modext(i)*sint*cosp
         ycoord=r_modext(i)*sint*sinp
         zcoord=r_modext(i)*cost
!find indices for 3d interpolation
         call find_index(xcoord, x, ndxmax, iim2, iim1, ii, iip1)
         call find_index(ycoord, y, ndymax, jjm2, jjm1, jj, jjp1)
         call find_index(zcoord, z, ndzmax, kkm2, kkm1, kk, kkp1)
!
!----------------------line source function-----------------------------
!
!if interpolated in log-space, take care about zeros
         sline = interpol3d_8p_lin(sline3d_cac(iim1,jjm1,kkm1), sline3d_cac(ii,jjm1,kkm1), &
                                   sline3d_cac(iim1,jj,kkm1),   sline3d_cac(ii,jj,kkm1), &
                                   sline3d_cac(iim1,jjm1,kk),   sline3d_cac(ii,jjm1,kk), &
                                   sline3d_cac(iim1,jj,kk),     sline3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont = interpol3d_8p_lin(scont3d_cac(iim1,jjm1,kkm1), scont3d_cac(ii,jjm1,kkm1), &
                                   scont3d_cac(iim1,jj,kkm1),   scont3d_cac(ii,jj,kkm1), &
                                   scont3d_cac(iim1,jjm1,kk),   scont3d_cac(ii,jjm1,kk), &
                                   scont3d_cac(iim1,jj,kk),     scont3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k)=velr3d_modext(i,j,k)*sint*cosp + velth3d_modext(i,j,k)*cost*cosp-velphi3d_modext(i,j,k)*sinp
         vely3d(i,j,k)=velr3d_modext(i,j,k)*sint*sinp + velth3d_modext(i,j,k)*cost*sinp+velphi3d_modext(i,j,k)*cosp
         velz3d(i,j,k)=velr3d_modext(i,j,k)*cost - velth3d_modext(i,j,k)*sint
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = t3d_modext(i,j,k)
!
!-----------------------------opacity-----------------------------------
!

!         opt_opal=1
!         kappa0=1.d0
!         alpha=0.d0
!         write(*,*) '1', vmax
!         vmax=3200.d5
         !         write(*,*) vmax
!         vth_fiducial = 30.d5
!         mdot = 3.3d-7 * xmsu / (3600.*24.*365.)
         rho = rho3d_modext(i,j,k)
         if(opt_opal.eq.0) then 
            opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
         elseif(opt_opal.eq.1) then
            opalbar = opalbar_model_hamann(sr, vmax, mdot, kappa0, alpha, vth_fiducial, r(i)*sr, rho)*sr
         else
            stop 'error in calc_model3d_spc3d: opt_opal not specified'
         endif
!
         if(t3d(i,j,k).gt.1.d5) opalbar=1.d-10
!set opalbar to cgs (actually, will be set in output3d_spc)
         opalbar3d(i,j,k) = opalbar!/sr
         opac3d(i,j,k) = 0.d0/sr
!         write(*,*) r(i), rho, opalbar
!
      enddo
   enddo
enddo
!
!vmax=5.d8
!test
!opalbar3d=opalbar3d*1.d4
!
!vmax=5.d8
!
end subroutine calc_model3d_spc3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_spc3db
!
!read in standard 3d input from sc3d.eo (cartesian coordinates)
!read in any other input file in 1d (radial stratification)
!calculate the opacity of 1d input file (radial stratification)
!create 3d spherical coordinates (with 1d radial stratification)
!interpolate the source functions from cartesian grid onto the spherical grid
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use hdf5
use mod_interp1d, only: find_index, interpol_yp
use mod_interp3d, only: interpol3d_8p_lin
use mod_grid, only: grid_equi
use mod_opacities, only: opalbar_model_kline, opalbar_model_hamann
use mod_iline, only: alpha, kappa0, kline, xnue0, na
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
integer(i4b) :: ndxmax, ndymax, ndzmax, nr_modext, ntheta_modext, nphi_modext
real(dp) :: xcoord, ycoord, zcoord, fdum, rad, scont, sline, rho, opalbar, mdot
real(dp) :: sint, cost, sinp, cosp
real(dp) :: kline_model, alpha_model, kappa0_model
!
! ... local arrays
real(dp), dimension(:), allocatable :: x, y, z
real(dp), dimension(:,:,:), allocatable :: velx3d_cac, vely3d_cac, velz3d_cac, &
                                           opac3d_cac, opalbar3d_cac, sline3d_cac, &
                                           scont3d_cac, t3d_cac
integer, dimension(:,:,:), allocatable :: imask3d_cac
real(dp), dimension(:), allocatable :: r_modext, rho1d_modext, velr1d_modext, t1d_modext
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_rad
integer(hsize_t), dimension(3) :: dims3d_cac
!
!
write(*,*) '----------------reading standard 3d input: cartesian coordinates---------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'options', group_id, err)
         call h5aopen_f(group_id, 'opt_opal', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kline_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(kline.ne.kline_model.and.opt_opal.eq.0) then
            write(*,*) 'warning in calc_model3d_spc3db: kline ne kline_model'
            write(*,*) 'kline used for upcoming calculations (from indat file)', kline
            write(*,*) 'kline used within sc3d:', kline_model
            write(*,*)
         endif
         call h5aopen_f(group_id, 'alpha', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, alpha_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(alpha.ne.alpha_model.and.opt_opal.eq.1) then
            write(*,*) 'warning in calc_model3d_spc3db: alpha ne alpha_model'
            write(*,*) 'alpha used for upcoming calculations (from indat file)', alpha
            write(*,*) 'alpha used within sc3d:', alpha_model
            write(*,*)
         endif            
         call h5aopen_f(group_id, 'kappa0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kappa0_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         if(kappa0.ne.kappa0_model.and.opt_opal.eq.1) then
            write(*,*) 'warning in calc_model3d_spc3db: kappa0 ne kappa0_model'
            write(*,*) 'kappa0 used for upcoming calculations (from indat file)', kappa0
            write(*,*) 'kappa0 used within sc3d:', kappa0_model
            write(*,*)
         endif
         call h5aopen_f(group_id, 'eps_line', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)         
         vmax=vmax*1.d5
         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'yhe', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'hei', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, hei, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'mdot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xmloss, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         mdot = xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         sr=rstar*rsu
         call h5aopen_f(group_id, 'lstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xlogg', attr_id, err)
            if(err.eq.-1) then
               write(*,*) 'error reading h5 file: use previously read xlogg'
             else
               call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            endif
         call h5aopen_f(group_id, 'vrot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         vrot=vrot*1.d5 !in cm/s
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'bcondition', group_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'ndxmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_x = (/ ndxmax /)
         call h5aopen_f(group_id, 'ndymax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_y = (/ ndymax /)
         call h5aopen_f(group_id, 'ndzmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_z = (/ ndzmax /)
      call h5gclose_f(group_id, err)

      dims3d_cac = (/ ndxmax, ndymax, ndzmax /)

      allocate(x(ndxmax), stat=err)
      allocate(y(ndymax), stat=err)
      allocate(z(ndzmax), stat=err)
      allocate(imask3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(t3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opac3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opalbar3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velx3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(vely3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velz3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(sline3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(scont3d_cac(ndxmax,ndymax,ndzmax), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 'mask3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_integer, imask3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opac3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opac3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
!convert opalbar and opac to 1/cm
         opac3d_cac=opac3d_cac/sr
         opalbar3d_cac=opalbar3d_cac/sr
         call h5dopen_f(group_id, 'velx3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)

      scont3d_cac=0.d0
      call h5gopen_f(file_id, 'solution3d', group_id, err)
         call h5dopen_f(group_id, 'sline3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, sline3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!calculate new vmax (adapted to maximum absolute velocities occuring in the atmosphere
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         fdum = sqrt(velx3d_cac(i,j,k)**2+vely3d_cac(i,j,k)+velz3d_cac(i,j,k)**2)
         if(fdum.gt.vmax) vmax=fdum
      enddo
   enddo
enddo
!
!-----------------------------------------------------------------------
!
write(*,*) '------------reading 1d atmospheric structure with radial stratification--------'
write(*,*) 'input_file: ', trim(input_file2)
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3db: allocation'
allocate(rho1d_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3db: allocation'
allocate(t1d_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3db: allocation'
allocate(velr1d_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3db: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
!normalize radius to radius 1
      r_modext=r_modext/r_modext(1) 
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho1d_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr1d_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t1d_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!do i=1, nr_modext
!   write(*,'(10es20.8)') r_modext(i), rho1d_modext(i), velr1d_modext(i), t1d_modext(i)
!enddo
!stop 'bla1'
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------interpolating source function and setting up opacities------------'
write(*,*)
!
nr=nr_modext
ntheta=51
nphi=2*ntheta-1
!
!-----------------------------------------------------------------------
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3db'

r=r_modext
call grid_equi(0.d0,2.d0*pi,nphi,phi)
!call grid_equi(1.d0,-1.d0,ntheta,theta)
!theta=acos(theta)
call grid_equi(0.d0,pi,ntheta,theta)
!
!
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!        
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
!        
         xcoord=r(i)*sint*cosp
         ycoord=r(i)*sint*sinp
         zcoord=r(i)*cost
!find indices for 3d interpolation
         call find_index(xcoord, x, ndxmax, iim2, iim1, ii, iip1)
         call find_index(ycoord, y, ndymax, jjm2, jjm1, jj, jjp1)
         call find_index(zcoord, z, ndzmax, kkm2, kkm1, kk, kkp1)
!
!----------------------line source function-----------------------------
!
!if interpolated in log-space, take care about zeros
         sline = interpol3d_8p_lin(sline3d_cac(iim1,jjm1,kkm1), sline3d_cac(ii,jjm1,kkm1), &
                                   sline3d_cac(iim1,jj,kkm1),   sline3d_cac(ii,jj,kkm1), &
                                   sline3d_cac(iim1,jjm1,kk),   sline3d_cac(ii,jjm1,kk), &
                                   sline3d_cac(iim1,jj,kk),     sline3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont = interpol3d_8p_lin(scont3d_cac(iim1,jjm1,kkm1), scont3d_cac(ii,jjm1,kkm1), &
                                   scont3d_cac(iim1,jj,kkm1),   scont3d_cac(ii,jj,kkm1), &
                                   scont3d_cac(iim1,jjm1,kk),   scont3d_cac(ii,jjm1,kk), &
                                   scont3d_cac(iim1,jj,kk),     scont3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k)=velr1d_modext(i)*sint*cosp
         vely3d(i,j,k)=velr1d_modext(i)*sint*sinp
         velz3d(i,j,k)=velr1d_modext(i)*cost
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = t1d_modext(i)
!
!-----------------------------opacity-----------------------------------
!
         rho = rho1d_modext(i)
         if(opt_opal.eq.0) then 
            opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
         elseif(opt_opal.eq.1) then
            opalbar = opalbar_model_hamann(sr, vmax, mdot, kappa0, alpha, vth_fiducial, r(i)*sr, rho)*sr
         else
            stop 'error in calc_model3d_spc3d: opt_opal not specified'
         endif
!
         if(t3d(i,j,k).gt.1.d5) opalbar=1.d-10
!set opalbar to cgs (actually, will be set in output3d_spc)
         opalbar3d(i,j,k) = opalbar!/sr
         opac3d(i,j,k) = 0.d0/sr
!         write(*,*) r(i), rho, opalbar
!
      enddo
   enddo
enddo
!
!
!
end subroutine calc_model3d_spc3db
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_spc3dc
!
!read in standard 3d input from sc3d.eo (cartesian coordinates)
!read in any other input file in 2d (spherical coordinates)
!calculate the opacity of 2d input file (spherical coordinates)
!interpolate the source functions from cartesian grid onto the spherical grid
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi, zero, two
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use hdf5
use mod_interp1d, only: find_index, interpol_yp
use mod_interp3d, only: interpol3d_8p_lin
use mod_grid, only: grid_equi
use mod_opacities, only: opalbar_model_kline, opalbar_model_hamann
use mod_iline, only: alpha, kappa0, kline, xnue0, na
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
integer(i4b) :: ndxmax, ndymax, ndzmax, nr_modext, ntheta_modext, nphi_modext
real(dp) :: xcoord, ycoord, zcoord, fdum, rad, scont, sline, rho, opalbar, mdot
real(dp) :: sint, cost, sinp, cosp
!
! ... local arrays
real(dp), dimension(:), allocatable :: x, y, z
real(dp), dimension(:,:,:), allocatable :: velx3d_cac, vely3d_cac, velz3d_cac, &
                                           opac3d_cac, opalbar3d_cac, sline3d_cac, &
                                           scont3d_cac, t3d_cac
integer, dimension(:,:,:), allocatable :: imask3d_cac
real(dp), dimension(:), allocatable :: r_modext, theta_modext, phi_modext
real(dp), dimension(:,:), allocatable :: velr2d_modext, velth2d_modext, velphi2d_modext, &
                                         rho2d_modext, t2d_modext
real(dp), dimension(:,:), allocatable :: velr2d_flip, velth2d_flip, velphi2d_flip, rho2d_flip, t2d_flip
real(dp), dimension(:,:,:), allocatable :: scont3d_flip, sline3d_flip
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_rad, dims_theta, dims_phi
integer(hsize_t), dimension(2) :: dims2d_spc
integer(hsize_t), dimension(3) :: dims3d_cac
!
!
write(*,*) '----------------reading standard 3d input: cartesian coordinates---------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'options', group_id, err)
         call h5aopen_f(group_id, 'opt_opal', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kline, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'alpha', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'kappa0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'eps_line', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)         
         vmax=vmax*1.d5
         call h5aopen_f(group_id, 'vmicro', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'yhe', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'hei', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, hei, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'mdot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xmloss, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         mdot = xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         sr=rstar*rsu
         call h5aopen_f(group_id, 'lstar', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'xlogg', attr_id, err)
            if(err.eq.-1) then
               write(*,*) 'error reading h5 file: use previously read xlogg'
             else
               call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            endif
         call h5aopen_f(group_id, 'vrot', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         vrot=vrot*1.d5 !in cm/s
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'bcondition', group_id, err)
         call h5aopen_f(group_id, 'xic1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'ndxmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_x = (/ ndxmax /)
         call h5aopen_f(group_id, 'ndymax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_y = (/ ndymax /)
         call h5aopen_f(group_id, 'ndzmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_z = (/ ndzmax /)
      call h5gclose_f(group_id, err)

      dims3d_cac = (/ ndxmax, ndymax, ndzmax /)

      allocate(x(ndxmax), stat=err)
      allocate(y(ndymax), stat=err)
      allocate(z(ndzmax), stat=err)
      allocate(imask3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(t3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opac3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(opalbar3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velx3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(vely3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(velz3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(sline3d_cac(ndxmax,ndymax,ndzmax), stat=err)
      allocate(scont3d_cac(ndxmax,ndymax,ndzmax), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 'mask3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_integer, imask3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opac3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opac3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, opalbar3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
!convert opalbar and opac to 1/cm
         opac3d_cac=opac3d_cac/sr
         opalbar3d_cac=opalbar3d_cac/sr
         call h5dopen_f(group_id, 'velx3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)

      scont3d_cac=zero
      call h5gopen_f(file_id, 'solution3d', group_id, err)
         call h5dopen_f(group_id, 'sline3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, sline3d_cac, dims3d_cac, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!
!-----------------------------------------------------------------------
!
write(*,*) '------------reading 2d atmospheric structure in spherical coordinates----------'
write(*,*) 'input_file: ', trim(input_file2)
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)

call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims_phi=(/nphi_modext/)
dims2d_spc=(/ nr_modext, ntheta_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3dc: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3dc: allocation'
allocate(rho2d_modext(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3dc: allocation'
allocate(t2d_modext(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velr2d_modext(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velth2d_modext(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velphi2d_modext(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'

!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
!normalize radius to radius 1
      r_modext=r_modext/r_modext(1) 
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho2d_modext, dims2d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr2d_modext, dims2d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth2d_modext, dims2d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi2d_modext, dims2d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t2d_modext, dims2d_spc, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)!
!------------------flip evereything (for debug reasons)-----------------
!
!allocate(velr2d_flip(nr_modext,ntheta_modext))
!allocate(velth2d_flip(nr_modext,ntheta_modext))
!allocate(velphi2d_flip(nr_modext,ntheta_modext))
!allocate(rho2d_flip(nr_modext,ntheta_modext))
!allocate(t2d_flip(nr_modext,ntheta_modext))
!allocate(scont3d_flip(ndxmax,ndymax,ndzmax))
!allocate(sline3d_flip(ndxmax,ndymax,ndzmax))
!!
!!write(*,*) 'test1'
!do i=1, ntheta_modext
!   velr2d_flip(:,i)=velr2d_modext(:,ntheta_modext+1-i)
!   velth2d_flip(:,i)= -velth2d_modext(:,ntheta_modext+1-i)
!   velphi2d_flip(:,i)=velphi2d_modext(:,ntheta_modext+1-i)
!   rho2d_flip(:,i)=rho2d_modext(:,ntheta_modext+1-i)
!   t2d_flip(:,i)=t2d_modext(:,ntheta_modext+1-i)
!enddo
!!write(*,*) 'test2'
!!
!do i=1, ndzmax
!   scont3d_flip(:,:,i)=scont3d_cac(:,:,ndzmax-i+1)
!   sline3d_flip(:,:,i)=sline3d_cac(:,:,ndzmax-i+1)
!enddo
!!write(*,*) 'test3'
!velr2d_modext=velr2d_flip
!velth2d_modext=velth2d_flip
!velphi2d_modext=velphi2d_flip
!rho2d_modext=rho2d_flip
!t2d_modext=t2d_flip
!
!write(*,*) 'test4'
!scont3d_cac=scont3d_flip
!sline3d_cac=sline3d_flip
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------interpolating source function and setting up opacities------------'
write(*,*)
!
nr=nr_modext
ntheta=ntheta_modext
!nphi=2*ntheta-1
nphi=101
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3dc'

r=r_modext
theta=theta_modext
call grid_equi(0.d0,two*pi,nphi,phi)
!
vmax = maxval(sqrt(velr2d_modext**2+velth2d_modext**2+velphi2d_modext**2))
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!        
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
!        
         xcoord=r(i)*sint*cosp
         ycoord=r(i)*sint*sinp
         zcoord=r(i)*cost
!find indices for 3d interpolation
         call find_index(xcoord, x, ndxmax, iim2, iim1, ii, iip1)
         call find_index(ycoord, y, ndymax, jjm2, jjm1, jj, jjp1)
         call find_index(zcoord, z, ndzmax, kkm2, kkm1, kk, kkp1)
!
!----------------------line source function-----------------------------
!
!if interpolated in log-space, take care about zeros
         sline = interpol3d_8p_lin(sline3d_cac(iim1,jjm1,kkm1), sline3d_cac(ii,jjm1,kkm1), &
                                   sline3d_cac(iim1,jj,kkm1),   sline3d_cac(ii,jj,kkm1), &
                                   sline3d_cac(iim1,jjm1,kk),   sline3d_cac(ii,jjm1,kk), &
                                   sline3d_cac(iim1,jj,kk),     sline3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont = interpol3d_8p_lin(scont3d_cac(iim1,jjm1,kkm1), scont3d_cac(ii,jjm1,kkm1), &
                                   scont3d_cac(iim1,jj,kkm1),   scont3d_cac(ii,jj,kkm1), &
                                   scont3d_cac(iim1,jjm1,kk),   scont3d_cac(ii,jjm1,kk), &
                                   scont3d_cac(iim1,jj,kk),     scont3d_cac(ii,jj,kk), &
                                   x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), &
                                   xcoord, ycoord, zcoord)
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k)=velr2d_modext(i,j)*sint*cosp + velth2d_modext(i,j)*cost*cosp-velphi2d_modext(i,j)*sinp
         vely3d(i,j,k)=velr2d_modext(i,j)*sint*sinp + velth2d_modext(i,j)*cost*sinp+velphi2d_modext(i,j)*cosp
         velz3d(i,j,k)=velr2d_modext(i,j)*cost - velth2d_modext(i,j)*sint
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = t2d_modext(i,j)
!
!-----------------------------opacity-----------------------------------
!
         rho = rho2d_modext(i,j)
         if(opt_opal.eq.0) then 
            opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
         elseif(opt_opal.eq.1) then
            opalbar = opalbar_model_hamann(sr, vmax, mdot, kappa0, alpha, vth_fiducial, r(i)*sr, rho)*sr
         else
            stop 'error in calc_model3d_spc3dc: opt_opal not specified'
         endif
!
         if(t3d(i,j,k).gt.1.d5) opalbar=1.d-10
!set opalbar to cgs (actually, will be set in output3d_spc)
         opalbar3d(i,j,k) = opalbar!/sr
         opac3d(i,j,k) = zero/sr
!         write(*,*) r(i), rho, opalbar
!
      enddo
   enddo
enddo
!
!write(*,*) maxval(sqrt(velx3d**2+vely3d**2+velz3d**2))
!stop
!
!
end subroutine calc_model3d_spc3dc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_be
!
!star has a Be-type disc from Dylans initial conditions
!
use prog_type
use fund_const
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi, zero, two
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: rmin, rmax, teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
     opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: alpha, kappa0, kline, xnue0, na, iline
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nz=301
integer(i4b) :: i, j, k, err, idum1, idum2, nr1, nr2, ii, opt_temperature, np
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, temp, vth
real(dp) :: sint, cost, sinp, cosp, rad
real(dp) :: zmin, zmax, tauz, dtau

real(dp) :: hi, mmw, mstar_cgs, tdisc, mdisc, dtheta_disc, theta_discl, theta_discu, rho0_disc, rdisc_max, slope
real(dp) :: csound, velphi, b2, b3, kcont, tdisc1, tdisc2, tacoeff, tbcoeff, vrot_scalefac
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
real(dp), dimension(:), allocatable :: p1d, tauz1d
real(dp), dimension(nz) :: z1d, opac1d
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff, dilfac
!
! ... for hdf5 file
!
opt_temperature = 1
!=0 for constant temperature at tdisc
!=1 for linear temperature law
!=2 for 1/r temperature law
!=3 for lucy optically thin temperature law
!=4 for carciofi temperature law
!
vrot_scalefac=1.d0
!
kcont=one
lstar=4.d0*pi*sr**2*cgs_sb*teff**4/xlsu
xic1=bnue(xnue0,trad)
!
!define disc 1
hi = 1.d0   !number free electrons for each hydrogen atom
mmw = mean_molecular_weight(hi,hei,yhe)  !mean molecular weight
mstar_cgs = sr**2 * ten**xlogg/cgs_grav!70.d0*xmsu !
vrot = sqrt(cgs_grav*mstar_cgs/sr)    !breakup velocity (consistent with v_phi of disc model)
!write(*,*) mstar_cgs/xmsu, vrot/1.d5
!stop

!either isothermal temperature or linear decay with radius
tdisc = 10.d3            !isothermal temperature of the disc
tdisc1 = 10.d3
tdisc2 = 10.d3
!
mdisc = 5.d-14 * xmsu     !*1.d-10 set to almost zero
dtheta_disc = 45.d0*pi/180.d0    !opening angle of the disc (plus/minus 45 degree here)
theta_discl = pi/two - dtheta_disc
theta_discu = pi/two + dtheta_disc
!maximum radius of the disc (in rstar)
rdisc_max=13.5d0*5.4d0/rstar
!write(*,*) rdisc_max
!slope of the disc
slope=1.5d0
!
!slope of temperature
tacoeff = tdisc1-(tdisc2-tdisc1)*rmin/(rdisc_max-rmin)
tbcoeff = (tdisc2-tdisc1)/(rdisc_max-rmin)
!
write(*,*) '----------------calc_model3d_be: Be-type disc for star 1 (Dylan)--------------'
write(*,*)
!
!
nr1=86
nr2=15
nr=nr1+nr2
idum1=12
idum2=20
ntheta=2*idum1+2*idum2+1
nphi=2*ntheta-1
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_be'
!
!radius in units of sr
allocate(fdum1_arr(nr1))
call grid_log(rmin, rdisc_max-1.d-5, nr1, fdum1_arr)
allocate(fdum2_arr(nr2))
call grid_log(rdisc_max+1.d-5, rmax, nr2, fdum2_arr)
ii=1
do i=1, nr1
   r(ii) = fdum1_arr(i)
!   write(*,*) ii, r(ii)   
   ii=ii+1
enddo
do i=1, nr2
   r(ii) = fdum2_arr(i)
!   write(*,*) ii, r(ii)   
   ii=ii+1
enddo
deallocate(fdum1_arr)
deallocate(fdum2_arr)
!
if(r(nr-1).le.rdisc_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
theta(1)=zero
theta(ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, nphi, phi)
!
!
!-------calculate temperature structure as function of r*sin(theta)-----
!
np=nr
allocate(p1d(np))
allocate(tauz1d(np))
p1d=r
!
csound = vsound(teff,mmw)            
rho0_disc = mdisc*sqrt(cgs_grav*mstar_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr**3.5d0
do i=1, np
   zmin = p1d(i)*cos(theta_discu)   
   zmax = p1d(i)*cos(theta_discl)
   call grid_equi(zmin, zmax, nz, z1d)
   do j=1, nz
      rad = sqrt(z1d(j)**2 + p1d(i)**2)
      rho = rho0_disc*(p1d(i))**(-slope) * exp(cgs_grav*mstar_cgs/csound**2 * (one/sr/rad-one/sr/p1d(i)))
      opac1d(j) = opac_thomson(yhe, hei, rho, kcont)*sr
   enddo
   tauz = 0.d0
   do j=2, nz
      dtau = 0.5d0*(opac1d(j-1)+opac1d(j))*(z1d(j)-z1d(j-1))
      tauz = tauz + dtau
   enddo
   tauz1d(i) = tauz
!   write(*,*) p1d(i), zmin, zmax, tauz
enddo
!stop
!
!-----------------------------------------------------------------------
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary         
         rad=r(i)
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
         
         if(rad.le.rdisc_max.and. &
            theta(j).gt.theta_discl.and.&
            theta(j).lt.theta_discu) then
!
            velphi = vrot_scalefac*sqrt(cgs_grav*mstar_cgs/rad/sr/sint)

            if(opt_temperature.eq.0) then
               temp=tdisc
            elseif(opt_temperature.eq.1) then
               temp=tacoeff + tbcoeff*rad
            elseif(opt_temperature.eq.2) then
               temp = max(teff/rad,5d3)
            elseif(opt_temperature.eq.3) then
               temp = teff * dilfac(rad)**0.25d0
            elseif(opt_temperature.eq.4) then
!neglect sint term
!               temp = teff/pi**0.25d0 * (asin(one/rad/sint) - sqrt(one-one/(rad*sint)**2)/rad/sint)**0.25d0
               temp = teff/pi**0.25d0 * (asin(one/rad) - sqrt(one-one/(rad)**2)/rad)**0.25d0               
!               if(temp.ne.temp) temp = 0.8d0*teff
               if(tauz1d(i).lt.0.1d0) temp = max(temp,0.6d0*teff)
            else
               stop 'opt_temperature not properly specified'
            endif
            vth = vthermal(vmicro, temp, na)
            csound = vsound(temp,mmw)            
            rho0_disc = mdisc*sqrt(cgs_grav*mstar_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr**3.5d0            
            rho = rho0_disc*(rad*sint)**(-slope) * exp(cgs_grav*mstar_cgs/csound**2 * (one/sr/rad-one/sr/rad/sint))
            if(j.eq.(ntheta/2+1).and.k.eq.1) then
               write(*,*) theta(j), r(i), rho, rho0_disc
            endif
!
!departure coeffcients in LTE
            b2=one!0.9!one
            b3=one!1.1!one
            opalbar = get_opalbar(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)/sr   !in cgs
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)
            opac = opac_thomson(yhe, hei, rho, kcont)
!            write(*,*) rho, opalbar
!            write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
         else
            velphi = zero
            temp = 0.8d0*teff
!            temp = teff * dilfac(rad)**0.25d0
            vth = vthermal(vmicro, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif   
!        
!----------------------line source function-----------------------------
!
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k) = -velphi*sinp
         vely3d(i,j,k) = velphi*cosp
         velz3d(i,j,k) = zero
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         opalbar3d(i,j,k) = opalbar*sr   !in units of sr
         opac3d(i,j,k) = opac*sr
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
         vel=sqrt((velx3d(i,j,k))**2+(vely3d(i,j,k))**2+(velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!vrot=0.d0*vrot
!
end subroutine calc_model3d_be
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_ileyk
!
!star a Be-type disc from ileyks simulation
!
use prog_type
use fund_const
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi, zero, two
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: rmin, rmax, teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use hdf5
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: alpha, kappa0, kline, xnue0, na, iline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont, velr, velth, velphi
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, temp, vth
real(dp) :: sint, cost, sinp, cosp, rad

real(dp) :: hi, mmw, mstar_cgs, tdisc, mdisc, dtheta_disc, theta_discl, theta_discu, rho0_disc, rdisc_max, slope
real(dp) :: csound, b2, b3, kcont, tdisc1, tdisc2, tacoeff, tbcoeff
real(dp), parameter :: rho_scalefac=.1d0, vrot_scalefac=1.3d0
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
real(dp), dimension(:,:,:), allocatable :: rho3d, velr3d, velth3d, velphi3d
!
! ... local characters
character(len=19), parameter :: fname='models/ileyk/f10.h5'
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff
!
! ... for hdf5 file
integer(hid_t) :: file_id, dset_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
kcont=one
lstar=4.d0*pi*sr**2*cgs_sb*teff**4/xlsu
xic1=bnue(xnue0,trad)
!
!define disc 1
hi = 1.d0   !number free electrons for each hydrogen atom
mmw = mean_molecular_weight(hi,hei,yhe)  !mean molecular weight
mstar_cgs = sr**2 * ten**xlogg/cgs_grav
vrot = sqrt(cgs_grav*mstar_cgs/sr)    !breakup velocity (consistent with v_phi of disc model)


!either isothermal temperature or linear decay with radius
tdisc = 10.d3            !isothermal temperature of the disc
tdisc1 = 10.d3
tdisc2 = 10.d3
!
!maximum radius of the disc (in rstar)
rdisc_max=10.d0
!
!slope of temperature
tacoeff = tdisc1-(tdisc2-tdisc1)*rmin/(rdisc_max-rmin)
tbcoeff = (tdisc2-tdisc1)/(rdisc_max-rmin)
!
write(*,*) '----------------calc_model3d_ileyk: Be-type disc from hydro-------------------'
write(*,*)
!

call h5open_f (err)
   call h5fopen_f(fname, h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nr, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ntheta, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'nphi', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nphi, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      dims_r = (/ nr /)
      dims_theta = (/ ntheta /)
      dims_phi = (/ nphi /)      
      dims_3d = (/ nr, ntheta, nphi /)

!allocate arrays
      allocate(r(nr), stat=err)
      allocate(theta(ntheta), stat=err)
      allocate(phi(nphi), stat=err)
      allocate(rho3d(nr,ntheta,nphi), stat=err)
      allocate(velr3d(nr,ntheta,nphi), stat=err)
      allocate(velth3d(nr,ntheta,nphi), stat=err)
      allocate(velphi3d(nr,ntheta,nphi), stat=err)      
      allocate(t3d(nr,ntheta,nphi), stat=err)

!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, theta, dims_theta, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'phi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, phi, dims_phi, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'rho3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         rho3d=rho3d*rho_scalefac   
         call h5dopen_f(group_id, 'velr3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velth3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velphi3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         velphi3d=velphi3d*vrot_scalefac   
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!scale to rstar
r=r/r(1)*sr

!radius in units of sr1
r=r/sr

if(abs(rmin-minval(r)).gt.small_number) then
   write(*,*) 'error in calc_model3d_ileyk: rmin not matching'
   write(*,*) 'set rmin to', minval(r)
   stop
endif
if(abs(rmax-maxval(r)).gt.small_number) then
   write(*,*) 'error in calc_model3d_ileyk: rmax not matching'
   write(*,*) 'set rmax to', maxval(r)
   stop
endif

allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_ileyk'

!
!
!
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
         rad=r(i)
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
         
!manipulate temperature
         t3d(i,j,k) = tacoeff + tbcoeff*rad
!
         velr = velr3d(i,j,k)
         velth = velth3d(i,j,k)
         velphi = velphi3d(i,j,k)
         temp = t3d(i,j,k)
         rho = rho3d(i,j,k)
!
!departure coeffcients in LTE
         b2=one
         b3=one
         opalbar = get_opalbar(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)/sr   !in cgs
         sline = sline_depcoeff(xnue0, temp, b2, b3)
         scont = bnue(xnue0,temp)
         opac = opac_thomson(yhe, hei, rho, kcont)
!         write(*,*) rho, opalbar
!         write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
!        
!----------------------line source function-----------------------------
!
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k) = -velphi*sinp
         vely3d(i,j,k) = velphi*cosp
         velz3d(i,j,k) = zero
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         opalbar3d(i,j,k) = opalbar*sr   !in units of sr
         opac3d(i,j,k) = opac*sr
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
         vel=sqrt((velx3d(i,j,k))**2+(vely3d(i,j,k))**2+(velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
end subroutine calc_model3d_ileyk

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_windlte
!
!beta-velocity wind in lte
!
use prog_type
use fund_const
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi, zero, two
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d
use params_modspec, only: rmin, rmax, teff, trad, xic1, vth_fiducial, sr, rstar, beta, vmin, vmax, vmicro, &
     opt_opal, eps_line, lstar, xlogg, vrot, yhe, hei, xmloss
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: alpha, kappa0, kline, xnue0, na, iline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, velr, velx, vely, velz, opac, opalbar, rho, sline, scont
real(dp) :: mdot_cgs, temp, vth
real(dp) :: sint, cost, sinp, cosp, rad

real(dp) :: mstar_cgs, twind
real(dp) :: csound, b2, b3, kcont, bconst
!
! ... local arrays
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff
!
! ... for hdf5 file
!
kcont=one
kline=one
lstar=4.d0*pi*sr**2*cgs_sb*teff**4/xlsu
xic1=bnue(xnue0,trad)
!
!define disc 1
mstar_cgs = sr**2 * ten**xlogg/cgs_grav
vrot = sqrt(cgs_grav*mstar_cgs/sr)    !breakup velocity (consistent with v_phi of disc model)
vrot = 0.d0
!
twind=0.8d0*teff
mdot_cgs = xmloss*xmsu/365.25d0/24.d0/3600.d0
bconst = (1.d0-vmin/vmax)**(1.d0/beta)
!
write(*,*) '----------------calc_model3d_windlte: beta-vel wind in lte--------------------'
write(*,*)
!
!
nr=101
ntheta=65
nphi=2*ntheta-1
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_windlte'
!
!radius in units of sr
call grid_log(rmin, rmax, nr,  r)
!
call grid_equi(zero, pi, ntheta, theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, nphi, phi)
!
!
!
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary         
         rad=r(i)

         velr = bvel(rad, vmax, bconst, beta)
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))         
         cosp=cos(phi(k))
!
         velx = velr*sint*cosp
         vely = velr*sint*sinp
         velz = velr*cost
!             
         temp = twind
         vth = vthermal(vmicro, temp, na)
         rho = mdot_cgs/4.d0/pi/(rad*sr)**2/velr
!departure coeffcients in LTE
         b2=one
         b3=one
         opalbar = get_opalbar(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)/sr   !in cgs
         sline = sline_depcoeff(xnue0, temp, b2, b3)
         scont = bnue(xnue0,twind)
         opac = opac_thomson(yhe, hei, rho, kcont)
!        
!----------------------line source function-----------------------------
!
         sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         velx3d(i,j,k) = velx
         vely3d(i,j,k) = vely
         velz3d(i,j,k) = velz
!
!---------------------------temperature---------------------------------
!
         t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         opalbar3d(i,j,k) = opalbar*sr   !in units of sr
         opac3d(i,j,k) = opac*sr
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
         vel=sqrt((velx3d(i,j,k))**2+(vely3d(i,j,k))**2+(velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
end subroutine calc_model3d_windlte
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model1d_sobolev
!
use prog_type
use fund_const, only: rsu, cgs_clight, xmsu, pi, one
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d, sline1d, scont1d
use params_modspec, only: vmin, vmax, beta, teff, trad, tmin, yhe, hei, eps_line, xmloss, &
                          rmax, sr, vmicro, vth_fiducial, xic1, sr, rstar
use options_modspec, only: input_mod
use mod_opacities, only: opalbar_model_kline, opalbar_model_hamann, get_opalbar
use mod_iline, only: alpha, kappa0, kline, xnue0, na, iline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: del, bconst, xmloss_cgs, rho, gradv
!
! ... local characters
!
! ... local functions
real(dp) :: bvel, vthermal, sobo1d, bnue
!
write(*,*) '---------calculating sobolev solution for line-strength parameterization-------'
write(*,*)
write(*,*) 'vmin', vmin
write(*,*) 'vmax', vmax
write(*,*) 'beta', beta
write(*,*) 'teff', teff
write(*,*) 'trad', trad
write(*,*) 'tmin', tmin
write(*,*) 'yhe', yhe
write(*,*) 'hei', hei
write(*,*) 'kline', kline
write(*,*) 'alpha', alpha
write(*,*) 'kappa0', kappa0
write(*,*) 'eps_line', eps_line
!
!-----------------------------------------------------------------------
!
n1d=41
!
!allocate arrays
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(sline1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error calc_model1d_kline'
allocate(velr1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(vth1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(opac1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(scont1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
allocate(t1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_kline'
!
!---------------calculate radial stratification (in log)----------------
!
r1d(1)=1.d0
r1d(2)=1.d0+1.d-4
del=log10(rmax/r1d(2))/(n1d-2)
!
do i=3, n1d
   r1d(i) = r1d(i-1)*10.d0**del
enddo
!write(*,*) r1d
!
!-----------------------------------------------------------------------
!
bconst=1.d0-(vmin/vmax)**(1.d0/beta)
xmloss_cgs = xmloss*xmsu/(365.25d0*24.d0*3600.d0)

xic1=bnue(xnue0, trad)
!
do i=1, n1d
   velr1d(i) = bvel(r1d(i), vmax, bconst, beta)
!
!set temperature to trad or to tmin
   t1d(i) = max(trad,tmin)
!thermal velocity
   vth1d(i) = vthermal(vmicro, t1d(i), na)
!density
   rho = xmloss_cgs/(4.d0*pi*r1d(i)**2 * sr**2 * velr1d(i))
!
!opacity
   if(input_mod.eq.1) then
      opalbar1d(i) = get_opalbar(iline, kline, sr, yhe, hei, t1d(i), vth_fiducial, xnue0, one, one, rho)/sr  !in cgs      
!      opalbar1d(i) = opalbar_model_kline(yhe, hei, rho, kline)
!      opalbar1d(i) = opalbar_model_hamann(sr, vmax, xmloss_cgs, kappa0, alpha, vth_fiducial, r1d(i), rho)
   elseif(input_mod.eq.2) then
      opalbar1d(i) = opalbar_model_hamann(sr, vmax, xmloss_cgs, kappa0, alpha, vth_fiducial, r1d(i), rho)
   else
      stop 'error in calc_model1d_sobolev: input_mod not appropriatelly set'
   endif
!
!sobolev line source function
   gradv = velr1d(i) * bconst*beta/r1d(i)**2/(1.d0-bconst/r1d(i))
   sline1d(i) = sobo1d(r1d(i), velr1d(i)/vth1d(i), gradv/vth1d(i), opalbar1d(i)*sr, t1d(i), xic1, xnue0, eps_line)
   
enddo
!


!write(*,*) iline, kline, rstar, sr, yhe, hei, t1d(1), vth_fiducial, xnue0
!stop 'go on in sobo'

scont1d=0.d0
opac1d=0.d0
!
!
!
end subroutine calc_model1d_sobolev
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model1d_luka
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, one
use options_modspec, only: input_file
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d, sline1d, scont1d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, yhe, hei
use mod_opacities, only: opalbar_halpha
use mod_iline, only: xnue0
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: lambda, fdum, b2, b3
!
! ... local arrays
real(dp), dimension(:), allocatable :: rho1d, kappa1d
!
! ... local characters
character(len=50) :: cdum
!
! ... local logcials
logical :: ldum
!
! ... local functions
real(dp) :: bnue, sline_depcoeff, sline_petrenz
!
!
write(*,*) '-------------------model for Lukas Wolf-Rayet wind-----------------------------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
!
n1d=2000

allocate(r1d(n1d))
allocate(velr1d(n1d))
allocate(t1d(n1d))
allocate(rho1d(n1d))
allocate(kappa1d(n1d))
!
open(1, file=trim(input_file), form='formatted')
   read(1,*) cdum
!
   do i=1, n1d
      read(1,'(5es13.5)') r1d(i), velr1d(i), rho1d(i), t1d(i), kappa1d(i)
!      write(*,'(5es20.8)') r1d(i), velr1d(i), rho1d(i), t1d(i), kappa1d(i)
   enddo
!
close(1)
!
sr=r1d(1)
rstar=sr/rsu
!
allocate(opac1d(n1d))
allocate(opalbar1d(n1d))
allocate(sline1d(n1d))
allocate(scont1d(n1d))
!
do i=1, n1d
   opac1d(i) = rho1d(i)*kappa1d(i)
   scont1d(i) = bnue(xnue0,t1d(i))

!departure coefficients (thus far in LTE)
   b2=one
   b3=one
   opalbar1d(i) = opalbar_halpha(sr, yhe, hei, t1d(i), vth_fiducial, xnue0, b2, b3, rho1d(i))
   sline1d(i) = sline_depcoeff(xnue0, t1d(i), b2, b3)
!
enddo
!
!since final output shall be in cgs to be properly read in in spec.f90
opalbar1d=opalbar1d/sr
velr1d=velr1d!/vth_fiducial
!
xic1=bnue(xnue0,trad)
!
!final output for radius in rstar
r1d=r1d/sr
!
end subroutine calc_model1d_luka
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model1d_petrenz
!
use prog_type
use fund_const, only: rsu, xmsu, pi, cgs_clight
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d, sline1d, scont1d
use params_modspec, only: vmin, vmax, beta, teff, trad, tmin, yhe, hei, eps_line, xmloss, &
                          rmax, sr, vmicro, vth_fiducial, xic1
use mod_opacities, only: opalbar_petrenz, depcoeff_petrenz
use mod_iline, only: alpha, kappa0, kline, xnue0, na
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err, idum
real(dp) :: deldop_fiducial, xmloss_cgs, b2, b3, t_el, rho, temp, &
            fdum1, fdum2, fdum3, delr, velr
!
! ... local characters
character (len=500) :: header
!
! ... local functions
real(dp) :: sline_petrenz
real(dp) :: bnue
!
!

write(*,*) '------------calculating 1d H-alpha model (from Petrenz&Puls 1995)--------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
n1d=41
!
!allocate arrays
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(sline1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error calc_model1d_petrenz'
allocate(velr1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(vth1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(opac1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(scont1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
allocate(t1d(n1d), stat=err)
   if(err.ne.0) stop 'allocation error calc_model1d_petrenz'
!
!----------------create radial grid equidistant in log-log--------------
!
r1d(1) = 1.d0
r1d(2) = r1d(1) + 1.d-3
delr = log10(rmax)/log10(r1d(2))
delr = log10(delr)/(n1d-2)
!
do i=3, n1d
   r1d(i) = 10.d0**(log10(r1d(i-1))*10.d0**delr)
enddo
!
!----------------------calculate halpha-model---------------------------
!
!for debug
!open(1, file='TRASH/forsol/model_test.dat')
!
do i=1, n1d
!velr1d in cgs so far
   call bvel3d(vmin, vmax, beta, r1d(i), 0.d0, 0.d0, velr, fdum1, fdum2, fdum3)
   velr1d(i) = velr
   velr=velr/vmax
!
!departure coefficients
   call depcoeff_petrenz(velr, b2, b3)
!in lte
!   b2=1.d0
!   b3=1.d0
!
!integrated line opacity
   temp=0.75d0*teff
   xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
   rho=xmloss_cgs/(4.d0*pi*r1d(i)*r1d(i)*sr*sr*velr1d(i))
   opalbar1d(i) = opalbar_petrenz(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
!
!line source function
   sline1d(i) = sline_petrenz(temp, b2, b3)/bnue(xnue0,trad)
!
!temperature
   t1d(i) = temp

!   write(1,'(i5,7es20.8)') i, r1d(i), velr1d(i), rho, sline1d(i), 0.d0, opalbar1d(i)*sr, 0.d0
!
enddo
!
scont1d=0.d0
opac1d=0.d0
!
!since final output shall be in cgs to be properly read in in spec.f90
opalbar1d=opalbar1d/sr
velr1d=velr1d!/vth_fiducial
!
xic1=1.d0
!
!close(1)

!stop 'go on from here'

!
!
end subroutine calc_model1d_petrenz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model2d_petrenz
!
use prog_type
use options_modspec, only: input_file
use fund_const, only: rsu, xmsu, pi, cgs_clight
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, sline3d, scont3d, opalbar3d, opac3d, t3d, velx3d, vely3d, velz3d
use params_modspec, only: rmax, sr, yhe, hei, vth_fiducial, trad, xic1
use hdf5
use mod_interp1d, only: find_index, interpol_yp
use mod_interp2d, only: interpol2d_4p_lin
use mod_grid, only: grid_loglog, grid_equi
use mod_opacities, only: opalbar_petrenz, depcoeff_petrenz
use mod_iline, only: xnue0
!
!***debug start
use params_modspec, only: vmin, vmax, teff, beta, xmloss
!***debug end
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, &
                nr_modext, ntheta_modext
real(dp) :: delr, rho, velr, velth, velphi, opalbar, sline, temp, b2, b3, &
            sint, cost, sinp, cosp, vinft
!
!***debug start
real(dp) :: fdum1, fdum2, fdum3, xmloss_cgs

!***debug end
!
! ... local arrays
real(dp), dimension(:), allocatable :: r_modext, theta_modext
real(dp), dimension(:,:), allocatable :: rho2d_modext, velr2d_modext, velth2d_modext, velphi2d_modext, &
                                         t2d_modext, vth2d_modext
!
! ... local characters
character (len=500) :: header
!
! ... local functions
real(dp) :: sline_petrenz
real(dp) :: bnue
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r, dims_theta
integer(hsize_t), dimension(2) :: dims_2d
!
!
write(*,*) '--------calculating/reading 2d H-alpha model (from Petrenz&Puls 1995)----------'
write(*,*) 'input_file: ', trim(input_file)
write(*,*)
 
call h5open_f (err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_r = (/ nr_modext /)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         dims_theta = (/ ntheta_modext /)
      call h5gclose_f(group_id, err)

      dims_2d = (/ nr_modext, ntheta_modext /)

      allocate(r_modext(nr_modext), stat=err)
      allocate(theta_modext(ntheta_modext), stat=err)
      allocate(t2d_modext(nr_modext,ntheta_modext), stat=err)
      allocate(rho2d_modext(nr_modext,ntheta_modext), stat=err)
      allocate(velr2d_modext(nr_modext,ntheta_modext), stat=err)
      allocate(velth2d_modext(nr_modext,ntheta_modext), stat=err)
      allocate(velphi2d_modext(nr_modext,ntheta_modext), stat=err)
      allocate(vth2d_modext(nr_modext,ntheta_modext), stat=err)
!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r_modext, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model', group_id, err)
         call h5dopen_f(group_id, 'rho', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'temperature', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velr', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velth2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velphi2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vth2d_modext, dims_2d, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
!
nr=81
ntheta=41
nphi=81
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model2d_petrenz'
!
!---------------------create 3d spherical grid--------------------------
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(1.d0, 1.d0+1.d-4, rmax, nr, r)
!
!theta-grid equidistant
call grid_equi(0.d0, pi, ntheta, theta)
!
!phi-grid equidistant
call grid_equi(0.d0, 2.d0*pi-1.d-3, nphi, phi)
!
!----------------------calculate halpha-model---------------------------
!
do i=1, nr
   do j=1, ntheta
!find indices for 2d interpolation
      call find_index(r(i), r_modext, nr_modext, iim2, iim1, ii, iip1)
      call find_index(theta(j), theta_modext, ntheta_modext, jjm2, jjm1, jj, jjp1)
!
!------------------------velocity components----------------------------
!
      velr = interpol2d_4p_lin(velr2d_modext(iim1,jjm1), velr2d_modext(ii,jjm1), &
                               velr2d_modext(iim1,jj),  velr2d_modext(ii,jj), &
                               r_modext(iim1), r_modext(ii), theta_modext(jjm1), theta_modext(jj), &
                               r(i), theta(j))
      velth = interpol2d_4p_lin(velth2d_modext(iim1,jjm1), velth2d_modext(ii,jjm1), &
                                velth2d_modext(iim1,jj),  velth2d_modext(ii,jj), &
                                r_modext(iim1), r_modext(ii), theta_modext(jjm1), theta_modext(jj), &
                                r(i), theta(j))
!to be consistent with petrenz paper: set polar velocity to zero for line spectrum
!      velth = 0.d0
      velphi = interpol2d_4p_lin(velphi2d_modext(iim1,jjm1), velphi2d_modext(ii,jjm1), &
                                velphi2d_modext(iim1,jj),  velphi2d_modext(ii,jj), &
                                r_modext(iim1), r_modext(ii), theta_modext(jjm1), theta_modext(jj), &
                                r(i), theta(j))
!terminal velocity as a function of theta
      vinft = interpol_yp(theta_modext(jjm1),theta_modext(jj),velr2d_modext(nr_modext,jjm1),velr2d_modext(nr_modext,jj), theta(j))
!
!---------------------------temperature---------------------------------
!
      temp = interpol2d_4p_lin(t2d_modext(iim1,jjm1), t2d_modext(ii,jjm1), &
                               t2d_modext(iim1,jj),  t2d_modext(ii,jj), &
                               r_modext(iim1), r_modext(ii), theta_modext(jjm1), theta_modext(jj), &
                               r(i), theta(j))
!
!-----------------------------opacity-----------------------------------
!
!interpolate density (in log-log for radial and log-lin for theta)
      rho = interpol2d_4p_lin(log10(rho2d_modext(iim1,jjm1)), log10(rho2d_modext(ii,jjm1)), &
                              log10(rho2d_modext(iim1,jj)),   log10(rho2d_modext(ii,jj)), &
                              log10(r_modext(iim1)), log10(r_modext(ii)), theta_modext(jjm1), theta_modext(jj), &
                              log10(r(i)), theta(j))
      rho = 10.d0**rho
      call depcoeff_petrenz(velr/vinft, b2, b3)
!in lte
!     b2=1.d0
!     b3=1.d0
!
      opalbar = opalbar_petrenz(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
      sline = sline_petrenz(temp,b2,b3)/bnue(xnue0,trad) !normalized to trad

      do k=1, nphi
         opalbar3d(i,j,k) = opalbar
         opac3d(i,j,k) = 0.d0
         scont3d(i,j,k) = 0.d0
         sline3d(i,j,k) = sline
         t3d(i,j,k) = temp

         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
         velx3d(i,j,k)=velr*sint*cosp + velth*cost*cosp-velphi*sinp
         vely3d(i,j,k)=velr*sint*sinp + velth*cost*sinp+velphi*cosp
         velz3d(i,j,k)=velr*cost - velth*sint
!
      enddo
   enddo
enddo
!
xic1=1.d0   !since sline normalized to bnue(xnue0,trad)
!
!
!****debug start
!write(*,*) vmin, xnue0, trad
!do i=1, nr
!!velr1d in cgs so far
!   call bvel3d(vmin, vmax, beta, r(i), 0.d0, 0.d0, velr, fdum1, fdum2, fdum3)
!
!   velr=velr/vmax
!!
!departure coefficients
!   call depcoeff_petrenz(velr, b2, b3)
!!in lte
!!   b2=1.d0
!!   b3=1.d0
!!
!!integrated line opacity
!   temp=0.75d0*teff
!   xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!   rho=xmloss_cgs/(4.d0*pi*r(i)*r(i)*sr*sr*velr*vmax)
!   opalbar = opalbar_petrenz(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
!!
!!line source function
!   sline = sline_petrenz(temp, b2, b3)/bnue(xnue0,trad)
!
!   write(*,'(20es20.8)') r(i), velr*vmax, sqrt(velx3d(i,1,1)**2 + vely3d(i,1,1)**2 + velz3d(i,1,1)**2), sline, sline3d(i,1,1), &
!                        opalbar, opalbar3d(i,1,1), sline_petrenz(temp,b2,b3), bnue(xnue0,trad), b2, b3
!
!enddo
!stop
!
!***debug end
!!
!!
!!
!!
end subroutine calc_model2d_petrenz
!
!***************************TO IMPLEMENT PROPERLY***********************
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!subroutine read_model1d_halpha
!!
!use prog_type
!use fund_const, only: rsu, cgs_clight
!use dime_model1d, only: n1d, r1d, opalbar1d, velr1d, sline1d, &
!                         acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d, &
!                         acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d, &
!                         acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d
!!use options_spec, only: inputdir_fs
!use params_spec, only: teff, trad, tmin, xnue0, sr, beta, vmax, vmin, xic1, rstar, xmloss, &
!                     xlogg, yhe, beta, vth_fiducial, vmicro
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i, err, idum
!real(dp) :: xnue1, fdum1, fdum2, fdum3, deldop_fiducial
!!
!! ... local characters
!character (len=500) :: header
!!
!! ... local functions
!real(dp) :: bnue
!!
!!
!write(*,*) '------------reading 1d h-alpha-model-------------------'
!write(*,*) '   sline from jos program'
!write(*,*) '   line opacity from jos program'
!write(*,*) '   radial velocity from jos program'
!write(*,*)
!!
!!open(1, file=trim(inputdir_fs)//'/model1d.dat', form='formatted')
!!
!!read physical parameter
!   read(1,*) header
!   read(1,'(i5, 12es20.8)')  n1d, teff, trad, tmin, xlogg, yhe, xmloss, vmax, vmin, beta, xnue0, xnue1, rstar
!
!   if(xnue0.ne.xnue1) stop 'error read_model1d: xnue0, xnue1 are not equal'
!   xic1=bnue(xnue0,trad)
!   xic1=1.d0
!   vmax=vmax*1.d-5
!   vmin=vmin*1.d-5
!   sr = rstar*rsu
!!
!!allocate arrays
!   allocate(r1d(n1d), stat=err)
!      if(err.ne.0) stop 'error read_model1d: r1d'
!   allocate(sline1d(n1d), stat=err)
!      if(err.ne.0) stop 'error read_model1d: sline1d'
!   allocate(opalbar1d(n1d), stat=err)
!      if(err.ne.0) stop 'error read_model1d: opalbar1d'
!   allocate(velr1d(n1d), stat=err)
!      if(err.ne.0) stop 'error read_model1d: velr1d'
!
!
!   read(1,*) header
!   do i=1, n1d
!      read(1,'(i5, 7es20.8)') idum, r1d(n1d+1-i), velr1d(n1d+1-i), fdum1, sline1d(n1d+1-i), fdum2, opalbar1d(n1d+1-i), fdum3
!   enddo
!!
!close(1)
!!
!!convert line opacity to own units
!deldop_fiducial=xnue0*vth_fiducial/cgs_clight
!opalbar1d=opalbar1d/deldop_fiducial
!!
!!convert velocity to vth_fiducial
!velr1d=velr1d/vth_fiducial
!!
!!--------prepare coefficients for cubic spline interpolation------------
!!
!write(*,*) '-------preparing cubic spline coefficients-------------'
!write(*,*)
!
!allocate(acoeff_sline1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: acoeff_sline1d'
!allocate(acoeff_opalbar1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: acoeff_opalbar1d'
!allocate(acoeff_velr1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: acoeff_velr1d'
!allocate(bcoeff_sline1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: bcoeff_sline1d'
!allocate(bcoeff_opalbar1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: bcoeff_opalbar1d'
!allocate(bcoeff_velr1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: bcoeff_velr1d'
!allocate(ccoeff_sline1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: ccoeff_sline1d'
!allocate(ccoeff_opalbar1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: ccoeff_opalbar1d'
!allocate(ccoeff_velr1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: ccoeff_velr1d'
!allocate(dcoeff_sline1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: dcoeff_sline1d'
!allocate(dcoeff_opalbar1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: dcoeff_opalbar1d'
!allocate(dcoeff_velr1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model1d: dcoeff_velr1d'
!!
!!monotonic spline interpolation
!call cube_mono_coeff(n1d, r1d, sline1d, acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d)
!call cube_mono_coeff(n1d, r1d, opalbar1d, acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d)
!call cube_mono_coeff(n1d, r1d, velr1d, acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d)
!!
!!or normal splines:
!!call spline_coeff(n1d, r1d, sline1d, acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d)
!!call spline_coeff(n1d, r1d, opalbar1d, acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d)
!!call spline_coeff(n1d, r1d, velr1d, acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d)
!!
!!
!end subroutine read_model1d_halpha
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!subroutine read_model2d_adm
!!
!use prog_type
!use fund_const, only: rsu, cgs_clight, pi
!use dime_model1d, only: n1d, r1d
!use dime_model2d, only: ntheta, latitude, sline2d, opalbar2d, &
!                         vel_r2d, vel_theta2d, temperature2d
!use params_spec, only: teff, trad, tmin, xnue0, sr, beta, vmax, vmin, xic1, rstar, xmloss, &
!                     xlogg, yhe, beta, sr, vth_fiducial, vmicro, rmax, chi_inf, &
!                     rhoc_star, delta, t_inf, v_esc, v_inf, vmin, ralfven
!!use options_spec, only: inputdir_fs
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i, j, err, tindx, n1d_jo
!real(dp) :: mu, mu_star, mu_min, mu_lower, mu_upper, mu_shock, rad, r_apex, r_shock
!real(dp) :: del, max, min, rho, temp, temp1d, vel, vel_r, vel_theta, vel_r2, vel_theta2
!real(dp) :: fdum1, fdum2, fdum3, b2, b3, opalbar, sline
!!
!! ... local arrays
!real(dp), dimension(:), allocatable :: r1d_jo, t1d_jo, r1d_rev, t1d_rev
!! ... local characters
!!
!! ... local functions
!real(dp) :: opalbar_petrenz
!real(dp) :: sline_petrenz
!real(dp) :: bnue
!real(dp) :: interpol_yp
!!
!!
!write(*,*) '-----------reading 2d adm model (h-alpha)--------------'
!write(*,*) '   sline is calculated using departure coefficients (puls et al 1996)'
!write(*,*) '   line opacity is calculated using departure coefficients (puls et al 1996)'
!write(*,*) '   velocity field and density are calculated from adm-model'
!write(*,*) '   temperature stratification is from a 1d-input-file'
!write(*,*)
!!
!!--------------------read in temperature--------------------------------
!!
!n1d_jo=n1d
!allocate(r1d_jo(n1d_jo), stat=err)
!allocate(t1d_jo(n1d_jo), stat=err)
!allocate(r1d_rev(n1d_jo), stat=err)
!allocate(t1d_rev(n1d_jo), stat=err)
!
!open(1, file='temp_jo.dat', form='unformatted')
!   read(1), t1d_rev
!close(1)
!!
!open(1, file='rad_jo.dat', form='unformatted')
!!   read(1), r1d_rev
!close(1)
!!
!!reverse arrays
!do i=1, n1d_jo
!   r1d_jo(i)=r1d_rev(n1d_jo+1-i)
!   t1d_jo(i)=t1d_rev(n1d_jo+1-i)
!enddo
!!
!deallocate(r1d_rev)
!deallocate(t1d_rev)
!!
!!-----------------------alllocate arrays--------------------------------
!!
!allocate(r1d(n1d), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: r1d'
!allocate(latitude(ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: latitude'
!allocate(sline2d(n1d, ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: sline2d'
!allocate(opalbar2d(n1d, ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: opalbar2d'
!allocate(vel_r2d(n1d, ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: vel_r2d'
!allocate(vel_theta2d(n1d, ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: vel_theta2d'
!allocate(temperature2d(n1d, ntheta), stat=err)
!   if(err.ne.0) stop 'error read_model2d_adm: temperature2d'
!!
!!----------------------create spatial grids-----------------------------
!!
!!radial grid: avoid grid point r=1 to avoid division by zero in e.g. wind upflow component
!r1d(1)=1.d0 + 1.d-8
!r1d(2)=1.d0 + 1.d-4
!del=log10(rmax)/(n1d-2)
!do i=3, n1d
!   r1d(i) = r1d(i-1)*10.d0**del
!enddo
!r1d(n1d)=rmax
!!
!!latitude grid: avoid values [0,pi/2,pi] to avoid divison by zero or negative roots
!min=1.d-4
!max=pi-1.d-5
!do i=1, ntheta
!   latitude(i) = min + (i-1)*(max-min)/(ntheta-1)
!enddo
!!
!!------------------------calculate 2d-grids-----------------------------
!!
!sline2d=0.d0
!opalbar2d=0.d0
!vel_r2d=0.d0
!vel_theta2d=0.d0
!!
!mu_min=1.d0-1.d-6
!!
!do i=1, n1d
!!
!   rad=r1d(i)
!!interpolate temperature: depends only on radius (from jo's program)
!   call find_index(rad, r1d_jo, n1d_jo, tindx)
!   temp1d = interpol_yp(r1d_jo(tindx-1), r1d_jo(tindx), t1d_jo(tindx-1), t1d_jo(tindx), rad)
!!
!   do j=1, ntheta
!!
!      mu=cos(latitude(j))
!!intersection point of photosphere with closed loop field line
!      mu_star=sqrt(1.d0-(1.d0-mu*mu)/rad)
!      mu_star=minval((/mu_star,mu_min/))
!!apex-radius
!      r_apex=1.d0/(1.d0-mu_star*mu_star)
!!shock-radius
!      mu_lower=0.d0
!      mu_upper=mu_star
!      call get_mu_shock(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
!      r_shock=r_apex*(1.d0-mu_shock*mu_shock)
!!
!!calculate wind upflow component everywhere
!      call component_w(rad, mu, rho, temp, vel, vel_r, vel_theta, &
!!                      0.75d0*teff, v_esc, v_inf, vmin*1.d5, rhoc_star)
!      rho=0.d0
!      vel=0.d0
!      vel_r=0.d0
!      vel_theta=0.d0 
!!
!      if(r_apex.lt.ralfven) then
!!in that regime: calculate cool downflow component and overwrite variables
!         call component_c(rad, mu, r_apex, rho, temp, vel, vel_r, vel_theta, &
!                          teff, v_esc, rhoc_star, delta) !instead of 0.75d0*teff
!         rho=rho/2.d0
!!
!         if(rad.gt.r_shock) then
!!in that regime: calculate post-shock component and overwrite variables
!!            call component_s(rad, mu, mu_star, mu_shock, rho, temp, vel, vel_r, vel_theta, &
!!                             0.75*teff, t_inf, v_esc, v_inf, vmin*1.d5, rhoc_star, delta)
!!             rho=0.d0
!!         endif
!      endif
!!
!!temperatures
!      temp=teff !0.75d0*teff
!!      temp=temp1d
!!
!!velocities
!      vel_r2d(i,j)=vel_r/vth_fiducial
!      vel_theta2d(i,j)=vel_theta/vth_fiducial
!!
!!use spherical symmetric departure coefficients
 !     call bvel3d(vmin*1.d5, vmax*1.d5, beta, r1d(i), 0.d0, 0.d0, vel_r, fdum1, fdum2, fdum3)
 !     vel_r=vel_r/vmax/1.d5
!!departure coefficients
 !     call depcoeff_petrenz(vel_r, b2, b3)
!!in lte
 !     b2=1.d0
 !     b3=1.d0
!!
!!integrated line opacity
!      opalbar2d(i,j) = opalbar_petrenz(rho, temp, b2, b3)
!!line source function
!      sline2d(i,j) = sline_petrenz(temp, b2, b3)/bnue(xnue0,teff)  !instead of 0.77d0*teff
!!temperature
!      temperature2d(i,j) = temp!
!
!!      write(*,'(8es20.8)') rad, mu, rho, temp, opalbar2d(i,j), sline2d(i,j), b2, b3
!!
!   enddo
!
!enddo
!!
!xic1=1.d0
!!
!!stop
!!
!!
!end subroutine read_model2d_adm
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output1d
!
use prog_type
use fund_const
use dime_modspec, only: n1d, r1d, velr1d, opalbar1d, opac1d, t1d, sline1d, scont1d
use params_modspec, only: teff, trad, xic1, xic2, vth_fiducial, rstar, vmicro, vrot, vmax, yhe, xlogg, vmin, lstar
use options_modspec, only: output_file
use hdf5
use mod_iline, only: xnue0, na
!
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
!
! ... local arrays
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
dims_r = (/ n1d /)
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
      call h5gcreate_f(file_id, 'dimensions', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'n1d', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, n1d, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
      call h5gcreate_f(file_id, 'input_parameters', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'teff', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'trad', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'xnue0', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'lstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vth_fiducial', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmicro', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vrot', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmax', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmin', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'yhe', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'logg', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------boundary condition-----------------------------
!
      call h5gcreate_f(file_id, 'bcondition', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'xic1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
      call h5gcreate_f(file_id, 'coordinates', group_id, err)
         call h5screate_simple_f(1, dims_r, dspace_id, err)
            call h5dcreate_f(group_id, 'r1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, r1d, dims_r, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------1d solution--------------------------------
!
     call h5gcreate_f(file_id, 'solution1d', group_id, err)
         call h5screate_simple_f(1, dims_r, dspace_id, err)
            call h5dcreate_f(group_id, 'sline1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, sline1d, dims_r, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'scont1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, scont1d, dims_r, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!--------------------------------1d model-------------------------------
!
     call h5gcreate_f(file_id, 'model1d', group_id, err)
         call h5screate_simple_f(1, dims_r, dspace_id, err)
            call h5dcreate_f(group_id, 't1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, t1d, dims_r, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'opac1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opac1d, dims_r, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'opalbar1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opalbar1d, dims_r, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'velr1d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velr1d, dims_r, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
write(*,'(a20,i20)') 'n1d', n1d
write(*,'(a20,es20.8)') 'teff', teff
write(*,'(a20,es20.8)') 'trad', trad
write(*,'(a20,es20.8)') 'xnue0', xnue0
write(*,'(a20,es20.8)') 'rstar', rstar
write(*,'(a20,es20.8)') 'vth_fiducial', vth_fiducial
write(*,'(a20,es20.8)') 'vmicro', vmicro
write(*,'(a20,es20.8)') 'vinf', vmax
write(*,'(a20,es20.8)') 'vmin', vmin
write(*,'(a20,es20.8)') 'vrot', vrot
write(*,'(a20,es20.8)') 'yhe', yhe
write(*,'(a20,es20.8)') 'logg', xlogg
write(*,'(a20,i20)') 'na', na
write(*,'(a20,es20.8)') 'xic1', xic1
write(*,'(a20,es20.8)') 'xic2', xic2
write(*,*)
write(*,*)
write(*,'(7a20)') 'r', 'T', 'opac', 'opalbar', 'velr', 'sline', 'scont'
do i=1, n1d
   write(*,'(7es20.8)') r1d(i), t1d(i), opac1d(i), opalbar1d(i), velr1d(i), sline1d(i), scont1d(i)
enddo
!
!
!
end subroutine output1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output3d
!
use prog_type
use fund_const
use dime_modspec, only: ndxmax, ndymax, ndzmax, x, y, z, t3d, opac3d, opalbar3d, &
                        velx3d, vely3d, velz3d, sline3d, scont3d, imask3d
use params_modspec, only: teff, trad, xic1, xic2, vth_fiducial, rstar, vmicro, vmax, xlogg, vrot, lstar
use options_modspec, only: output_file
use mod_iline, only: xnue0, na
use hdf5
!
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k
!
! ... local arrays
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
integer(hsize_t), dimension(3) :: dims_3d
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
dims_x = (/ ndxmax /)
dims_y = (/ ndymax /)
dims_z = (/ ndzmax /)
dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
      call h5gcreate_f(file_id, 'dimensions', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'ndymax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
      call h5gcreate_f(file_id, 'input_parameters', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'teff', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'trad', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'xnue0', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'lstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vth_fiducial', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmicro', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmax', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vrot', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'logg', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------boundary condition-----------------------------
!
      call h5gcreate_f(file_id, 'bcondition', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'xic1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
      call h5gcreate_f(file_id, 'coordinates', group_id, err)
         call h5screate_simple_f(1, dims_x, dspace_id, err)
            call h5dcreate_f(group_id, 'x', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_y, dspace_id, err)
            call h5dcreate_f(group_id, 'y', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_x, dspace_id, err)
            call h5dcreate_f(group_id, 'z', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------1d solution--------------------------------
!
     call h5gcreate_f(file_id, 'solution3d', group_id, err)
         call h5screate_simple_f(3, dims_3d, dspace_id, err)
            call h5dcreate_f(group_id, 'sline3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'scont3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!--------------------------------1d model-------------------------------
!
     call h5gcreate_f(file_id, 'model3d', group_id, err)
         call h5screate_simple_f(3, dims_3d, dspace_id, err)
            call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_integer, imask3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 't3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'opac3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'opalbar3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'velx3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'vely3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'velz3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
write(*,'(a20,i20)') 'ndxmax', ndxmax
write(*,'(a20,i20)') 'ndymax', ndymax
write(*,'(a20,i20)') 'ndzmax', ndzmax
write(*,'(a20,es20.8)') 'teff', teff
write(*,'(a20,es20.8)') 'trad', trad
write(*,'(a20,es20.8)') 'xnue0', xnue0
write(*,'(a20,es20.8)') 'rstar', rstar
write(*,'(a20,es20.8)') 'lstar', lstar
write(*,'(a20,es20.8)') 'logg', xlogg
write(*,'(a20,es20.8)') 'vth_fiducial', vth_fiducial
write(*,'(a20,es20.8)') 'vmicro', vmicro
write(*,'(a20,es20.8)') 'vrot', vrot
write(*,'(a20,i20)') 'na', na
write(*,'(a20,es20.8)') 'xic1', xic1
write(*,'(a20,es20.8)') 'xic2', xic2
write(*,*)
write(*,*)
write(*,'(10a16)') 'x(y=0,z=0)', 'T(x,y=0,z=0)', 'opac(x,y=0,z=0)', 'opalbar(x,y=0,z=0)', 'kappa', &
                  'velx(x,y=0,z=0)', 'vely(x,y=0,z=0)', 'velz(x,y=0,z=0)', 'sline(x,y=0,z=0)', 'scont(x,y=0,z=0)'
j=ndymax/2+1
k=ndzmax/2+1
do i=1, ndxmax
   write(*,'(10es16.6)') x(i), t3d(i,j,k), opac3d(i,j,k), opalbar3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k), &
                        vely3d(i,j,k), velz3d(i,j,k), sline3d(i,j,k), scont3d(i,j,k)
enddo
!
!
!
end subroutine output3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output3d_spc
!
use prog_type
use fund_const
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, t3d, opac3d, opalbar3d, &
                        velx3d, vely3d, velz3d, sline3d, scont3d, imask3d
use params_modspec, only: teff, trad, xic1, xic2, vth_fiducial, rstar, vmicro, vmax, sr, xlogg, yhe, vmin, lstar, vrot
use options_modspec, only: output_file
use mod_iline, only: xnue0, na
use hdf5
!
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k
!
! ... local arrays
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
dims_r = (/ nr /)
dims_theta = (/ ntheta /)
dims_phi = (/ nphi /)
dims_3d = (/ nr, ntheta, nphi /)
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
      call h5gcreate_f(file_id, 'dimensions', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nr, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'ntheta', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, ntheta, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'nphi', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nphi, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
      call h5gcreate_f(file_id, 'input_parameters', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'teff', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'logg', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'yhe', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'trad', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'xnue0', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'lstar', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vth_fiducial', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmicro', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmax', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vrot', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmin', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------boundary condition-----------------------------
!
      call h5gcreate_f(file_id, 'bcondition', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'xic1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'xic2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xic2, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
      call h5gcreate_f(file_id, 'coordinates', group_id, err)
         call h5screate_simple_f(1, dims_r, dspace_id, err)
            call h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, r, dims_r, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_theta, dspace_id, err)
            call h5dcreate_f(group_id, 'theta', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, theta, dims_theta, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_phi, dspace_id, err)
            call h5dcreate_f(group_id, 'phi', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, phi, dims_phi, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------1d solution--------------------------------
!
     call h5gcreate_f(file_id, 'solution3d', group_id, err)
         call h5screate_simple_f(3, dims_3d, dspace_id, err)
            call h5dcreate_f(group_id, 'sline3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'scont3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!     
     call h5gcreate_f(file_id, 'model3d', group_id, err)
         call h5screate_simple_f(3, dims_3d, dspace_id, err)
            call h5dcreate_f(group_id, 't3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            opac3d=opac3d/sr   !output in cgs
            call h5dcreate_f(group_id, 'opac3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            opac3d=opac3d*sr   
            opalbar3d=opalbar3d/sr   !output in cgs
            call h5dcreate_f(group_id, 'opalbar3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            opalbar3d=opalbar3d*sr
            call h5dcreate_f(group_id, 'velx3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'vely3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'velz3d', h5t_native_double, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
write(*,'(a20,i20)') 'nr', nr
write(*,'(a20,i20)') 'ntheta', ntheta
write(*,'(a20,i20)') 'nphi', nphi
write(*,'(a20,es20.8)') 'teff', teff
write(*,'(a20,es20.8)') 'trad', trad
write(*,'(a20,es20.8)') 'logg', xlogg
write(*,'(a20,es20.8)') 'yhe', yhe
write(*,'(a20,es20.8)') 'xnue0', xnue0
write(*,'(a20,es20.8)') 'rstar', rstar
write(*,'(a20,es20.8)') 'lstar', lstar
write(*,'(a20,es20.8)') 'vth_fiducial', vth_fiducial
write(*,'(a20,es20.8)') 'vmicro', vmicro
write(*,'(a20,es20.8)') 'vrot', vrot
write(*,'(a20,i20)') 'na', na
write(*,'(a20,es20.8)') 'xic1', xic1
write(*,'(a20,es20.8)') 'xic2', xic2
write(*,*)
write(*,*)
!
j=ntheta/2+1
k=1
!
write(*,'(a20,2f8.4)') 'at (theta,phi)= ', theta(j)*180./pi, phi(k)*180./pi
write(*,'(9a20)') 'r', 'T(r)', 'opac(r)', 'opalbar(r)', &
                  'velx(r)', 'vely(r)', 'velz(r)', 'sline(r)', 'scont(r)'
write(*,'(9a20)') '[rstar]', '[K]', '[1/sr]', '[1/s/sr]', &
                  '[cm/s]', '[cm/s]', '[cm/s]', '[erg/s/cm^2/Hz]', '[erg/s/cm^2/Hz]'
do i=1, nr
   write(*,'(9es20.8)') r(i), t3d(i,j,k), opac3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k), &
                        vely3d(i,j,k), velz3d(i,j,k), sline3d(i,j,k), scont3d(i,j,k)
enddo
!
!
!
end subroutine output3d_spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_nicowr3d
!
!model3d.h5 file in 3d (spherical coordinates)
!calculate the opacity for input file (spherical coordinates)
!calculate the source function for input file (spherical coordinates)
!
use prog_type
use fund_const, only: cgs_planck, cgs_clight, cgs_kb, rsu, xmsu, pi, zero, one, half
use options_modspec, only: indat_file, input_file, input_file2, input_mod
use dime_modspec, only: nr, ntheta, nphi, r, theta, phi, velx3d, vely3d, velz3d, opac3d, &
                        opalbar3d, sline3d, scont3d, imask3d, t3d, trad3d
use params_modspec, only: teff, trad, xic1, vth_fiducial, sr, rstar, vmax, vmicro, &
                          opt_opal, eps_line, lstar, xlogg, vrot, yhe, yhe_mass, hei, xmloss, unit_length
use mod_opal
use mod_lte
use hdf5
use mod_opacities, only: opac_thomson, opac_opal, get_opalbar, opalbar_lte_table
use mod_iline, only: alpha, kappa0, kline, xnue0, na, gl, flu, na, iline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: opt_scont, opt_sline, opt_opac, opt_vlat
integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
integer(i4b) :: ndxmax, ndymax, ndzmax, nr_modext, ntheta_modext, nphi_modext
real(dp) :: xcoord, ycoord, zcoord, fdum, rad, scont, sline, rho, opalbar, mdot, b2, b3, vel
real(dp) :: sint, cost, sinp, cosp
real(dp) :: kline_model, alpha_model, kappa0_model
real(dp) :: kcont, tfloor
real(dp) :: fdum1, fdum2, fdum3, fdum4
!
! ... local arrays
real(dp), dimension(:), allocatable :: r_modext, theta_modext, phi_modext
real(dp), dimension(:,:,:), allocatable :: velr3d_modext, velth3d_modext, velphi3d_modext, &
                                           rho3d_modext, t3d_modext, trad3d_modext
!
! ... local characters
character(len=500) :: cdum
!
! ... local functions
real(dp) :: dilfac, bnue, sline_depcoeff
!
! ... for hdf5 file
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_rad, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims3d_spc
!
! ... namelists
namelist / input_usr/ opt_sline, opt_opac, opt_scont, opt_vlat, kcont, tfloor
!
!
!only for testing if opacity subroutines are reasonable
!temperature
!fdum1 = 40.d3
!do i=1, 20
!   !density
!   fdum2 = 10.**(-15. + (i-1)*(-10.+15.)/19.)
!   fdum3 = get_opalbar(1, 1.d0, unit_length, yhe, hei, fdum1, vth_fiducial, xnue0, 1.d0, 1.d0, fdum2) !in 1/unit_length
!   fdum4 = get_opalbar(0, 1.d0, unit_length, yhe, hei, fdum1, vth_fiducial, xnue0, 1.d0, 1.d0, fdum2) !in 1/unit_length
!   write(*,*) fdum2, fdum3, fdum4
!enddo
!stop 'go on here'
!
write(*,*) '--------------------reading structure from input file--------------------------'
write(*,*) 'reading usr input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!-----------------------------------------------------------------------
!
write(*,*) '------------reading 3d atmospheric structure in spherical coordinates----------'
write(*,*) 'input_file: ', trim(input_file2)
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'nphi', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nphi_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)

call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims_phi=(/nphi_modext/)
dims3d_spc=(/ nr_modext, ntheta_modext, nphi_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(phi_modext(nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(rho3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(t3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(trad3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velr3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velth3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velphi3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
!normalize radius to radius 1
      r_modext=r_modext/r_modext(1) 
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'phi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, phi_modext, dims_phi, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'trad', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, trad3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
write(*,*) '-----------------calculating source functions and opacities--------------------'
write(*,*)
!
nr=nr_modext
ntheta=ntheta_modext
nphi=nphi_modext
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(trad3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'

r=r_modext
theta=theta_modext
phi=phi_modext
!
!
!do j=1, ntheta
!   write(*,*) j, theta(j)*180./pi
!enddo


do i=1, nr
   rad = r(i)   
   do j=1, ntheta
      do k=1, nphi
!        
         sint=sin(theta(j))
         cost=cos(theta(j))
         sinp=sin(phi(k))
         cosp=cos(phi(k))
!
!--------------------------temperatures---------------------------------
!
         t3d(i,j,k) = max(t3d_modext(i,j,k),tfloor)
         trad3d(i,j,k) = max(trad3d_modext(i,j,k),tfloor)
!
!------------------------velocity components----------------------------
!
         select case(opt_vlat)
         case(0)
            velx3d(i,j,k) = velr3d_modext(i,j,k)*sint*cosp
            vely3d(i,j,k) = velr3d_modext(i,j,k)*sint*sinp
            velz3d(i,j,k) = velr3d_modext(i,j,k)*cost
         case(1)         
            velx3d(i,j,k) = velr3d_modext(i,j,k)*sint*cosp + velth3d_modext(i,j,k)*cost*cosp-velphi3d_modext(i,j,k)*sinp
            vely3d(i,j,k) = velr3d_modext(i,j,k)*sint*sinp + velth3d_modext(i,j,k)*cost*sinp+velphi3d_modext(i,j,k)*cosp
            velz3d(i,j,k) = velr3d_modext(i,j,k)*cost - velth3d_modext(i,j,k)*sint
         case default
            stop 'error in read_model3d: opt_vlat not specified'
         end select
!
!---------------------calculate continuum opacity-----------------------
!         
         select case(opt_opac)
            case(0)
               opac3d(i,j,k) = zero
            case(1)
               !in 1/unit_length               
               opac3d(i,j,k) = opac_opal(kcont, yhe_mass, hei, log10(rho3d_modext(i,j,k)),log10(t3d_modext(i,j,k)), nrho_opal, ntemp_opal, rho_opal, temp_opal, kappa_opal)*unit_length
            case(2)
!               opac3d(i,j,k) = opac_thomson(yhe,hei,rho3d_modext(i,j,k),kcont) !in cgs
               opac3d(i,j,k) = opac_thomson(yhe,hei,rho3d_modext(i,j,k),kcont)*unit_length
            case default
               stop 'error in read_model3d: opt_opac not specified'
            end select

!         write(*,*) opac3d(i,j,k)/rho3d_modext(i,j,k)/unit_length, unit_length, yhe2
!         
!---------------------calculate frequency integrated line opacity-------
!
         b2 = one
         b3 = one
         opalbar3d(i,j,k) = get_opalbar(iline, kline, unit_length, yhe, hei, t3d(i,j,k), vth_fiducial, xnue0, b2, b3, rho3d_modext(i,j,k)) !in 1/unit_length
!         
!---------------------calculate continuum source function---------------
!         
         select case(opt_scont)
            case(0)
               scont3d(i,j,k) = zero
            case(1)
               scont3d(i,j,k) = bnue(xnue0,trad3d(i,j,k))
            case(2)
               if(rad.lt.one) rad=one
               dilfac = 0.5d0*(one-sqrt(one-one/rad**2))
               scont3d(i,j,k) = xic1*dilfac
!            case(3)
!               scont3d(i,j,k) = bnue(xnue0,trad3d(i,j,k))*bandwidth
            case default
              stop 'error in read_model3d: opt_scont not specified'
         end select
!         
!---------------------calculate line source function--------------------
!         
         select case(opt_sline)
            case(0)
               sline3d(i,j,k) = zero
            case(1)
               sline3d(i,j,k) = bnue(xnue0,trad3d(i,j,k))
            case(2)
               if(rad.lt.one) rad=one
               dilfac = half*(one-sqrt(one-one/rad**2))               
               sline3d(i,j,k) = xic1*dilfac
            case(3)
               b2 = one
               b3 = one
               sline3d(i,j,k) = sline_depcoeff(xnue0, t3d(i,j,k), b2, b3)
            case default
               stop 'error in read_model3d: opt_sline not specified'
         end select
!
!
      enddo
   enddo
enddo
!
!stop 'go on here'
!-----------------------------set velocities------------------------------
!
xic1=bnue(xnue0,trad)
!
!maximum velocity in global system
vmax=zero
do i=1, nr
   do j=1, ntheta
      do k=1, nphi
         vel=sqrt((velx3d(i,j,k))**2+(vely3d(i,j,k))**2+(velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

!stop 'go on in here'
!
end subroutine calc_model3d_nicowr3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_test_model_js_lh

  !First read in model created in the model-step

  use options_modspec
  use dime_modspec
  use params_modspec
  use hdf5
  use mod_opacities
  use mod_iline
  implicit none
  
  real(dp), dimension(:), allocatable :: r_modext, theta_modext, phi_modext
  real(dp), dimension(:,:,:), allocatable :: velr3d_modext, velth3d_modext, velphi3d_modext,rho3d_modext, t3d_modext  
  integer(i4b) :: nr_modext, ntheta_modext, nphi_modext, err, i,j,k
  real(dp) :: sint,cost,sinp,cosp
  
  !The below is for hdf5 format 
  integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
  integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
  integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_rad, dims_theta, dims_phi
  integer(hsize_t), dimension(3) :: dims3d_cac, dims3d_spc

  real(dp) :: bnue,dilfac,b2,b3 
    
!-----------------------------------------------------------------------
!
write(*,*) '------------reading 3d atmospheric structure in spherical coordinates----------'
write(*,*) 'input_file: ', trim(input_file2)
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'nphi', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nphi_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)

call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims_phi=(/nphi_modext/)
dims3d_spc=(/ nr_modext, ntheta_modext, nphi_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(phi_modext(nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(rho3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(t3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velr3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velth3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
allocate(velphi3d_modext(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in calc_model3d_spc3d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
!normalize radius to radius 1
      r_modext=r_modext/r_modext(1) 
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'phi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, phi_modext, dims_phi, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t3d_modext, dims3d_spc, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!  

!What happens here below? aloocate r, th, etc arrays used for LINE calcs
nr=nr_modext
ntheta=ntheta_modext
nphi=nphi_modext
!
!allocate arrays
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opac3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(opalbar3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(scont3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(sline3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velx3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(vely3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(velz3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
allocate(t3d(nr,ntheta,nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_spc3d'
r=r_modext
theta=theta_modext
phi=phi_modext

!JS-NOTE: opac, are EXTINCITION coefficients (1/cm) 
!BUT opalbar is NOT: it is FREQUENCY-INTEGRATED LINE EXTINCTION
!COEFFICIENT, I.E.  Hz/cm  
do k=1,nphi
   do j=1,ntheta 
      do i=1,nr
         opac3d(i,j,k) = 0.d0 !sr*0.34d0*rho3d_modext(i,j,k) 
         !opalbar(i,j,k) =
         opalbar3d(i,j,k) = sr*opalbar_model_kline(0.1d0, 2.d0, rho3d_modext(i,j,k), kline)
         !sscall opalbar_model_kline(yhe, hei, rho, kline)
         !iline is the line-identifier -- in src/mode_iline.f90 (or something similar)  
         !b2 = 1.d0
         !b3 = 1.d0  !these are departire coefficients (something dummy, depends on iline)  
         !opalbar3d(i,j,k) = get_opalbar(iline, kline, sr, yhe, hei, t3d_modext(i,j,k), &
         !     vth_fiducial, xnue0, b2, b3, rho3d_modext(i,j,k))/sr   !in cgs
         !This is a Levin-function for parameterised extinctions
         scont3d(i,j,k) = 0.d0
         sline3d(i,j,k) = bnue(xnue0,t3d_modext(i,j,k))*dilfac(r(i))
         t3d(i,j,k) = t3d_modext(i,j,k)
         !Express velocity vectors in x,y,z, NOT in input r,th,phi 
         sint = sin(theta(j))
         cost = cos(theta(j))
         sinp = sin(phi(k)) 
         cosp = cos(phi(k)) 
         velx3d(i,j,k)=velr3d_modext(i,j,k)*sint*cosp + velth3d_modext(i,j,k)*cost*cosp-velphi3d_modext(i,j,k)*sinp
         vely3d(i,j,k)=velr3d_modext(i,j,k)*sint*sinp + velth3d_modext(i,j,k)*cost*sinp+velphi3d_modext(i,j,k)*cosp
         velz3d(i,j,k)=velr3d_modext(i,j,k)*cost - velth3d_modext(i,j,k)*sint
!         print*,sr/rsu,r(i)
      enddo      
!      stop 'test...'
   enddo
enddo

!First lower boundary condition for intensity
!xic1 = bnue(xnue0,40000.d0)
!t3d(i,j,k) = t3d_modext(i,j,k)
!call diffus(xnue0, t_iim1, t_ii, t_iip1, r_iim1, r_ii, r_iip1, xic1, xic2)
call diffus(xnue0, t3d(1,1,1), t3d(1,1,1), t3d(2,1,1), r(1), r(1), r(2), xic1, xic2)  
xic2 = -xic1 !xic2/opac3d(1,1,1)
!print*,xic1,xic2,t3d(1,1,1),t3d(2,1,1),r(1),r(2),bnue(xnue0,40000.d0)
!stop 'test...'
!in spec.f90 there is option that places xic1 at an array in space and frequency
!(no angle included yet)
!there is a variable xic2, not yet really used though (or debugged etc)...  


end subroutine calc_test_model_js_lh
