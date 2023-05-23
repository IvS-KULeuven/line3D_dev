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
   use omp_lib
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
      call read_model1db
      call output1d
    case(19)
      call calc_model3d_amrvac
      call output3d_spc
    case(20)
      call calc_model3d_costum
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
   use mod_iline, only: iline, get_iline, na, kline, alpha, kappa0
!
!
   implicit none
!
! ... local scalars
   real(dp) :: fdum
!
! ... local scalars
! LP
   integer(i4b) :: iotstat
!
! ... local characters
! LP
   character(len=300) :: input_arg
!
! ... namelist
   namelist / input_options / input_file, input_file2, output_file, input_mod
   namelist / input_model / teff, trad, xlogg, rstar, rmax, tmin, xmloss, vmin, vmax, vrot, &
      vmicro, vth_fiducial, beta, yhe, hei
   namelist / input_line / iline, eps_line, kline, alpha, kappa0
!
!-----------------------------------------------------------------------
!
! changes : LP updating the io for the innput nml file reading
   call get_command_argument(number=1, value=input_arg, status=iotstat)
! print*,inputstat

   IF (iotstat.EQ.0) then
      ! do somethind
      write(*,*) 'Priceeding with input argument: ', trim(input_arg)
      indat_file = trim(input_arg)
      print*, indat_file
   else
      write(*,*) '----------------------------read input-----------------------------------------'
      write(*,*) 'input file name (*.nml) to define model-atmosphere'
      read(*,*) indat_file
      write(*,*) 'reading input from file: ', trim(indat_file)
      write(*,*)
      ! do somtheing else
   endif
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
!
!calculate mass fraction where yhe = n_he/n_h = yhe/(1-yhe), where yhe=n_he/(n_h+n_he)
   yhe_mass = one / (one + one/four/yhe)
   call get_iline(iline, yhe_mass)
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
!opal tables
   call get_opal_table()
!
!lte tables
   if(iline.eq.0) then
      call get_lte_table()
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
subroutine calc_model3d_amrvac
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
      real(dp) :: dilfac, bnue, sline_depcoeff, cpu_time_start, cpu_time_finish
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
      write(*,'(4(A, I4))')  'Nr = ', nr, "  Nt = ", ntheta, '  Np = ', nphi, '  Tot = ', nr*ntheta*nphi
      call cpu_time(cpu_time_finish)
      !$OMP DO SCHEDULE(GUIDED)
      do i=1, nr
         rad = r(i)
         call cpu_time(cpu_time_start)
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
         call cpu_time(cpu_time_finish)
         write(*,"(F6.2, A, F10.2, A, F6.3)") (100.0d0*i*j*k)/(nr*ntheta*nphi), &
         "%   ETC(h) =", (nr*ntheta*nphi - i*j*k)*(cpu_time_finish-cpu_time_start)/60/60, &
                     "   ET(s) =", cpu_time_finish-cpu_time_start
      enddo
      !$OMP END DO
   !
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
   end subroutine calc_model3d_amrvac
   !
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !
   subroutine calc_model3d_costum
   
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
   
   end subroutine calc_model3d_costum
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
