!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3d_spc
!
   use prog_type
   use fund_const
   use options_spec, only: input_file, opt_incl_gdark, opt_incl_sdist
   use dime_model3d, only: nr_spc, ntheta_spc, nphi_spc, r_spc, theta_spc, phi_spc, t3d, vth3d, velx3d, vely3d, velz3d, &
      opac3d, opalbar3d, scont3d, sline3d, imask3d
   use params_spec, only: teff, trad, tmin, xnue0, rstar, xic1, xic2, vth_fiducial, sr, na, vmax, xlogg, yhe, vmin, lstar, vrot
   use params_model, only: vth_fiducial_model, vmicro_model, vmax_model
   use hdf5
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   integer(i4b) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2
   real(dp) :: dum, xmloss_cgs
   real(dp) :: r_scalar, theta_scalar, phi_scalar, sint, cost, sinp, cosp
   real(dp) :: x1, xp, x2, y1, yp, y2, z1, zp, z2, &
      rada, radb, radc, radd, rade, radf, radg, radh, radp, &
      vala, valb, valc, vald, vale, valf, valg, valh, valp
   real(dp) :: velr, velth, velphi, rho
!
! ... local logicals
   logical :: linfo_phot, linfo_max, expol, llogx, llogy, llogz, llogf, lr2, lcheck
!
! ... for hdf5
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1), parameter :: dims_scalars = (/ 1 /)
   integer(hsize_t), dimension(1) :: dims_radius, dims_theta, dims_phi
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local function
   real(dp) :: opalbar_model_hamann
!
   write(*,*) '------reading 3d model atmosphere (spc) from-----------'
   write(*,*) 'file-name: ', input_file
!
!---------------------------read h5-file--------------------------------
!
!open hdf5-interface
   call h5open_f(err)
   call h5fopen_f(trim(input_file), h5f_acc_rdonly_f, file_id, err)
!
!--------------------------read dimensions------------------------------
!
   call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', dset_id, err)
   call h5aread_f(dset_id, h5t_native_integer, nr_spc, dims_scalars, err)
   call h5aclose_f(dset_id, err)
   call h5aopen_f(group_id, 'ntheta', dset_id, err)
   call h5aread_f(dset_id, h5t_native_integer, ntheta_spc, dims_scalars, err)
   call h5aclose_f(dset_id, err)
   call h5aopen_f(group_id, 'nphi', dset_id, err)
   call h5aread_f(dset_id, h5t_native_integer, nphi_spc, dims_scalars, err)
   call h5aclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!
!set dimension-arrays
   dims_radius=(/ nr_spc /)
   dims_theta=(/ ntheta_spc /)
   dims_phi=(/ nphi_spc /)
   dims_3d=(/nr_spc, ntheta_spc, nphi_spc /)
!
!----------------------------input parameters---------------------------
!
   call h5gopen_f(file_id, 'input_parameters', group_id, err)
   call h5aopen_f(group_id, 'teff', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'yhe', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'trad', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmin', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'xnue0', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'rstar', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'lstar', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark) then
      write(*,*) 'error in read_model3d_spc: opt_incl_gdark only possible with given lstar'
      write(*,*) '    but lstar not specified in h5-file'
      write(*,*) '    => change option opt_incl_gdark, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'logg', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark.or.err.eq.-1.and.opt_incl_sdist) then
      write(*,*) 'error in read_model3d_spc: opt_incl_gdark/opt_incl_sdist only possible with given logg'
      write(*,*) '    but loggr not specified in h5-file'
      write(*,*) '    => change options opt_incl_gdark/opt_incl_sdist, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vrot', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vth_fiducial_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmicro', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmicro_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'na', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
   sr=rstar*rsu
   vrot=vrot/1.d5 !in km/s
!
!------------------------boundary condition-----------------------------
!
   call h5gopen_f(file_id, 'bcondition', group_id, err)
   call h5aopen_f(group_id, 'xic1', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aexists_f(group_id, 'xic2', lcheck, err)
   if(lcheck) then
      call h5aopen_f(group_id, 'xic2', attr_id, err)
      call h5aread_f(attr_id, h5t_native_double, xic2, dims_scalars, err)
      call h5aclose_f(attr_id, err)
   else
      xic2=zero
   endif
   call h5gclose_f(group_id, err)
!
!-----------------------coordinates-------------------------------------
!
!allocate arrays
   allocate(r_spc(nr_spc), stat=err)
   allocate(theta_spc(ntheta_spc), stat=err)
   allocate(phi_spc(nphi_spc), stat=err)

   allocate(t3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(sline3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(scont3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(velx3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(vely3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(velz3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(opac3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(opalbar3d(nr_spc,ntheta_spc,nphi_spc), stat=err)
   allocate(vth3d(nr_spc,ntheta_spc,nphi_spc), stat=err)

   call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, r_spc, dims_radius, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, theta_spc, dims_theta, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'phi', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, phi_spc, dims_phi, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------3d solution--------------------------------
!
   call h5gopen_f(file_id, 'solution3d', group_id, err)
   call h5dopen_f(group_id, 'sline3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'scont3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
   call h5gopen_f(file_id, 'model3d', group_id, err)
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
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!-----------------------------------------------------------------------
!
!transform opacity to (eventually) different vth_fiducial, and to 1/rstar
   opalbar3d=opalbar3d*vth_fiducial_model/vth_fiducial * sr
   opac3d=opac3d * sr
!
!measure velocities in vth_fiducial
   velx3d=velx3d/vth_fiducial
   vely3d=vely3d/vth_fiducial
   velz3d=velz3d/vth_fiducial
!
!minimum temperature
   tmin=minval(t3d)
!
!maximum terminal velocity (to extend frequency range)
   vmax_model=vmax!*1.d5
   i=nr_spc
   do j=1, ntheta_spc
      do k=1, nphi_spc
         velr = sqrt(velx3d(i,j,k)**2 + vely3d(i,j,k)**2 + velz3d(i,j,k)**2)
!       write(*,*) velr, vth_fiducial, velr*vth_fiducial, vmax_model
         if(velr*vth_fiducial.gt.vmax_model) then
            vmax_model = velr*vth_fiducial
         endif
      enddo
   enddo

!write(*,*) vmax
!write(*,*) vmax_model
!stop 'go on here'
!
!open(1, file='TRASH/test2.dat', form='formatted')
!do i=1, nr_spc
!   write(1,'(10es20.8)') r_spc(i), opalbar3d(i,1,1), sline3d(i,1,1), velz3d(i,1,1)
!enddo
!close(1)
!stop 'go on in read_model3d_spc'
!
!
end subroutine read_model3d_spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine transform_model3d
!
!transform 3d-arrays to magnetic system, and then back to
!   formal solver 3d-arrays
   use prog_type
   use fund_const, only: cgs_clight, pi
   use params_mod3d, only: obliquity_mod3d, vthfiducial_mod3d, xnue0_mod3d
   use params_spec, only: obliquity, sr, vth_fiducial
   use dime_model1d, only: n1d
   use dime_model3d, only: opac3d, opalbar3d, velx3d, vely3d, velz3d, &
      t3d, scont3d, sline3d, x, y, z, &
      ndxmax, ndymax, ndzmax
   use mod_spectrum, only: rmin, rmax
   use mod_interp3d, only: get_xyz_indx, get_xyz_values1, get_xyz_values2, trilin_complete
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   integer(i4b) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2
   real(dp) :: xm, ym, zm, x1, x2, y1, y2, z1, z2
   real(dp) :: rada, radb, radc, radd, rade, radf, radg, radh, radp
   real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh
   real(dp) :: velx_m, vely_m, velz_m, velx, vely, velz, opal, opac, sline, scont, temp
   real(dp) :: del_obliquity
   real(dp) :: deldop, deldop_mod3d
!
! ... local logicals
   logical :: linfo_phot, linfo_max, expol
   logical :: llogx, llogy, llogz, llogf, lr2
!
! ... local arrays
   real(dp), dimension(:,:,:), allocatable :: opac3d_tmp, opalbar3d_tmp, velx3d_tmp, vely3d_tmp, &
      velz3d_tmp, t3d_tmp, scont3d_tmp, sline3d_tmp
!
!
!
   write(*,*) '-----------transforming 3d-coordinate system-----------'
   write(*,*)
   write(*,*) 'obliquity used in 3d-grid: ', obliquity_mod3d
   write(*,*) 'obliquity used here: ', obliquity
!
!-----------------------------------------------------------------------
!
   allocate(opac3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: opac3d_tmp'
!
   allocate(opalbar3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: opalbar3d_tmp'
!
   allocate(velx3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: velx3d_tmp'
!
   allocate(vely3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: vely3d_tmp'
!
   allocate(velz3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: velz3d_tmp'
!
   allocate(t3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: t3d_tmp'
!
   allocate(scont3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: scont3d_tmp'
!
   allocate(sline3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in transform_model3d: sline3d_tmp'
!
!-----------------------------------------------------------------------
!
   del_obliquity = obliquity - obliquity_mod3d

   if(del_obliquity.gt.pi) then
      stop 'error in transform_model3d: obliquities probably not in rad'
   endif
!
!-----------------------------------------------------------------------
!
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
!transform coordinates to magnetic system (del_oblquity in rad)
            call transform_sigr_sigm(x(i), y(j), z(k), xm, ym, zm, del_obliquity)
!
            call info_region(xm, ym, zm, rmin, rmax, linfo_phot, linfo_max)

            if(linfo_phot.and.linfo_max) then
!
!search for indices of a cube surrounding the point of interest (or neighbouring cubes for extrapolation)
               call get_xyz_indx(xm, ym, zm, x, y, z, ndxmax, ndymax, ndzmax, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol, rmin, rmax)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
               call get_xyz_values1(xm, ym, zm, x, y, z, ndxmax, ndymax, ndzmax, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  x1, x2, y1, y2, z1, z2, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  llogx, llogy, llogz, expol)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, velx3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
               llogx=.false.
               llogy=.false.
               llogz=.false.
               llogf=.false.
               lr2=.false.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, velx_m)
!
!get vely on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, vely3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
               llogx=.false.
               llogy=.false.
               llogz=.false.
               llogf=.false.
               lr2=.false.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, vely_m)
!
!get velz on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, velz3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
               llogx=.false.
               llogy=.false.
               llogz=.false.
               llogf=.false.
               lr2=.false.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, velz_m)
!
!backtransformation
               call transform_sigm_sigr(velx_m, vely_m, velz_m, velx, vely, velz, del_obliquity)
!
!-------------------interpolation of continuum opacity------------------
!
!get line-opacities on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, opac3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!here, can apply only lin-lin interpolation by setting llogx, llogy, llogz, llogf=.false.
!            llogx=.false.
!            llogy=.false.
!            llogz=.false.
!            llogf=.false.
!here, can decide if values shall be interpolated by function*r^2
               lr2=.true.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, opac)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, opalbar3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!here, can apply only lin-lin interpolation by setting llogx, llogy, llogz, llogf=.false.
!            llogx=.false.
!            llogy=.false.
!            llogz=.false.
!            llogf=.false.
!here, can decide if values shall be interpolated by function*r^2
               lr2=.true.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, opal)
!
!---------------------continuum source function-------------------------
!
!get sline on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, scont3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!            llogx=.false.
!            llogy=.false.
!            llogz=.false.
!            llogf=.false.
               lr2=.true.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, scont)
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, sline3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!            llogx=.false.
!            llogy=.false.
!            llogz=.false.
!            llogf=.false.
               lr2=.true.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, sline)
!
!
!------------------------interpolation of temperature-------------------
!
!get sline on cube vertices
               call get_xyz_values2(ndxmax, ndymax, ndzmax, t3d, &
                  indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!            llogx=.false.
!            llogy=.false.
!            llogz=.false.
!            llogf=.false.
               lr2=.true.
!actual interpolation
               call trilin_complete(xm, ym, zm, x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, .true., llogx, llogy, llogz, llogf, lr2, temp)
!
            else
               opac = 0.d0
               opal = 0.d0
               sline = 0.d0
               scont = 0.d0
               velx = 0.d0
               vely = 0.d0
               velz = 0.d0
               temp = 0.d0
            endif
!

!         write(*,*) opalbar3d(i,j,k), opal
            opac3d_tmp(i,j,k)=opac
            opalbar3d_tmp(i,j,k)=opal
            sline3d_tmp(i,j,k)=sline
            scont3d_tmp(i,j,k)=scont
            velx3d_tmp(i,j,k)=velx
            vely3d_tmp(i,j,k)=vely
            velz3d_tmp(i,j,k)=velz
            t3d_tmp(i,j,k)=temp
!
         enddo
      enddo
   enddo
!
!****************************debug outputs******************************
!
!open(1, file='trash/dime.dat', form='formatted')
!   write(1,*) ndxmax, ndymax, ndzmax
!close(1)
!
!for same units as original output
!deldop_mod3d = vth_fiducial_mod3d*xnue0_mod3d/cgs_clight
!deldop = vth_fiducial*xnue0_mod3d/cgs_clight
!opalbar3d_tmp = opalbar3d_tmp/deldop_mod3d
!opalbar3d_tmp = opalbar3d_tmp*deldop
!opalbar3d_tmp = opalbar3d_tmp/sr
!
!velx3d_tmp=velx3d_tmp*vth_fiducial
!vely3d_tmp=vely3d_tmp*vth_fiducial
!velz3d_tmp=velz3d_tmp*vth_fiducial
!write(*,*) velz3d_tmp
!!
!open(1, file='trash/model_opalbar3d.dat', form='unformatted')
!   write(1) opalbar3d_tmp
!close(1)
!open(1, file='trash/model_velx3d.dat', form='unformatted')
!   write(1) velx3d_tmp
!close(1)
!open(1, file='trash/model_vely3d.dat', form='unformatted')
!   write(1) vely3d_tmp
!close(1)
!open(1, file='trash/model_velz3d.dat', form='unformatted')
!   write(1) velz3d_tmp
!close(1)
!
!open(1, file='trash/opac3d.dat', form='unformatted')
!   write(1) opac3d
!close(1)
!
!***********************************************************************
!
   opac3d=opac3d_tmp
   opalbar3d=opalbar3d_tmp
   sline3d=sline3d_tmp
   scont3d=scont3d_tmp
   velx3d=velx3d_tmp
   vely3d=vely3d_tmp
   velz3d=velz3d_tmp
   t3d=t3d_tmp
!
!
!
end subroutine transform_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_params
!
   use prog_type
   use fund_const, only: pi
   use params_mod3d, only: sr_mod3d, teff_mod3d, xnue0_mod3d, beta_mod3d, vmax_mod3d, vmin_mod3d, &
      xlogg_mod3d, yhe_mod3d, hei_mod3d, ralfven_mod3d, chiinf_mod3d, obliquity_mod3d
   use params_spec, only: sr, teff, xnue0, beta, vmax, vmin, xlogg, yhe, hei, ralfven, chi_inf
!
   implicit none
!
   if(xnue0.ne.xnue0_mod3d) stop 'error check_params: xnue0_mod3d ne xnue0'
   if(abs(sr-sr_mod3d)/sr_mod3d.gt.1.d-7) stop 'error check_params: sr_mod3d ne sr'
   if(teff.ne.teff_mod3d) stop 'error check_params: teff_mod3d ne teff'
   if(beta.ne.beta_mod3d) stop 'error check_params: beta_mod3d ne beta'
   if(vmax.ne.vmax_mod3d) stop 'error check_params: vmax_mod3d ne vmax'
   if(vmin.ne.vmin_mod3d) stop 'error check_params: vmin_mod3d ne vmin'
!
   if(xlogg.ne.xlogg_mod3d) stop 'error check_params: xlogg_mod3d ne xlogg'
   if(yhe.ne.yhe_mod3d) stop 'error check_params: yhe_mod3d ne yhe'
   if(hei.ne.hei_mod3d) stop 'error check_params: hei_mod3d ne hei'
   if(chi_inf.ne.chiinf_mod3d) stop 'error check_params: chiinf_mod3d ne chi_inf'
   if(ralfven.ne.ralfven_mod3d) stop 'error check_params: ralfven_mod3d ne ralfven'
!
   if(obliquity_mod3d.gt.pi) stop 'error check_params: obliquity from 3d-input not in radians'
!
end subroutine check_params
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read in input parameter-----------------------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: zero
   use options_spec, only: input_file, output_dir, input_mod, opt_photprof, opt_obsdir_read, opt_surface, opt_int2d, opt_incl_gdark, opt_incl_sdist
   use params_spec, only: vth_fiducial, vmicro, vrot
   use mod_spectrum, only: rmin, rmax
   use mod_surfb, only: nsurfb, xobss_surface, alphas_surface, gammas_surface
   use dime_spec, only: nalpha, ngamma
!
   implicit none
!
! ... local scalars
   integer(i4b) :: err
!
! ... local scalars
! LP
   integer(i4b) :: iotstat
!
! ... local characters
!LP
   character(len=300) :: input_arg
   character(len=100) :: fname
!
! ... local arrays
   real(dp), dimension(100) :: alpha_surface=zero, gamma_surface=zero, xobs_surface=zero
!
! ... namelist
   namelist / input_options / input_mod, input_file, output_dir, opt_photprof, opt_obsdir_read, opt_surface, opt_int2d, opt_incl_gdark, opt_incl_sdist, nalpha, ngamma
   namelist / input_model / vrot, vth_fiducial, vmicro, rmin, rmax
   namelist / input_surface / nsurfb, alpha_surface, gamma_surface, xobs_surface
!
! ... local functions
!
!
!----------------------read input/output directories--------------------
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
!-----------------------------defaults----------------------------------
!
   output_dir = './outputFILES'
   opt_obsdir_read = .true.
   opt_surface = .false.
   opt_int2d = .false.
   opt_incl_gdark = .false.
   opt_incl_sdist = .false.
   nalpha = 1
   ngamma = 1
!
   vrot = zero
   vth_fiducial = 1.d2
   vmicro = 1.d2
   rmin = 1.d0
   rmax = 10.d0
!
   nsurfb = 1
   alpha_surface = zero
   gamma_surface = zero
   xobs_surface = zero
!
!-----------------------------------------------------------------------
!
   open(1, file=trim(fname), status='old', form='formatted')
!
   read(1, nml=input_options)
!
   read(1, nml=input_model)
   vth_fiducial = vth_fiducial * 1.d5
   vmicro = vmicro * 1.d5
   !
   read(1, nml=input_surface)
   !
   allocate(xobss_surface(nsurfb), stat=err)
   allocate(alphas_surface(nsurfb), stat=err)
   allocate(gammas_surface(nsurfb), stat=err)
   xobss_surface = xobs_surface(1:nsurfb)
   alphas_surface = alpha_surface(1:nsurfb)
   gammas_surface = gamma_surface(1:nsurfb)
!
   close(1)
!
!
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model1d
!
   use prog_type
   use fund_const
   use options_spec, only: input_file, opt_incl_gdark, opt_incl_sdist
   use dime_model1d, only: n1d, r1d, vth1d, &
      velr1d, acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d, &
      opalbar1d, acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d, &
      sline1d, acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d, &
      opac1d, acoeff_opac1d, bcoeff_opac1d, ccoeff_opac1d, dcoeff_opac1d, &
      scont1d, acoeff_scont1d, bcoeff_scont1d, ccoeff_scont1d, dcoeff_scont1d, &
      t1d, acoeff_t1d, bcoeff_t1d, ccoeff_t1d, dcoeff_t1d
   use params_spec, only: teff, trad, tmin, xnue0, rstar, xic1, vth_fiducial, sr, na, vmax, yhe, xlogg, vmin, lstar, vrot
   use params_model, only: vth_fiducial_model, vmicro_model
   use hdf5
   use mod_interp1d, only: cube_mono_coeff
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
! ... for hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r
!
!-----------------------------------------------------------------------
!
   write(*,*) '-------------------------reading model atmosphere from-------------------------'
   write(*,*) 'file ', trim(input_file)
   write(*,*)
!
!-----------------------------------------------------------------------
!
   call h5open_f(err)
   call h5fopen_f(trim(input_file), h5f_acc_rdwr_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
   call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'n1d', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, n1d, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
   dims_r = (/ n1d /)
!
!----------------------------input parameters---------------------------
!
   call h5gopen_f(file_id, 'input_parameters', group_id, err)
   call h5aopen_f(group_id, 'teff', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'trad', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmin', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'yhe', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'logg', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark.or.err.eq.-1.and.opt_incl_sdist) then
      write(*,*) 'error in read_model3d_spc: opt_incl_gdark/opt_incl_sdist only possible with given logg'
      write(*,*) '    but loggr not specified in h5-file'
      write(*,*) '    => change options opt_incl_gdark/opt_incl_sdist, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'xnue0', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'rstar', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'lstar', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark) then
      write(*,*) 'error in read_model1d: opt_incl_gdark only possible with given lstar'
      write(*,*) '    but lstar not specified in h5-file'
      write(*,*) '    => change option opt_incl_gdark, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vth_fiducial_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmicro', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmicro_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vrot', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'na', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
   sr=rstar*rsu
   vrot=vrot/1.d5 !in km/s
!
!------------------------boundary condition-----------------------------
!
   call h5gopen_f(file_id, 'bcondition', group_id, err)
   call h5aopen_f(group_id, 'xic1', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
   allocate(r1d(n1d))
   allocate(sline1d(n1d))
   allocate(scont1d(n1d))
   allocate(velr1d(n1d))
   allocate(t1d(n1d))
   allocate(opac1d(n1d))
   allocate(opalbar1d(n1d))
   allocate(vth1d(n1d))

   call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, r1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------1d solution--------------------------------
!

   call h5gopen_f(file_id, 'solution1d', group_id, err)
   call h5dopen_f(group_id, 'sline1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, sline1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'scont1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, scont1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!--------------------------------1d model-------------------------------
!
   call h5gopen_f(file_id, 'model1d', group_id, err)
   call h5dopen_f(group_id, 't1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, t1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'opac1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, opac1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'opalbar1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, opalbar1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr1d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, velr1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!transform opacity to (eventually) different vth_fiducial, and to 1/rstar
   opalbar1d = opalbar1d*vth_fiducial_model/vth_fiducial * sr
   opac1d = opac1d*sr
!
!measure velocities in vth_fiducial
   velr1d=velr1d/vth_fiducial
!
!for thin continuum models
!opac1d=0.d0
!scont1d=0.d0
!
   tmin=minval(t1d)
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!--------prepare coefficients for cubic spline interpolation------------
!
   write(*,*) '-------preparing cubic spline coefficients-------------'
   write(*,*)

   allocate(acoeff_sline1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_sline1d'
   allocate(acoeff_scont1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_scont1d'
   allocate(acoeff_opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_opalbar1d'
   allocate(acoeff_opac1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_opac1d'
   allocate(acoeff_velr1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_velr1d'
   allocate(acoeff_t1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: acoeff_t1d'
   allocate(bcoeff_sline1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_sline1d'
   allocate(bcoeff_scont1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_scont1d'
   allocate(bcoeff_opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_opalbar1d'
   allocate(bcoeff_opac1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_opac1d'
   allocate(bcoeff_velr1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_velr1d'
   allocate(bcoeff_t1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: bcoeff_t1d'
   allocate(ccoeff_sline1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_sline1d'
   allocate(ccoeff_scont1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_scont1d'
   allocate(ccoeff_opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_opalbar1d'
   allocate(ccoeff_opac1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_opac1d'
   allocate(ccoeff_velr1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_velr1d'
   allocate(ccoeff_t1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: ccoeff_t1d'
   allocate(dcoeff_sline1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_sline1d'
   allocate(dcoeff_scont1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_scont1d'
   allocate(dcoeff_opalbar1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_opalbar1d'
   allocate(dcoeff_opac1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_opac1d'
   allocate(dcoeff_velr1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_velr1d'
   allocate(dcoeff_t1d(n1d), stat=err)
   if(err.ne.0) stop 'error read_model1d: dcoeff_t1d'
!
!monotonic spline interpolation
   call cube_mono_coeff(n1d, r1d, sline1d, acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d)
   call cube_mono_coeff(n1d, r1d, scont1d, acoeff_scont1d, bcoeff_scont1d, ccoeff_scont1d, dcoeff_scont1d)
   call cube_mono_coeff(n1d, r1d, opalbar1d, acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d)
   call cube_mono_coeff(n1d, r1d, opac1d, acoeff_opac1d, bcoeff_opac1d, ccoeff_opac1d, dcoeff_opac1d)
   call cube_mono_coeff(n1d, r1d, velr1d, acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d)
   call cube_mono_coeff(n1d, r1d, t1d, acoeff_t1d, bcoeff_t1d, ccoeff_t1d, dcoeff_t1d)
!
!or normal splines:
!call spline_coeff(n1d, r1d, sline1d, acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d)
!call spline_coeff(n1d, r1d, scont1d, acoeff_scont1d, bcoeff_scont1d, ccoeff_scont1d, dcoeff_scont1d)
!call spline_coeff(n1d, r1d, opalbar1d, acoeff_opalbar1d, bcoeff_opalbar1d, ccoeff_opalbar1d, dcoeff_opalbar1d)
!call spline_coeff(n1d, r1d, opac1d, acoeff_opac1d, bcoeff_opac1d, ccoeff_opac1d, dcoeff_opac1d)
!call spline_coeff(n1d, r1d, velr1d, acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d)
!call spline_coeff(n1d, r1d, t1d, acoeff_t1d, bcoeff_t1d, ccoeff_t1d, dcoeff_t1d)
!
!write(*,*) sr/rsu
!write(*,*) opalbar1d/sr
!stop

end subroutine read_model1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3d
!
   use prog_type
   use fund_const
   use options_spec, only: input_file, opt_incl_gdark, opt_incl_sdist
   use dime_model3d, only: ndxmax, ndymax, ndzmax, x, y, z, t3d, vth3d, velx3d, vely3d, velz3d, &
      opac3d, opalbar3d, scont3d, sline3d, imask3d
   use params_spec, only: teff, trad, tmin, xnue0, rstar, xlogg, xic1, vth_fiducial, sr, na, vmax, lstar, vrot
   use params_model, only: vth_fiducial_model, vmicro_model
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
! ... for hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d
!
!-----------------------------------------------------------------------
!
   write(*,*) '-------------------------reading model atmosphere from-------------------------'
   write(*,*) 'file ', trim(input_file)
   write(*,*)
!
!-----------------------------------------------------------------------
!
   call h5open_f(err)
   call h5fopen_f(trim(input_file), h5f_acc_rdwr_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
   call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'ndxmax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ndymax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, ndymax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ndzmax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
!----------------------------input parameters---------------------------
!
   call h5gopen_f(file_id, 'input_parameters', group_id, err)
   call h5aopen_f(group_id, 'teff', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'trad', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmax', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'xnue0', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'rstar', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'lstar', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark) then
      write(*,*) 'error in read_model3d: opt_incl_gdark only possible with given lstar'
      write(*,*) '    but lstar not specified in h5-file'
      write(*,*) '    => change option opt_incl_gdark, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vth_fiducial_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vmicro', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vmicro_model, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'logg', attr_id, err)
   if(err.eq.-1.and.opt_incl_gdark.or.err.eq.-1.and.opt_incl_sdist) then
      write(*,*) 'error in read_model3d_spc: opt_incl_gdark/opt_incl_sdist only possible with given logg'
      write(*,*) '    but loggr not specified in h5-file'
      write(*,*) '    => change options opt_incl_gdark/opt_incl_sdist, or update modelspec-routine'
      stop
   endif
   call h5aread_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'vrot', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'na', attr_id, err)
   call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
   sr=rstar*rsu
   vrot=vrot/1.d5 !in km/s
!
!------------------------boundary condition-----------------------------
!
   call h5gopen_f(file_id, 'bcondition', group_id, err)
   call h5aopen_f(group_id, 'xic1', attr_id, err)
   call h5aread_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
   allocate(x(ndxmax), stat=err)
   allocate(y(ndymax), stat=err)
   allocate(z(ndzmax), stat=err)

   allocate(t3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(imask3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(sline3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(scont3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(velx3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(vely3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(velz3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(opac3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(opalbar3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(vth3d(ndxmax,ndymax,ndzmax), stat=err)

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
!----------------------------3d solution--------------------------------
!
   call h5gopen_f(file_id, 'solution3d', group_id, err)
   call h5dopen_f(group_id, 'sline3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'scont3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
!--------------------------------1d model-------------------------------
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
!
!transform opacity to (eventually) different vth_fiducial, and to 1/rstar
   opalbar3d=opalbar3d*vth_fiducial_model/vth_fiducial * sr
!
!measure velocities in vth_fiducial
   velx3d=velx3d/vth_fiducial
   vely3d=vely3d/vth_fiducial
   velz3d=velz3d/vth_fiducial
!
!for thin continuum models
   opac3d=0.d0
   scont3d=0.d0
!
!minimum temperature
   tmin=maxval(t3d)
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            if(imask3d(i,j,k).ne.0) then
               if(t3d(i,j,k).lt.tmin) tmin=t3d(i,j,k)
            endif
         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!open(1, file='TRASH/test1.dat', form='formatted')
!i=ndxmax/2+1
!j=ndymax/2+1
!do k=1, ndzmax
!   write(1,'(10es20.8)') z(k), opalbar3d(i,j,k), sline3d(i,j,k), velz3d(i,j,k)
!enddo
!close(1)
!stop 'go on in read_model3d'
!
end subroutine read_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_model1d
!
   use prog_type
   use dime_model1d, only: n1d, r1d, opac1d, t1d, vth1d, scont1d, &
      velr1d, opalbar1d, sline1d
   use params_spec, only: xic1, vmicro, xnue0, vth_fiducial, na, tmin, vmax, vrot, rstar, lstar, xlogg
!
   implicit none
!
! ... arguments
   integer(i4b) :: i
!
   write(*,*) '----------------------1d (radial) atmospheric structure------------------------'
   write(*,*)
!
   write(*,'(a20, es20.8)') 'I_c', xic1
   write(*,'(a20, es20.8)') 'xnue0', xnue0
   write(*,'(a20, es20.8)') 'tmin [K]', tmin
   write(*,'(a20, i20)') 'na', na
   write(*,'(a20, es20.8)') 'vmax [km/s]', vmax/1.d5
   write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro/1.d5
   write(*,'(a20, es20.8)') 'vth_fiducial [km/s]', vth_fiducial/1.d5
   write(*,'(a20, es20.8)') 'vrot [km/s]', vrot
   write(*,'(a20, es20.8)') 'rstar [r_sun]', rstar
   write(*,'(a20, es20.8)') 'lstar [l_sun]', lstar
   write(*,'(a20, es20.8)') 'logg ', xlogg
   write(*,*)

   write(*,'(a4, 7(a20))') '#', 'r [r_star]', 'opac [1/sr]', 'opalbar [1/sr]', 'velr [vth*]', 't [K]', 's_cont', 's_line'
   do i=1, n1d
      write(*,'(i4, 7(es20.8))')  i, r1d(i), opac1d(i), opalbar1d(i), velr1d(i), t1d(i), scont1d(i), sline1d(i)
   enddo
   write(*,*)
!
end subroutine print_model1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_model3d
!
   use prog_type
   use dime_model3d, only: ndxmax, ndymax, ndzmax, x, opac3d, opalbar3d, &
      velx3d, vely3d, velz3d, sline3d, scont3d, t3d
   use params_spec, only: xic1, vmicro, xnue0, vth_fiducial, na, tmin, vmax, vrot, lstar, xlogg, rstar, vrot, rstar
!
   implicit none
!
! ... arguments
   integer(i4b) :: i, j, k
!
   write(*,*) '--------------------3d atmospheric structure (along x-axis)--------------------'
   write(*,*)
!
   write(*,'(a20, es20.8)') 'I_c', xic1
   write(*,'(a20, es20.8)') 'xnue0', xnue0
   write(*,'(a20, es20.8)') 'tmin [K]', tmin
   write(*,'(a20, i20)') 'na', na
   write(*,'(a20, es20.8)') 'vmax [km/s]', vmax/1.d5
   write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro/1.d5
   write(*,'(a20, es20.8)') 'vth_fiducial [km/s]', vth_fiducial/1.d5
   write(*,'(a20, es20.8)') 'vrot [km/s]', vrot
   write(*,'(a20, es20.8)') 'rstar [r_sun]', rstar
   write(*,'(a20, es20.8)') 'lstar [l_sun]', lstar
   write(*,'(a20, es20.8)') 'logg ', xlogg
   write(*,*)
!
   j=ndymax/2+1
   k=ndzmax/2+1
!
   write(*,'(a4, 9(a20))') '#', 'x [r_star]', 'opac [1/cm]', 'opalbar [1/sr]', 'velx [vth*]', &
      'vely [vth*]', 'velz [vth*]', 't [K]', 's_cont', 's_line'
   do i=1, ndxmax
      write(*,'(i4, 9(es20.8))')  i, x(i), opac3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k), vely3d(i,j,k), &
         velz3d(i,j,k), t3d(i,j,k), scont3d(i,j,k), sline3d(i,j,k)
   enddo
   write(*,*)
!
end subroutine print_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_model3d_spc
!
   use prog_type
   use dime_model3d, only: nr_spc, ntheta_spc, nphi_spc, r_spc, theta_spc, phi_spc, &
      opac3d, opalbar3d, velx3d, vely3d, velz3d, sline3d, scont3d, t3d
   use params_spec, only: xic1, vmicro, xnue0, vth_fiducial, na, tmin, vmax, lstar, vrot, xlogg, rstar
!
   implicit none
!
! ... arguments
   integer(i4b) :: i, j, k
!
   write(*,*) '------3d atmospheric structure (spc) as fct of radius at theta, phi=0, 0-------'
   write(*,*)
!
   write(*,'(a20, es20.8)') 'I_c', xic1
   write(*,'(a20, es20.8)') 'xnue0', xnue0
   write(*,'(a20, es20.8)') 'tmin [K]', tmin
   write(*,'(a20, i20)') 'na', na
   write(*,'(a20, es20.8)') 'vmax [km/s]', vmax/1.d5
   write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro/1.d5
   write(*,'(a20, es20.8)') 'vth_fiducial [km/s]', vth_fiducial/1.d5
   write(*,'(a20, es20.8)') 'vrot [km/s]', vrot
   write(*,'(a20, es20.8)') 'rstar [r_sun]', rstar
   write(*,'(a20, es20.8)') 'lstar [l_sun]', lstar
   write(*,'(a20, es20.8)') 'logg ', xlogg
   write(*,*)
!
   j=1!ntheta_spc/2+1!1
   k=1
!
   write(*,*) 'at theta, phi', theta_spc(j), phi_spc(k)
   write(*,'(a4, 9(a20))') '#', 'r [r_star]', 'opac [1/sr]', 'opalbar [1/sr]', 'velx [vth*]', &
      'vely [vth*]', 'velz [vth*]', 't [K]', 's_cont', 's_line'
   do i=1, nr_spc
      write(*,'(i4, 9(es20.8))')  i, r_spc(i), opac3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k), vely3d(i,j,k), &
         velz3d(i,j,k), t3d(i,j,k), scont3d(i,j,k), sline3d(i,j,k)
   enddo
   write(*,*)
!
end subroutine print_model3d_spc
