!
!***********************************************************************
!--------------------------OUTPUT---------------------------------------
!***********************************************************************
!
subroutine output
!
   use prog_type
   use fund_const
   use angles, only: dim_mu, dim_phi, dim_omega, nodes_mu, n_x, n_y, n_z
   use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, &
      scont3d, mint3d, sline3d, mintbar3d, ssobo3d, &
      fcontx3d, fconty3d, fcontz3d, &
      kcontxx3d, kcontyy3d, kcontzz3d, kcontxy3d, kcontxz3d, kcontyz3d, &
      imask3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d, opac3d, opalbar3d, &
      velx3d, vely3d, velz3d, t3d, eps_cont3d
   use freq, only: nxobs, nodes_xobs, xnue0
   use iter, only: itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
   use mod_directories, only: output_dir, output_file
   use params_input, only: teff, trad, rstar, vmin, vmax, beta, vth_fiducial, vmicro, yhe, hei, xmloss, &
      kcont, kline, eps_line, na, alpha, kappa0, lstar, vrot, xlogg
   use bcondition, only: xic1, ntheta_gdark, theta_gdark, xic1_gdark, teff_gdark
   use options
   use hdf5
!
!
   implicit none
!
! ... arguments
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: iopt_sol2d, iopt_incl_cont, iopt_start_cont, iopt_ng_cont, iopt_ait_cont, &
      iopt_incl_line, iopt_start_line, iopt_ng_line, iopt_ait_line, &
      iopt_incl_gdark, iopt_incl_sdist
   real(dp) :: eps_cont
!
! ... local arrays
   integer(i4b), dimension(:,:,:), allocatable :: mask_totreg3d, mask_innreg3d, mask_bpoint3d, mask3d
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z, dims_mu, dims_phi, dims_xobs, &
      dims_itc, dims_itl, dims_omega, dims_gdark
   integer(hsize_t), dimension(3) :: dims_3d
!
!-----------------------------------------------------------------------
!
   write(*,*) '-------------------------output to directory-----------------------------------'
   write(*,*) 'output to file ', output_dir//'/'//trim(output_file)
   write(*,*)
!
!-----------------------------------------------------------------------
!
   dims_gdark = (/ ntheta_gdark /)
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_mu = (/ dim_mu /)
   dims_phi = (/ dim_phi /)
   dims_omega = (/ dim_omega /)
   dims_xobs = (/ nxobs /)
   dims_itc = (/ itmaxc /)
   dims_itl = (/ itmaxl /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
!in order to have correct (4-byte) integer in hdf5-file
   allocate(mask3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(mask_totreg3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(mask_innreg3d(ndxmax,ndymax,ndzmax), stat=err)
   allocate(mask_bpoint3d(ndxmax,ndymax,ndzmax), stat=err)
   mask_totreg3d = imask_totreg3d
   mask_innreg3d = imask_innreg3d
   mask_bpoint3d = imask_bpoint3d
   mask3d = imask3d
!
!--------convert all logicals to integers (to be read in idl)-----------
!
!convert options
   iopt_sol2d=0
   iopt_incl_cont=0
   iopt_incl_line=0
   iopt_start_cont=0
   iopt_ng_cont=0
   iopt_ait_cont=0
   iopt_start_line=0
   iopt_ng_line=0
   iopt_ait_line=0
   iopt_incl_gdark=0
   iopt_incl_sdist=0
   if(opt_sol2d) iopt_sol2d=1
   if(opt_incl_cont) iopt_incl_cont=1
   if(opt_start_cont) iopt_start_cont=1
   if(opt_ng_cont) iopt_ng_cont=1
   if(opt_ait_cont) iopt_ait_cont=1
   if(opt_incl_line) iopt_incl_line=1
   if(opt_start_line) iopt_start_line=1
   if(opt_ng_line) iopt_ng_line=1
   if(opt_ait_line) iopt_ait_line=1
   if(opt_incl_gdark) iopt_incl_gdark=1
   if(opt_incl_sdist) iopt_incl_sdist=1
!
!
!calculate mean eps_cont
   call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
!-----------------------------------------------------------------------
!
   call h5open_f(err)
   call h5fcreate_f(output_dir//'/'//trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!------------------------------options----------------------------------
!
   call h5gcreate_f(file_id, 'options', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'input_mod_dim', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, input_mod_dim, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'spatial_grid1d', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, spatial_grid1d, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'spatial_grid3d', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, spatial_grid3d, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_method', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_method, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_opal', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_opac', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_opac, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_angint_method', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_angint_method, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_ltec', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_ltec, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_sol2d', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_sol2d, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_incl_cont', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_incl_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_start_cont', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_start_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_ng_cont', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_ng_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_ait_cont', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_ait_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_incl_line', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_incl_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_start_line', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_start_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_ng_line', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_ng_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_ait_line', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_ait_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_alo_cont', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_alo_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_alo_line', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_alo_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_incl_gdark', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_incl_gdark, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'opt_incl_sdist', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, iopt_incl_sdist, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
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
   call h5acreate_f(group_id, 'dim_mu', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_mu, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'dim_phi', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_phi, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'dim_omega', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_omega, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nxobs', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nxobs, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'kcont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kcont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kline', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kline, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kappa0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_line', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_cont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'teff', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xlogg', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xlogg, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'trad', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'rstar', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'lstar', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, lstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmicro', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'yhe', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'hei', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, hei, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'mdot', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xmloss, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)

   call h5acreate_f(group_id, 'vmin', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmax', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'beta', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, beta, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vrot', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vrot, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
!------------------------boundary condition-----------------------------
!
   call h5gcreate_f(file_id, 'bcondition', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'ntheta_gdark', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ntheta_gdark, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_gdark, dspace_id, err)
   call h5dcreate_f(group_id, 'theta_gdark', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, theta_gdark, dims_gdark, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'teff_gdark', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, teff_gdark, dims_gdark, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xic1_gdark', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xic1_gdark, dims_gdark, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!---------------------angular grids-------------------------------------
!
   call h5gcreate_f(file_id, 'angles', group_id, err)
   call h5screate_simple_f(1, dims_mu, dspace_id, err)
   call h5dcreate_f(group_id, 'nodes_mu', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, nodes_mu, dims_mu, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_omega, dspace_id, err)
   call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'n_y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_y, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!-----------------frequency grid----------------------------------------
!
   call h5gcreate_f(file_id, 'frequencies', group_id, err)
   call h5screate_simple_f(1, dims_xobs, dspace_id, err)
   call h5dcreate_f(group_id, 'nodes_xobs', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, nodes_xobs, dims_xobs, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------convergence behaviour----------------------------
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxc', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'itmaxl', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxl, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'devmaxc', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, devmaxc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'devmaxl', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, devmaxl, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
!
   call h5screate_simple_f(1, dims_itc, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_arr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_arr, dims_itc, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
!
   call h5screate_simple_f(1, dims_itl, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_arr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_arr, dims_itl, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------3d solution--------------------------------
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'scont3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontx3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontx3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fconty3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fconty3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontxx3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontxx3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontyy3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontyy3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontzz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontzz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontxy3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontxy3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontxz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontxz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontyz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontyz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mintbar3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mintbar3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)


!
!--------------------------------3d model-------------------------------
!
   call h5gcreate_f(file_id, 'model3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mask_totreg3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask_totreg3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mask_innreg3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask_innreg3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mask_bpoint3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask_bpoint3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
!
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'eps_cont3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, eps_cont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'velx3d', h5t_native_real, dspace_id, &
      dset_id, err)
   velx3d=velx3d*vth_fiducial
   call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
   velx3d=velx3d/vth_fiducial
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vely3d', h5t_native_real, dspace_id, &
      dset_id, err)
   vely3d=vely3d*vth_fiducial
   call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
   vely3d=vely3d/vth_fiducial
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'velz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   velz3d=velz3d*vth_fiducial
   call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
   velz3d=velz3d/vth_fiducial
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!------------------------------debugging--------------------------------
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!!
!!------------------------output warnings--------------------------------
!!
!write(*,*) 'warnings.dat'
!open (1,file=trim(dirout)//'/warnings.dat', form='formatted')
!   write(1,*) 'warn_grid1', warn_grid1
!   write(1,*) 'warn_grid2', warn_grid2
!   write(1,*) 'warn_grid3', warn_grid3
!   write(1,*) 'warn_xobs1', warn_xobs1
!   write(1,*) 'warn_xobs2', warn_xobs2
!   write(1,*) 'warn_mu', warn_mu
!   write(1,*) 'warn_phi', warn_phi
!   write(1,*) 'warn_angles', warn_angles
!   write(1,*) 'warn_sym3d', warn_sym3d
!   write(1,*) 'warn_itmaxl', warn_itmaxl
!   write(1,*) 'warn_itmaxc', warn_itmaxc
!close(1)
!!
!!-----------------------output timing properties------------------------
!!
!write(*,*) 'timing.dat'
!open(1, file=trim(dirout)//'/timing.dat', form='formatted')
!   write(1,*) '--------------timing searchlight test----------------'
!   write(1,*) 'total time for all angles:, ', ttot_sl
!   write(1,*) 'number of calculated angles: ', dim_mu*dim_phi
!   write(1,*) 'time per angle (diffmeth): ', ttot_sl/(dim_mu*dim_phi)
!   write(1,*)
!!
!   write(1,*) '---continuum: timing for inverting alo---------------'
!   write(1,*) 'total time for all iterations:, ', ttot_alo
!   write(1,*) 'number of caluculating alo: ', it_tot-1
!   write(1,*) 'average time for inversion: ', ttot_alo/(it_tot-1)
!   write(1,*)
!!
!   write(1,*) '------------continuum: timing------------------------'
!   write(1,*) 'total time for all iterations', ttot_it
!   write(1,*) 'number of iterations: ', it_tot
!   write(1,*) 'average time per iteration: ', ttot_it/it_tot
!   write(1,*) 'total compuation time: ', te_tot-ts_tot
!   write(1,*)
!!
!   write(1,*) '--------line: timing for inverting alo---------------'
!   write(1,*) 'total time for all iterations (cpu): ', ttotl_alo
!   write(1,*) 'number of caluculating alo: ', it_totl-1
!   write(1,*) 'average time for inversion: ', ttotl_alo/(it_totl-1)
!   write(1,*)
!!
!   write(1,*) '-----------------line: timing------------------------'
!   write(1,*) 'total time for all iterations (cpu): ', ttotl_it
!   write(1,*) 'total time for all iterations (sys): ', ttotl_it_sys
!   write(1,*) 'number of iterations: ', it_totl
!   write(1,*) 'average time per iteration (cpu): ', ttotl_it/it_totl
!   write(1,*) 'average time per iteration (sys): ', ttotl_it_sys/it_totl
!   write(1,*) 'total compuation time (cpu): ', te_tot-ts_tot
!   write(1,*)
!close(1)
!!
!!
!!----------------------output 1-d grids in cgs--------------------------
!!
!write(*,*) 'grids1d.dat'
!open(1, file=trim(dirout)//'/grids1d.dat', form='formatted')
!   write(1,104) 'radius [cm]', 'vel [cm/s]', 'opath [1/cm]', 'opal [1/cm]', 't [k]'
!   do i=1, n1d_cr
!      write(1,105) r1d_cr(i), vel1d_cr(i), opath1d_cr(i), opalbar1d_cr(i), t1d_cr(i)
!   enddo
!close(1)
!104 format(5(a20))
!105 format(5(e20.8))
!!
!!-----------------output solution on central ray------------------------
!!
!write(*,*) 'solution1d.dat'
!open(1, file=trim(dirout)//'/solution1d.dat', form='formatted')
!   write(1,106) 'radius', 'j', 's_cont', 'jbar', 's_line', 'ssobo', 'ssoboc'
!   do i=1, n1d_cr
!      write(1,107) r1d_cr(i), mint1d_cr(i), scont1d_cr(i), mintbar1d_cr(i), sline1d_cr(i), ssobo1d_cr(i), ssoboc1d_cr(i)
!   enddo
!close(1)
!106 format(7(a20))
!107 format(7(e20.8))
!
!-----------------output 3-d grids in cgs (model)-----------------------
!
!write(*,*) 'model_r3d.dat'
!open(1, file=trim(dirout)//'/model_r3d.dat', form='unformatted')
!   write(1) r3d
!close(1)
!!
!!
!!
!!
!!
!!
!!
!write(*,*) 'solution_mintbar3d.dat'
!open(1, file=trim(dirout)//'/solution_mintbar3d.dat', form='unformatted')
!   write(1) mintbar3d
!close(1)
!!
!!
!
!
end subroutine output
!
!***********************************************************************
!***********************************************************************
!
!        SUBROUTINES FOR INPUT/OUTPUT OF ALO OPERATORS
!
!***********************************************************************
!***********************************************************************
!
subroutine output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, b_vec, nd, nnz)
!
!for linear system alocont_data*x_vec+b_vec=0.d0
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only: ndxmax, ndymax, ndzmax
   use hdf5
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: nd, nnz
   real(dp), dimension(nd), intent(in) :: alocont_data_diag, b_vec
   real(dp), dimension(nnz), intent(in) :: alocont_data
   integer(i4b), dimension(nnz), intent(in) :: alocont_colindx, alocont_rowindx
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_alo, dims_diag
!
! ... local characters
   character(len=14) :: fname='alo_cont.h5'
!
!--------------------------output grid----------------------------------
!
   dims_alo = (/ nnz /)
   dims_diag = (/ nd /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nd_alo', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nnz, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nd_diag', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd, dims_scalars, err)
   call h5aclose_f(attr_id, err)
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
!----------------------------coordinates--------------------------------
!
   call h5gcreate_f(file_id, 'alo', group_id, err)
   call h5screate_simple_f(1, dims_alo, dspace_id, err)
   call h5dcreate_f(group_id, 'alocont_data', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, alocont_data, dims_alo, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'alocont_rowindx', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, alocont_rowindx, dims_alo, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'alocont_colindx', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, alocont_colindx, dims_alo, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_diag, dspace_id, err)
   call h5dcreate_f(group_id, 'alocont_data_diag', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, alocont_data_diag, dims_diag, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'b_vec', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, b_vec, dims_diag, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
end subroutine output_alocont_coo
!
!***********************************************************************
!***********************************************************************
!
!   SUBROUTINES FOR INPUT/OUTPUT AFTER EACH ITERATION STEP
!
!***********************************************************************
!***********************************************************************
!
subroutine output_itstep_cont(itnr)
!
!------------------output after each iteration step---------------------
!input: itnr: current iteration step
!
   use prog_type
   use fund_const
   use dime3d, only: scont3d, mint3d
   use iter, only: epsmaxc_arr, itmaxc, devmaxc
   use mod_directories, only: output_dir_temp
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: itnr
!
! ... local scalars
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
   write(*,*) '-----------------------storing current iterate in------------------------------'
   write(*,*) output_dir_temp//'/iterparam_cont.dat'
   write(*,*) output_dir_temp//'/iterepsmax_cont.dat'
   write(*,*) output_dir_temp//'/solution_scont3d.dat'
   write(*,*) output_dir_temp//'/solution_mint3d.dat'
   write(*,*)
!
!----------------------output iteration parameter-----------------------
!
   open(1, file=trim(output_dir_temp)//'/iterparam_cont.dat', form='formatted')
   write(1,'(3(a20))') 'itmaxc', 'current it', 'devmax'
   write(1,'(2i20, e20.8)') itmaxc, itnr, devmaxc
   close(1)

   open(1, file=trim(output_dir_temp)//'/iterepsmax_cont.dat', form='unformatted')
   write(1) epsmaxc_arr
   close(1)
!
!-----------------output 3-d solution grids-----------------------------
!
   open(1, file=trim(output_dir_temp)//'/solution_scont3d.dat', form='unformatted')
   write(1) scont3d
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/solution_mint3d.dat', form='unformatted')
   write(1) mint3d
   close(1)
!
!
!
end subroutine output_itstep_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_itstep_line(itnr)
!
!------------------output after each iteration step---------------------
!input: itnr: current iteration step
!
   use prog_type
   use fund_const
   use dime3d, only: sline3d, mintbar3d
   use iter, only: epsmaxl_arr, itmaxl, devmaxl
   use mod_directories, only: output_dir_temp
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: itnr
!
! ... local scalars
!
! ... local arrays
!-----------------------------------------------------------------------
!
   write(*,*) '-----------------------storing current iterate in------------------------------'
   write(*,*) output_dir_temp//'/iterparam_line.dat'
   write(*,*) output_dir_temp//'/iterepsmax_line.dat'
   write(*,*) output_dir_temp//'/solution_sline3d.dat'
   write(*,*) output_dir_temp//'/solution_mintbar3d.dat'
   write(*,*)
!
!----------------------output iteration parameter-----------------------
!
   open(1, file=trim(output_dir_temp)//'/iterparam_line.dat', form='formatted')
   write(1,'(3(a20))') 'itmaxl', 'current it', 'devmax'
   write(1,'(2i20, e20.8)') itmaxl, itnr, devmaxl
   close(1)

   open(1, file=trim(output_dir_temp)//'/iterepsmax_line.dat', form='unformatted')
   write(1) epsmaxl_arr
   close(1)
!
!-----------------output 3-d solution grids-----------------------------
!
   open(1, file=trim(output_dir_temp)//'/solution_sline3d.dat', form='unformatted')
   write(1) sline3d
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/solution_mintbar3d.dat', form='unformatted')
   write(1) mintbar3d
   close(1)
!
!
!
end subroutine output_itstep_line
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine input_itstep_cont(it_start)
!
!
!
   use prog_type
   use fund_const
   use dime3d, only: scont3d, mint3d
   use iter, only: epsmaxc_arr
   use mod_directories, only: output_dir_temp
!
   implicit none
!
! ... arguments
   integer(i4b) :: it_start
!
! ... local scalars
   integer(i4b) :: dumi
   real(dp) :: dumr
!
! ... local character
   character(len=1000) :: header
!
! ... local arrays
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   write(*,*) '-----------------------reading current iterate from----------------------------'
   write(*,*) output_dir_temp//'/iterparam_cont.dat'
   write(*,*) output_dir_temp//'/iterepsmax_cont.dat'
   write(*,*) output_dir_temp//'/solution_scont3d.dat'
   write(*,*) output_dir_temp//'/solution_mint3d.dat'
   write(*,*)
!
!----------------------output iteration parameter-----------------------

   open(1, file=trim(output_dir_temp)//'/iterparam_cont.dat', form='formatted')
   read(1,*) header
   read(1,'(2i20, e20.8)') dumi, it_start, dumr
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/iterepsmax_cont.dat', form='unformatted')
   read(1) epsmaxc_arr
   close(1)
!
!-----------------output 3-d solution grids-----------------------------
!
   open(1, file=trim(output_dir_temp)//'/solution_scont3d.dat', form='unformatted')
   read(1) scont3d
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/solution_mint3d.dat', form='unformatted')
   read(1) mint3d
   close(1)
!
!
!
end subroutine input_itstep_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine input_itstep_line(it_start)
!
!
!
   use prog_type
   use fund_const
   use dime3d, only: sline3d, mintbar3d
   use iter, only: epsmaxl_arr
   use mod_directories, only: output_dir_temp
!
   implicit none
!
! ... arguments
   integer(i4b) :: it_start
!
! ... local scalars
   integer(i4b) :: dumi
   real(dp) :: dumr
!
! ... local character
   character(len=1000) :: header
!y
! ... local arrays
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   write(*,*) '-----------------------reading current iterate from----------------------------'
   write(*,*) output_dir_temp//'/iterparam_line.dat'
   write(*,*) output_dir_temp//'/iterepsmax_line.dat'
   write(*,*) output_dir_temp//'/solution_sline3d.dat'
   write(*,*) output_dir_temp//'/solution_mintbar3d.dat'
   write(*,*)
!
!----------------------output iteration parameter-----------------------

   open(1, file=trim(output_dir_temp)//'/iterparam_line.dat', form='formatted')
   read(1,*) header
   read(1,'(2i20, e20.8)') dumi, it_start, dumr
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/iterepsmax_line.dat', form='unformatted')
   read(1) epsmaxl_arr
   close(1)
!
!-----------------output 3-d solution grids-----------------------------
!
   open(1, file=trim(output_dir_temp)//'/solution_sline3d.dat', form='unformatted')
   read(1) sline3d
   close(1)
!
   open(1, file=trim(output_dir_temp)//'/solution_mintbar3d.dat', form='unformatted')
   read(1) mintbar3d
   close(1)
!
!
!
end subroutine input_itstep_line
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!
!SUBROUTINE OUTPUT_ITSTEP_LINE(ITNR, dirOUT)
!!
!!------------------OUTPUT AFTER EACH ITERATION STEP---------------------
!!INPUT: ITNR: CURRENT ITERATION STEP
!!
!USE prog_type
!USE fund_const
!USE DIME3D, ONLY: SLINE3D, MINTBAR3D
!USE ITER, ONLY: EPSMAXL_ARR, ITMAXL, DEVMAXL
!!
!IMPLICIT NONE
!!
!! ... arguments
!INTEGER(I4B), INTENT(IN) :: ITNR
!CHARACTER(len=100), INTENT(IN) :: dirOUT
!!
!! ... local scalars
!!
!! ... local arrays
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!WRITE(*,*) '--------------OUTPUT TO FILE------------------'
!!
!!----------------------OUTPUT ITERATION PARAMETER-----------------------
!
!OPEN(1, FILE=TRIM(dirOUT)//'/iterParam_LINE.dat', FORM='FORMATTED')
!   WRITE(1,'(3(A20))') 'ITMAXL', 'CURRENT IT', 'DEVMAX'
!   WRITE(1,'(2I20, E20.8)') ITMAXL, ITNR, DEVMAXL
!CLOSE(1)
!!
!OPEN(1, FILE=TRIM(dirOUT)//'/iterEPSMAX_LINE.dat', FORM='UNFORMATTED')
!   WRITE(1) EPSMAXL_ARR
!CLOSE(1)
!!
!!-----------------OUTPUT 3-D SOLUTION GRIDS-----------------------------
!!
!OPEN(1, FILE=TRIM(dirOUT)//'/SOLUTION_SLINE3D.dat', FORM='UNFORMATTED')
!   WRITE(1) SLINE3D
!CLOSE(1)
!!
!OPEN(1, FILE=TRIM(dirOUT)//'/SOLUTION_MINTBAR3D.dat', FORM='UNFORMATTED')
!   WRITE(1) MINTBAR3D
!CLOSE(1)
!!
!!
!!
!END SUBROUTINE OUTPUT_ITSTEP_LINE
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!
!SUBROUTINE INPUT_ITSTEP_LINE(IT_START, dirIN)
!!
!!
!!
!USE prog_type
!USE fund_const
!USE DIME3D, ONLY: SLINE3D, MINTBAR3D
!USE ITER, ONLY: EPSMAXL_ARR
!!
!IMPLICIT NONE
!!
!! ... arguments
!INTEGER(I4B) :: IT_START
!CHARACTER(len=100), INTENT(IN) :: dirIN
!!
!! ... local scalars
!INTEGER(I4B) :: DUMI
!REAL(DP) :: DUMR
!!
!! ... local character
!CHARACTER(len=1000) :: HEADER
!!
!! ... local arrays
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!
!WRITE(*,*) '--------------INPUT FROM FILE-----------------'
!!
!!----------------------INPUT ITERATION PARAMETER------------------------
!!
!OPEN(1, FILE=TRIM(dirIN)//'/iterParam_LINE.dat', FORM='FORMATTED')
!   READ(1,*) HEADER
!   READ(1,'(2I20, E20.8)') DUMI, IT_START, DUMR
!CLOSE(1)
!!
!OPEN(1, FILE=TRIM(dirIN)//'/iterEPSMAX_LINE.dat', FORM='UNFORMATTED')
!   READ(1) EPSMAXL_ARR
!CLOSE(1)
!!
!!-----------------INPUT 3-D SOLUTION GRIDS------------------------------
!!
!OPEN(1, FILE=TRIM(dirIN)//'/SOLUTION_SLINE3D.dat', FORM='UNFORMATTED')
!   READ(1) SLINE3D
!CLOSE(1)
!!
!OPEN(1, FILE=TRIM(dirIN)//'/SOLUTION_MINTBAR3D.dat', FORM='UNFORMATTED')
!   READ(1) MINTBAR3D
!CLOSE(1)
!!
!!
!!
!END SUBROUTINE INPUT_ITSTEP_LINE
!
!***********************************************************************
!***********************************************************************
!                 OUTPUT FOR BENCHMARKS
!***********************************************************************
!***********************************************************************
!
subroutine output_benchmark01
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only: x, z, ndxmax, ndzmax
   use mod_benchmark, only: nr, r1d, x_u2d, z_u2d, x_d2d, z_d2d, int2d_sc, int2d_fvm, &
      xu_interp_opac2d, xd_interp_opac2d, &
      zu_interp_opac2d, zd_interp_opac2d, &
      xu_interp_int2d, zu_interp_int2d, &
      opac2d, scont2d, &
      opac_u2d, scont_u2d, opac_d2d, scont_d2d, int_u2d, &
      int1d, abs1d, contr1d, int1d_sc, abs1d_sc, contr1d_sc, opac1d, scont1d, &
      int1d_scray, int1d_fvmray, abs2d_sc, abs2d_fvm, contr2d_sc, contr2d_fvm, int1d_fvm, abs1d_fvm, contr1d_fvm
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x , dims_z, dims_r
   integer(hsize_t), dimension(2) :: dims_2d
   integer(hsize_t), dimension(3) :: dims_interp2d
!
! ... local characters
   character(len=14) :: fname='benchmark01.h5'
!
!--------------------------output grid----------------------------------
!
   dims_r = (/ nr /)
   dims_x = (/ ndxmax /)
   dims_z = (/ ndzmax /)
   dims_2d = (/ ndxmax, ndzmax /)
   dims_interp2d = (/ ndxmax, ndzmax, 8 /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------solution 1d--------------------------------
!
   call h5gcreate_f(file_id, 'solution1d', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'opac1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1d, dims_r, err)
   call h5dclose_f(dset_id, err)

   call h5dcreate_f(group_id, 'int1d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'abs1d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, abs1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'contr1d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, contr1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)

   call h5dcreate_f(group_id, 'int1d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'abs1d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, abs1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'contr1d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, contr1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)

   call h5dcreate_f(group_id, 'int1d_scray', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_scray, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_fvmray', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_fvmray, dims_r, err)
   call h5dclose_f(dset_id, err)
!
   call h5dcreate_f(group_id, 'int1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'abs1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, abs1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'contr1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, contr1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------solution 2d--------------------------------
!
   call h5gcreate_f(file_id, 'solution2d', group_id, err)
   call h5screate_simple_f(2, dims_2d, dspace_id, err)
   call h5dcreate_f(group_id, 'int2d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2d_sc, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'abs2d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, abs2d_sc, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'contr2d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, contr2d_sc, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int2d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2d_fvm, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'abs2d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, abs2d_fvm, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'contr2d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, contr2d_fvm, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!---------------------model: oapcity and source functions---------------
!
   call h5gcreate_f(file_id, 'model2d', group_id, err)
   call h5screate_simple_f(2, dims_2d, dspace_id, err)
   call h5dcreate_f(group_id, 'opac2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!------------------------------debugging--------------------------------
!
   call h5gcreate_f(file_id, 'debug', group_id, err)
   call h5screate_simple_f(2, dims_2d, dspace_id, err)
   call h5dcreate_f(group_id, 'x_upwind', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x_u2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'x_downwind', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x_d2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'z_upwind', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z_u2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'z_downwind', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z_d2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
!
   call h5screate_simple_f(3, dims_interp2d, dspace_id, err)
   call h5dcreate_f(group_id, 'xu_int2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xu_interp_int2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'zu_int2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, zu_interp_int2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xu_opac2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xu_interp_opac2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'zu_opac2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, zu_interp_opac2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xd_opac2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xd_interp_opac2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'zd_opac2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, zd_interp_opac2d, dims_interp2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
!
   call h5screate_simple_f(2, dims_2d, dspace_id, err)
   call h5dcreate_f(group_id, 'int_u2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int_u2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac_u2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac_u2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont_u2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont_u2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac_d2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac_d2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont_d2d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont_d2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)

   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
end subroutine output_benchmark01
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark02
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only: x, z, ndxmax, ndzmax
   use angles, only: dim_mu, nodes_mu
   use mod_benchmark, only: int2d_sc, int2d_fvm, n_z
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x , dims_z
   integer(hsize_t), dimension(2) :: dims_2d
!
! ... local characters
   character(len=14) :: fname='benchmark02.h5'
!
!
!
   dims_x = (/ ndxmax /)
   dims_z = (/ ndzmax /)
   dims_2d = (/ ndxmax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'direction', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'mu', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, n_z, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'searchlight2d', group_id, err)
   call h5screate_simple_f(2, dims_2d, dspace_id, err)
   call h5dcreate_f(group_id, 'int2d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2d_sc, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int2d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2d_fvm, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
!
end subroutine output_benchmark02
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark03
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only: x, z, ndxmax, ndzmax
   use angles, only: dim_mu, nodes_mu
   use mod_benchmark, only: int2dsc_angdep, int2dfvm_angdep
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x , dims_z, dims_mu
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark03.h5'
!
!
!
   dims_x = (/ ndxmax /)
   dims_z = (/ ndzmax /)
   dims_mu = (/ dim_mu /)
!
   dims_3d = (/ ndxmax, ndzmax, dim_mu /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndxmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndzmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'dim_mu', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_mu, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_mu, dspace_id, err)
   call h5dcreate_f(group_id, 'nodes_mu', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, nodes_mu, dims_mu, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'searchlight2d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'int2dsc_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2dsc_angdep, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int2dfvm_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int2dfvm_angdep, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
   deallocate(int2dsc_angdep, stat=err)
   if(err.ne.0) stop 'error searchlight_test: deallocation'
   deallocate(int2dfvm_angdep, stat=err)
   if(err.ne.0) stop 'error searchlight_test: deallocation'
!
!
!
end subroutine output_benchmark03
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark04
!
   use prog_type
   use mod_directories, only: output_dir_test
   use angles, only: dim_mu, nodes_mu
   use mod_benchmark, only: mint1d_theo, mint1d_sc, mint1d_fvm, intsc_angdep, intfvm_angdep, n1d_angdep, r1d_angdep
   use dimecr, only: n1d_cr, r1d_cr
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_mu
   integer(hsize_t), dimension(2) :: dims_angdep
!
! ... local characters
   character(len=14) :: fname='benchmark04.h5'
!
!
!
   dims_r = (/ n1d_cr /)
   dims_mu = (/ dim_mu /)
   dims_angdep = (/ n1d_angdep, dim_mu /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'dim_mu', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_mu, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'n1d_angdep', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_angdep, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'angles', group_id, err)
   call h5screate_simple_f(1, dims_mu, dspace_id, err)
   call h5dcreate_f(group_id, 'nodes_mu', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, nodes_mu, dims_mu, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_angdep, dspace_id, err)
   call h5dcreate_f(group_id, 'r1d_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_angdep, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'mint_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_theo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_angdep', group_id, err)
   call h5screate_simple_f(2, dims_angdep, dspace_id, err)
   call h5dcreate_f(group_id, 'intsc_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, intsc_angdep, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'intfvm_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, intfvm_angdep, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark04
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark05
!
   use prog_type
   use mod_directories, only: output_dir_test
   use mod_benchmark, only: mint1d_sc, mint1d_fvm, mint1d_joray, mint1d_jomom, epsmaxc_sc, epsmaxc_fvm, &
      t1d_cr, t1d_jo, opac1d_cr, opac1d_jo
   use dimecr, only: n1d_cr, r1d_cr
   use iter, only: itmaxc
   use params_input, only: kcont
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_iter
!
! ... local characters
   character(len=14) :: fname='benchmark05.h5'
!
!
!
   dims_r = (/ n1d_cr /)
   dims_iter = (/ itmaxc /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'kcont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kcont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxc', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_iter, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_sc, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_fvm, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 't1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'mint_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_joray', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_joray, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_jomom', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_jomom, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark05
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark06
!
   use prog_type
   use mod_directories, only: output_dir_test
   use mod_benchmark, only: nd_method01, nd_method02, nd_method03, nd_method04, nd_fine, &
      s1d_method01, s1d_method02, s1d_method03, s1d_method04, s1d_fine, &
      tau1d_method01, tau1d_method02, tau1d_method03, tau1d_method04, tau1d_fine, &
      vth1d_method01, vth1d_method02, vth1d_method03, vth1d_method04, vth1d_fine, &
      vel1d_method01, vel1d_method02, vel1d_method03, vel1d_method04, vel1d_fine, &
      xcmf1d_method01, xcmf1d_method02, xcmf1d_method03, xcmf1d_method04, xcmf1d_fine, &
      profile1d_method01, profile1d_method02, profile1d_method03, profile1d_method04, profile1d_fine, &
      opalbar1d_method01, opalbar1d_method02, opalbar1d_method03, opalbar1d_method04, opalbar1d_fine, &
      opal1d_method01, opal1d_method02, opal1d_method03, opal1d_method04, opal1d_fine, &
      sline1d_method01, sline1d_method02, sline1d_method03, sline1d_method04, sline1d_fine, &
      int1d_method01, int1d_method02, int1d_method03, int1d_method04, int1d_fine, &
      ndc_method01, ndc_method02, ndc_method03, ndc_method04, ndc_fine, &
      s1dc_method01, s1dc_method02, s1dc_method03, s1dc_method04, s1dc_fine, &
      tau1dc_method01, tau1dc_method02, tau1dc_method03, tau1dc_method04, tau1dc_fine, &
      vth1dc_method01, vth1dc_method02, vth1dc_method03, vth1dc_method04, vth1dc_fine, &
      vel1dc_method01, vel1dc_method02, vel1dc_method03, vel1dc_method04, vel1dc_fine, &
      xcmf1dc_method01, xcmf1dc_method02, xcmf1dc_method03, xcmf1dc_method04, xcmf1dc_fine, &
      profile1dc_method01, profile1dc_method02, profile1dc_method03, profile1dc_method04, profile1dc_fine, &
      opalbar1dc_method01, opalbar1dc_method02, opalbar1dc_method03, opalbar1dc_method04, opalbar1dc_fine, &
      opal1dc_method01, opal1dc_method02, opal1dc_method03, opal1dc_method04, opal1dc_fine, &
      sline1dc_method01, sline1dc_method02, sline1dc_method03, sline1dc_method04, sline1dc_fine, &
      scont1dc_method01, scont1dc_method02, scont1dc_method03, scont1dc_method04, scont1dc_fine, &
      stot1dc_method01, stot1dc_method02, stot1dc_method03, stot1dc_method04, stot1dc_fine, &
      opac1dc_method01, opac1dc_method02, opac1dc_method03, opac1dc_method04, opac1dc_fine, &
      opatot1dc_method01, opatot1dc_method02, opatot1dc_method03, opatot1dc_method04, opatot1dc_fine, &
      int1dc_method01, int1dc_method02, int1dc_method03, int1dc_method04, int1dc_fine
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_fine, dims_method01, dims_method02, dims_method03, dims_method04, &
      dimsc_fine, dimsc_method01, dimsc_method02, dimsc_method03, dimsc_method04
!
! ... local characters
   character(len=14) :: fname='benchmark06.h5'
!
!
!
   dims_fine = (/ nd_fine /)
   dims_method01 = (/ nd_method01 /)
   dims_method02 = (/ nd_method02 /)
   dims_method03 = (/ nd_method03 /)
   dims_method04 = (/ nd_method04 /)

   dimsc_fine = (/ ndc_fine /)
   dimsc_method01 = (/ ndc_method01 /)
   dimsc_method02 = (/ ndc_method02 /)
   dimsc_method03 = (/ ndc_method03 /)
   dimsc_method04 = (/ ndc_method04 /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!------------------------------------------------------------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nd_method01', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd_method01, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nd_method02', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd_method02, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nd_method03', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd_method03, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nd_method04', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd_method04, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'nd_fine', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, nd_fine, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndc_method01', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndc_method01, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndc_method02', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndc_method02, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndc_method03', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndc_method03, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndc_method04', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndc_method04, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'ndc_fine', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, ndc_fine, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method00', group_id, err)
   call h5screate_simple_f(1, dims_fine, dspace_id, err)
   call h5dcreate_f(group_id, 's1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_fine, dims_fine, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method01', group_id, err)
   call h5screate_simple_f(1, dims_method01, dspace_id, err)
   call h5dcreate_f(group_id, 's1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_method01, dims_method01, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method02', group_id, err)
   call h5screate_simple_f(1, dims_method02, dspace_id, err)
   call h5dcreate_f(group_id, 's1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_method02, dims_method02, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method03', group_id, err)
   call h5screate_simple_f(1, dims_method03, dspace_id, err)
   call h5dcreate_f(group_id, 's1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_method03, dims_method03, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method04', group_id, err)
   call h5screate_simple_f(1, dims_method04, dspace_id, err)
   call h5dcreate_f(group_id, 's1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1d_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1d_method04, dims_method04, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method00c', group_id, err)
   call h5screate_simple_f(1, dimsc_fine, dspace_id, err)
   call h5dcreate_f(group_id, 's1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opatot1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opatot1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'stot1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, stot1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1dc_fine', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1dc_fine, dimsc_fine, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method01c', group_id, err)
   call h5screate_simple_f(1, dimsc_method01, dspace_id, err)
   call h5dcreate_f(group_id, 's1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opatot1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opatot1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'stot1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, stot1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1dc_method01', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1dc_method01, dimsc_method01, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method02c', group_id, err)
   call h5screate_simple_f(1, dimsc_method02, dspace_id, err)
   call h5dcreate_f(group_id, 's1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opatot1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opatot1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'stot1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, stot1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1dc_method02', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1dc_method02, dimsc_method02, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method03c', group_id, err)
   call h5screate_simple_f(1, dimsc_method03, dspace_id, err)
   call h5dcreate_f(group_id, 's1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opatot1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opatot1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'stot1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, stot1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1dc_method03', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1dc_method03, dimsc_method03, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'method04c', group_id, err)
   call h5screate_simple_f(1, dimsc_method04, dspace_id, err)
   call h5dcreate_f(group_id, 's1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'tau1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, tau1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vth1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vth1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'vel1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vel1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'xcmf1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcmf1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opal1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opal1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opatot1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opatot1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'profile1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, profile1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'stot1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, stot1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int1dc_method04', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int1dc_method04, dimsc_method04, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)


   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark06
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark07
!
   use prog_type
   use mod_directories, only: output_dir_test
   use mod_benchmark, only: ssobo1d_crx, ssobo1d_cry, ssobo1d_crz, velr1d_cr, opalbar1d_cr, t1d_cr, &
      ssobo1d
   use dimecr, only: n1d_cr, r1d_cr
   use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, ssobo3d
   use options, only: opt_opal
   use params_input, only: kline, eps_line
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_x, dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d

!
! ... local characters
   character(len=14) :: fname='benchmark07.h5'
!
!
!
   dims_r = (/ n1d_cr /)
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'kline', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kline, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_line', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
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
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'velr1d_cr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velr1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_cr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't1d_cr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'solution1d', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_crx', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_crx, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_cry', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_cry, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_crz', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_crz, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'ssobo3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark07
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark08
!
   use prog_type
   use mod_directories, only: output_dir_test
   use mod_benchmark, only: sline1d_sc, sline1d_fvm, sline1d_jo, ssobo1d_jo, ssobo1d_cr, &
      epsmaxl_sc, epsmaxl_fvm, &
      t1d_cr, t1d_jo, opalbar1d_cr, opalbar1d_jo, velr1d_cr
   use dimecr, only: n1d_cr, r1d_cr
   use iter, only: itmaxl
   use options, only: opt_opal
   use params_input, only: kline, alpha, kappa0, eps_line, vmax, rstar, trad, na, vmicro, teff, vth_fiducial
   use params_stellar, only: sr
   use freq, only: xnue0
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_iter
!
! ... local characters
   character(len=14) :: fname='benchmark08.h5'
!
!
!
   dims_r = (/ n1d_cr /)
   dims_iter = (/ itmaxl /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'opt_opal', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kline', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kline, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kappa0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_line', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmax', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmicro', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'trad', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'teff', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'rstar', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxl', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxl, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_iter, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_sc, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_fvm, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 't1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   opalbar1d_cr=opalbar1d_cr/sr
   call h5dcreate_f(group_id, 'opalbar1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   opalbar1d_cr=opalbar1d_cr*sr
   call h5dcreate_f(group_id, 'velr1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velr1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'sline_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_sc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_fvm, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_cr', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark08
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark10
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only: x, y, z, ndxmax, ndymax, ndzmax
   use angles, only: dim_mu, nodes_mu
   use mod_benchmark, only: int3d_sc, int3d_fvm, n_y, n_z
   use hdf5
!
   implicit none
!
! ... local scalars
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark10.h5'
!
!
!
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
   call h5gcreate_f(file_id, 'direction', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'n_x', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, sqrt(1.d0-n_y**2-n_z**2), dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'n_y', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, n_y, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'n_z', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, n_z, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

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
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'searchlight3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'int3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'int3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, int3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
!
end subroutine output_benchmark10
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark11
!
   use prog_type
   use mod_directories, only: output_dir_test
   use angles, only: dim_omega, n_x, n_y, n_z, weight_omega
   use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask3d
   use mod_benchmark, only: n1d_angdep, mint3d_theo, mint3d_sc, mint3d_fvm, &
      intsc_angdep2, intfvm_angdep2, xcoord_angdep2, ycoord_angdep2, zcoord_angdep2, &
      fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm, &
      fcontr3d_theo, fcontth3d_theo, fcontphi3d_theo
   use options, only: opt_angint_method
   use hdf5
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id, type_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_omega
   integer(hsize_t), dimension(2) :: dims_angdep
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark11.h5'
!
   integer, dimension(:,:,:), allocatable :: mask3d
   integer, dimension(:,:), allocatable :: phi_mask_int
!
!in order to have correct (4-byte) integer in hdf5-file
   allocate(mask3d(ndxmax,ndymax,ndzmax))
   mask3d=imask3d
!
!
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
!
   dims_omega = (/ dim_omega /)
   dims_angdep = (/ n1d_angdep, dim_omega /)
!
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
   call h5gcreate_f(file_id, 'options', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'opt_angint_method', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_angint_method, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
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
   call h5acreate_f(group_id, 'dim_omega', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, dim_omega, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'n1d_angdep', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_angdep, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'angles', group_id, err)
   call h5screate_simple_f(1, dims_omega, dspace_id, err)
   call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'n_y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_y, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'weight_omega', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, weight_omega, dims_omega, err)
   call h5dclose_f(dset_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_angdep, dspace_id, err)
   call h5dcreate_f(group_id, 'xcoord_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, xcoord_angdep2, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ycoord_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ycoord_angdep2, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'zcoord_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, zcoord_angdep2, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mint3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint3d_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_theo, dims_3d, err)
   call h5dclose_f(dset_id, err)


   call h5dcreate_f(group_id, 'fcontr3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontr3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontr3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontr3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontr3d_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontr3d_theo, dims_3d, err)
   call h5dclose_f(dset_id, err)


   call h5dcreate_f(group_id, 'fcontth3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontth3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontth3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontth3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontth3d_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontth3d_theo, dims_3d, err)
   call h5dclose_f(dset_id, err)

   call h5dcreate_f(group_id, 'fcontphi3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontphi3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontphi3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontphi3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontphi3d_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontphi3d_theo, dims_3d, err)
   call h5dclose_f(dset_id, err)



   call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_angdep', group_id, err)
   call h5screate_simple_f(2, dims_angdep, dspace_id, err)
   call h5dcreate_f(group_id, 'intsc_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, intsc_angdep2, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'intfvm_angdep', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, intfvm_angdep2, dims_angdep, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark11
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark12
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only:  t3d, eps_cont3d, opac3d, imask3d, x, y, z, ndxmax, ndymax, ndzmax, imask_totreg3d
   use mod_benchmark, only: mint3d_sc, mint3d_fvm, mint1d_joray, mint1d_jomom, fcont1d_joray, fcont1d_jomom, &
      epsmaxc_sc, epsmaxc_fvm, &
      fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm, &
      kcontrr3d_sc, kcontrth3d_sc, kcontrphi3d_sc, kcontthth3d_sc, kcontthphi3d_sc, kcontphiphi3d_sc, &
      kcontrr3d_fvm, kcontrth3d_fvm, kcontrphi3d_fvm, kcontthth3d_fvm, kcontthphi3d_fvm, kcontphiphi3d_fvm, &
      t1d_jo, opac1d_jo
   use dimecr, only: n1d_cr, r1d_cr
   use iter, only: itmaxc
   use params_input, only: kcont
   use params_stellar, only: sr
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars
   real(dp) :: eps_cont
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_iter
   integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark12.h5'
   integer, dimension(:,:,:), allocatable :: mask3d
!
!in order to have correct (4-byte) integer in hdf5-file
   allocate(mask3d(ndxmax,ndymax,ndzmax))
   mask3d=imask3d
!
!-----------------------------------------------------------------------
!
!calculate a mean of eps_cont
   call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
   dims_r = (/ n1d_cr /)
   dims_iter = (/ itmaxc /)
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'kcont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kcont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_cont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

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
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxc', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_iter, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_sc, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_fvm, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 't1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'eps_cont3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, eps_cont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d/sr
   call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d*sr
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'mint_joray', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_joray, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint_jomom', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint1d_jomom, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcont_joray', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcont1d_joray, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcont_jomom', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcont1d_jomom, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mint3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontr3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontr3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontr3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontr3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontth3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontth3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontth3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontth3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontphi3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontphi3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'fcontphi3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, fcontphi3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrr3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrr3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontthth3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontthth3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontphiphi3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontphiphi3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrth3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrth3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrphi3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrphi3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontthphi3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontthphi3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrr3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrr3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontthth3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontthth3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontphiphi3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontphiphi3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrth3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrth3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontrphi3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontrphi3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'kcontthphi3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, kcontthphi3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark12
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark13
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only:  t3d, opac3d, opalbar3d, ssobo3d, velx3d, vely3d, velz3d, imask3d, x, y, z, ndxmax, ndymax, ndzmax
   use mod_benchmark, only: mintbar3d_sc, mintbar3d_fvm, sline3d_sc, sline3d_fvm, &
      epsmaxl_sc, epsmaxl_fvm, velr1d_cr, &
      t1d_jo, opac1d_jo, opalbar1d_jo, sline1d_jo, ssobo1d_jo
   use dimecr, only: n1d_cr, r1d_cr
   use iter, only: itmaxl
   use options, only: opt_opal
   use params_input, only: kline, eps_line, na, alpha, kappa0, vmax, vmicro, trad, teff, &
      rstar, vth_fiducial
   use freq, only: xnue0
   use params_stellar, only: sr
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars

! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_iter
   integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark13.h5'
   integer, dimension(:,:,:), allocatable :: mask3d
!
!in order to have correct (4-byte) integer in hdf5-file
   allocate(mask3d(ndxmax,ndymax,ndzmax))
   mask3d=imask3d
!
!-----------------------------------------------------------------------

   dims_r = (/ n1d_cr /)
   dims_iter = (/ itmaxl /)
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
!
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'opt_opal', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kline', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kline, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kappa0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_line', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmax', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmicro', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'trad', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'teff', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'rstar', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

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
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxl', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxl, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_iter, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_sc, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_fvm, dims_iter, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 't1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'velr1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velr1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d/sr
   call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d*sr
   opalbar3d=opalbar3d/sr
   call h5dcreate_f(group_id, 'opalbar3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opalbar3d=opalbar3d*sr
   velx3d=velx3d*vth_fiducial
   call h5dcreate_f(group_id, 'velx3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   velx3d=velx3d/vth_fiducial
   vely3d=vely3d*vth_fiducial
   call h5dcreate_f(group_id, 'vely3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   vely3d=vely3d/vth_fiducial
   velz3d=velz3d*vth_fiducial
   call h5dcreate_f(group_id, 'velz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   velz3d=velz3d/vth_fiducial
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'sline1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mintbar3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mintbar3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mintbar3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mintbar3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark13
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark14
!
   use prog_type
   use mod_directories, only: output_dir_test
   use dime3d, only:  t3d, eps_cont3d, opac3d, opalbar3d, ssobo3d, velx3d, vely3d, velz3d, imask3d, imask_totreg3d, x, y, z, ndxmax, ndymax, ndzmax
   use mod_benchmark, only: mint3d_sc, mint3d_fvm, mintbar3d_sc, mintbar3d_fvm, sline3d_sc, sline3d_fvm, scont3d_sc, scont3d_fvm, &
      epsmaxl_sc, epsmaxl_fvm, epsmaxc_sc, epsmaxc_fvm, velr1d_cr, &
      t1d_jo, opac1d_jo, opalbar1d_jo, sline1d_jo, scont1d_jo, ssobo1d_jo
   use dimecr, only: n1d_cr, r1d_cr
   use iter, only: itmaxl, itmaxc
   use options, only: opt_opal
   use params_input, only: kline, eps_line, kcont, na, alpha, kappa0, vmax, vmicro, trad, teff, &
      rstar, vth_fiducial
   use freq, only: xnue0
   use params_stellar, only: sr
   use bcondition, only: xic1
   use hdf5
!
   implicit none
!
! ... local scalars
   real(dp) :: eps_cont
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_iterl, dims_iterc
   integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z
   integer(hsize_t), dimension(3) :: dims_3d
!
! ... local characters
   character(len=14) :: fname='benchmark14.h5'
   integer, dimension(:,:,:), allocatable :: mask3d
!
!in order to have correct (4-byte) integer in hdf5-file
   allocate(mask3d(ndxmax,ndymax,ndzmax))
   mask3d=imask3d
!
!calculate a mean of eps_cont
   call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
!-----------------------------------------------------------------------

   dims_r = (/ n1d_cr /)
   dims_iterl = (/ itmaxl /)
   dims_iterc = (/ itmaxc /)
   dims_x = (/ ndxmax /)
   dims_y = (/ ndymax /)
   dims_z = (/ ndzmax /)
   dims_3d = (/ ndxmax, ndymax, ndzmax /)
!
   write(*, *) 'output to file ', output_dir_test//'/'//fname

   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)

   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'opt_opal', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, opt_opal, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kcont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kcont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kline', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kline, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'kappa0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, kappa0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_line', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_line, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'eps_cont', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, eps_cont, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmax', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vmicro', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vmicro, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'trad', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'teff', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'rstar', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)

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
   call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n1d_cr, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_x, dspace_id, err)
   call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_y, dspace_id, err)
   call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_z, dspace_id, err)
   call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'r', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, r1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'itmaxc', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'itmaxl', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, itmaxl, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5screate_simple_f(1, dims_iterc, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_sc, dims_iterc, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxc_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_fvm, dims_iterc, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_iterl, dspace_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_sc, dims_iterl, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'epsmaxl_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, epsmaxl_fvm, dims_iterl, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'model_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 't1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opac1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'opalbar1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'velr1d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velr1d_cr, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)

   call h5gcreate_f(file_id, 'model3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_integer, mask3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'eps_cont3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, eps_cont3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 't3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, t3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d/sr
   call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opac3d=opac3d*sr
   opalbar3d=opalbar3d/sr
   call h5dcreate_f(group_id, 'opalbar3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, opalbar3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   opalbar3d=opalbar3d*sr
   velx3d=velx3d*vth_fiducial
   call h5dcreate_f(group_id, 'velx3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   velx3d=velx3d/vth_fiducial
   vely3d=vely3d*vth_fiducial
   call h5dcreate_f(group_id, 'vely3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   vely3d=vely3d/vth_fiducial
   velz3d=velz3d*vth_fiducial
   call h5dcreate_f(group_id, 'velz3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   velz3d=velz3d/vth_fiducial
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution_cr', group_id, err)
   call h5screate_simple_f(1, dims_r, dspace_id, err)
   call h5dcreate_f(group_id, 'scont1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo1d_jo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo1d_jo, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'solution3d', group_id, err)
   call h5screate_simple_f(3, dims_3d, dspace_id, err)
   call h5dcreate_f(group_id, 'mint3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mint3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mint3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mintbar3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mintbar3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'mintbar3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, mintbar3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'ssobo3d', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, ssobo3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'scont3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, scont3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline3d_sc', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline3d_sc, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'sline3d_fvm', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, sline3d_fvm, dims_3d, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine output_benchmark14

