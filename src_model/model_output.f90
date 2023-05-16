!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod1d
!
!-----------------save 1d model atmosphere as h5-file--------------------
!
   use prog_type
   use mod_directories, only: model_dir, model1d_file
   use dime_modext, only: nr_modext
   use model1d, only: r_modext1d, velr_modext1d, rho_modext1d, t_modext1d, vth_modext1d, eps_cont_modext1d
   use params_stellar, only: sr
   use hdf5
!
   implicit none
!
! ... local scalars
   integer(i4b) :: err
!
! ... output to hdf5
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
   integer(hsize_t), dimension(1) :: dims_r
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
   dims_r = (/ nr_modext /)
!
   write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model1d_file)
   write(*,*)

   call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model1d_file), h5f_acc_trunc_f, file_id, err)
   CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
   CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
   CALL h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)
   CALL h5sclose_f(aspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
   CALL h5screate_simple_f(1, dims_r, dspace_id, err)
   CALL h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, r_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'model', group_id, err)
   CALL h5screate_simple_f(1, dims_r, dspace_id, err)
   CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velr', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velr_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'temperature', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, t_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'vth', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, vth_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'eps_cont', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, eps_cont_modext1d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5fclose_f(file_id, err)
   CALL h5close_f(err)
!
!
end subroutine output_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod2d
!
!-----------------save 2d model atmosphere as h5-file--------------------
!
   use prog_type
   use mod_directories, only: model_dir, model2d_file
   use dime_modext, only: nr_modext, ntheta_modext
   use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
      velphi_modext2d, rho_modext2d, t_modext2d, vth_modext2d, eps_cont_modext2d
   use hdf5
!
   implicit none
!
! ... local scalars
   integer(i4b) :: err
!
! ... output to hdf5
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_theta
   integer(hsize_t), dimension(2) :: dims
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
   dims_r = (/ nr_modext /)
   dims_theta = (/ ntheta_modext /)
   dims = (/ nr_modext, ntheta_modext /)
!
   write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model2d_file)
   write(*,*)

   call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model2d_file), h5f_acc_trunc_f, file_id, err)
   CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
   CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
   CALL h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)
   CALL h5acreate_f(group_id, 'ntheta', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)

   CALL h5sclose_f(aspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
   CALL h5screate_simple_f(1, dims_r, dspace_id, err)
   CALL h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, r_modext2d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5screate_simple_f(1, dims_theta, dspace_id, err)
   CALL h5dcreate_f(group_id, 'theta', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, theta_modext2d, dims_theta, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'model', group_id, err)
   CALL h5screate_simple_f(2, dims, dspace_id, err)
   CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velr', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velr_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velth', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velth_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velphi', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velphi_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'temperature', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, t_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'vth', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, vth_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'eps_cont', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, eps_cont_modext2d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)

   CALL h5gclose_f(group_id, err)
   CALL h5fclose_f(file_id, err)
   CALL h5close_f(err)
!
!
end subroutine output_mod2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod3d
!
!-----------------save 3d model atmosphere as h5-file--------------------
!
   use prog_type
   use mod_directories, only: model_dir, model3d_file
   use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
   use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, &
      velr_modext3d, velth_modext3d, velphi_modext3d, &
      rho_modext3d, t_modext3d, trad_modext3d, vth_modext3d, eps_cont_modext3d
   use hdf5
!
   implicit none
!
! ... local scalars
   integer(i4b) :: err
!
! ... output to hdf5
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_r, dims_theta, dims_phi
   integer(hsize_t), dimension(3) :: dims
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
   dims_r = (/ nr_modext /)
   dims_theta = (/ ntheta_modext /)
   dims_phi = (/ nphi_modext /)
   dims = (/ nr_modext, ntheta_modext, nphi_modext /)
!
   write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model3d_file)
   write(*,*)

   call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model3d_file), h5f_acc_trunc_f, file_id, err)
   CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
   CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
   CALL h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)
   CALL h5acreate_f(group_id, 'ntheta', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)
   CALL h5acreate_f(group_id, 'nphi', h5t_native_integer, aspace_id, &
      attr_id, err)
   CALL h5awrite_f(attr_id, h5t_native_integer, nphi_modext, dims_scalars, err)
   CALL h5aclose_f(attr_id, err)
   CALL h5sclose_f(aspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
   CALL h5screate_simple_f(1, dims_r, dspace_id, err)
   CALL h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, r_modext3d, dims_r, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5screate_simple_f(1, dims_theta, dspace_id, err)
   CALL h5dcreate_f(group_id, 'theta', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, theta_modext3d, dims_theta, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5screate_simple_f(1, dims_phi, dspace_id, err)
   CALL h5dcreate_f(group_id, 'phi', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, phi_modext3d, dims_phi, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)
   CALL h5gclose_f(group_id, err)
!
   CALL h5gcreate_f(file_id, 'model', group_id, err)
   CALL h5screate_simple_f(3, dims, dspace_id, err)
   CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velr', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velr_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velth', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velth_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'velphi', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, velphi_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'temperature', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, t_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'trad', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, trad_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'vth', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, vth_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5dcreate_f(group_id, 'eps_cont', h5t_native_double, dspace_id, &
      dset_id, err)
   CALL h5dwrite_f(dset_id, h5t_native_double, eps_cont_modext3d, dims, err)
   CALL h5dclose_f(dset_id, err)
   CALL h5sclose_f(dspace_id, err)

   CALL h5gclose_f(group_id, err)
   CALL h5fclose_f(file_id, err)
   CALL h5close_f(err)
!
!
end subroutine output_mod3d
