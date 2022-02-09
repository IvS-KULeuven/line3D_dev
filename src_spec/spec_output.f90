!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_surface(indx, xobs, alpha, gamma)
!
!-------------output of surface emergent intensity for given------------
!--------------xobs and direction specified by alpha, gamma-------------
!
use prog_type
use options_spec, only: output_dir
use dime_spec, only: np, nzeta
use mod_spectrum, only: p, zeta
use mod_surfb, only: iem_surface, iemi_surface, iabs_surface, icont_surface
use params_spec, only: xic1, xnue0, vth_fiducial
use hdf5
!
implicit none
!
! ... arguments
integer(i4b) :: indx
real(dp), intent(in) :: xobs, alpha, gamma
!
! ... local scalars
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, attr_id
!
! ... local arrays
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_p , dims_zeta
integer(hsize_t), dimension(2) :: dims_2d
!
! ... local characters
character(len=8) :: name
character(len=14) :: fname_base='/spec_surface_'
character(len=512) :: fname
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
fname = trim(output_dir)//fname_base//trim(name)//'.h5'
!
!
!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(fname)
!
!-----------------------------------------------------------------------
!
dims_p = (/ np /)
dims_zeta = (/ nzeta /)
dims_2d = (/ np, nzeta /)
!
!-----------------------------------------------------------------------
!
!initialize fortran interface
call h5open_f (err)
!
!create hdf5 file
call h5fcreate_f(trim(fname), h5f_acc_trunc_f, file_id, err)
!
!---------------------------first group: dimensions---------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'np', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, np, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'nzeta', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nzeta, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!---------------------------second group: parameter---------------------
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
!
!create the dataspace to store parameter (each parameter stored as a single scalar)
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!
         call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)            
!            
         call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'xobs', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xobs, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'gamma', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, gamma, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
      call h5sclose_f(aspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!------------------third group: coordinates-----------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
!
      call h5screate_simple_f(1, dims_p, dspace_id, err)
         call h5dcreate_f(group_id, 'p', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, p, dims_p, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)
!
      call h5screate_simple_f(1, dims_zeta, dspace_id, err)
         call h5dcreate_f(group_id, 'zeta', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, zeta, dims_zeta, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)

   call h5gclose_f(group_id, err)
!
!----------------fourth group: intensity on surface---------------------
!
   call h5gcreate_f(file_id, 'surfb', group_id, err)
!
      call h5screate_simple_f(2, dims_2d, dspace_id, err)
!
         call h5dcreate_f(group_id, 'iem_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iem_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)

         call h5dcreate_f(group_id, 'iemi_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iemi_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)

         call h5dcreate_f(group_id, 'iabs_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iabs_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)

         call h5dcreate_f(group_id, 'icont_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, icont_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)            
!
      call h5sclose_f(dspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!
!
call h5fclose_f(file_id, err)
!
call h5close_f(err)
!

end subroutine output_surface
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_int2d(indx, xobs, alpha, gamma)
!
!-------------output of surface emergent intensity for given------------
!--------------xobs and direction specified by alpha, gamma-------------
!
use prog_type
use options_spec, only: output_dir
use dime_spec, only: np, nzeta
use mod_spectrum, only: p, zeta
use mod_int2d, only: nz_ray_max, int_2d, tau_2d, zcoord_2d, xcoord_2d, iemi_2d, iabs_2d, vn_2d
use params_spec, only: xic1, obliquity, vth_fiducial, xnue0
use hdf5
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx
real(dp), intent(in) :: xobs, alpha, gamma
!
! ... local scalars
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, attr_id
!
! ... local arrays
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(2) :: dims_2d


! ... local characters
character(len=8) :: name
character(len=12) :: fname_base='/spec_int2d_'
character(len=512) :: fname
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
fname = trim(output_dir)//fname_base//trim(name)//'.h5'
!
!
!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(fname)
!
!-----------------------------------------------------------------------
!
dims_2d = (/ np, nz_ray_max /)
!
!-----------------------------------------------------------------------
!
!initialize fortran interface
call h5open_f (err)
!
!create hdf5 file
call h5fcreate_f(trim(fname), h5f_acc_trunc_f, file_id, err)
!
!---------------------------first group: dimensions---------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!
         call h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, np, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nz_ray_max, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
      call h5sclose_f(aspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!---------------------------second group: parameter---------------------
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!
         call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'xnue0', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'vth_fiducial', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
         call h5aclose_f(attr_id, err)            
!
         call h5acreate_f(group_id, 'xobs', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xobs, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5acreate_f(group_id, 'gamma', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, gamma, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
      !close dataspace
      call h5sclose_f(aspace_id, err)
!
   !close second group
   call h5gclose_f(group_id, err)
!
!------------------third group: coordinates-----------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
!
      call h5screate_simple_f(2, dims_2d, dspace_id, err)
         call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xcoord_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
         call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, zcoord_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)

   call h5gclose_f(group_id, err)
!
!----------------fourth group: intensity and optical depth--------------
!
   call h5gcreate_f(file_id, 'along_ray', group_id, err)
!
      call h5screate_simple_f(2, dims_2d, dspace_id, err)
!
         call h5dcreate_f(group_id, 'int2d', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, int_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
         call h5dcreate_f(group_id, 'iemi2d', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iemi_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
         call h5dcreate_f(group_id, 'iabs2d', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iabs_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
         call h5dcreate_f(group_id, 'tau2d', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, tau_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
         call h5dcreate_f(group_id, 'vn2d', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, vn_2d, dims_2d, err)
         call h5dclose_f(dset_id, err)
!
      call h5sclose_f(dspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!
!
call h5fclose_f(file_id, err)
!
call h5close_f(err)
!

end subroutine output_int2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_fluxem(indx)
!
!-----------------output of emergent flux profile-----------------------
!
use prog_type
use options_spec, only: output_dir
use fund_const, only: pi
use dime_spec, only: nxobs_fs, nalpha, ngamma
use mod_spectrum, only: flux_tot, flux_cont, xobs, alpha, gamma
use params_spec, only: vmax, vth_fiducial
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx
!
! ... local scalars
integer(i4b) :: i
!
! ... local characters
character(len=8) :: name
!
!

write(*,*) '------------------------------output to file-----------------------------------'
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
write(*,*) trim(output_dir)//'/FLUXEM_'//trim(name)//'.dat'
!
!xobs in units of vmax, instead of vthfiducial
!
open(1, file=trim(output_dir)//'/FLUXEM_'//trim(name)//'.dat', form='formatted')
   write(1,'(i5)') nxobs_fs 
   write(1,'(a6, e20.8)') 'alpha=', alpha
   write(1,'(a6, e20.8)') 'gamma=', gamma
   write(1,'(4(a20))') 'xobs[vth*]', 'flux_tot', 'flux_cont', 'f_tot/f_cont'
   do i=1, nxobs_fs
      write(1,'(4(e20.8))')  xobs(i), flux_tot(i), flux_cont(i), flux_tot(i)/flux_cont(i)
   enddo
close(1)
!
end subroutine output_fluxem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_fluxem_debug(indx)
!
!-----------------output of emergent flux profile-----------------------
!
use prog_type
use options_spec, only: output_dir
use fund_const, only: pi
use dime_spec, only: nxobs_fs, nalpha, ngamma
use mod_spectrum, only: flux_tot, flux_cont, xobs, alpha, gamma, flux_emi, flux_abs
use params_spec, only: vmax, vth_fiducial
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx
!
! ... local scalars
integer(i4b) :: i
real(dp) :: ftot, fcont, femi, fabs
!
! ... local characters
character(len=8) :: name
!
!
write(*,*) '------------------------------output to file-----------------------------------'
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
write(*,*) trim(output_dir)//'/FLUXEM_DEBUG_'//trim(name)//'.dat'
!
!xobs in units of vmax, instead of vthfiducial
!
open(1, file=trim(output_dir)//'/FLUXEM_DEBUG_'//trim(name)//'.dat', form='formatted')
   write(1,'(i5)') nxobs_fs 
   write(1,'(a6, e20.8)') 'alpha=', alpha
   write(1,'(a6, e20.8)') 'gamma=', gamma
   write(1,'(5(a20))') 'xobs[vth*]', 'flux_tot', 'flux_cont', 'flux_emi', 'flux_abs'
   do i=1, nxobs_fs
      ftot=flux_tot(i)
      fcont=flux_cont(i)
      femi=flux_emi(i)
      fabs=flux_abs(i)
      if(ftot.lt.1.d-15) ftot=0.d0
      if(fcont.lt.1.d-15) fcont=0.d0
      if(femi.lt.1.d-80) femi=0.d0
      if(fabs.lt.1.d-80) fabs=0.d0
      write(1,'(5(e20.8))')  xobs(i), ftot, fcont, femi, fabs
   enddo
close(1)
!
end subroutine output_fluxem_debug
