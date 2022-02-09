!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_surface(xobs, alpha, gamma)
!
!-------------output of surface emergent intensity for given------------
!--------------xobs and direction specified by alpha, gamma-------------
!
use prog_type
use options_spec, only: output_dir
use mod_surfb, only: iem_surface, iemi_surface, iabs_surface
use mod_triangles, only: npoints, points_xcoord, points_ycoord
use params_spec, only: xic1, xic2
use hdf5
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, alpha, gamma
!
! ... local scalars
integer(i4b) :: err, i
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, attr_id
!
! ... local characters
character(len=17), parameter :: fname='spec_surface.h5'
!
! ... local arrays
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_points
!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(output_dir)//'/'//fname
!
!-----------------------------------------------------------------------
!
dims_points = (/ npoints /)
!
!do i=4018, 4018
!   write(*,*) points_xcoord(i), points_ycoord(i)
!enddo


!-----------------------------------------------------------------------
!
!initialize fortran interface
call h5open_f (err)
!
!create hdf5 file
call h5fcreate_f(trim(output_dir)//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!---------------------------first group: dimensions---------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'npoints', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, npoints, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!---------------------------second group: parameter---------------------
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'xic1', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'xic2', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xic2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'xobs', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xobs, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'alpha', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'gamma', h5t_native_real, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, gamma, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
!
   call h5gclose_f(group_id, err)
!
!------------------third group: coordinates-----------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
!
      call h5screate_simple_f(1, dims_points, dspace_id, err)
         call h5dcreate_f(group_id, 'points_xcoord', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, points_xcoord, dims_points, err)
            call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id, 'points_ycoord', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, points_ycoord, dims_points, err)
         call h5dclose_f(dset_id, err)            
      call h5sclose_f(dspace_id, err)

      call h5gclose_f(group_id, err)

!      write(*,*) points_xcoord(2843), points_ycoord(2843)
!
!----------------fourth group: intensity on surface---------------------
!
   call h5gcreate_f(file_id, 'surfaceb', group_id, err)
!
      call h5screate_simple_f(1, dims_points, dspace_id, err)
         call h5dcreate_f(group_id, 'iem_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iem_surface, dims_points, err)
         call h5dclose_f(dset_id, err)

         call h5dcreate_f(group_id, 'iemi_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iemi_surface, dims_points, err)
         call h5dclose_f(dset_id, err)

         call h5dcreate_f(group_id, 'iabs_surface', h5t_native_real, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iabs_surface, dims_points, err)
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
subroutine output_int2d(xobs, alpha, gamma)
!
!-------------output of surface emergent intensity for given------------
!--------------xobs and direction specified by alpha, gamma-------------
!
use prog_type
use options_spec, only: output_dir
!use dime_spec, only: np, nzeta
!use mod_spectrum, only: p, zeta
!use mod_int2d, only: nz_ray_max, int_2d, tau_2d, zcoord_2d, xcoord_2d, iemi_2d, iabs_2d, vn_2d
!!use params_spec, only: xic1, obliquity
!use hdf5_fstructure
use hdf5
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, alpha, gamma
!
! ... local scalars
real(dp) :: xic1, obliquity
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, attr_id
!
! ... local characters
character(len=15), parameter :: fname='spec_int2d.h5'

! ... local arrays
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(2) :: dims_2d
!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(output_dir)//'/'//fname
!
!-----------------------------------------------------------------------
!
!dims_2d = (/ np, nz_ray_max /)
!
!-----------------------------------------------------------------------
!
!initialize fortran interface
call h5open_f (err)
!
!create hdf5 file
!call h5fcreate_f(trim(output_dir)//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!---------------------------first group: dimensions---------------------
!
!   call h5gcreate_f(file_id, gname_dim, group_id, err)
!
!create the dataspace to store dimensions (each dimension-parameter stored as a single scalar)
!      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!!
!         !create attribute for n1d attached to group
!         call h5acreate_f(group_id, aname_nx, h5t_native_integer, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_integer, np, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!         !create attribute for ntheta attached to group (same dataspace as n1d)
!         call h5acreate_f(group_id, aname_nz, h5t_native_integer, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_integer, nz_ray_max, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!      !close dataspace
!      call h5sclose_f(aspace_id, err)
!!
!   !close first group
!   call h5gclose_f(group_id, err)
!!
!!---------------------------second group: parameter---------------------
!!
!   call h5gcreate_f(file_id, gname_param, group_id, err)
!!
!!create the dataspace to store parameter (each parameter stored as a single scalar)
!      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!!
!         !create attribute for xic1 attached to group
!         call h5acreate_f(group_id, aname_xic, h5t_native_real, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!         !create attribute for xobs
!         call h5acreate_f(group_id, aname_xobs, h5t_native_real, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_double, xobs, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!         !create attribute for obliquity
!         call h5acreate_f(group_id, aname_obliquity, h5t_native_real, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_double, obliquity, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!         !create attribute for alpha
!         call h5acreate_f(group_id, aname_alpha, h5t_native_real, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_double, alpha, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
!         !create attribute for obliquity
!         call h5acreate_f(group_id, aname_gamma, h5t_native_real, aspace_id, &
!         attr_id, err)
!            !write attribute to file
!            call h5awrite_f(attr_id, h5t_native_double, gamma, dims_scalars, err)
!         !close attribute
!         call h5aclose_f(attr_id, err)
!!
      !close dataspace
!      call h5sclose_f(aspace_id, err)
!!
!   !close second group
!   call h5gclose_f(group_id, err)
!!
!!------------------third group: coordinates-----------------------------
!!
!   call h5gcreate_f(file_id, gname_coord, group_id, err)
!!
!      !p-coordinate
!      call h5screate_simple_f(2, dims_2d, dspace_id, err)
!         call h5dcreate_f(group_id, dname_xgrid, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, xcoord_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!      call h5sclose_f(dspace_id, err)
!!
!      !zeta-oordinate
!      call h5screate_simple_f(2, dims_2d, dspace_id, err)
!         call h5dcreate_f(group_id, dname_zgrid, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, zcoord_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!      call h5sclose_f(dspace_id, err)
!
!   !close third group
!   call h5gclose_f(group_id, err)
!!
!!----------------fourth group: intensity and optical depth--------------
!!
!   call h5gcreate_f(file_id, gname_aray, group_id, err)
!!
!      call h5screate_simple_f(2, dims_2d, dspace_id, err)
!!
!         !intensity along ray
!         call h5dcreate_f(group_id, dname_int2d, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, int_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!!
!         !intensity along ray (emission part)
!         call h5dcreate_f(group_id, dname_iemi2d, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, iemi_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!!
!         !intensity along ray (absorption part)
!         call h5dcreate_f(group_id, dname_iabs2d, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, iabs_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!!
!         !optical depth along ray
!         call h5dcreate_f(group_id, dname_tau2d, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, tau_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!!
!         !velocity along ray
!         call h5dcreate_f(group_id, dname_vn2d, h5t_native_real, dspace_id, &
!         dset_id, err)
!            !write the dataset to file
!            call h5dwrite_f(dset_id, h5t_native_double, vn_2d, dims_2d, err)
!         call h5dclose_f(dset_id, err)
!!
!      call h5sclose_f(dspace_id, err)
!!
!   call h5gclose_f(group_id, err)
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
subroutine output_triangles(indx)
!
!-----------------output of triangulation for the given direction-------
!
use prog_type
use options_spec, only: output_dir
use mod_triangles, only: npoints, ntriangles, points_xcoord, points_ycoord, points_indx, triangles_ip
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx
!
! ... local scalars
integer(i4b) :: i, ip1, ip2, ip3, findx
!
! ... local characters
character(len=8) :: name
character(len=200) :: fname
!
!
write(*,*) '-------------------------output triangles to files-----------------------------'
findx=200*indx
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
fname=trim(output_dir)//'/spec_points_'//trim(name)//'.dat'
write(*,*) trim(fname)
!
open(findx, file=trim(fname), form='formatted')
write(findx,'(2a10,2a20)') 'index', 'pindex', 'xcoord', 'ycoord'
do i=1, npoints
!   write(*,'(2i10,2es20.8)') i, points_indx(i), points_xcoord(i), points_ycoord(i)   
   write(findx,'(2i10,2es20.8)') i, points_indx(i), points_xcoord(i), points_ycoord(i)
enddo
close(findx)
!

!do i=4000, 4100
!   write(*,*) i, points_indx(i), points_xcoord(i), points_ycoord(i)
!enddo
!stop
!
!
fname=trim(output_dir)//'/spec_triangles_'//trim(name)//'.dat'
write(*,*) trim(fname)
!
open(findx, file=trim(fname), form='formatted')
write(findx,'(3a10,6a20)') 'ip1', 'ip2', 'ip3', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3'
do i=1, ntriangles
   ip1 = triangles_ip(1,i)
   ip2 = triangles_ip(2,i)
   ip3 = triangles_ip(3,i)   
   write(findx,'(3i10,6es20.8)') ip1, ip2, ip3, points_xcoord(ip1), points_ycoord(ip1), &
                                           points_xcoord(ip2), points_ycoord(ip2), &
                                           points_xcoord(ip3), points_ycoord(ip3)
enddo
close(findx)
!
!
end subroutine output_triangles
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
use fund_const, only: pi, zero
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
integer(i4b) :: i, findx
!
! ... local characters
character(len=8) :: name
!
findx=200*indx
!
write(*,*) '------------------------------output to file-----------------------------------'
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
write(*,*) trim(output_dir)//'/FLUXEM_'//trim(name)//'.dat'
!
!xobs in units of vmax, instead of vthfiducial
!
open(findx, file=trim(output_dir)//'/FLUXEM_'//trim(name)//'.dat', form='formatted')
   write(findx,'(i5)') nxobs_fs 
   write(findx,'(a6, e20.8)') 'alpha=', alpha
   write(findx,'(a6, e20.8)') 'gamma=', gamma
   write(findx,'(4(a20))') 'xobs[vth*]', 'flux_tot', 'flux_cont', 'f_tot/f_cont'
   do i=1, nxobs_fs
      if(flux_tot(i).lt.1.d-15) flux_tot(i)=zero
      if(flux_cont(i).lt.1.d-15) flux_cont(i)=zero
      write(findx,'(4(e20.8))')  xobs(i), flux_tot(i), flux_cont(i), flux_tot(i)/flux_cont(i)
   enddo
close(findx)
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
integer(i4b) :: i, findx
real(dp) :: ftot, fcont, femi, fabs
!
! ... local characters
character(len=8) :: name
!
findx=200*indx
!
write(*,*) '------------------------------output to file-----------------------------------'
!
!store indx in string-variable 'name' with 5 digits and 5 zeros to the left
write(name,'(i5.5)') indx
write(*,*) trim(output_dir)//'/FLUXEM_DEBUG_'//trim(name)//'.dat'
!
!xobs in units of vmax, instead of vthfiducial
!
open(findx, file=trim(output_dir)//'/FLUXEM_DEBUG_'//trim(name)//'.dat', form='formatted')
   write(findx,'(i5)') nxobs_fs 
   write(findx,'(a6, e20.8)') 'alpha=', alpha
   write(findx,'(a6, e20.8)') 'gamma=', gamma
   write(findx,'(5(a20))') 'xobs[vth*]', 'flux_tot', 'flux_cont', 'flux_emi', 'flux_abs'
   do i=1, nxobs_fs
      ftot=flux_tot(i)
      fcont=flux_cont(i)
      femi=flux_emi(i)
      fabs=flux_abs(i)
      if(ftot.lt.1.d-15) ftot=0.d0
      if(fcont.lt.1.d-15) fcont=0.d0
      if(femi.lt.1.d-80) femi=0.d0
      if(fabs.lt.1.d-80) fabs=0.d0
      write(findx,'(5(e20.8))')  xobs(i), ftot, fcont, femi, fabs
   enddo
close(findx)
!
end subroutine output_fluxem_debug
