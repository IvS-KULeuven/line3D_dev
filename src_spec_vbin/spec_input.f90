!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3d_spc
!
use prog_type
use fund_const
use options_spec, only: input_file
use dime_model3d, only: cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_r_spc, cs1_theta_spc, cs1_phi_spc, &
                        cs1_t3d, cs1_vth3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, &
                        cs1_opac3d, cs1_opalbar3d, cs1_scont3d, cs1_sline3d, cs1_imask3d, &
                        cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_r_spc, cs2_theta_spc, cs2_phi_spc, &
                        cs2_t3d, cs2_vth3d, cs2_velx3d, cs2_vely3d, cs2_velz3d, &
                        cs2_opac3d, cs2_opalbar3d, cs2_scont3d, cs2_sline3d, cs2_imask3d
use params_spec, only: vx01, vy01, vz01, x01, y01, z01, teff1, trad1, lstar1, logg1, yhe1, fehe1, aenh1, rstar1, rmin1, rmax1, sr1, xic1, vrot1, vmicro1, &
                       vx02, vy02, vz02, x02, y02, z02, teff2, trad2, lstar2, logg2, yhe2, fehe2, aenh2, rstar2, rmin2, rmax2, sr2, xic2, vrot2, vmicro2, &
                       ex01, ey01, ez01, rot_axis01, &
                       ex02, ey02, ez02, rot_axis02, &
                       vth_fiducial, na, vmax, tmin, xnue0, iline
use params_model, only: vth_fiducial_model, vmax_model
use mod_spectrum, only: unit_length
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: tmin1, tmin2
!
! ... local logicals
!
! ... for hdf5
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id, group_id0
integer(hsize_t), dimension(1), parameter :: dims_scalars = (/ 1 /), dims_vectors = (/ 3 /)
integer(hsize_t), dimension(1) :: dims_radius, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
! ... local function
real(dp) :: bnue
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
!----------------------------input parameters---------------------------
!
!x01, y01, z01, x02, y02, z02 in unit_length   
      call h5gopen_f(file_id, 'input_parameters', group_id, err)
         call h5aopen_f(group_id, 'x01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, x01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'y01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, y01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'z01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, z01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'x02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, x02, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'y02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, y02, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'z02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, z02, dims_scalars, err)
         call h5aclose_f(attr_id, err)

!vx01, vy01, vz01, vx02, vy02, vz02 (cm/s)
         call h5aopen_f(group_id, 'vx01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vx01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vy01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vy01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vz01', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vz01, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vx02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vx02, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vy02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vy02, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vz02', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vz02, dims_scalars, err)
         call h5aclose_f(attr_id, err)            

!temperatures in kelvin
         call h5aopen_f(group_id, 'teff1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'teff2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, teff2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'trad2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, trad2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'logg1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, logg1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'logg2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, logg2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'lstar1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'lstar2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, lstar2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'yhe1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, yhe1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'yhe2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, yhe2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'fehe1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, fehe1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'fehe2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, fehe2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'aenh1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, aenh1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'aenh2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, aenh2, dims_scalars, err)
         call h5aclose_f(attr_id, err)

!rstar1, rstar2 in r_sun
!rmin1, rmax1 in rstar1
!rmin2, rmax2 in rstar2 
         call h5aopen_f(group_id, 'rstar1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rstar2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rstar2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rmin1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rmin1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rmin2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rmin2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rmax1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rmax1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'rmax2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, rmax2, dims_scalars, err)
         call h5aclose_f(attr_id, err)

!velocities in  in cm/s
         call h5aopen_f(group_id, 'vrot1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vrot2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vrot2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
         call h5aopen_f(group_id, 'vth_fiducial', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vth_fiducial_model, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmicro1', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'vmicro2', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmicro2, dims_scalars, err)
         call h5aclose_f(attr_id, err)
            

!vmax is the maximum found velocity in the global coordinate system
         call h5aopen_f(group_id, 'vmax', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
!unit length of global coordinate system [rsun]
         call h5aopen_f(group_id, 'unit_length', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, unit_length, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         
         call h5aopen_f(group_id, 'xnue0', attr_id, err)
            call h5aread_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'na', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, na, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'iline', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, iline, dims_scalars, err)
         call h5aclose_f(attr_id, err)
!
!coordinate basis of individual coordinate systems (w.r.t global system)
         call h5dopen_f(group_id, 'ex01', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ex01, dims_vectors, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'ey01', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ey01, dims_vectors, err)
         call h5dclose_f(dset_id, err)                     
         call h5dopen_f(group_id, 'ez01', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ez01, dims_vectors, err)
         call h5dclose_f(dset_id, err)         
         call h5dopen_f(group_id, 'ex02', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ex02, dims_vectors, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'ey02', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ey02, dims_vectors, err)
         call h5dclose_f(dset_id, err)                     
         call h5dopen_f(group_id, 'ez02', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, ez02, dims_vectors, err)
         call h5dclose_f(dset_id, err)         

!rotation axis of individual stars (w.r.t global system)
         call h5dopen_f(group_id, 'rot_axis01', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rot_axis01, dims_vectors, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'rot_axis02', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rot_axis02, dims_vectors, err)
         call h5dclose_f(dset_id, err)                     
      call h5gclose_f(group_id, err)
!
      if(rot_axis01(3).ne.one) stop 'error in read_model3d_spc: rot_axis01 ne (0,0,1) not implemented yet'
      if(rot_axis02(3).ne.one) stop 'error in read_model3d_spc: rot_axis02 ne (0,0,1) not implemented yet'      
!
      sr1=rstar1*rsu
      sr2=rstar2*rsu      
!
!--------------------------star 1---------------------------------------
!
!dimensions
      call h5gopen_f(file_id, 'star1', group_id0, err)      
         call h5gopen_f(group_id0, 'dimensions', group_id, err)
            call h5aopen_f(group_id, 'nr', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs1_nr_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
            call h5aopen_f(group_id, 'ntheta', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs1_ntheta_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
            call h5aopen_f(group_id, 'nphi', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs1_nphi_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!set dimension-arrays
         dims_radius=(/ cs1_nr_spc /)
         dims_theta=(/ cs1_ntheta_spc /)
         dims_phi=(/ cs1_nphi_spc /)
         dims_3d=(/cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc /)
!
!
!
!coordinates
!allocate arrays
         allocate(cs1_r_spc(cs1_nr_spc), stat=err)
         allocate(cs1_theta_spc(cs1_ntheta_spc), stat=err)
         allocate(cs1_phi_spc(cs1_nphi_spc), stat=err)
   
         allocate(cs1_t3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_sline3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_scont3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_velx3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_vely3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_velz3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_opac3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_opalbar3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
         allocate(cs1_vth3d(cs1_nr_spc,cs1_ntheta_spc,cs1_nphi_spc), stat=err)
   
         call h5gopen_f(group_id0, 'coordinates', group_id, err)
            call h5dopen_f(group_id, 'r', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_r_spc, dims_radius, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'theta', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_theta_spc, dims_theta, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'phi', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_phi_spc, dims_phi, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!3d solution
         call h5gopen_f(group_id0, 'solution3d', group_id, err)
            call h5dopen_f(group_id, 'sline3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_sline3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'scont3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_scont3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!3d model
         call h5gopen_f(group_id0, 'model3d', group_id, err)
            call h5dopen_f(group_id, 't3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_t3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'opac3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_opac3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_opalbar3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'velx3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_velx3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'vely3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_vely3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'velz3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs1_velz3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!
      call h5gclose_f(group_id0, err)
!
!--------------------------star 2---------------------------------------
!
!dimensions
      call h5gopen_f(file_id, 'star2', group_id0, err)      
         call h5gopen_f(group_id0, 'dimensions', group_id, err)
            call h5aopen_f(group_id, 'nr', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs2_nr_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
            call h5aopen_f(group_id, 'ntheta', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs2_ntheta_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
            call h5aopen_f(group_id, 'nphi', dset_id, err)
               call h5aread_f(dset_id, h5t_native_integer, cs2_nphi_spc, dims_scalars, err)
            call h5aclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!set dimension-arrays
         dims_radius=(/ cs2_nr_spc /)
         dims_theta=(/ cs2_ntheta_spc /)
         dims_phi=(/ cs2_nphi_spc /)
         dims_3d=(/cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc /)
!
!
!
!coordinates
!allocate arrays
         allocate(cs2_r_spc(cs2_nr_spc), stat=err)
         allocate(cs2_theta_spc(cs2_ntheta_spc), stat=err)
         allocate(cs2_phi_spc(cs2_nphi_spc), stat=err)
   
         allocate(cs2_t3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_sline3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_scont3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_velx3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_vely3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_velz3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_opac3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_opalbar3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
         allocate(cs2_vth3d(cs2_nr_spc,cs2_ntheta_spc,cs2_nphi_spc), stat=err)
   
         call h5gopen_f(group_id0, 'coordinates', group_id, err)
            call h5dopen_f(group_id, 'r', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_r_spc, dims_radius, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'theta', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_theta_spc, dims_theta, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'phi', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_phi_spc, dims_phi, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!3d solution
         call h5gopen_f(group_id0, 'solution3d', group_id, err)
            call h5dopen_f(group_id, 'sline3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_sline3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'scont3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_scont3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!3d model
         call h5gopen_f(group_id0, 'model3d', group_id, err)
            call h5dopen_f(group_id, 't3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_t3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'opac3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_opac3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'opalbar3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_opalbar3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'velx3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_velx3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'vely3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_vely3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
            call h5dopen_f(group_id, 'velz3d', dset_id, err)
               call h5dread_f(dset_id, h5t_native_double, cs2_velz3d, dims_3d, err)
            call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
!
!
      call h5gclose_f(group_id0, err)
!      
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
!transform opacity to (eventually) different vth_fiducial, and to 1/unit_length
!(on input in 1/cm)
cs1_opalbar3d = cs1_opalbar3d*vth_fiducial_model/vth_fiducial * unit_length*rsu
cs2_opalbar3d = cs2_opalbar3d*vth_fiducial_model/vth_fiducial * unit_length*rsu
cs1_opac3d = cs1_opac3d * unit_length*rsu
cs2_opac3d = cs2_opac3d * unit_length*rsu
!
!measure velocities in vth_fiducial (except vmax, vmicro)
cs1_velx3d=cs1_velx3d/vth_fiducial
cs1_vely3d=cs1_vely3d/vth_fiducial
cs1_velz3d=cs1_velz3d/vth_fiducial
cs2_velx3d=cs2_velx3d/vth_fiducial
cs2_vely3d=cs2_vely3d/vth_fiducial
cs2_velz3d=cs2_velz3d/vth_fiducial
vx01=vx01/vth_fiducial
vy01=vy01/vth_fiducial
vz01=vz01/vth_fiducial
vx02=vx02/vth_fiducial
vy02=vy02/vth_fiducial
vz02=vz02/vth_fiducial
vrot1=vrot1/vth_fiducial
vrot2=vrot2/vth_fiducial
!
!
!minimum temperature
tmin1=minval(cs1_t3d)
tmin2=minval(cs2_t3d)
tmin=minval((/tmin1,tmin2/))
!
xic1=bnue(xnue0,trad1)
xic2=bnue(xnue0,trad2)
!
!xic1=one
!xic2=one/two
!test
!cs1_opac3d=0.d0
!cs1_scont3d=0.d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0


!write(*,*) cs2_r_spc

!stop 'warning: dont forget to put on again'
!
end subroutine read_model3d_spc
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
use options_spec, only: input_file, output_dir, input_mod, opt_photprof1, opt_photprof2, &
                        opt_obsdir_read, opt_surface, opt_int2d, opt_incl_gdark1, opt_incl_sdist1, &
                        opt_incl_gdark2, opt_incl_sdist2, &
                        opt_pgrid01, opt_pgrid02, opt_rgrid01, opt_rgrid02
use params_spec, only: vth_fiducial
use mod_surfb, only: xobs_surface, alpha_surface, gamma_surface
use dime_spec, only: nalpha, ngamma
!
implicit none
!
! ... local scalars
!
! ... local characters
character(len=100) :: fname
!
! ... namelist
namelist / input_options / input_mod, input_file, output_dir, opt_photprof1, opt_photprof2, &
                           opt_obsdir_read, opt_surface, opt_int2d, opt_incl_gdark1, opt_incl_sdist1, &
                           opt_incl_gdark2, opt_incl_sdist2, opt_pgrid01, opt_pgrid02, &
                           opt_rgrid01, opt_rgrid02, nalpha, ngamma
namelist / input_model / vth_fiducial
namelist / input_surface / alpha_surface, gamma_surface, xobs_surface
!
! ... local functions
!
!----------------------read input/output directories--------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere and options'
read(*,*) fname
write(*,*) 'reading input from file: ', trim(fname)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(fname), status='old', form='formatted')
!
   read(1, nml=input_options)
!
   read(1, nml=input_model)
   vth_fiducial = vth_fiducial * 1.d5
!
   read(1, nml=input_surface)
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
subroutine print_model3d_spc
!
use prog_type
use dime_model3d, only: cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_r_spc, cs1_theta_spc, cs1_phi_spc, &
                        cs1_opac3d, cs1_opalbar3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_sline3d, cs1_scont3d, cs1_t3d, &
                        cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_r_spc, cs2_theta_spc, cs2_phi_spc, &
                        cs2_opac3d, cs2_opalbar3d, cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_sline3d, cs2_scont3d, cs2_t3d
use params_spec, only: xnue0, vth_fiducial, na, iline, tmin, vmax, &
                       vx01, vy01, vz01, x01, y01, z01, lstar1, logg1, yhe1, fehe1, aenh1, teff1, trad1, rstar1, vrot1, xic1, vmicro1, &
                       vx02, vy02, vz02, x02, y02, z02, lstar2, logg2, yhe2, fehe2, aenh2, teff2, trad2, rstar2, vrot2, xic2, vmicro2
!
implicit none
!
! ... arguments
integer(i4b) :: i, j, k
!
write(*,*) '------3d atmospheric structure (spc) as fct of radius at theta, phi=0, 0-------'
write(*,*)

write(*,'(a20, i20)') 'iline', iline
write(*,'(a20, es20.8)') 'xnue0', xnue0
write(*,'(a20, es20.8)') 'tmin [K]', tmin
write(*,'(a20, i20)') 'na', na
write(*,'(a20, es20.8)') 'vmax [km/s]', vmax/1.d5
write(*,'(a20, es20.8)') 'vth_fiducial [km/s]', vth_fiducial/1.d5
write(*,*)
!
write(*,*) '------------------star 1--------------------------'
write(*,'(a20, es20.8)') 'I_c', xic1
write(*,'(a20, es20.8)') 'vrot [km/s]', vrot1*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro1/1.d5
write(*,'(a20, es20.8)') 'rstar [r_sun]', rstar1
write(*,'(a20, es20.8)') 'lstar [l_sun]', lstar1
write(*,'(a20, es20.8)') 'teff', teff1
write(*,'(a20, es20.8)') 'trad', trad1
write(*,'(a20, es20.8)') 'logg', logg1
write(*,'(a20, es20.8)') 'yhe', yhe1
write(*,'(a20, es20.8)') 'fehe', fehe1
write(*,'(a20, es20.8)') 'aenh', aenh1
write(*,'(a20, es20.8)') 'x01 [unit_l]', x01
write(*,'(a20, es20.8)') 'y01 [unit_l]', y01
write(*,'(a20, es20.8)') 'z01 [unit_l]', z01
write(*,'(a20, es20.8)') 'vx01 [km/s]', vx01*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'vy01 [km/s]', vy01*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'vz01 [km/s]', vz01*vth_fiducial/1.d5
write(*,*)
!
j=1
k=1
!
write(*,'(a4, 9(a20))') '#', 'r [r_star1]', 'opac [1/unit_l]', 'opalbar [1/unit_l]', 'velx [vth*]', &
                        'vely [vth*]', 'velz [vth*]', 't [K]', 's_cont', 's_line'
do i=1, cs1_nr_spc
   write(*,'(i4, 9(es20.8))')  i, cs1_r_spc(i), cs1_opac3d(i,j,k), cs1_opalbar3d(i,j,k), cs1_velx3d(i,j,k), cs1_vely3d(i,j,k), &
                              cs1_velz3d(i,j,k), cs1_t3d(i,j,k), cs1_scont3d(i,j,k), cs1_sline3d(i,j,k)
enddo
write(*,*)


write(*,*) '------------------star 2--------------------------'
write(*,'(a20, es20.8)') 'I_c', xic2
write(*,'(a20, es20.8)') 'vrot [km/s]', vrot2*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro2/1.d5
write(*,'(a20, es20.8)') 'rstar [r_sun]', rstar2
write(*,'(a20, es20.8)') 'lstar [l_sun]', lstar2
write(*,'(a20, es20.8)') 'teff', teff2
write(*,'(a20, es20.8)') 'trad', trad2
write(*,'(a20, es20.8)') 'logg', logg2
write(*,'(a20, es20.8)') 'yhe', yhe2
write(*,'(a20, es20.8)') 'fehe', fehe2
write(*,'(a20, es20.8)') 'aenh', aenh2
write(*,'(a20, es20.8)') 'x02 [unit_l]', x02
write(*,'(a20, es20.8)') 'y02 [unit_l]', y02
write(*,'(a20, es20.8)') 'z02 [unit_l]', z02
write(*,'(a20, es20.8)') 'vx01 [km/s]', vx02*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'vy01 [km/s]', vy02*vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'vz01 [km/s]', vz02*vth_fiducial/1.d5
write(*,*)
!
j=1
k=1
!
write(*,'(a4, 9(a20))') '#', 'r [r_star2]', 'opac [1/unit_l]', 'opalbar [1/unit_l]', 'velx [vth*]', &
                        'vely [vth*]', 'velz [vth*]', 't [K]', 's_cont', 's_line'
do i=1, cs2_nr_spc
   write(*,'(i4, 9(es20.8))')  i, cs2_r_spc(i), cs2_opac3d(i,j,k), cs2_opalbar3d(i,j,k), cs2_velx3d(i,j,k), cs2_vely3d(i,j,k), &
                              cs2_velz3d(i,j,k), cs2_t3d(i,j,k), cs2_scont3d(i,j,k), cs2_sline3d(i,j,k)
enddo
write(*,*)


!
end subroutine print_model3d_spc
