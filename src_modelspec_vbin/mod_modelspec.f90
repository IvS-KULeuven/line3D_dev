!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options_modspec
!
use prog_type
!
implicit none
!
character(len=100) :: indat_file, input_file, input_file2, output_file
!
integer(i4b) :: input_mod
!input_mod=0  3D model: calculate stuff
!             output to two spherical coordinate systems
!                       
!
!
end module options_modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_modspec
!
use prog_type
!
implicit none
!
!-----------------dimensions of spherical 3d grids-------------------
!------------------------(one for each star)-------------------------
!
integer(i4b) :: cs1_nr, cs1_ntheta, cs1_nphi
real(dp), dimension(:), allocatable :: cs1_r, cs1_theta, cs1_phi
real(dp), dimension(:,:,:), allocatable :: cs1_opalbar3d, cs1_opac3d, cs1_sline3d, cs1_scont3d, &
                                           cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d, cs1_rho3d
integer(i4b), dimension(:,:,:), allocatable :: cs1_imask3d
!
!
integer(i4b) :: cs2_nr, cs2_ntheta, cs2_nphi
real(dp), dimension(:), allocatable :: cs2_r, cs2_theta, cs2_phi
real(dp), dimension(:,:,:), allocatable :: cs2_opalbar3d, cs2_opac3d, cs2_sline3d, cs2_scont3d, &
                                           cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d, cs2_rho3d
integer(i4b), dimension(:,:,:), allocatable :: cs2_imask3d
!
!
!
end module dime_modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_modspec
!
!
use prog_type
!
implicit none
!
real(dp) :: x01, y01, z01, rstar1, rmin1, rmax1, sr1, teff1, trad1, vmicro1, &
     vx01, vy01, vz01, vmin1, vmax1, vrot1, beta1, mdot1, logg1, lstar1, yhe1, yhe1_mass, hei1, &
     fehe1, aenh1
real(dp), dimension(3) :: ex01, ey01, ez01, rot_axis01
real(dp) :: x02, y02, z02, rstar2, rmin2, rmax2, sr2, teff2, trad2, vmicro2, &
     vx02, vy02, vz02, vmin2, vmax2, vrot2, beta2, mdot2, logg2, lstar2, yhe2, yhe2_mass, hei2, &
     fehe2, aenh2
real(dp), dimension(3) :: ex02, ey02, ez02, rot_axis02
!
real(dp) :: eps_line
!
real(dp) :: vth_fiducial, vmax, unit_length
!
!
!
end module params_modspec
