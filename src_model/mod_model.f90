!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_directories
!
implicit none
!
character(len=100) :: indat_file
!
character(len=300) :: model_dir
!
character(len=10), parameter :: model1d_file='model1d.h5'
character(len=10), parameter :: model2d_file='model2d.h5'
character(len=10), parameter :: model3d_file='model3d.h5'
!modelxd_file:   file where model-atmosphere is stored
!
end module mod_directories
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_modext
!
!--------------------dimensions of external atmoshpere------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: nr_modext, ntheta_modext, nphi_modext
!
!
end module dime_modext
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model1d
!
!----------------------1d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: r_modext1d, velr_modext1d, &
                                       rho_modext1d, t_modext1d, vth_modext1d, &
                                       eps_cont_modext1d
!
!
end module model1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model2d
!
!----------------------2d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: r_modext2d, theta_modext2d
real(dp), dimension(:,:), allocatable :: velr_modext2d, velth_modext2d, velphi_modext2d, &
                                         rho_modext2d, t_modext2d, vth_modext2d, &
                                         eps_cont_modext2d
!
end module model2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model3d
!
!----------------------3d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: r_modext3d, theta_modext3d, phi_modext3d
real(dp), dimension(:,:,:), allocatable :: velr_modext3d, velth_modext3d, velphi_modext3d, &
                                           rho_modext3d, t_modext3d, trad_modext3d, vth_modext3d, &
                                           eps_cont_modext3d
!
end module model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_input
!
!--------------------------input parameters-----------------------------
!
use prog_type
!
implicit none
!
logical :: opt_incl_line
!
integer(i4b) :: input_mod
!
real(dp) :: teff, trad, xlogg, rstar, lstar, rmin, rmax, rlim, tmin, xmloss, vmin, vmax, vmicro, &
            vth_fiducial, beta, yhe, hei, xnue0, eps_cont, kcont, vrot, mstar
!
integer(i4b) :: na
!
end module params_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_stellar
!
!--------------------------stellar parameters---------------------------
!
use prog_type
!
implicit none
!
real(dp) :: sr
!
end module params_stellar
