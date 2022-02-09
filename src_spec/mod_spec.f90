!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options_spec
!
implicit none
!
character(len=100) :: input_file, output_dir
!
integer :: opt_photprof
!!opt_photprof = 0    if no photospheric profile is used (=> diffusion approximation on inner boundary)
!!opt_photprof = 1    if photospheric profile is read in (specified by fnamephotprof)
!
logical :: opt_obsdir_read
!opt_obsdir_read=t   if angles alpha and gamma shall be read in
!opt_obsdir_read=f   if angles alpha and gamma are calculated equidistantly
!
logical :: opt_surface
!opt_surface=t   if (only) surface brightness shall be calculated for given xobs, and direction
!opt_surface=f   if standard formal solver is applied
!
logical :: opt_int2d
!opt_surface=t   if (only) intensity and optical depth are calculated for a given xobs, and direction
!opt_surface=f   if standard formal solver is applied
!
logical :: opt_incl_gdark, opt_incl_sdist
!opt_incl_gdark=t   if gravity darkening shall be included for rotating stars
!opt_incl_sdist=t   if surface distortion shall be included for rotating stars
!
integer, parameter :: interp_photprof=0
!interp_photprof=0   if linear interpolation is performed for photospheric profile
!interp_photprof=1   if monotonic cubic spline interpolation is performed for photospheric profile
!
!
integer :: input_mod=0
!input_mod = 0 if 1d-input is used
!input_mod = 1 if 3d-input is used
!input_mod = 2 if 3d-input is used (in spherical coordinates)
!
end module options_spec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_spec
!
use prog_type
!
implicit none
!
!--------------------dimension of p, zeta and radial grid---------------
!
integer(i4b) :: nxobs_fs
!standard resolution
integer(i4b), parameter :: np1_nc=40, np1_c=10, ni=4, np_c=np1_c*ni, np_nc=np1_nc*ni, np=np_nc+np_c-(ni-1), nzeta1=20, nzeta=ni*nzeta1 - (ni-1)
integer(i4b), parameter :: nr=400
!double resolution
!integer(i4b), parameter :: np1_nc=80, np1_c=20, ni=4, np_c=np1_c*ni, np_nc=np1_nc*ni, np=np_nc+np_c-(ni-1), nzeta1=40, nzeta=ni*nzeta1 - (ni-1)
!integer(i4b), parameter :: nr=800
!quad resolution
!integer(i4b), parameter :: np1_nc=160, np1_c=40, ni=4, np_c=np1_c*ni, np_nc=np1_nc*ni, np=np_nc+np_c-(ni-1), nzeta1=80, nzeta=ni*nzeta1 - (ni-1)
!integer(i4b), parameter :: nr=1600
integer(i4b) :: nalpha, ngamma
!
!nxobs_fs: number of used frequency points
!np1_nc: number of non-core rays (are allowed to have unequal steps)
!np1_c:  number of     core rays (are allowed to have unequal steps)
!ni:     number of equidistant intervals added between each (unequal) p-grid-points
!np_c:   total number of core p-rays
!np_nc:  total number of non-core rays
!np:     total number of p-rays used for flux-calculation
!
!nzeta1: number of zeta angles (allowed to have unequal steps)
!nzeta:  total number of zeta-angles used for flux-calculation
!        (with ni equidistant subintervals between each two zeta-grid points)
!nr:     number of radial grid points
!
!nalpha:   number of directions alpha, that shall be calculated
!ngamma:   number of directions gamma, that shall be calculated
!           alpha and gamma angles specify direction to observer
!
end module dime_spec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_model1d
!
!--------------------------1d-model atmosphere--------------------------
!
use prog_type
!
implicit none
!
!
integer(i4b) :: n1d
real(dp), dimension(:), allocatable :: sline1d, scont1d, r1d, velr1d, opalbar1d, opac1d, t1d, vth1d
!
!for cubic spline interpolation
real(dp), dimension(:), allocatable :: acoeff_sline1d, bcoeff_sline1d, &
                                       ccoeff_sline1d, dcoeff_sline1d, &
                                       acoeff_scont1d, bcoeff_scont1d, &
                                       ccoeff_scont1d, dcoeff_scont1d, &
                                       acoeff_velr1d, bcoeff_velr1d, &
                                       ccoeff_velr1d, dcoeff_velr1d, &
                                       acoeff_opalbar1d, bcoeff_opalbar1d, &
                                       dcoeff_opalbar1d, ccoeff_opalbar1d, &
                                       acoeff_opac1d, bcoeff_opac1d, &
                                       dcoeff_opac1d, ccoeff_opac1d, &
                                       acoeff_t1d, bcoeff_t1d, &
                                       dcoeff_t1d, ccoeff_t1d
!
!
!
end module dime_model1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_model3d
!
!--------------------------3d-model atmosphere--------------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: ndxmax, ndymax, ndzmax         !cartesian coordinates
integer(i4b) :: nr_spc, ntheta_spc, nphi_spc   !spherical coordinates
!
real(dp), dimension(:), allocatable :: x, y, z
real(dp), dimension(:), allocatable :: r_spc, theta_spc, phi_spc
real(dp), dimension(:,:,:), allocatable :: opac3d, opalbar3d, t3d, r3d, &
                                           velx3d, vely3d, velz3d, &
                                           vth3d, sline3d, scont3d
integer, dimension(:,:,:), allocatable :: imask3d
!
!
end module dime_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_spectrum
!
use prog_type
use fund_const
!
implicit none
!
integer(i4b) :: nz_ray
real(dp), dimension(:), allocatable :: z_ray, opalbar_ray, opac_ray, &
                                       sline_ray, scont_ray, &
                                       profile_ray, velz_ray, vth_ray, &
                                       temp_ray
!$omp threadprivate(nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, velz_ray, &
!$omp               vth_ray, temp_ray)
!nz_ray:      number of data points along an arbitrary ray
!z_ray:       line of sight coordinates of an arbitrary ray
!opalbar_ray: frequency integrated line opacity along ray
!sline_ray:   line-source-function along ray
!profile_ray: profile function along ray
!velz_ray:    line of sight velocity along ray
!vth_ray:     thermal velocity along ray
!temp_ray:    temperature along ray
!
real(dp), dimension(:), allocatable :: p, pw, pw1, pw_err, r
real(dp), dimension(:), allocatable :: zeta, zetaw, zetaw1, zetaw_err
!p:    p-ray grid
!r:    radial grid
!zeta:  grid of angles used for flux integration
!zetaw: corresponding integration weights
!
real(dp), dimension(:), allocatable :: gamma_arr
real(dp), dimension(:), allocatable :: alpha_arr
real(dp) :: alpha, gamma

real(dp), dimension(3) :: nhat
real(dp), dimension(3,3) :: transmat, rotmat, rotmat_trans
!alpha(_arr):     angle of observer wrt z-axis (either rotation or magnetic axis)
!gamma(_arr):     angle of observer wrt x-z-plane
!nhat:      direction vector to observer in x,y,z-coordinate system
!transmat:  transformation matrix from cylindrical coordinates to
!              x,y,z coordinate system
!rotmat: transformation matrix from carthesian coordinates in formal solver
!             to carthesian coordinates in 3d-grid
!rotmat_trans: transformation matrix from carhtesian coordinates in 3d-grid
!                 to carthesian coordinates in formal solver 
!                 (get velocity components correct!)
!
real(dp) :: rmax, rmin
!
real(dp), parameter :: del_xobs=1.d0/3.d0
!del_xobs:  maximum allowed frequency steps 
!           (shifted frequency from line center in fiducial doppler widths)
!
real(dp), dimension(:), allocatable :: xobs, xnue, xic_nue, xicc_nue, flux_tot, flux_cont, &
                                       flux_emi, flux_abs, femi_p, fabs_p, &
                                       ftot_p, fcont_p, ftot_errp, fcont_errp, &
                                       ftot_err, fcont_err, normt_p, normt, &
                                       acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic, &
                                       flux_tot_tmp, flux_cont_tmp, normt_tmp, ftot_err_tmp, &
                                       fcont_err_tmp, flux_emi_tmp, flux_abs_tmp
!
!$omp threadprivate(flux_tot_tmp, flux_cont_tmp, normt_tmp, ftot_err_tmp, fcont_err_tmp, &
!$omp               ftot_p, fcont_p, ftot_errp, fcont_errp, normt_p, flux_emi_tmp, flux_abs_tmp, &
!$omp               femi_p, fabs_p)
!
!xobs:      grid of frequency points
!flux_tot:  total emergent flux as function of frequency
!flux_cont: continuum emergent flux as function of frequency
!flux_emi:  only emission part of profile
!flux_abs:  only absorption part of profile
!normt:     test normalization of p and zeta integration
!ftot_errp:  error of total flux in p-integration (for each xobs and each phi)
!fcont_errp: error of continuum flux in p-integration (for each xobs and each phi)
!ftot_err:   error of total flux in zeta-integration (for each xobs)
!fcont_err:  error of continuum flux in zeta-integration (for each xobs)
!*_tmp:      temporary arrays for integration in parallelized code
!
real(dp) :: relerr_contp, relerr_totp, relerr_cont, relerr_tot, hmax
!$omp threadprivate(relerr_contp, relerr_totp)
!
!relerr_contp: maximum error of continuum flux in p-integration
!relerr_totp:  maximum error of total flux in p-integration
!relerr_cont:  maximumn error of continuum flux in zeta-integration
!relerr_tot:   maximum error of total flux in p-integration
!hmax:         maximum step size of p-grid
!
real(dp), parameter :: del_vel=1.d0/3.d0
real(dp) :: del_vel2
!del_vel: maximum allowed velocity steps along ray in thermal velocities (in order to resolve resonance zones)
!del_vel2: maximum allowed velocity steps along ray in fiducial thermal velocities
!
real(dp) :: iin, iin_c, iem, iem_c, iemi, iabs
!$omp threadprivate(iin, iin_c, iem, iem_c, iemi, iabs)
!iin:   core intensity (read in planck-function and set iin to that value
!                       or use photospheric profile)
!iin_c: core intensity in the absence of line
!iem:   emergent intensity at outer boundaries (from line and continuum)
!iem_c: emergent intensity at outer boundaries (from continuum alone)
!iemi:  only emission part of profile
!iabs:  only absorption part of profile
!
logical :: lcore
!$omp threadprivate(lcore)
!lcore: logical to describe core and non-core rays
!
real(dp) :: vphot_proj, phinorm
!$omp threadprivate(vphot_proj, phinorm)
!vphot_proj: photospheric velocity projected onto ray direction (only rotational velocity)
!phinorm:    normalized profile
!
!
end module mod_spectrum
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module timing_spec
!
use prog_type
!
real(dp) :: ts_tot, te_tot, t_tot, ts_obsdir, te_obsdir, ttot_obsdir
real(dp) :: ts_setup, te_setup, t_setup_tmp, t_setup
real(dp) :: ts_hunt, te_hunt, t_hunt_tmp, t_hunt
real(dp) :: ts_trilin, te_trilin, t_trilin_tmp, t_trilin
!
!$omp threadprivate(ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, &
!$omp               t_hunt_tmp, ts_trilin, te_trilin, t_trilin_tmp)
!
!for system time (not cpu-time)
real(dp) :: ttot_obsdir_sys
integer(i4b) :: ticks_obsdir, ticks_initial_obsdir, ticks_final_obsdir
!
!ts_tot: start of program
!te_tot: end of program
!t_tot:  total computation time
!
!ts_obsdir:  start calculation for each direction
!te_obsdir:  end of calculation for each direction
!ttot_obsdir: end of calculation for each direction
!
!ts_setup: start of subroutine setup_ray
!te_setup: end of subroutine setup_ray
!t_setup: total compuation time (for setting up ray)
!
!ts_hunt: start of searching routine to find indices
!te_hunt: end of searching routine to find indices
!t_hunt:  total computation time (for finding indices)
!ts_trilin: start of trilinear interpolation
!te_trilin: end of trilinear interpolation
!t_trilin:  total computation time for all trilinear interpolations
!

!
end module timing_spec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_surfb
!
use prog_type
!
implicit none
!
integer(i4b) :: nsurfb
real(dp), dimension(:), allocatable :: xobss_surface, alphas_surface, gammas_surface
!
real(dp) :: xobs_surface, alpha_surface, gamma_surface, xic1_surface
real(dp), dimension(:,:), allocatable :: iem_surface, iemi_surface, iabs_surface, icont_surface
!
!iem_surface: emergent intensity on a p-zeta-surface
!             calculated at xobs_surface, alpha_surface, gamma_surface
!iemi_surface: only emission part
!iabs_surface: only absorption part
!icont_surface: only continuum part
!
!
end module mod_surfb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_int2d
!
use prog_type
!
implicit none
!
integer(i4b) :: nz_ray_max
!
real(dp), dimension(:,:), allocatable :: int_2d, tau_2d, zcoord_2d, xcoord_2d, iemi_2d, iabs_2d, vn_2d
!
!zcoord_2d: z-coordinate along ray for a given impact-parameter p
!xcoord_2d: x-coordinate perpendicular to ray (basically = impact-parameter)
!int_2d:    intensity on given coordinates
!tau_2d:    optical depths on given coordinates
!vn_2d:     velocity along ray in 2d
!
end module mod_int2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_gdark
!
!module for gravity darkening
!
use prog_type
!
implicit none
!
integer(i4b), parameter :: ntheta_gdark=51
real(dp), dimension(ntheta_gdark) :: theta_gdark, xic1_gdark, teff_gdark
real(dp) :: smajorax_a, smajorax_b, smajorax_c
!
real(dp) :: xic1_factor
!$omp threadprivate(xic1_factor)
end module mod_gdark
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_model
!
use prog_type
!
implicit none
!
real(dp) :: vth_fiducial_model, vmicro_model, vmax_model
!
!need vmax_model because may be larger than input v_inf, however, need large
!   enough frequency range
!
end module params_model
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_spec
!
use prog_type
!
implicit none
!
!input file
!
!input parameter
real(dp) :: xlogg, rstar, lstar, rmax, xmloss, &
            yhe, hei, sr, mstar

real(dp) :: teff, tmin, t_inf, xic1, trad

real(dp) :: vmin, vmax, vrot, vth_fiducial, vmicro, beta, b, vth_min, v_esc, v_inf

real(dp) :: eps_line, kline, alpha, kappa0
!
!real(dp) :: teff, trad, xic1
!real(dp) :: vth_fiducial, vmicro, vrot
real(dp) :: xnue0
!transition frequency of considered line-transition
!
integer(i4b) :: na
!atomic mass number of considered atom (needed for correct calculation of
!   thermal velocity)
!
real(dp), parameter :: vth_lim=5.d5
!vth_lim is minimum adopted thermal velocity if vmicro is small and temperature is small (for computational reasons)
!
real(dp), parameter :: xlim=3.d0
!xlim is the maximum allowed shift of (nue_obs-nue_0)/vmax
!
!for adm-model
real(dp) :: ralfven, delta, chi_inf, obliquity, rhoc_star, rhow_star, rmax_sfct
!
!
!
end module params_spec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_mod3d
!
use prog_type
!
implicit none
!
!input parameter (from 3d-model)
!
real(dp) :: sr_mod3d, xlogg_mod3d, rstar_mod3d
!
real(dp) :: teff_mod3d, tmin_mod3d, xmloss_mod3d, xic1_mod3d
!
real(dp) :: vmin_mod3d, vmax_mod3d, vthfiducial_mod3d, beta_mod3d, vmicro_mod3d
!
real(dp) :: eps_line_mod3d
!
real(dp) :: xnue0_mod3d
!
real(dp) :: obliquity_mod3d, ralfven_mod3d, delta_mod3d, chiinf_mod3d
!
real(dp) :: theta_d_mod3d, dtheta_abl_mod3d, beta_accr_mod3d, tau_d_mod3d
!
real(dp) :: yhe_mod3d, hei_mod3d
!
integer(i4b) :: na_mod3d
!
integer(i4b) :: inputmod_mod3d

end module params_mod3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_spectests
!
use prog_type
!
integer(i4b), parameter :: opt_test_fs=2
!
!opt_test_fs = 0   for sline=0.d0, profile=1.d0, opac=opac_const (constant continuum absorption)
!opt_test_fs = 1   for sline=0.d0, profile=1.d0, opac = opac/(z+1)^2 (decaying opacity, only absorption)
!opt_test_fs = 2   for sline=0.d0, profile=1.d0, opac = opac/(z+1)^3 (decaying opacity, only absorption)
!opt_test_fs = 3   for sline=const, profile=1.d0, opac=opac_const (constant continnum absorption + const. source terms)
!
!--------------testing continuum transport along ray--------------------
!
real(dp), parameter :: zmin=0.d0, zmax=10.d0
real(dp), parameter :: opac_const=1.d0, scont_const=1.d0
real(dp), parameter :: opac_max=10.d0
real(dp), parameter :: iin=1.d0
!
!-----------------------testing line flux-------------------------------
!
real(dp) :: kappa0, chi_alpha
!
end module mod_spectests
