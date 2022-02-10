!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options_spec
!
implicit none 
!
character(len=500) :: input_file, output_dir
!
integer :: opt_photprof1
integer :: opt_photprof2
!opt_photprof = 0    if no photospheric profile is used (=> diffusion approximation on inner boundary)
!opt_photprof = 1    if photospheric profile is read in (specified by fnamephotprof)
!opt_photprof1, opt_photprof2 refers to each star
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
logical :: opt_incl_gdark1, opt_incl_sdist1, opt_incl_gdark2, opt_incl_sdist2
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
character(len=5) :: opt_pgrid01, opt_pgrid02, opt_rgrid01, opt_rgrid02
!opt_pgrid = 'lin', 'log', 'llog' if p-grid is calculated with linear, log or log-log spacing
!opt_rgrid = 'lin', 'log', 'llog' if r-grid is calculated with linear, log or log-log spacing
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
!------------------------------for each star----------------------------
!
integer(i4b) :: nxobs_fs
!
!for jet model
integer(i4b), parameter :: cs1_np_c=100, cs1_np_nc=20, cs1_np=cs1_np_nc+cs1_np_c-1, cs1_nzeta=81
integer(i4b), parameter :: cs1_nr=20
integer(i4b), parameter :: cs2_np_c=40, cs2_np_nc=160, cs2_np=cs2_np_nc+cs2_np_c-1, cs2_nzeta=81
integer(i4b), parameter :: cs2_nr=400
!
!standard resolution star 1 and star 2
!integer(i4b), parameter :: cs1_np_c=40, cs1_np_nc=160, cs1_np=cs1_np_nc+cs1_np_c-1, cs1_nzeta=81
!integer(i4b), parameter :: cs1_nr=400
!integer(i4b), parameter :: cs2_np_c=40, cs2_np_nc=160, cs2_np=cs2_np_nc+cs2_np_c-1, cs2_nzeta=81
!integer(i4b), parameter :: cs2_nr=400
!
!double resolution star 1 and star 2
!integer(i4b), parameter :: cs1_np_c=40, cs1_np_nc=320, cs1_np=cs1_np_nc+cs1_np_c-1, cs1_nzeta=181
!integer(i4b), parameter :: cs1_nr=800
!integer(i4b), parameter :: cs2_np_c=40, cs2_np_nc=320, cs2_np=cs2_np_nc+cs2_np_c-1, cs2_nzeta=181
!integer(i4b), parameter :: cs2_nr=800
!
!half resolution star 1 and star 2
!integer(i4b), parameter :: cs1_np_c=20, cs1_np_nc=80, cs1_np=cs1_np_nc+cs1_np_c-1, cs1_nzeta=81
!integer(i4b), parameter :: cs1_nr=100
!integer(i4b), parameter :: cs2_np_c=10, cs2_np_nc=40, cs2_np=cs2_np_nc+cs2_np_c-1, cs2_nzeta=81
!integer(i4b), parameter :: cs2_nr=50


integer(i4b) :: nalpha, ngamma, ndirs
!
!nxobs_fs: number of used frequency points
!np_c:   total number of core p-rays
!np_nc:  total number of non-core rays
!np:     total number of p-rays used for flux-calculation
!
!nzeta:  total number of zeta-angles used for flux-calculation
!        (with ni equidistant subintervals between each two zeta-grid points)
!nr:     number of radial grid points
!
!nalpha:   number of directions alpha, that shall be calculated
!ngamma:   number of directions gamma, that shall be calculated
!           alpha and gamma angles specify direction to observer
!ndirs:    total number of directions
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
!star 1
integer(i4b) :: cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc   !spherical coordinates
real(dp), dimension(:), allocatable :: cs1_r_spc, cs1_theta_spc, cs1_phi_spc
real(dp), dimension(:,:,:), allocatable :: cs1_opac3d, cs1_opalbar3d, cs1_t3d, cs1_r3d, &
                                           cs1_velx3d, cs1_vely3d, cs1_velz3d, &
                                           cs1_vth3d, cs1_sline3d, cs1_scont3d
integer, dimension(:,:,:), allocatable :: cs1_imask3d
!
!
!
integer(i4b) :: cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc   !spherical coordinates
real(dp), dimension(:), allocatable :: cs2_r_spc, cs2_theta_spc, cs2_phi_spc
real(dp), dimension(:,:,:), allocatable :: cs2_opac3d, cs2_opalbar3d, cs2_t3d, cs2_r3d, &
                                           cs2_velx3d, cs2_vely3d, cs2_velz3d, &
                                           cs2_vth3d, cs2_sline3d, cs2_scont3d
integer, dimension(:,:,:), allocatable :: cs2_imask3d
!
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
integer(i1b), dimension(:), allocatable :: imask_ray
real(dp), dimension(:), allocatable :: z_ray, opalbar_ray, opac_ray, &
                                       sline_ray, scont_ray, &
                                       profile_ray, velz_ray, vth_ray, &
                                       temp_ray
!$omp threadprivate(nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, velz_ray, &
!$omp               vth_ray, temp_ray, imask_ray)
!
integer(i4b) :: cs1_nz_ray
integer(i1b), dimension(:), allocatable :: cs1_imask_ray
real(dp), dimension(:), allocatable :: cs1_z_ray, cs1_opalbar_ray, cs1_opac_ray, &
                                       cs1_sline_ray, cs1_scont_ray, &
                                       cs1_velz_ray, cs1_vth_ray, &
                                       cs1_temp_ray
!$omp threadprivate(cs1_nz_ray, cs1_z_ray, cs1_opalbar_ray, cs1_opac_ray, cs1_sline_ray, &
!$omp               cs1_scont_ray, cs1_velz_ray, cs1_vth_ray, cs1_temp_ray, cs1_imask_ray)

integer(i4b) :: cs2_nz_ray
integer(i1b), dimension(:), allocatable :: cs2_imask_ray
real(dp), dimension(:), allocatable :: cs2_z_ray, cs2_opalbar_ray, cs2_opac_ray, &
                                       cs2_sline_ray, cs2_scont_ray, &
                                       cs2_velz_ray, cs2_vth_ray, &
                                       cs2_temp_ray
!$omp threadprivate(cs2_nz_ray, cs2_z_ray, cs2_opalbar_ray, cs2_opac_ray, cs2_sline_ray, &
!$omp               cs2_scont_ray, cs2_velz_ray, cs2_vth_ray, cs2_temp_ray, cs2_imask_ray)
!
!nz_ray:      number of data points along an arbitrary ray
!z_ray:       line of sight coordinates of an arbitrary ray
!opalbar_ray: frequency integrated line opacity along ray
!sline_ray:   line-source-function along ray
!profile_ray: profile function along ray
!velz_ray:    line of sight velocity along ray
!vth_ray:     thermal velocity along ray
!temp_ray:    temperature along ray
!
real(dp), dimension(:), allocatable :: gamma_arr
real(dp), dimension(:), allocatable :: alpha_arr
real(dp) :: alpha, gamma
!$omp threadprivate(alpha,gamma)

real(dp) :: unit_length
real(dp), dimension(3) :: nhat, translvec1, translvec2
real(dp), dimension(3,3) :: transmat, transmat_inv, rotmat1, rotmat2, rotmat1_inv, rotmat2_inv
!$omp threadprivate(nhat, transmat, transmat_inv, &
!$omp               translvec1, translvec2)

!unit_length:  length scale of global coordinate system
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
!the maximum radius of global coordinate system such that each system is included
!real(dp) :: rmax
!t
!
real(dp), parameter :: del_xobs=1.d0/3.d0
!del_xobs:  maximum allowed frequency steps 
!           (shifted frequency from line center in fiducial doppler widths)
!



!binary version     
real(dp), dimension(:), allocatable :: xobs, xnue, &
                                       xic1_nue, xicc1_nue, xic2_nue, xicc2_nue, &
                                       acoeff_xic1, bcoeff_xic1, ccoeff_xic1, dcoeff_xic1, &     
                                       acoeff_xic2, bcoeff_xic2, ccoeff_xic2, dcoeff_xic2
real(dp), dimension(:), allocatable :: flux_tot, flux_cont, flux_emi, flux_abs, normt
!$omp threadprivate(flux_tot, flux_cont, flux_emi, flux_abs, normt)




real(dp), dimension(:), allocatable :: cs1_p, cs2_p, cs1_r, cs2_r
real(dp), dimension(:), allocatable :: cs1_zeta, cs2_zeta

!xobs:      grid of frequency points
!flux_tot:  total emergent flux as function of frequency
!flux_cont: continuum emergent flux as function of frequency
!flux_emi:  only emission part of profile
!flux_abs:  only absorption part of profile
!normt:     test normalization of p and zeta integration
!

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
logical :: lcore1, lcore2
!$omp threadprivate(lcore1, lcore2)
!lcore: logical to describe core and non-core rays
!
real(dp) :: vphot_proj, vphot1_proj, vphot2_proj, phinorm
!$omp threadprivate(vphot_proj, vphot1_proj, vphot2_proj, phinorm)
!vphot_proj: photospheric velocity projected onto ray direction (only rotational velocity+velocity in global coordinate system)
!phinorm:    normalized profile
!
!
end module mod_spectrum
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_triangles
!
use prog_type 
!
implicit none
!
integer(i4b) :: ntriangles, npoints
!$omp threadprivate(npoints, ntriangles)
!
!indices for each triangle and the index of adjacent triangles, index of vertices
integer(i4b), dimension(:,:), allocatable :: triangles_ip, triangles_if
!$omp threadprivate(triangles_ip, triangles_if)
!
!all points
integer(i4b), dimension(:), allocatable :: points_indx
real(dp), dimension(:), allocatable :: points_xcoord, points_ycoord, points_weight
real(dp), dimension(:,:), allocatable :: points_coords
!$omp threadprivate(points_indx, points_xcoord, points_ycoord, points_weight, points_coords)


end module mod_triangles
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module timing_spec
!
use prog_type
!
real(dp) :: t_tot
!
real(dp) :: ttot_obsdir, t_setup1, t_setup2, t_setup3, t_triangles
!$omp threadprivate(ttot_obsdir, t_setup1, t_setup2, t_setup3, t_triangles)
!
!t_tot:  total computation time
!
!ttot_obsdir: calculation time of each direction
!
!t_setup1: total compuation time (for setting up ray of star1)
!t_setup2: total compuation time (for setting up ray of star2)
!t_setup3: total compuation time (for setting up ray in global coordinates)
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
real(dp) :: xobs_surface, alpha_surface, gamma_surface
!
real(dp), dimension(:), allocatable :: iem_surface, iemi_surface, iabs_surface
!
!iem_surface: emergent intensity on a p-zeta-surface
!             calculated at xobs_surface, alpha_surface, gamma_surface
!iemi_surface: only emission part
!iabs_surface: only absorption part
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
integer(i4b), parameter :: ntheta1_gdark=51, ntheta2_gdark=51
real(dp), dimension(ntheta1_gdark) :: theta1_gdark, xic1_gdark, teff1_gdark
real(dp), dimension(ntheta2_gdark) :: theta2_gdark, xic2_gdark, teff2_gdark
real(dp) :: smajorax1_a, smajorax1_b, smajorax1_c, &
            smajorax2_a, smajorax2_b, smajorax2_c
!
real(dp) :: xic1_factor, xic2_factor
!$omp threadprivate(xic1_factor, xic2_factor)
!
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
real(dp) :: vth_fiducial_model, vmax_model
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
real(dp) :: rmax0
!$omp threadprivate(rmax0)
real(dp) :: x01, y01, z01, vx01, vy01, vz01, rstar1, rmin1, rmax1, sr1, teff1, trad1, xic1, vrot1, logg1, lstar1, yhe1, fehe1, aenh1, vmicro1
real(dp) :: x02, y02, z02, vx02, vy02, vz02, rstar2, rmin2, rmax2, sr2, teff2, trad2, xic2, vrot2, logg2, lstar2, yhe2, fehe2, aenh2, vmicro2
!
real(dp), dimension(3) :: ex01, ey01, ez01, rot_axis01
real(dp), dimension(3) :: ex02, ey02, ez02, rot_axis02
!
real(dp) :: kline, eps_line
!
real(dp) :: vth_fiducial, vmax, vth_min
!
real(dp) :: tmin
!
real(dp) :: xnue0
!transition frequency of considered line-transition
!
integer(i4b) :: na
!
integer(i4b) :: iline
!
!
real(dp), parameter :: xlim=20.d0!135.d0
!xlim is the maximum allowed shift of (nue_obs-nue_0)/vmax
!note: sometimes up to 135 (fastwind photospheric profile...)
!
!
end module params_spec
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
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
