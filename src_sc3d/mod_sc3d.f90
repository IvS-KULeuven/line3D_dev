!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_directories
!
!---------------------module for directories----------------------------
!
implicit none
!
character(len=100) :: input_file
character(len=100) :: output_file
character(len=300) :: model_dir
character(len=11), parameter :: output_dir='outputFILES'
character(len=16), parameter :: output_dir_temp='outputFILES_TEMP'
character(len=16), parameter :: output_dir_test='outputFILES_TEST'
!input_file:      namelist-file with all input-data
!output_file:     name of the final output file (has to end with .h5)
!input_dir:       directory of input models
!output_dir:      directory of final output
!output_dir_temp: directory where solution after each iteration step is stored
!all strings are only the directory-name (without '/')
!
character(len=10), parameter :: model1d_file='model1d.h5'
character(len=10), parameter :: model2d_file='model2d.h5'
character(len=10), parameter :: model3d_file='model3d.h5'
!modelxd_file:   file where model-atmosphere is stored
!
character(len=16), parameter :: model1d_phot_file='model1d_phot.dat'
!model1d_phot_file: file where photospheric model is stored
!
end module mod_directories
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options
!
!-----------------------------------------------------------------------
!------------------all options (e.g. grid-options, etc)-----------------
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: input_mod, input_mod_dim
!input_mod: which input-model is used (see *.nml file for more information)
!           however, not needed, because only external file will be read in
!input_mod_dim: which external file shall be read: 1d-model, 2d-model or 3d-model
!
integer(i4b) :: spatial_grid1d
!spatial_grid1d=0 if equidistant radial grid is used (subroutine grid1d_r_equi)
!spatial_grid1d=1 if equidistant velocity grid is used (subroutine grid1d_vel_equi)
!spatial_grid1d=2 if equidistant tau_thomson grid is used (subroutine grid1d_tau_equi)
!spatial_grid1d=3 if equidistant log(tau_thomson) grid is used (subroutine grid1d_tau_log)
!spatial_grid1d=4 if combination is used (see subroutine grid1d_final for details)
!spatial_grid1d=5 if combination is used (see subroutine grid1d_final_2 for details)
!spatial_grid1d=6 if grid is calucalted equidistant in log-space (subroutine grid1d_r_log)

integer(i4b) :: spatial_grid3d
!spatial_grid3d=0 if 3d grid is calculated from 1d grid with equidistant core
!                 points 
!spatial_grid3d=1 if 3d grid is calculated from a mean-value approach
!                 (minimizing distance of subsequent coordinates from 1d-grid)
!spatial_grid3d=2 if 3d grid is calculated from a mean-value approach 
!                 (minimizing distance of subsequent coordinates from original
!                 input-grid)
!spatial_grid3d=3 if 3d grid is calculated completely equidistant

integer(i4b) :: opt_opac
!opt_opac = 0 if continuum opacity is parameterized by kcont*chi_thomson
!
integer(i4b) :: opt_opal
!opt_opal = 0 if line is parameterized as in hennicker et al 2017
!opt_opal = 1 if line is parameterized as in hamann 1981
!
integer(i4b) :: opt_angint_method
!opt_angint_method=0 if angular integration is used with trapezoidal rule (nodes equidistant in theta and phi)
!opt_angint_method=1 if angular integration is used with trapezoidal rule (nodes from Lobell & Blomme 2008)
!opt_angint_method=2 if angular integration is used with simpsons rule (nodes equidistant in theta and phi)
!                       (note: mu-grid and phi-grid will be made equidistant for three subsequent points)
!opt_angint_method=3 if angular integration is used with simpson rule corrected for the error
!                       from a grid with half resolution (also known as boole's rule)
!opt_angint_method=4 if angular integration is used with cubic splines (catmull-rom-spline, nodes equidistant in theta and phi)
!opt_angint_method=5 if angular integration is used with gauss-legendre-integration (for each octant)
!opt_angint_method=6 if angular integration is used with gauss-chebyshev-integration (for each octant)
!opt_angint_method=7 if angular integration is used with triangulation (linear integrals)
!opt_angint_method=8 if angular integration is used with triangulation ('pseudo'-gauss integrals per triangle)
!opt_angint_method=9 if angular integration is used with lebedev interpolation (optimized nodes on the sphere)
!
integer(i4b) :: opt_method
!opt_method=0  if finite volume method shall be used
!opt_method=1  if linear short characteristics method shall be used
!opt_method=2  if quadratic bezier short characteristics method shall be used
!
logical :: opt_sol2d
!opt_sol2d = true if pseudo-2d solution shall be performed
!
logical :: opt_incl_cont
!opt_incl_cont = true if continuum shall be included in calculations
!
logical :: opt_start_cont
!opt_start_cont = true if continuum is calculated from beginning on
!opt_start_cont = false if current iterate is read in and calculation continues from that
!
logical :: opt_incl_line
!opt_incl_line = true if line shall be included in calculations
!
logical :: opt_start_line
!opt_start_line = true if line is calculated from beginning on
!opt_start_line = false if current iterate is read in and calculation continues from that
!
logical :: opt_ng_cont, opt_ait_cont
!opt_ng_cont = true if ng-extrapolation shall be performed for continuum transfer
!opt_ait_cont = true if aitken extrapolation shall be performed for continuum transfer
!
logical :: opt_ng_line, opt_ait_line
!opt_ng_line = true if ng-extrapolation shall be performed for line transfer
!opt_ait_line = true if aitken extrapolation shall be performed for line transfer
!
integer(i4b) :: opt_alo_cont   !(affects only 3d solution scheme)
!opt_alo_cont = 0 if classical lambda iteration shall be used for continuum transfer
!opt_alo_cont = 1 if diagonal alo is used for ALI scheme for continuum transfer
!opt_alo_cont = 2 if direct neighbour alo is used for ALI scheme for continuum transfer (6 neighbours)
!opt_alo_cont = 3 if nearest neighbour alo is used for ALI scheme for continuum transfer (26 neighbours)
!
integer(i4b) :: opt_alo_line   !(affects only 3d solution scheme)
!opt_alo_line = 0 if classical lambda iteration shall be used for line transfer
!opt_alo_line = 1 if diagonal alo is used for ALI scheme for line transfer
!opt_alo_line = 2 if nearest neighbour alo is used for ALI scheme for line transfer
!
logical :: opt_incl_gdark, opt_incl_sdist
!opt_incl_gdark = true if gravity darkening shall be included
!opt_incl_sdist = true if surface distortion of rotating star shall be included
!
integer(i4b) :: opt_ltec
!opt_ltec = 0 if only one frequency point is considered
!opt_ltec = 1 if grey temperature stratification is calculated
!
end module options
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module warnings
!
logical :: warn_mu, warn_phi
!warn_mu:  warning if mu-grid not symmetric
!warn_phi: warning if phi-grid not symmetric
!
logical :: warn_angles, warn_angles2
!warn_angles:  warning if n_x, n_y, n_z not symmetric
!warn_angles2: warning if required resolution is not met
!
logical :: warn_itmaxc, warn_itmaxl
!
end module warnings
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module modext
!
!-----------------------------------------------------------------------
!----------------------external model atmosphere------------------------
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: nr_modext, ntheta_modext, nphi_modext, nr_modphot
!
!for photospheric model
real(dp), dimension(:), allocatable :: r_modphot, t_modphot1d, rho_modphot1d
!
!for wind model
real(dp), dimension(:), allocatable :: r_modext, theta_modext, phi_modext
real(dp), dimension(:), allocatable :: rho_modext1d, t_modext1d, velr_modext1d, vth_modext1d, &
                                       eps_cont_modext1d
real(dp), dimension(:,:), allocatable :: rho_modext2d, t_modext2d, velr_modext2d, &
                                         velth_modext2d, velphi_modext2d, vth_modext2d, &
                                         eps_cont_modext2d
real(dp), dimension(:,:,:), allocatable :: rho_modext3d, t_modext3d, velr_modext3d, &
                                           velth_modext3d, velphi_modext3d, vth_modext3d, &
                                           eps_cont_modext3d

!
end module modext
!
!-----------------------------------------------------------------------
!-----------------------diffusion approximation-------------------------
!-----------------------------------------------------------------------
!
module bcondition
!
use prog_type
!
implicit none
!
real(dp) :: xic1, xic2, corrfc, dbdtau, opathboundary
!
!xic1: planck-function b
!xic2: db/dz
!
!number of points at inner boundary, and cooresponding arrays
!describing the x, y, z position of the points, and the intensity
integer(i4b) :: n_inner
integer, dimension(:), allocatable :: indx_xinner, indx_yinner, indx_zinner
real(dp), dimension(:), allocatable :: int_inner
!
!for gravity darkening
integer(i4b), parameter :: ntheta_gdark=51
real(dp), dimension(ntheta_gdark) :: theta_gdark, xic1_gdark, teff_gdark
!
end module bcondition
!
!-----------------------------------------------------------------------
!-------------------------------1d-grids--------------------------------
!-----------------------------------------------------------------------
!
module dime1d
!
use prog_type
!
integer(i4b) :: n1d     != 60!69 !60 !14 !36
integer(i4b) :: n1d_t   !=81 !61 !101
integer(i4b) :: n1d_r   !=22 !42 !2
integer(i4b) :: n1d_dum
!
real(dp) :: delv
!
real(dp), dimension(:), allocatable :: r1d, r1d_dum
!
!general grid construction:
!1. create dummy radial grid with n1d_dum grid points
!       n1d_dum=(vmax-vmin)/delv  if grid shall be equidistant in velocity space
!       n1d_dum=n1d_t+n1d_r       else   
!2. create actually used radial grid with n1d grid points from dummy grid
!
!if dummy radial grid is calculated equidistant in velocity-space
!delv: steps of velocity grid in units of vth_fiducial
!
!if dummy radial grid is calculated equidistant in tau:
!n1d_t: number of grid points equidistant in tau
!n1d_r: number of grid points added to have good radial resolution
!n1d_dum: total number of grid points
!
!
end module dime1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dimecr
!
!for grids on central ray (along z-axis)
use prog_type
!
integer(i4b) :: n1d_cr
real(dp), dimension(:), allocatable :: r1d_cr, int1d_cr, mint1d_cr, mintbar1d_cr, &
                                       scont1d_cr, sline1d_cr, norm1d_cr
real(dp), dimension(:), allocatable :: alocont1d_diag_cr, aloline1d_diag_cr

!relative deviation after each iteration step for continuum
real(dp), dimension(:), allocatable :: epsc1d_cr, epsl1d_cr
real(dp), dimension(:,:), allocatable :: epshistoryc1d_cr, epshistoryl1d_cr
!
end module dimecr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module angles
!
! angular integration parameters and arrays
!
use prog_type
use fund_const
!
implicit none
!
integer(i4b) :: n_theta
integer(i4b) :: dim_mu, dim_phi, dim_omega
!
real(dp), dimension(:), allocatable :: nodes_mu, weight_mu
real(dp), dimension(:), allocatable :: n_x, n_y, n_z, weight_omega
!
integer, dimension(:,:), allocatable :: q_alo
!
!
!dim_mu: total number of mu-integration-nodes
!nodes_mu: mu-integration-nodes
!weight_mu: mu-integration-weight
!
!dim_phi: total number of phi-integration-nodes
!weight_phi: phi-integration-weights
!weight_omega: total solid-angle integration weight for each direction
!n_x, n_y, n_z: direction vector components (n_x,n_y,n_z)
!
!n_theta: number of integration nodes for theta=[0,pi/2]
!
!q_alo: indices for alo-calculations for each direction
!       to store nearest neighbours correctly
!
end module angles
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module freq
!
!   frequencies
!
use prog_type
!
implicit none
!
! ... scalars
integer(i4b) :: nxobs
real(dp) :: deltax, xcmf_max, xcmf_min, xmax
real(dp) :: xnue0
!
! ... logicals
logical :: lfreqint
!
! ... arrays
real(dp), dimension(:), allocatable :: nodes_xobs, nodes_xnue, weight_xobs, xic_nue
!
!lfreqint: if frequency integrated values are used
!xnue0: transition frequency of the line
!nxobs: number of frequency integration nodes
!nodes_xobs: frequency integration nodes in observers frame
!            (shift from line center in fiducial doppler units)
!nodes_xnue: corresponding frequency points
!deltax: allowed step-size of x-obs-grid
!xmax: range of frequency grid
!xcmf_min, xcmf_max: range of cmf frequencies accounted for profile calculations
!xic_nue: photospheric profile
!
end module freq
!
!-----------------------------------------------------------------------
!----------------------------3d - grids---------------------------------
!-----------------------------------------------------------------------
!
module dime3d
!
use prog_type
use omp_lib
!
implicit none
!
integer(i4b) :: ncx, ncy, ncz, ndx, ndy, ndz, ndxmax, ndymax, ndzmax
real(dp) :: delx_max, dely_max, delz_max
real(dp), dimension(:), allocatable :: x, y, z
!
!x,y,z:    x,y,z-coordinates, dimension: ndxmax, ndymax, ndzmax
!note: ndx, ndy, ndz have additional point (+1) because phantom-point is included
!delx_max, dely_max, delz_max: maximum allowed increment in x,y,z-coordinates

real(dp), dimension(:,:,:), allocatable :: rho3d, opac3d, opalbar3d, t3d, velx3d, vely3d, velz3d, &
                                           vth3d, int3d, scont3d, sline3d, ssobo3d, eps_cont3d, &
                                           normalization3d
real(dp), dimension(:,:,:), allocatable :: mint3d, mintbar3d, fcontx3d, fconty3d, fcontz3d, kcontxx3d, kcontyy3d, kcontzz3d, kcontxy3d, kcontxz3d, kcontyz3d
integer(i1b), dimension(:,:,:), allocatable :: imask3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d
!logical, dimension(:,:,:), allocatable :: mask_totreg3d , mask_innreg3d, mask_bpoint3d
!imask_totreg3d:  is 1 if within computational domain, 0 else
!imask_innreg3d:  is 1 if inside central star, 0 else
!imask_bpoint3d:  is 1 if boundary point, 0 else

!for finite volume method
real(dp), dimension(:,:,:,:), allocatable :: delx_arr, dely_arr, delz_arr
!
!for all alo terms
real(dp), dimension(:,:,:,:), allocatable :: alocont_nn3d, aloline_nn3d, alocont_o_nn3d, aloline_on_nn3d
real(dp), dimension(:), allocatable :: alocont_data, alocont_data_diag, aloline_data, aloline_data_diag
integer(i4b), dimension(:), allocatable :: alocont_colindx, alocont_rowindx, aloline_colindx, aloline_rowindx
!
!for paralleliztion
real(dp), dimension(:,:,:,:), allocatable :: aloline_nn3d_tmp, alocont_nn3d_tmp
real(dp), dimension(:,:,:), allocatable :: mint3d_tmp, mintbar3d_tmp, normalization3d_tmp, fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
                                           kcontxx3d_tmp, kcontyy3d_tmp, kcontzz3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, kcontyz3d_tmp     
!$omp threadprivate(int3d, aloline_on_nn3d, aloline_nn3d_tmp, mintbar3d_tmp, normalization3d_tmp)
!$omp threadprivate(alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp)
!$omp threadprivate(kcontxx3d_tmp, kcontyy3d_tmp, kcontzz3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, kcontyz3d_tmp)
!
!
!
!alocont_o_nn3d: nearest neighbout alo dependent on omega
!aloline_on_nn3d: nearest neighbour alo dependent on omega and nue
!
!int3d: intensity on 3d-grid
!scont3d, sline3d: source-functions on 3d-grid
!
!mask_totreg3d: mask to specify where intensity shall be calculated
!                  (skip grid points with r larger rlim and inside core)
!mask_innreg3d: mast to specify where intensity shall be calculated
!                  (only inner region: skip grid-points inside star)
!bmask3d: true where grid-points lie exactly on inner boundary
!
!alocont_nn3d: nearest neighbour alo for continuum (per angle)
!alocont_nn3d_aint: angle-integrated nearest neighbour alo for continuum
!aloline_on_nn3d: nearest neighbour alo for line (per solid angle omega and frequency nu)
!aloline_nn3d: nearest neighbour alo for line (integrated over solid angle and frequency)
!
!all _tmp-arrays are neaded for parallelization (thread-private arrays)
!
!
end module dime3d
!
!-----------------------------------------------------------------------
!------------------------input (from *.nml file)------------------------
!-----------------------------------------------------------------------
!
module params_input
!
use prog_type
!
implicit none
!
!for model atmosphere
real(dp) :: teff, trad, tmin, xlogg, rstar, rmax, vmin, vmax, beta, &
            vth_fiducial, vmicro, yhe, hei, xmloss, vrot, lstar
!
!for modelling the continuum
real(dp) :: kcont
!
!for modelling the line
real(dp) :: kline, eps_line, alpha, kappa0
!
!atomic number
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
use prog_type
!
implicit none
!
real(dp) :: sr, taux, tauy, tauz, smajorax_a, smajorax_b, smajorax_c
!
!
end module params_stellar
!
!-----------------------------------------------------------------------
!-----------------------iteration-variables-----------------------------
!-----------------------------------------------------------------------
!
module iter
!
use prog_type
implicit none
!
integer(i4b), parameter :: itmaxc=150, itmaxl=150
integer(i4b) :: it_start_cont, it_start_line
real(dp), parameter :: devmaxc=1.d-3, devmaxl=1.d-3
real(dp), dimension(itmaxc) :: epsmaxc_arr
real(dp), dimension(itmaxl) :: epsmaxl_arr
!
!devmax(l): maximum required percentage-deviation (continuum, line)
!eps1d: error for central-ray mean intensities
!
end module iter
!
!-----------------------------------------------------------------------
!----------------------------ng-extrapolation---------------------------
!-----------------------------------------------------------------------
!
module ng_extra
!
use prog_type
!
implicit none
!
integer(i4b), parameter :: ng_const=6
!
!ng_const: ng (or aitken)-extrapolation is performed at each ng_const iterate
!
end module ng_extra
!
!-----------------------------------------------------------------------
!---------------------------profile function----------------------------
!-----------------------------------------------------------------------
!
module profile
!
use prog_type
!
implicit none
!
real(dp), dimension(:,:,:,:,:), allocatable :: renormalization
!
!renormalization, dim(ndxmax,ndymax,ndzmax,dim_mu,dim_phi,nxobs)
!   => will never be used!!!
!
end module profile
!
!-----------------------------------------------------------------------
!------------------------sobolev approximation--------------------------
!-----------------------------------------------------------------------
!
module soboapprox
!
use prog_type
!
implicit none
!
real(dp), dimension(9,11) :: ufunc09
real(dp), dimension(10,11) :: ufunc10
!
end module soboapprox
!
!-----------------------------------------------------------------------
!------------------------time profiling---------------------------------
!-----------------------------------------------------------------------
!
module timing
!
use prog_type
!
implicit none
!
real(dp) :: ts_tot, te_tot
real(dp) :: ts_it, te_it, ttot_it
real(dp) :: tsl_it, tel_it, ttotl_it, ttotl_it_sys
real(dp) :: ts_sl, te_sl, ttot_sl
real(dp) :: ts_alo, te_alo, ttot_alo
real(dp) :: tsl_alo, tel_alo, ttotl_alo
real(dp) :: elapsed_time_tot
integer(i4b) :: it_tot, it_totl
real(dp) :: ttot_it_fvm, ttot_it_sc
!
real(dp) :: ts_case1, te_case1, tt_case1, &
            ts_case2, te_case2, tt_case2, &
            ts_case3, te_case3, tt_case3, &
            ts_interpu, te_interpu, tt_interpu, &
            ts_interpd, te_interpd, tt_interpd, &
            ts_aloo, te_aloo, tt_aloo, ts_fs1d, te_fs1d, tt_fs1d, &
            ts_integ, te_integ, tt_integ
!
!elapsed_time*: system-time, not wallclock time
!
!ts_tot: start of complete program
!te_tot: end   of complete program
!
!ts_it:  start of iteration (of complete difference-method)
!te_it:  end   of iteration (of complete difference-method)
!ttot_it: total time needed for all iterations
!it_tot: total number of iterations needed
!
!tsl_it:  start of iteration (of complete difference-method) (line-transport)
!tel_it:  end   of iteration (of complete difference-method) (line-transport)
!ttotl_it: total time needed for all iterations (line-transport)
!it_totl: total number of iterations needed (line-transport)
!
!ts_sl: start of searchlight test
!te_sl: end   of searchlight test
!ttot_sl: total time needed for complete searchlight test
!
!ts_alo: start of inverting alo (in each iteration)
!te_alo: end   of inverting alo (in each iteration)
!ttot_alo: total time to invert alo (each iterations added)
!
!tsl_alo: start of inverting alo (in each iteration) (line-transport)
!tel_alo: end   of inverting alo (in each iteration) (line-transport)
!ttotl_alo: total time to invert alo (each iterations added) (line-transport)
!
end module timing
!
!-----------------------------------------------------------------------
!----------------limits for information regions-------------------------
!-----------------------------------------------------------------------
!
module inf_reg
!
use prog_type
!
implicit none
!
real(dp) :: rmin, rlim !=9.96d0 !13.2d0 !3.2d0 !4.6d0 !3.d0 !12.1d0
!
!rlim is the maximum radius up to which the radiative transport shall
!     be calculated
!     (should be chosten to be slightly(.5 to 1 %) larger than rmax
end module inf_reg
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_benchmark
!
use prog_type
!
implicit none
!
!input variables
integer(i4b) :: im_source, im_opacity, im_vel, benchmark_mod
real(dp) :: tau_min, tau_max, source_min, source_max, n_y, n_z, nn_y, nn_z
!
!reference 1d grid
integer(i4b), parameter :: nr=101
real(dp), dimension(nr) :: opac1d, scont1d, r1d, tau1d, int1d, int1d_sc, int1d_scray, int1d_fvmray, int1d_fvm, &
                           abs1d, abs1d_sc, abs1d_fvm, contr1d, contr1d_sc, contr1d_fvm
!
!for benchmark01 (and benchmark 02)
real(dp), dimension(:,:), allocatable :: int2d_fvm, abs2d_fvm, contr2d_fvm, int2d_sc, abs2d_sc, contr2d_sc, opac2d, scont2d
real(dp), dimension(:,:), allocatable :: x_u2d, z_u2d, x_d2d, z_d2d,  &
                                         scont_u2d, opac_u2d, scont_d2d, opac_d2d, int_u2d
real(dp), dimension(:,:,:), allocatable :: xu_interp_int2d, xu_interp_opac2d, zu_interp_int2d, zu_interp_opac2d
real(dp), dimension(:,:,:), allocatable :: xd_interp_opac2d, zd_interp_opac2d
!
!for benchmark03
real(dp), dimension(:,:,:), allocatable :: int2dsc_angdep, int2dfvm_angdep
!
!for benchmark04
integer(i4b) :: n1d_angdep
real(dp), dimension(:), allocatable :: mint1d_theo, mint1d_sc, mint1d_fvm
real(dp), dimension(:,:), allocatable :: intsc_angdep, intfvm_angdep
real(dp), dimension(:), allocatable :: r1d_angdep
!
!for benchmark05 (and benchmark12)
real(dp), dimension(:), allocatable :: mint1d_joray, mint1d_jomom, epsmaxc_sc, epsmaxc_fvm, &
                                       fcont1d_joray, fcont1d_jomom, &
                                       opac1d_cr, opac1d_jo, t1d_cr, t1d_jo
!
!for benchmark06
integer(i4b) :: nd, nd_fine, ndc_fine, nd_method01, nd_method02, nd_method03, nd_method04, &
                ndc_method01, ndc_method02, ndc_method03, ndc_method04
real(dp) :: xobs1, xobs2, xobs3
real(dp), dimension(:), allocatable :: s1d_method01, s1d_method02, s1d_method03, s1d_method04, s1d_fine, &
                                       tau1d_method01, tau1d_method02, tau1d_method03, tau1d_method04, tau1d_fine, &
                                       vth1d_method01, vth1d_method02, vth1d_method03, vth1d_method04, vth1d_fine, &
                                       vel1d_method01, vel1d_method02, vel1d_method03, vel1d_method04, vel1d_fine, &
                                       xcmf1d_method01, xcmf1d_method02, xcmf1d_method03, xcmf1d_method04, xcmf1d_fine, &
                                       profile1d_method01, profile1d_method02, profile1d_method03, profile1d_method04, profile1d_fine, &
                                       opalbar1d_method01, opalbar1d_method02, opalbar1d_method03, opalbar1d_method04, opalbar1d_fine, &
                                       opal1d_method01, opal1d_method02, opal1d_method03, opal1d_method04, opal1d_fine, &
                                       sline1d_method01, sline1d_method02, sline1d_method03, sline1d_method04, sline1d_fine, &
                                       int1d_method01, int1d_method02, int1d_method03, int1d_method04, int1d_fine, &
                                       s1dc_method01, s1dc_method02, s1dc_method03, s1dc_method04, s1dc_fine, &
                                       tau1dc_method01, tau1dc_method02, tau1dc_method03, tau1dc_method04, tau1dc_fine, &
                                       vth1dc_method01, vth1dc_method02, vth1dc_method03, vth1dc_method04, vth1dc_fine, &
                                       vel1dc_method01, vel1dc_method02, vel1dc_method03, vel1dc_method04, vel1dc_fine, &
                                       xcmf1dc_method01, xcmf1dc_method02, xcmf1dc_method03, xcmf1dc_method04, xcmf1dc_fine, &
                                       profile1dc_method01, profile1dc_method02, profile1dc_method03, profile1dc_method04, profile1dc_fine, &
                                       opalbar1dc_method01, opalbar1dc_method02, opalbar1dc_method03, opalbar1dc_method04, opalbar1dc_fine, &
                                       opal1dc_method01, opal1dc_method02, opal1dc_method03, opal1dc_method04, opal1dc_fine, &
                                       opac1dc_method01, opac1dc_method02, opac1dc_method03, opac1dc_method04, opac1dc_fine, &
                                       opatot1dc_method01, opatot1dc_method02, opatot1dc_method03, opatot1dc_method04, opatot1dc_fine, &
                                       scont1dc_method01, scont1dc_method02, scont1dc_method03, scont1dc_method04, scont1dc_fine, &
                                       stot1dc_method01, stot1dc_method02, stot1dc_method03, stot1dc_method04, stot1dc_fine, &
                                       sline1dc_method01, sline1dc_method02, sline1dc_method03, sline1dc_method04, sline1dc_fine, &
                                       int1dc_method01, int1dc_method02, int1dc_method03, int1dc_method04, int1dc_fine
!
!for benchmark07
real(dp), dimension(:), allocatable :: ssobo1d_crx, ssobo1d_cry, ssobo1d_crz, velr1d_cr, opalbar1d_cr, &
                                       ssobo1d
!
!for benchmark08 (and benchmark13 and benchmark14)
real(dp), dimension(:), allocatable :: scont1d_jo, sline1d_jo, ssobo1d_jo, sline1d_sc, sline1d_fvm, ssobo1d_cr, &
                                       epsmaxl_sc, epsmaxl_fvm, opalbar1d_jo, mintbar1d_sc, mintbar1d_fvm
!
!for benchmark10
real(dp), dimension(:,:,:), allocatable :: int3d_fvm, int3d_sc
real(dp), dimension(:,:), allocatable :: intsc_angdep2, intfvm_angdep2
real(dp), dimension(:), allocatable :: xcoord_angdep2, ycoord_angdep2, zcoord_angdep2
!
!for benchmark11
real(dp), dimension(:,:,:), allocatable :: mint3d_theo, fcontr3d_theo, fcontth3d_theo, fcontphi3d_theo
!
!for benchmark12
real(dp), dimension(:,:,:), allocatable :: mint3d_fvm, mint3d_sc, scont3d_fvm, scont3d_sc, &
                                           fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm, &
                                           fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, &
                                           kcontrr3d_fvm, kcontthth3d_fvm, kcontphiphi3d_fvm, &
                                           kcontrth3d_fvm, kcontrphi3d_fvm, kcontthphi3d_fvm, &
                                           kcontrr3d_sc, kcontthth3d_sc, kcontphiphi3d_sc, &
                                           kcontrth3d_sc, kcontrphi3d_sc, kcontthphi3d_sc
                                           
!
!for benchmark13
real(dp), dimension(:,:,:), allocatable :: mintbar3d_fvm, mintbar3d_sc, sline3d_fvm, sline3d_sc
!
!benchmark_mod=1 for 2d solution along a given direction specified by n_z
!benchmark_mod=2 for 2d searchlight beam test along a given direction specified by n_z
!benchmark_mod=3 for 2d searchlight beam test along all given mu-directions
!benchmark_mod=4 for 2d optically thin continuum (theoretical solution by dilution factor)
!benchmark_mod=5 for 2d spherically symmetric solutions of continuum (theoretical solution from JOs 1d program)
!benchmark_mod=6 for 1d line solution given different velocity fields
!benchmark_mod=7 for 3d sobolev source function for spherically symmetric problems, compared with 1d sobolev solution
!benchmark_mod=8 for 2d spherically symmetric solutions of the line (theoretical solution from JOs 1d program and
!                       from 3d sobolev approach)
!benchmark_mod=9 for 2d timing properties
!benchmark_mod=10 for 3d searchlight beam test along a given direction specified by n_y, n_z
!benchmark_mod=11 for 3d optically thin continuum (theoretical solution by dilution factor)
!benchmark_mod=12 for 3d spherically symmetric solutions of continuum (theoretical solution from JOs 1d program)
!benchmark_mod=13 for 3d spherically symmetric solutions of line transfer (theoretical solution from JOs 1d program)
!
!
end module mod_benchmark
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_debug
!
use prog_type
!
integer(i4b) :: indx1, indx2, indxx, indxy, indxz
real(dp) :: xu_debug, yu_debug, zu_debug, &
            nnx_debug, nny_debug, nnz_debug, &
            x_ii, y_ii, z_ii, &
            opalbar_ii, velx_ii, vely_ii, velz_ii, opac_ii
!$omp threadprivate(xu_debug, yu_debug, zu_debug, nnx_debug, nny_debug, nnz_debug, x_ii, y_ii, z_ii, opalbar_ii, velx_ii, vely_ii, velz_ii, opac_ii)
!
!
end module mod_debug
