!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------------PROGRAM----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!calculates radiative transfer by using short-characteristics method
!
!v0: including subroutines to read atmospheric strucutre from s
!       external files
!    including subroutines to make opacity-models
!    including subroutines to create angular grids
!    including subroutines to create frequency grid
!    including subroutine fsc_cont1d (1d formal solution using short characteristics
!                                     for continuum radiative transfer)
!    including subroutine fsc_cont2d (2d formal solution using short characteristics
!                                     for continuum radiative transfer)
!    including subroutine fsc_cont2d_debug (2d formal solution using short characteristics
!                                           for continuum radiative transfer
!                                           -> including debugging output)
!    including subroutine make_benchmark (for different benchmark models)
!    including benchmark01: 1d and 2d short characteristics solution with
!                           analytic solutions along a direction
!    including benchmark02: 2d searchlight beam test along a direction
!    including benchmark03: 2d searchlight beam test along all mu directions
!    including benchmark04: 2d optically thin continuum (compared to dilution factor)
!    including subroutine mint_sc2d to calculate mean intensity on z-axis for 2d short characteristics
!    including subroutines ng_expol1d, ait_expol1d to extrapolate subsequent solutions
!    including subroutine conttrans_sc2d to calculate complete continuum radiative transfer
!                           -> iteration scheme, so far using ng-extrapolation or aitkens-extrapolation
!    including subroutine ffvm_cont2d (2d formal solution using finite volume method
!                                      for continuum radiative transfer)
!    including subroutine mint_fvm2d to calculate mean intensity on z-axis for 2d fvm method
!                                                 (and pseudo alo)
!    including subroutine conttrans_fvm2d to calculate complete continuum radiative transfer
!                           -> iteration scheme, so far using ng-extrapolation or atikens-extrapolation
!
!    including benchmark05: 2d continuum solution for 1d spherically symmetric models
!                           using short characteristics and finite volume method
!                           (compared to input-solution from JOs 1d program)
!
!v1:  including subroutine fsc_line1d (1d formal solution using short characteristics
!                                      for line transfer)
!     including subroutine fsc_line1d_sspace (1d formal solution using short characteristics
!                                      for line transfer, integration in s-space (spatial))
!     including subroutine fsc_line1d_debug (1d formal soluion using short characteristics
!                                      for line transfer, with additional output for debugging reasons)
!     including subroutine benchmark06: 1d short characteristics solution to test the refinement
!                                       of resonance zones
!     including subroutine mintbar_sc2d to calculate scattering integral on z-axis for 2d sc method
!     including subroutine linetrans_sc2d to calculate complete line radiative transfer
!                           -> iteration scheme, so far using ng-extrapolation or aitkens-extrapolation
!     including subroutine ffvm_line2d (1d formal soluion using finite volume method
!                                      for line transfer)
!     including subroutine mintbar_fvm2d to calculate scattering integral on z-axis for 2d fvm method
!     including subroutine linetrans_fvm2d to calculate complete line radiative transfer
!                           -> iteration scheme, so far using ng-extrapolation or aitkens-extrapolation
!     including subroutine output_itestep_line: saving current iterate in file
!     including subroutine input_itstep_line: reading iterate from file (if option set)
!     including subroutine sobo3d to calculate 3d sobolev source functions (neglecting continuum and multiple resonances)
!     including subroutine calc_startval_line: reading iterate from file (if option set) or
!                                              setting start-value to sobolev source function
!     including subroutine benchmark07: benchmark for 3d sobolev solution (compared with 1d sobolev solution)
!     including subroutine benchmark08: 2d line solution for 1d spherically symmetric models using
!                                       sc, fvm, and sobolev method (compared to input-solution
!                                       from JOs 1d program)
!     including subroutine benchmark09: test timing for 2d line-solution
!
!v2: new versions of fsc_cont2d, fsc_line2d, fsc_cont1d, fsc_cont1db
!       -> including (pseudo) diagonal alo
!    including subroutine print_input: sumnmary of all input parameters
!    including subroutine fsc_cont3d: 3d formal solution using short characteristics
!                                     for continuum radiative transfer
!    including subroutine benchmark10: searchlight test for 3d routines
!    including subroutine setup_phot1d: interpolating photospheric model from
!                                       1d Kurucz data onto own grid
!    including new calculation of boundary condition
!
!v3: including subroutine interpolation_test2d: testing 2d upwind interpolation
!                                               to be used in fsc_cont3d
!    including subroutine mint_sc3d: angular loop for intensities to obtain mean intensity
!    including angular integration in subroutine fsc_cont3d
!    including subroutine calc_dev3d: calculating deviation of 3d arrays
!    including subroutine scont_new3d: new iterate of continuum source function
!    including subroutine scont_new3d_classic: new iterate via classical lambda iteration
!    including subroutine scont_new3d_diag: new iterate via ALI with diagonal ALO
!    including subroutines store_ng3d: stores 3d source function as 1d array in each iterations tep
!    including subroutines ng_expol3d, ait_expol3d: extrapolation of 3d source functions given
!                                                   at subsequent iteration steps and as 1d array!
!    including benchmark12: 3d continuum solution for 1d spherically symmetric models
!                           using short characteristics and finite volume method
!                           (compared to input-solution from JOs 1d program)
!    including subroutines matmul_coo, jsor_coo, newit_jor_coo, newit_sor_crs, coocsr, 
!                          scont_new3d_nn, calc_alocont_nn_coo for direct neighbour alo-calculations
!    including nearest neighbour alo for linear source integrals
!
!v4: several bug fixes
!
!v5: deleting options mu_grid, phi_grid, since automatically coupled to integration method
!    combining subroutines for calculating angular integration nodes (mu and phi together)
!    including subroutines for calculating lebedev integration nodes and triangulated
!         integration nodes
!
!v6: including nearest neighbour ALO for SC solution: extending to all 26 surrounding points;
!    updating subroutines scont_new3d_diag, scont_new3d_nn for new ALO structure;
!    including subroutines coeff2d_* to calculate interpolation coefficients;
!
!v7: including subroutine coeff2d_contu: performing interpolation and coefficient-calculations
!                                        for all physical properties used in continuum calculations
!                                        in one washup
!    including subroutine coeff2d_contd: performing interpolation and coefficient-calculations
!                                        for all physical properties used in continuum calculations
!                                        in one washup
!    including subroutine coeff3d_contu: performing interpolation and coefficient-calculations
!                                        for all physical properties used in continuum calculations
!                                        in one washup
!
!v8: including subroutine ffvm_line3d: finite volume formal solution for line transfer
!    including subroutine fsc_line3d_lin: 3d short characteristic solution for line transfer
!                                       using purely linear interpolations and including refinement
!    including subroutine fsc_line1d_lin: 1d short characteristic solution for line transfer
!                                         using purely linear interpolations and including refinement
!    including subroutine coeff2d_lineu_lin: performing interpolation and coefficient-calculations
!                                            for all physical properties used in line calculations
!                                            in one washup (linear interpolations)
!    including subroutine coeff2d_lined_lin: performing interpolation and coefficient-calculations
!                                            for all physical properties used in line calculations
!                                            in one washup (linear interpolations)
!    including subroutine coeff3d_lineu_lin: performing interpolation and coefficient-calculations
!                                            for all physical properties used in line calculations
!                                            in one washup (linear interpolations)
!    including benchmark13: 3d line transfer solution for 1d spherically symmetric models
!                           using short characteristics and finite volume method
!                           (compared to input-solution from JOs 1d program)
!
!v9: splitting some subroutines
!    including analytical expressions for velocity and frequency integrated
!        line opacity in subrouines coeff2d_lineu_lin and coeff2d_lined_lin
!    rewriting subroutine fsc_line1d and fsc_line1d_debug:
!        analytic calculation of del-tau steps
!    including subroutine fsc_linec1d_lin:  solution for line and continuum transfer 
!        using purely linear interpolations
!    including subroutine fsc_linec3d_lin: 3d short characteristic solution for line and continuum transfer
!         using purely linear interpolations
!    including subroutine coeff2d_linecu_lin: performing interpolation and coefficient-calculations
!         for all physical properties used in line + continuum calculations
!         in one washup (linear interpolations)
!    including subroutine coeff3d_linecu_lin: performing interpolation and coefficient-calculations
!         for all physical properties used in line + continuum calculations
!         in one washup (linear interpolations)
!    including benchmark14: 3d line + continnum transfer solution for 1d spherically symmetric models
!                         using short characteristics and finite volume method
!                         (compared to input-solution from JOs 1d program)
!    including subroutine fsc_linec1d:  solution for line and continuum transfer 
!        using bezier interpolations
!    including subroutine fsc_linec3d: 3d short characteristic solution for line and continuum transfer
!         using bezier interpolations
!    including subroutine coeff2d_linecu: performing interpolation and coefficient-calculations
!         for all physical properties used in line + continuum calculations
!         in one washup (bezier interpolations) on upwind face
!    including subroutine coeff2d_linecd: performing interpolation and coefficient-calculations
!         for all physical properties used in line + continuum calculations
!         in one washup (bezier interpolations) on downwind face
!    some debugs in fsc_line1d (deleting error-function approach, and using numerical approach
!          for delta-tau integrals)
!    including extrapolation of input-model outside of calculation volume
!          (required for short-characteristics interpolations)
!    including subroutine recalc_vgridxyz: recalculates xyz-grids from x=rmin to x=rmax
!          using equidistant velocity steps (if delv>1/3), and recalculates 
!          model on that new grid
!
!v10 completely new grid construction procedures: 
!    including subroutine gridxyz_optb: calculation of probability densities for x, y, z coordinates
!    including subroutine recalc_vgridxyz2: recalculates xyz-grids from x=rmin to x=rmax
!          using equidistant velocity steps if necessary (the 'true' velocity law shall be given by
!          a 1d beta velocity law, important for low beta!!!)
!
!v11   including subroutines to calculate dels_r for distorted stellar surface
!          (stellar surface is approximated by ellipsoid!!!)
!          -> dist_surface_rot
!          -> info_region_rot
!       including opt_ltec to calculate grey temperature stratification in radiative equilibrium
!       including opt_opac to define the opacity law
!
!v12    including opt_method=2 for quadratic bezier SC method within contiinuum transport
!       combining subroutines conttrans_sc3d, conttrans_fvm3d in conttrans3d
!       combining subroutines conttrans_sc2d, conttrans_fvm2d in conttrans2d
!       updating benchmark11 for parallel computing
!
!v13  Including a depth-dependent eps_cont -> eps_cont3d
!      sobo3d moved to src/
!      Replacing mask_totreg3d with imask_totreg3d
!                mask_innreg3d with imask_innreg3d
!                mask_bpoint3d with imask_bpoint3d
!      Inclouding calculation of the pressure tensor
!
!vNEW  Restructuring all interpolation modules (now in separate module files)
!      Including a 3D AMRVAC reader in modules mod_model.f90
!      debugging: setup_mod1d, setup_mod2d, setup_mod3d to enable calculation of opacities at the outer boundary
!                         (within r<= rmax, and not as before within r<=rlim)
!      including opt_method=2 for quadratic bezier SC method within 3d line transport
!      combining subroutines linetrans_sc3d, linetrans_fvm3d in linetrans3d
!      combining subroutines linetrans_sc2d, linetrans_fvm2d in linetrans2d (not tested yet)
!      all interpolation and integration routines modularized
!
!*************************NOTES/TODO's**********************************
!
!IMPORTANT NOTE: FVM-METHOD CAN NOT HANDLE ROTATIONAL DISTORTED SURFACE YET!!!!!
!                => TODO: CALCULATION OF DELX, DELY, DELZ STEPS + I_CORE ON STELLAR SURFACE
!
!IMPORTANT NOTE: ADAPT SOURCE-CONTRIBUTION INTEGRALS ETC. FOR OWN PURPOSE:
!                tau_u:=0.d0, s_u:=0.d0
!                => can accelerate code-performance!!!!
!
!IMPORTANT: update fsc_linec1d and fsc_line1d: when no refinement, and resonance zone in downwind interval, delt_d>>delt_u
!                => source contribution=?
!
!output_benchmark13: do not read in fvm-solution (only for tests)
!benchmark13: do not return before calculating fvm solution
!
!iteration scheme for line+(thick)continuum: solve problems
!          with overestimation of source corrections
!          (occurring only when bezier routines are used)
!
!interpolation of source functions on input-model grid!!!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------START-------------------------------------------
!-----------------------------------------------------------------------
!
program main
!
use prog_type
use options, only: spatial_grid3d, input_mod_dim, opt_opal, opt_incl_cont, opt_incl_line, opt_sol2d, opt_method, opt_opac, opt_ltec
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, velx3d, vely3d, velz3d, vth3d, t3d, opalbar3d, ssobo3d, imask_totreg3d
use mod_benchmark, only: benchmark_mod
use bcondition, only: xic1
use freq, only: xnue0
use params_input, only: vmin, vmax, beta, eps_line, vth_fiducial
use omp_lib
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: twall_stot, twall_etot
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
!
!
!
!call test_trash
!stop
!
twall_stot=omp_get_wtime()
!
call init_timing
!
call init_warnings
!
!
!----------------read all input (including stellar parameters,----------
!------------------------options, grid-sizes, etc.)---------------------
!
call read_input
!
!-----------------check if input options are compatible-----------------
!
call check_options
!
!--------------------calculating x, y, z coordinates--------------------
!
write(*,*) '-------------------------calculating x, y, z grids-----------------------------'
write(*,*)
!
select case(spatial_grid3d)
   case(0)
      call gridxyz
   case(1)
!      call gridxyz_opt
      call gridxyz_optb
   case(2)
!      call gridxyz_opt2d
      call gridxyz_opt2db
   case(3)
      call gridxyz_equi
   case(4)
      call gridxyz_opt1d
   case(5)
      call gridxyz_opt3d
   case default
      stop 'error in main: option spatial_grid3d not set'
end select
!
!recalculate 3d grids, such that radial steps are not too large
call recalc_gridxyz
!
!check 3d grid
call check_gridxyz
!
!---------------------interpolating model atmosphere--------------------
!
write(*,*) '----------------------interpolating model atmosphere---------------------------'
write(*,*)
!
call allocate_global3d
!
select case(input_mod_dim)
   case(1)
      call setup_mod1d
   case(2)
      call setup_mod2d
   case(3)
      call setup_mod3d
   case default
      stop 'error in main: option input_mod_dim not set'
end select
!
!-----------------------------------------------------------------------
!
if(opt_incl_line) then
!
!only required for sc-method when line transfer is included: start
!recalculate 3d grids, sucht that velocity steps are not too large
!
   write(*,*) '-------------recalculating xyz-grids to optimize velocity steps----------------'
   write(*,*)
!
!recalculate the x,y,z grids, optimized for equal velocity steps
   select case(spatial_grid3d)
      case(0,3,4,5)
         write(*,*) 'no recalculation for spatial_grid3d=0,3,4,5, since grid completely defined (how bad it ever may be)'
         write(*,*)
      case(1)
         call recalc_vgridxyz2  !recalculate grid by interpolating velocities from 1d radial grid (required especially if velocity law is steep!!!)
      case(2)
         call recalc_vgridxyz   !recalculate grid by interpolating velocities from old grid
      case default
         stop 'error in main: option spatial_grid3d not set'
   end select
!
!check 3d grid
   call check_gridxyz
!
!interpolationg model atmosphere on new grid
   write(*,*) '----------------------interpolating model atmosphere---------------------------'
   write(*,*)
!
   select case(input_mod_dim)
      case(1)
         call setup_mod1d
      case(2)
         call setup_mod2d
      case(3)
         call setup_mod3d
      case default
         stop 'error in main: option input_mod_dim not set'
   end select
!
endif
!
!only required for sc-method: end
!
!-----------------------------------------------------------------------
!
!calculating new mask (required for sc method)
call setup_mask3d
!
!----------------------calculating opacity model------------------------
!
write(*,*) '-----------------------setting up continuum opacity model----------------------'
write(*,*)
!
select case(opt_opac)
   case(0)
     call setup_opac1
   case default
      stop 'error in main: option opt_opac not set'
end select
!
write(*,*) '-------------------------setting up line opacity model-------------------------'
write(*,*)
!
select case(opt_opal)
   case(0)
      call setup_opalbar1
   case(1)
      call setup_opalbar2
   case default
      stop 'error in main: option opt_opal not set'
end select
!
!-----------------calculate phantom points on the axes------------------
!
call set_xyz_phantom
!
!--------------------interpolating photospheric model-------------------
!
write(*,*) '---------------------interpolating photospheric model--------------------------'
write(*,*)
!
call setup_phot1d
!
!-----------------------------------------------------------------------
!
!calculate delx, dely, delz as array (for FVM)
if(opt_method.eq.0) call grid_delxyz
!
!--------------------------check units----------------------------------
!
call check_units
!
!------------calculating angular and frequency grids--------------------
!
!calculate number of frequency integration nodes
call dime_xobs_nodes
!
!frequency grid
call calcnodes_xobs
call check_xobs
!
!theta/mu and phi grid
call calcnodes_omega
!
!------------------------print out input model--------------------------
!
!print out all dimensions
call print_dimensions
!
!print out input parameters
call print_input
!
!print out atmsopheric structure
call print_model
!
!-----------------------------------------------------------------------
!
call check_grid3d
!
!-----------------------------------------------------------------------
!
!calculate grid for central rays
call grid_cr
!
!-----------------------------------------------------------------------
!
!boundary condition from diffusion approximation
call calc_bcondition
!
!boundary conditions inside the photosphere
call calc_bcondition_phot
!
!calculate and store intensities for inner boundary points
!call calc_bintensity
!
!---------------------performing several tests--------------------------
!
!check different integration schemes
call integration_test
call integration_test2
!
!check different interpolation schemes
call interpolation_test
!call interpolation_test2
!
!check 2d interpolation
!call interpolation_test2d
!
!check integration of profile function
!call profile_test
!
!----------------------dummy subroutine for debugging-------------------
!
!call make_debug
!
!-----------------------------------------------------------------------
!
!call test_timing
!stop
!
call make_benchmark
!
!stop program, because opacities and source functions have been overwritten
if(benchmark_mod.ne.0) then 
   twall_etot=omp_get_wtime()
!
   write(*,*) 'total computation (wallclock) time', twall_etot-twall_stot
!
   write(*,*)
   stop '----------------------benchmark and main program done--------------------------'
endif
!
!-----------------------------------------------------------------------
!
!call test_timing
!call output
!
!--------------test calculation of boundary intensities-----------------
!
!CALL TEST_GET_IC
!
!----------------------print out all input parameter--------------------
!
!and print out continuum optical depths along axes
!CALL CALC_TAUC_XYZ
!
!-----------------------output atmospheric structure--------------------
!
!-----------------------------------------------------------------------
!
if (opt_incl_cont) then
!
!---------------calculate continuum radiative transport-----------------
!
!conttrans2d calculates
!   1. intensities in x-z plane for all angles mu
!   2. mean intensities by integration on central z-axis (phi-symmetry)
!   3. source function on central z-axis
!   4. source function on 3d-grid by interpolation
!
   if(opt_sol2d) then
      call conttrans2d
   else
      call conttrans3d
   endif
!
endif
!
!
write(*,*) '----------------------------setting up temperatures----------------------------'
write(*,*)
!
!update temperatures if option set
select case(opt_ltec)
   case(0)
   case(1)
      call setup_temperatures1
   case default
      stop 'error in main: opt_ltec not properly specified'
end select
!
!
if(.not.opt_incl_line) then
   call output
   stop
endif
!
!
!
!if(opt_alocont.eq.3) then
!   call output_sparsity(outputdir, ndxmax*ndymax*ndzmax, alocont_row_indx, alocont_col_indx, alocont_data)
!endif
!
!
!if(.not.opt_incl_line) stop 'end of program: no line-transfer calculated'
!
!
!----------------CALCULATE LAMBDA MATRIX--------------------------------
!
!!CALL CALC_OPERATORS
!
!!CALL SAVE_OPERATORS
!
!!CALL OUTPUT(outputDIR)
!
!!STOP
!
!!CALL READ_OPERATORS

!!CALL CHECK_ALO
!!CALL CONTTRANS_3D
!
!-----------------------TEST LAMBDA MATRIX------------------------------
!
!!CALL TESTLAMBS_ITER
!!CALL TESTLAMBS_DIRECT
!
!!CALL OUTPUT(outputDIR)

!!STOP
!
!
!----------------------calculate line transport-------------------------
!------------------------numerical solution-----------------------------
!
!-------------------calculate sobolev source-function-------------------
!
!
!without continuum
call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
!
!!with continuum (needs to be debugged, thus, use always without continuum so far)
!IF (OPT_INCL_CONT) THEN
!   CALL CALC_SOBOLEV_CONT
!ENDIF
!
if(opt_sol2d) then
   call linetrans2d
else
   call linetrans3d
endif
!
!-----------------------------------------------------------------------
!-----------------------OUTPUT------------------------------------------
!-----------------------------------------------------------------------
!
!CALL CPU_TIME(TE_TOT)
!CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
!nb_ticks = nb_ticks_final - nb_ticks_initial
!elapsed_time_tot   = dble(nb_ticks) / nb_ticks_sec
!CALL PRRTIME
!CALL OUTPUT(outputDIR)
!

!print out warnings
call print_warnings
!
call output
!
twall_etot=omp_get_wtime()
!
write(*,*) 'total computation (wallclock) time', twall_etot-twall_stot
!
write(*,*)
stop '-----------------------------main program done---------------------------------'
!
end program main
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!---------------------SUBROUTINES AND FUNCTIONS-------------------------
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine test_trash
!
  use prog_type
  use mod_interp1d
  use mod_integ1d, only: integral2, integral2b, integral3, integral3b
!
implicit none
!
integer(i4b) :: i, li
integer(i4b), parameter :: nd=101
real(dp) :: a, b, c, dfdx_im1, dfdx_i
real(dp) :: x, y2, y3, y2b, y3b, dx
real(dp) :: x_im1, x_i, y2_im1, y2_i, y3_im1, y3_i, integ2, integ3
!
real(dp), dimension(3) :: xp, yp
!
!
!-----------------------testing interpolation---------------------------
!
xp=(/ 1.d0, 1.5d0, 3.d0 /)
yp=(/ 0.d0, 1.d0, 0.5d0 /)
!
!yp=(/ 1.45766689d-5, 1.62739987d-5, 1.66346311d-5 /)
!xp=(/ 0.d0, 8.09150723d-2, 2.09954004d-1 /)
!
open(1, file='TRASH/interp.dat', form='formatted')
write(1,'(5a20)') 'x', 'y_quad2', 'y_quad3', 'y_quad2b', 'y_quad3b'
!
do i=1, nd
   x = xp(1) + (i-1)*(xp(3)-xp(1))/(nd-1)
!   call coeff_typ_quad3(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x(i), &
!                        a, b, c)
!   y1 = a*yp(1) + b*yp(2) + c*yp(3)

   y2 = interpol_typ_quad2(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x)
   y3 = interpol_typ_quad3(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x)

   y2b = interpol_typ_quad2b(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x)
   y3b = interpol_typ_quad3b(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x)

   li=0
   if(x.eq.xp(1).or.x.eq.xp(2).or.x.eq.xp(3)) li=1
   write(1, '(5es20.8, i5)') x, y2, y3, y2b, y3b, li
enddo
!
!-----------------------testing integration-----------------------------
!
!trapezoidal integration from xp(0) to xp(1)
integ2=0.d0
integ3=0.d0
dx=(xp(2)-xp(1))/(nd-1)
y2_im1 = yp(1)
y3_im1 = yp(1)
x_im1 = xp(1)
do i=2, nd
   x_i = x_im1 + dx
   y2_i = interpol_typ_quad2(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x_i)
   y3_i = interpol_typ_quad3(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x_i)
!
   integ2 = integ2 + (y2_i+y2_im1)*dx/2.d0
   integ3 = integ3 + (y3_i+y3_im1)*dx/2.d0
!
   x_im1 = x_i
   y2_im1 = y2_i
   y3_im1 = y3_i
enddo
!
write(*,*) 'interval [x1,x2]'
write(*,*) 'quad2', integ2, integral2(xp(1), xp(2), xp(3), yp(1), yp(2), yp(3))
write(*,*) 'quad3', integ3, integral3(xp(1), xp(2), xp(3), yp(1), yp(2), yp(3))
!
!trapezoidal integration from xp(1) to xp(0)
integ2=0.d0
integ3=0.d0
dx=(xp(3)-xp(2))/(nd-1)
y2_im1 = yp(2)
y3_im1 = yp(2)
x_im1 = xp(2)
do i=2, nd
   x_i = x_im1 + dx
   y2_i = interpol_typ_quad2b(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x_i)
   y3_i = interpol_typ_quad3b(yp(1), yp(2), yp(3), xp(1), xp(2), xp(3), x_i)
!
   integ2 = integ2 + (y2_i+y2_im1)*dx/2.d0
   integ3 = integ3 + (y3_i+y3_im1)*dx/2.d0
!
   x_im1 = x_i
   y2_im1 = y2_i
   y3_im1 = y3_i
enddo
!
write(*,*) 'interval [x2,x3]'
write(*,*) 'quad2', integ2, integral2b(xp(1), xp(2), xp(3), yp(1), yp(2), yp(3))
write(*,*) 'quad3', integ3, integral3b(xp(1), xp(2), xp(3), yp(1), yp(2), yp(3))
!


!
end subroutine test_trash
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine init_timing
!
use prog_type
!
implicit none
!
write(*,*) 'TODO: subroutine init_timing still to be implemented'
write(*,*)
!
end subroutine init_timing
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine init_warnings
!
use prog_type
use warnings, only: warn_mu, warn_phi, warn_angles
!
implicit none
!
warn_mu=.false.
warn_phi=.false.
warn_angles=.false.
!
end subroutine init_warnings
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_timing
!
!---------------print out runtime properties----------------------------
!
use prog_type
!
implicit none
!
write(*,*) 'TODO: subroutine print_timing still to be implemented'
write(*,*)
!
end subroutine print_timing
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_warnings
!
!---------------print out warnings--------------------------------------
!
use prog_type
use warnings
!
implicit none
!
write(*,*) '-------------------------------warnings----------------------------------------'
write(*,*)
if(warn_mu) write(*,*) 'warning: mu-grid is not symmetric'
if(warn_phi) write(*,*) 'warning: phi-grid is not symmetric'
if(warn_angles) write(*,*) 'warning: n_x, n_y, n_z are not symmetric'
write(*,*)
!
end subroutine print_warnings
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_dimensions
!
!------------print out applied all dimensions and applied grids---------
!
use prog_type
use dime1d, only: n1d
use dime3d, only: ncx, ncy, ncz, ndxmax, ndymax, ndzmax
use angles, only: dim_mu, dim_phi, dim_omega
use options, only: spatial_grid3d
use freq, only: nxobs
!
implicit none
!
! ... local characters
character(len=100) :: gridtxt
!
select case(spatial_grid3d)
   case(0)
      gridtxt = 'radial grid with equidistant core points'
   case(1)
      gridtxt = 'radial grid with minimizing subsequent coordinate distances'
   case(2)
      gridtxt = '2d-grid from model-file and minimizing subsequent coordinate distances'
   case(3)
      gridtxt = 'completely equidistant grid in x, y, z'
   case(4)
      gridtxt = '1d-grid from model-file and using corresponding radial pdf'
   case(5)
      gridtxt = '3d-grid from model-file and using corresponding radial and latitudinal pdf'
   case default
      stop 'error in print_dimensions: spatial_grid3d not set'
end select
!
!
write(*,*) '-------------------------------dimensions--------------------------------------'
write(*,*)
write(*,*) 'grid option: ', trim(gridtxt)
write(*,*) 'n1d=', n1d
write(*,*) 'ncx=', ncx
write(*,*) 'ncy=', ncy
write(*,*) 'ncz=', ncz
write(*,*) 'ndxmax=', ndxmax
write(*,*) 'ndymax=', ndymax
write(*,*) 'ndzmax=', ndzmax
write(*,*) 'n_theta=', dim_mu
write(*,*) 'n_phi=', dim_phi
write(*,*) 'n_xobs=', nxobs
write(*,*) 'number of actually used direction vectors: ', dim_omega
write(*,*)
write(*,*)
!
end subroutine print_dimensions
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_model
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, &
                  opac3d, opalbar3d, velx3d, vely3d, velz3d, vth3d, t3d, eps_cont3d
use params_input, only: vth_fiducial
!
implicit none
!
! ... locals calars
integer(i4b) :: i, j, k
!
write(*,*) '--------------------------model along x-axis-----------------------------------'
write(*,*)
!
j=ndymax/2+1
k=ndzmax/2+1
write(*,'(a20, 2es20.8)') 'along x at (y,z)= ', y(j), z(k)
write(*,'(9a20)') 'x [r_star]', 'opath [1/r_star]', 'opalbar [1/r_star]', &
           'velx [km/s]', 'vely [km/s]', 'velz [km/s]', 'vth [km/s]', 'temp [K]', 'eps_cont'
do i=1, ndxmax
   write(*,'(9es20.8)') x(i), opac3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k)*vth_fiducial/1.d5, &
              vely3d(i,j,k)*vth_fiducial/1.d5, velz3d(i,j,k)*vth_fiducial/1.d5, &
              vth3d(i,j,k)/1.d5, t3d(i,j,k), eps_cont3d(i,j,k)
enddo
write(*,*)
!
end subroutine print_model
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_input
!
use prog_type
use fund_const
use params_input, only: teff, trad, tmin, xlogg, rstar, lstar, rmax, vmin, vmax, beta, &
                        vth_fiducial, vmicro, yhe, hei, xmloss, kcont, kline, &
                        eps_line, alpha, kappa0, na
use params_stellar, only: sr, smajorax_a, smajorax_b, smajorax_c
use options, only: input_mod, opt_opal, opt_incl_line, opt_incl_cont, input_mod_dim, &
                   spatial_grid1d, spatial_grid3d, opt_angint_method, &
                   opt_method, opt_sol2d, opt_start_cont, opt_ng_cont, opt_ait_cont, &
                   opt_start_line, opt_ng_line, opt_ait_line, opt_incl_gdark, opt_incl_sdist
use freq, only: xnue0
use bcondition, only: xic1, xic2
use mod_benchmark, only: benchmark_mod, im_source, im_opacity, im_vel, tau_min, &
                         tau_max, source_min, source_max, n_z, n_y
!
implicit none

write(*,*) '------------------------summary of input parameter-----------------------------'
write(*,*)
!
write(*,*) 'general input'
write(*,'(a20, 2e20.8)') 'r_star [cm], [r_sun]', sr, rstar
write(*,'(a20, 2e20.8)') 'l_star [erg/s], [l_sun]', lstar*xlsu, lstar
write(*,'(a20, e20.8)') 'v_min [km/s]', vmin
write(*,'(a20, e20.8)') 'v_max [km/s]', vmax
write(*,'(a20, e20.8)') 'beta', beta
write(*,'(a20, e20.8)') 'xlogg', xlogg
write(*,'(a20, e20.8)') 'mdot', xmloss
write(*,'(a20, e20.8)') 'yhe', yhe
write(*,'(a20, e20.8)') 'hei', hei
write(*,'(a20, i20)') 'na', na
write(*,'(a20, e20.8)') 't_min [K]', tmin
write(*,'(a20, e20.8)') 't_eff [K]', teff
write(*,'(a20, e20.8)') 't_rad [K]', trad
write(*,'(a20, e20.8)') 'nue_0 [1/s]', xnue0
write(*,'(a20, e20.8)') 'lambda [a]', cgs_clight/xnue0*1.d8
write(*,'(a20, e20.8)') 'xic1', xic1
write(*,'(a20, e20.8)') 'xic2', xic2
!
write(*,*)
write(*,'(a20, e20.8)') 'vmicro [km/s]', vmicro
write(*,'(a20, e20.8)') 'vth_fiducial [km/s]', vth_fiducial
!
write(*,*)
write(*,'(a20, e20.8)') 'smajorax_a', smajorax_a
write(*,'(a20, e20.8)') 'smajorax_b', smajorax_b
write(*,'(a20, e20.8)') 'smajorax_c', smajorax_c
!
write(*,*)
if(.not.opt_incl_cont) write(*,*) 'continuum parameter not needed'
write(*,'(a20, e20.8)') 'k_cont', kcont

write(*,*)
if(.not.opt_incl_line) write(*,*) 'line parameter not needed'
write(*,'(a20, e20.8)') 'epsilon_l', eps_line
write(*,'(a20, e20.8)') 'k_line', kline
write(*,'(a20, e20.8)') 'alpha', alpha
write(*,'(a20, e20.8)') 'kappa0', kappa0
write(*,'(a20, l20)') 'opt_opal', opt_opal
!
write(*,*)
write(*,*) 'options'
write(*,'(a20, i20)') 'input_mod', input_mod
write(*,'(a20, i20)') 'input_mod_dim', input_mod_dim
write(*,'(a20, i20)') 'spatial_grid1d', spatial_grid1d
write(*,'(a20, i20)') 'spatial_grid3d', spatial_grid3d
write(*,'(a20, i20)') 'opt_angint_method', opt_angint_method
write(*,'(a20, i20)') 'opt_method', opt_method
write(*,'(a20, l20)') 'opt_sol2d', opt_sol2d
write(*,'(a20, l20)') 'opt_incl_cont', opt_incl_cont
write(*,'(a20, l20)') 'opt_start_cont', opt_start_cont
write(*,'(a20, l20)') 'opt_ng_cont', opt_ng_cont
write(*,'(a20, l20)') 'opt_ait_cont', opt_ait_cont
write(*,'(a20, l20)') 'opt_incl_line', opt_incl_line
write(*,'(a20, l20)') 'opt_start_line', opt_start_line
write(*,'(a20, l20)') 'opt_ng_line', opt_ng_line
write(*,'(a20, l20)') 'opt_ait_line', opt_ait_line
write(*,'(a20, l20)') 'opt_incl_gdark', opt_incl_gdark
write(*,'(a20, l20)') 'opt_incl_sdist', opt_incl_sdist

write(*,*)
write(*,*) 'benchmarks'
write(*,'(a20, i20)') 'benchmark_mod', benchmark_mod
write(*,'(a20, i20)') 'im_source', im_source
write(*,'(a20, i20)') 'im_opacity', im_source
write(*,'(a20, i20)') 'im_vel', im_vel
write(*,'(a20, es20.8)') 'tau_min', tau_min
write(*,'(a20, es20.8)') 'tau_max', tau_max
write(*,'(a20, es20.8)') 'source_min', source_min
write(*,'(a20, es20.8)') 'source_max', source_max
write(*,'(a20, es20.8)') 'n_y', n_z
write(*,'(a20, es20.8)') 'n_z', n_z
!
write(*,*)
!
end subroutine print_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!

use prog_type
use mod_directories, only: model_dir, input_file, output_file
use options, only: input_mod_dim, input_mod, spatial_grid1d, spatial_grid3d, opt_opal, &
                   opt_angint_method, opt_sol2d, opt_incl_cont, opt_start_cont, &
                   opt_ng_cont, opt_ait_cont, opt_method, opt_incl_line, opt_start_line, &
                   opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, opt_incl_gdark, &
                   opt_incl_sdist, opt_ltec, opt_opac
use params_input, only: teff, trad, tmin, xlogg, rstar, lstar, rmax, tmin, xmloss, vmin, &
                        vmax, vrot, beta, yhe, hei, na, vmicro, &
                        kcont, eps_line, kline, alpha, kappa0, vth_fiducial
use params_stellar, only: sr, smajorax_a, smajorax_b, smajorax_c
use freq, only: deltax, xcmf_max, xcmf_min, xnue0, lfreqint
use fund_const, only: rsu, pi, cgs_grav, xmsu
use dime1d, only: n1d, n1d_t, n1d_r, delv
use dime3d, only: ncx, ncy, ncz, delx_max, dely_max, delz_max, & 
                  ndx, ndy, ndz, ndxmax, ndymax, ndzmax
use angles, only: n_theta
use inf_reg, only: rmin, rlim
use mod_benchmark, only: benchmark_mod, im_source, im_opacity, im_vel, &
                         tau_min, tau_max, source_min, source_max, n_z, n_y, nn_z, nn_y
!
!!
implicit none
!
! ... local scalars
real(dp) :: mstar, eps_cont
!
! ... local characters
!
! ... local functions
real(dp) :: calc_req
!
! ... namelist
namelist / input_options / model_dir, output_file, input_mod, input_mod_dim, spatial_grid1d, spatial_grid3d, opt_opal, &
                           opt_angint_method, opt_method, opt_sol2d, opt_incl_cont, &
                           opt_start_cont, opt_ng_cont, opt_ait_cont, opt_incl_line, opt_start_line, &
                           opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, &
                           opt_incl_gdark, opt_incl_sdist, opt_ltec, opt_opac
namelist / input_mod_1d / teff, trad, xlogg, rstar, lstar, rmax, tmin, xmloss, vmin, vmax, &
                          vmicro, vth_fiducial, vrot, beta, yhe, hei, xnue0, na
namelist / input_cont / eps_cont, kcont
namelist / input_line / eps_line, kline, alpha, kappa0
namelist / dimensions_freq / deltax, xcmf_max
namelist / dimensions_angles / n_theta
namelist / dimensions_1d / n1d, n1d_t, n1d_r, delv
namelist / dimensions_3d / ncx, ncy, ncz, delx_max, dely_max, delz_max
namelist / input_infreg / rmin, rlim
namelist / benchmark / benchmark_mod, im_source, im_opacity, im_vel, tau_min, &
                       tau_max, source_min, source_max, n_y, n_z
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) input_file
write(*,*) 'reading input from file: ', trim(input_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(input_file), status='old', form='formatted')

!read options, i.e. which model shall be used
   rewind 1
   read(1, nml=input_options)

!read 1d model parameters
   rewind 1
   read(1, nml=input_mod_1d)
!
   sr=rstar*rsu
   tmin=tmin*teff
   rmin=1.d0
   vth_fiducial=vth_fiducial*1.d5
   vmicro=vmicro*1.d5
!
   if(opt_incl_sdist) then
!calculate distortion of stellar surface: equatorial radius
      mstar = sr**2 * 10.d0**xlogg/cgs_grav
      mstar = mstar/xmsu
      smajorax_a = calc_req(rstar, mstar, vrot)
      smajorax_b = smajorax_a
      smajorax_c = 1.d0
   else
!spherical star
      smajorax_a = 1.d0
      smajorax_b = 1.d0
      smajorax_c = 1.d0
   endif

!
!-------read dimensions for own 1d and 3d grids and corresponding-------
!-----------------------grid parameters---------------------------------
!
   rewind 1
   read(1, nml=dimensions_1d)
!
   rewind 1
   read(1, nml=dimensions_3d)
!
   ndx=n1d+ncx+1
   ndxmax=2*ndx-1
   ndy=n1d+ncy+1
   ndymax=2*ndy-1
   ndz=n1d+ncz+1
   ndzmax=2*ndz-1
!
!---------read dimensions/grid-parameters for frequency grid------------
!
   rewind 1
   read(1, nml=dimensions_freq)
!
!symmetric profiles
   xcmf_min=-xcmf_max
!
!-----------read dimensions/grid-parameters for angle grid--------------
!
   rewind 1
   read(1, nml=dimensions_angles)
!
!--------------------read information region----------------------------
!
   rewind 1
   read(1, nml=input_infreg)
!
!--------read model parameters for modelling continuum and line---------
!
   rewind 1
   read(1, nml=input_cont)
!
   rewind 1
   read(1, nml=input_line)
!
!-----------------------read benchmark parameters-----------------------
!
   rewind 1
   read(1, nml=benchmark)
   nn_y=n_y
   nn_z=n_z
!
close(1)
!
select case(opt_ltec)
   case(0)
      lfreqint = .false.
   case(1)
      lfreqint = .true.
   case default
      stop 'error in read_input: opt_ltec not properly specified'
end select
!
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine allocate_global3d
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, rho3d, opac3d, opalbar3d, t3d, &
                  velx3d, vely3d, velz3d, vth3d, imask_totreg3d, int3d, &
                  sline3d, ssobo3d, scont3d, imask_innreg3d, imask_bpoint3d, imask3d, &
                  mint3d, mintbar3d, normalization3d, eps_cont3d, &
                  alocont_nn3d, alocont_o_nn3d, aloline_nn3d, aloline_on_nn3d, fcontx3d, fconty3d, fcontz3d, &
                  kcontxx3d, kcontyy3d, kcontzz3d, kcontxy3d, kcontxz3d, kcontyz3d
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: err
!
!
!
allocate(t3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(eps_cont3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(rho3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(opac3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(opalbar3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(vth3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(velx3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(vely3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(velz3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(imask3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(imask_totreg3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(imask_innreg3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(imask_bpoint3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
!
!set everything to negative values or to NaN to produce error if not properly set during code execution
t3d = -1.d0
eps_cont3d=-1.d0
rho3d = -1.d0
opac3d = -1.d0
opalbar3d = -1.d0
vth3d = -1.d0
!velx3d = 1.d0/0.d0
!vely3d = 1.d0/0.d0
!velz3d = 1.d0/0.d0
imask3d = -1
imask_totreg3d = -1
imask_innreg3d = -1
imask_bpoint3d = -1
!
allocate(int3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(mint3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(fcontx3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(fconty3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(fcontz3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontxx3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontyy3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontzz3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontxy3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontxz3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(kcontyz3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(mintbar3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(sline3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(ssobo3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(scont3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
!
allocate(alocont_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(alocont_o_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(aloline_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
allocate(aloline_on_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
!
allocate(normalization3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'error in allocate_global3d: allocation'
!
end subroutine allocate_global3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***********************************************************************
!********************MATHEMATICAL ROUTINES******************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine dime_xobs_nodes
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, velx3d, vely3d, velz3d, vth3d, imask_totreg3d
use freq, only: deltax, xmax, xcmf_max, nxobs, nodes_xobs, nodes_xnue, weight_xobs, xic_nue
use params_input, only: vth_fiducial
!!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: vmax, vth_min, vel
!
! ... local functions
!
!-----------------------------------------------------------------------
!
write(*,*) '--------------------calculating number frequency points------------------------'
write(*,*)
!
!------find minimum thermal velocity and maximum absolute velocities----
!
vmax=0.d0
vth_min=1.d8
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(imask_totreg3d(i,j,k).ne.0) then
            vel=sqrt(velx3d(i,j,k)**2 + vely3d(i,j,k)**2 + velz3d(i,j,k)**2)
            if(vel.gt.vmax) vmax=vel
            if(vth3d(i,j,k).lt.vth_min) then
               vth_min=vth3d(i,j,k)
!               write(*,*) x(i), y(j), z(k), vth_min, vth3d(i,j,k), sqrt(x(i)**2 + y(j)**2 + z(k)**2)
            endif
         endif
      enddo
   enddo
enddo

!write(*,*) vth_min/1.d5, vmax, vth_fiducial/vth_min
!stop 'go on here'
!
!vmax in vth_min
vmax = vmax*vth_fiducial/vth_min

!write(*,*) vmax, vth_fiducial ,vth_min
!
!-----------------------------------------------------------------------
!
!if(vth_fiducial.gt.vth_min+1.d-7) then
!   write(*,*) 'error in dime_nodes_xobs: minimum thermal velocity gt fiducial thermal velocity'
!   write(*,*) 'choose vth_fiducial=', vth_min/1.d5
!   stop
!endif
!
!-----------------------------------------------------------------------
!
!calculate number of required frequency points in vth_min
xmax = xcmf_max+vmax
!
nxobs=2*ceiling(xmax/deltax) + 1

!write(*,*) xmax, xcmf_max, deltax, vmax
!stop
!
!transform xmax to vth_fiducial
xmax = xmax*vth_min/vth_fiducial

!write(*,*) nxobs, xmax, vth_min, vth_fiducial
!stop
!
write(*,*) 'number of frequency points calculated for a range v_max/vth_fid', -xmax, xmax
write(*,*)
!
!allocate arrays
allocate(nodes_xobs(nxobs), stat=err)
   if(err.ne.0) stop 'error in dime_xobs_nodes: allocation'
allocate(weight_xobs(nxobs), stat=err)
   if(err.ne.0) stop 'error in dime_xobs_nodes: allocation'
allocate(nodes_xnue(nxobs), stat=err)
   if(err.ne.0) stop 'error in dime_xobs_nodes: allocation'
allocate(xic_nue(nxobs), stat=err)
   if(err.ne.0) stop 'error in dime_xobs_nodes: allocation'
!
end subroutine dime_xobs_nodes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_xobs
!
!------------calculates frequency integration nodes as------------------
!-------shift from line center in fiducial doppler widths---------------
!
use prog_type
use fund_const, ONLY: cgs_clight
use freq, only: xnue0, nodes_xobs, nodes_xnue, weight_xobs, xic_nue, nxobs, xmax
use params_input, only: vmax, teff, tmin, na, vth_fiducial, vmicro
use mod_integ1d, only: precalc_weight_trapez
!!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
real(dp) :: deldop_fiducial
!
! ... local functions
real(dp) :: deldop
!
!
write(*,*) '------------------------calculating frequency grid-----------------------------'
write(*,*)
!
deldop_fiducial = deldop(xnue0, vth_fiducial)
!
!calculate frequency grids

do i=1, nxobs
   nodes_xobs(i) = -xmax + (i-1)*2.d0*xmax/float(nxobs-1)
   nodes_xnue(i) = deldop_fiducial * nodes_xobs(i) + xnue0
enddo
!
!calculate integration weights
!
call precalc_weight_trapez(nodes_xobs, nxobs, weight_xobs)
!
!write(*,*) nodes_xobs
!stop 'go on in calcnodes_xobs'
!
end subroutine calcnodes_xobs
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_xobs
!
use prog_type
use freq, only: nodes_xobs, weight_xobs, nxobs, deltax
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp), parameter :: tol=1.d-2
real(dp) :: del, rel, val1, val2
!
! ... local functions
!
!
!check if delta-x steps are okay
do i=2, nxobs
   del=abs(nodes_xobs(i)-nodes_xobs(i-1))
   rel=(del-deltax)/deltax
   if(rel.gt.tol) then
      write(*,*) del, deltax, rel
      stop 'error in check_xobs: wrong delta-x steps'
   endif
enddo
!
!check if xobs-grid is symmetric
do i=1, nxobs
   val1=abs(nodes_xobs(nxobs+1-i))
   val2=abs(nodes_xobs(i))
   if(abs(val2-val1).gt.1.d-13) then
      stop 'error in check_xobs: x-obs grid not symmetric'
   endif
enddo
!
!
!
end subroutine check_xobs
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_units
!
use prog_type
use inf_reg, only: rmin
use params_stellar, only: taux, tauy, tauz, sr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, opac3d
!
implicit none
!
integer(i4b) :: i, j, k
real(dp) :: dtau
!
!
write(*,*) '-----------------------------checking units------------------------------------'
write(*,*)
!
!continuum optical depth along x-axis
j=ndymax/2+1
k=ndzmax/2+1
taux=0.d0
do i=ndxmax-1, ndxmax/2+1, -1
   if(x(i).lt.rmin) exit
   dtau=0.5d0*(opac3d(i,j,k)+opac3d(i+1,j,k))*(x(i+1)-x(i))
   taux=taux+dtau
enddo
!
!continuum optical depth along y-axis
i=ndxmax/2+1
k=ndzmax/2+1
tauy=0.d0
do j=ndymax-1, ndymax/2+1, -1
   if(y(j).lt.rmin) exit
   dtau=0.5d0*(opac3d(i,j,k)+opac3d(i,j+1,k))*(y(j+1)-y(j))
   tauy=tauy+dtau
enddo
!
!continuum optiocal depth along z-axis
i=ndxmax/2+1
j=ndymax/2+1
tauz=0.d0
do k=ndzmax-1, ndzmax/2+1, -1
   if(z(k).lt.rmin) exit
   dtau=0.5d0*(opac3d(i,j,k)+opac3d(i,j,k+1))*(z(k+1)-z(k))
   tauz=tauz+dtau
enddo
!
write(*,'(a40, 3es20.8)') 'continuum optical depth along x, y, z', taux, tauy, tauz
write(*,*) 
!
end subroutine check_units
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_grid3d
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, t3d, eps_cont3d, rho3d, opac3d, opalbar3d, vth3d, &
                  velx3d, vely3d, velz3d, imask3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d
use params_input, only: vth_fiducial
use freq, only: deltax, nodes_xobs, nxobs
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: delv, del, rel, val1, val2
real(dp), parameter :: tol=1.d-2
!
! ... local functions
!
!
!check if all nans or negative flags have been overwritten
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(t3d(i,j,k).lt.0.) stop 'error in check_grid3d: t3d < 0'
         if(eps_cont3d(i,j,k).lt.0.) stop 'error in check_grid3d: eps_cont3d < 0'
         if(rho3d(i,j,k).lt.0.) stop 'error in check_grid3d: rho3d < 0'
         if(opac3d(i,j,k).lt.0.) stop 'error in check_grid3d: opac3d < 0'
         if(opalbar3d(i,j,k).lt.0.) stop 'error in check_grid3d: opalbar3d < 0'
         if(vth3d(i,j,k).lt.0.) stop 'error in check_grid3d: vth3d < 0'
         if(imask3d(i,j,k).lt.0) stop 'error in check_grid3d: imask3d < 0'
         if(imask_totreg3d(i,j,k).lt.0) stop 'error in check_grid3d: imask_totreg3d < 0'
         if(imask_innreg3d(i,j,k).lt.0) stop 'error in check_grid3d: imask_innreg3d < 0'
         if(imask_bpoint3d(i,j,k).lt.0) stop 'error in check_grid3d: imask_bpoint3d < 0'
         if(t3d(i,j,k).ne.t3d(i,j,k)) then
!            write(*,*) i, j, k, x(i), y(j), z(k), sqrt(x(i)**2+y(j)**2+z(k)**2), t3d(i,j,k)
            stop 'error in check_grid3d: t3d is NaN'
         endif
         
         if(eps_cont3d(i,j,k).ne.eps_cont3d(i,j,k)) stop 'error in check_grid3d: eps_cont3d is NaN'
         if(rho3d(i,j,k).ne.rho3d(i,j,k)) stop 'error in check_grid3d: rho3d is NaN'
         if(opac3d(i,j,k).ne.opac3d(i,j,k)) stop 'error in check_grid3d: opac3d is NaN'
         if(opalbar3d(i,j,k).ne.opalbar3d(i,j,k)) stop 'error in check_grid3d: opalbar3d is NaN'
         if(vth3d(i,j,k).ne.vth3d(i,j,k)) stop 'error in check_grid3d: vth3d is NaN'
         if(velx3d(i,j,k).ne.velx3d(i,j,k)) stop 'error in check_grid3d: velx3d is NaN'
         if(vely3d(i,j,k).ne.vely3d(i,j,k)) stop 'error in check_grid3d: vely3d is NaN'         
         if(velz3d(i,j,k).ne.velz3d(i,j,k)) stop 'error in check_grid3d: velz3d is NaN'         
      enddo
   enddo
enddo
!
!
!
!check if delta-x steps are okay
do i=2, nxobs
   del=abs(nodes_xobs(i)-nodes_xobs(i-1))
   rel=(del-deltax)/deltax
   if(rel.gt.tol) then
      write(*,*) del, deltax, rel
      stop 'error in check_xobs: wrong delta-x steps'
   endif
enddo
!
!check if xobs-grid is symmetric
do i=1, nxobs
   val1=abs(nodes_xobs(nxobs+1-i))
   val2=abs(nodes_xobs(i))
   if(abs(val2-val1).gt.1.d-13) then
      stop 'error in check_xobs: x-obs grid not symmetric'
   endif
enddo
!
!
!
end subroutine check_grid3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_bnue
!
!------calculates planck-functions for all frequency grid values--------
!---------------------check if bnue approx const------------------------
!
use prog_type
use fund_const, ONLY: cgs_clight
use freq, only: nxobs, nodes_xobs, xnue0, lfreqint
use params_input, only: teff, vth_fiducial
use mod_directories, only: output_dir_test
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: deldop_fiducial, nue

! ... local arrays
real(dp), dimension(nxobs) :: planck_fct
!
! ... local functions
real(dp) :: bnue2
!
!doppler width
deldop_fiducial = xnue0 * vth_fiducial/cgs_clight
!
!
!calculate planck function for teff and at each frequency
open(1, file=trim(output_dir_test)//'planck_fct.dat', form='formatted')
   do i=1, nxobs
      nue=deldop_fiducial * nodes_xobs(i) + xnue0
      planck_fct(i) = bnue2(xnue0, teff, lfreqint)
      write(1,*) nue, planck_fct(i)
   enddo
close(1)
!
!
end subroutine check_bnue
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_options
!
!-------------checks whether all options are compatible-----------------
!
USE prog_type
use mod_directories, only: model_dir, output_dir, output_dir_temp, output_dir_test
use options, only: spatial_grid3d, spatial_grid1d, opt_angint_method, &
                   input_mod, input_mod_dim, opt_ng_cont, opt_ait_cont, opt_ng_line, opt_ait_line, opt_sol2d, &
                   opt_incl_gdark, opt_incl_sdist, opt_method, opt_ltec, opt_opac, opt_incl_cont, opt_incl_line
use mod_benchmark, only: benchmark_mod
!
implicit none
!
! ... arguments
!
! ... local characters
!
! ... local logicals
logical :: check_dir1, check_dir2, check_dir3, check_dir4
!
!--------------------check if directories exist-------------------------
!
!with gfortran
inquire(file=trim(output_dir)//'/.', exist=check_dir1)
inquire(file=trim(output_dir_temp)//'/.', exist=check_dir2)
inquire(file=trim(output_dir_test)//'/.', exist=check_dir3)
inquire(file=trim(model_dir)//'/.', exist=check_dir4)
!
!with ifort
!inquire(directory=trim(output_dir), exist=check_dir1)
!inquire(directory=trim(output_dir_temp), exist=check_dir2)
!inquire(directory=trim(output_dir_test), exist=check_dir3)
!inquire(directory=trim(model_dir), exist=check_dir4)
!
if(.not.check_dir1) then
   write(*,*) 'error in check_options: directory does not exist'
   write(*,*) output_dir
   stop
endif
if(.not.check_dir2) then
   write(*,*) 'error in check_options: directory does not exist'
   write(*,*) output_dir_temp
   stop
endif
if(.not.check_dir3) then
   write(*,*) 'error in check_options: directory does not exist'
   write(*,*) output_dir_test
   stop
endif
if(.not.check_dir4) then
   write(*,*) 'error in check_options: directory does not exist'
   write(*,*) model_dir
   stop
endif
!
write(*,*) 'TODO: subroutine check_options to be extended'
write(*,*)
!
!---IF LINE-TRANSPORT IS CALCULATED FROM A CERTAIN ITERATION STEP-------
!--------CONTINUUM TRANSPORT NEEDS TO BE READ IN AS WELL----------------
!-------------SINCE IT WAS PREVIOUSLY CALCULATED------------------------
!
!if(.not.opt_start_line.and.opt_start_cont) then
!   write(*,*) 'opt_start_cont=true not allowed if opt_start_line=false'
!   write(*,*) '   => set opt_start_cont to false'
!   stop
!endif
!
!-----EITHER NG-EXTRAPOLATION OR AITKEN-EXTRAPOLATION IS TO BE USED-----
!
if(opt_ng_cont.and.opt_ait_cont) then
   write(*,*) 'either ng-extrapol or aitken-extrapol or none of both is to be used'
   write(*,*) '   => set opt_ng_cont and/or opt_ait_cont to false'
   stop
endif
!
if(opt_ng_line.and.opt_ait_line) then
   write(*,*) 'either ng-extrapol or aitken-extrapol or none of both is to be used'
   write(*,*) '   => set opt_ng_line and/or opt_ait_line to false'
   stop
endif
!
!-----------------------------------------------------------------------
!
select case(opt_method)
   case(-1) !no RT, only initial conditions
   case(0) !fvm
   case(1) !linear sc
   case(2) !quadratic bezier
   case default
      stop 'error in check_options: opt_method not properly specified'
end select
!
!-----------------------------------------------------------------------
!
select case(opt_angint_method)
   case(7,8,9)
      if(opt_sol2d) then 
         write(*,*) 'error check_options: triangular solid angle integration'
         write(*,*) '   not allowed in 2d solution method'
         stop
      endif
   case default
end select
!
!---check if 3d-grid structure options are compatible with input-model--

if(input_mod_dim.ne.1.and.input_mod_dim.ne.2) then
   if(spatial_grid3d.eq.2) then
      write(*,*) 'spatial_grid3d=2 and input_mod ne 1 or 2 does not make sense:'
      write(*,*) '   3d-grid from dylan kees model atmosphere, although that model'
      write(*,*) '   shall not be calculated'
      stop
   endif
endif 
!
!-----------------------------------------------------------------------


!
!------------check if angular integration method is allowed-------------
!---------------when 2d benchmarks are performed------------------------
!
select case(opt_angint_method)
   case(7,8,9)
      select case(benchmark_mod)
         case(1,2,3,4,5,8,9)
            write(*,*) 'error check_options: triangular solid angle integration'
            write(*,*) '   not allowed in 2d solution method'
            stop
          case default
      end select
   case default
end select
!
!--------gravity darkening and surface distortion not implemented-------
!----------------------for 3d fvm method yet----------------------------
!
if(opt_method.eq.0) then
   if(opt_incl_gdark) stop 'error check_options: gravity darkening not implemented for fvm method yet'
   if(opt_incl_sdist) stop 'error check_options: surface distortion not implemented for fvm method yet'
endif
!
!--------opt_ltec only makes sense if continuum is included-------------
!
if(.not.opt_incl_cont) then
   if(opt_ltec.ge.1) stop 'error check_options: opt_ltec makes no sense if continnuum is not included'
endif
!
if(opt_incl_line) then
   if(opt_ltec.ge.1) stop 'error check_options: opt_ltec together with line transport needs to be implemented'
endif
!
!
end subroutine check_options
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_bcondition
!
!calculation of boundary condition: diffusion approximation with
!
!   xic1=bnue(xnue, trad)
!   xic2=dB/dtau(xnue, teff)
!or grey
!
!   xic1 = sigmab/pi * trad^4
!   xic2 = 0.d0   (temperature stratification will be calculated)
!
use prog_type
use freq, only: xnue0, lfreqint
use params_input, only: trad, vrot, lstar, xlogg, rstar
use params_stellar, only: sr
use bcondition, only: xic1, xic2, corrfc, ntheta_gdark, theta_gdark, xic1_gdark, teff_gdark
use fund_const, only: pi, rsu, xmsu, xlsu, cgs_sb, cgs_grav
use options, only: opt_incl_gdark, opt_ltec
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b), parameter :: ntheta_fine=101
integer(i4b) :: i, k, iim2, iim1, ii, iip1
real(dp) :: c1, c2
real(dp) :: t_iim1, t_ii, t_iip1, r_iim1, r_ii, r_iip1
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, l_star_cgs, w0, omega, sigma, sigma_fit
!
! ... for testing
!real(dp) :: xp, yp, zp, rad, theta_i, indx, calc_icore_gdark
!
! ... local arrays
real(dp), dimension(ntheta_fine) :: theta_fine, rsurface_fine, gr_fine, gtheta_fine, gperp_fine, integrand_fine, teff_fine
!
! ... local functions
real(dp) :: bnue2
!
!for now: calculate xic2=db/dr=0 through dummy values of r, and constant temperature
t_iim1=trad
t_ii=trad
t_iip1=trad
r_iim1=0.99d0
r_ii=1.d0
r_iip1=1.01d0
!
select case(opt_ltec)
   case(0) 
      call diffus(xnue0, t_iim1, t_ii, t_iip1, r_iim1, r_ii, r_iip1, xic1, xic2)
!corrfc from jo's 1d program
      corrfc=0.883513122620739d0
      corrfc=0.959147197062505d0
      corrfc=0.996799721102628d0
      corrfc=1.d0
      if(corrfc.eq.1.d0) xic2=0.d0
   case(1)
      xic1 = bnue2(xnue0, trad, lfreqint)
      xic2 = 0.d0
   case default
      stop 'error in calc_bcondition: opt_ltec not properly specified'
end select
!
!
!!calculate photospheric profile (set to constant value)
!call calc_photprof(nxobs, xic_nue, xic1)
!!read photospheric profile on external grid and interpolate
!call get_photprof(nxobs, nodes_xnue, xic_nue, xnue0, fnamephotprof)
!
!----------------------calculate gravity darkening law------------------
!-------------following petrenz/puls 1996, cranmer/owocki 1995----------
!---------------------(neglecting eddington factor!!!)------------------
!
!create theta-grid
!c1 = -pi/2.d0/0.9d0
!c2 = pi/2.d0/0.9d0
do i=1, ntheta_fine
   theta_fine(i) = 1.d-5 + (i-1)*(pi/2.d0-1.d-5)/(ntheta_fine-1)
!   theta_fine(i) = c1*10.d0**(-float(i-1)/float(ntheta_fine-1)) + c2
enddo
!theta_fine(1) = 1.d-5
!
!
!
r_pole_cgs = sr
m_star_cgs = sr**2 * 10.d0**xlogg/cgs_grav
v_rot_cgs = vrot*1.d5
l_star_cgs = lstar*xlsu
!
w0 = v_rot_cgs**2*r_pole_cgs/2.d0/cgs_grav/m_star_cgs
omega = sqrt(27.d0/4.d0*w0*(1.d0-w0)**2)
!
!
!
do i=1, ntheta_fine
!radius as function of theta in units of rstar
   if(abs(omega).lt.1.d-10) then
      rsurface_fine(i) = 1.d0
   else
      rsurface_fine(i)=3.d0/omega/sin(theta_fine(i))*cos((pi+acos(omega*sin(theta_fine(i))))/3.d0)
   endif
!
!gravity components
   gr_fine(i) = cgs_grav*m_star_cgs/r_pole_cgs**2 * (-1.d0/rsurface_fine(i)**2 + 8.d0/27.d0*rsurface_fine(i)*omega**2*sin(theta_fine(i))**2)
   gtheta_fine(i) = cgs_grav*m_star_cgs/r_pole_cgs**2 * 8.d0/27.d0*rsurface_fine(i)*omega**2*sin(theta_fine(i))*cos(theta_fine(i))
   gperp_fine(i) = sqrt(gr_fine(i)**2+gtheta_fine(i)**2)
!
!integrand to determine sigma
   integrand_fine(i) = 4.d0*pi*gperp_fine(i)*(rsurface_fine(i)*r_pole_cgs)**2*sin(theta_fine(i)) / (-gr_fine(i)/gperp_fine(i))
!   write(*,'(10es20.8)') theta_fine(i), rsurface_fine(i), gr_fine(i), gtheta_fine(i), gperp_fine(i), integrand_fine(i)
enddo
!
!calculate sigma
sigma=0.d0
do i=2, ntheta_fine
   sigma = sigma + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
enddo
sigma_fit = 4.d0*pi*cgs_grav*m_star_cgs*(1.d0-0.1969*omega**2 - 0.094292*omega**4 + &
            0.33812*omega**6 - 1.3066*omega**8 + 1.8286*omega**10 - 0.92714*omega**12)
!
!calculate effective temperature
teff_fine = (l_star_cgs*gperp_fine/cgs_sb/sigma)**0.25d0
!!
!write(*,*) 'sigma', sigma, sigma_fit
!write(*,*) 'r_pole', r_pole_cgs/rsu
!write(*,*) 'm_star', m_star_cgs/xmsu
!write(*,*) 'v_rot', v_rot_cgs/1.d5
!write(*,*) 'l_star', l_star_cgs/xlsu
!write(*,*) 'w0', w0
!write(*,*) 'omega', omega
!do i=1, ntheta_fine
!   write(*,'(10es20.8)') theta_fine(i)*180.d0/pi, rsurface_fine(i), log10(gperp_fine(i)), integrand_fine(i), teff_fine(i)
!enddo
!stop
!
!create final theta grid
!note: if theta-grid is changed, need also to change index-search in
!      subroutine calc_icore_gdark!!!!
c1 = -pi/2.d0/0.9d0
c2 = pi/2.d0/0.9d0
do i=1, ntheta_gdark
   theta_gdark(i) = c1*10.d0**(-float(i-1)/float(ntheta_gdark-1)) + c2
   call find_index(theta_gdark(i), theta_fine, ntheta_fine, iim2, iim1, ii, iip1)
   teff_gdark(i) = interpol_yp(theta_fine(iim1), theta_fine(ii), teff_fine(iim1), teff_fine(ii), theta_gdark(i))
!
!if grey atmosphere is calculated
   xic1_gdark(i) = bnue2(xnue0, teff_gdark(i), lfreqint)
!
enddo
!
!if gravity darkening is not included, overwrite arrays
if(opt_incl_gdark) then
   trad = teff_gdark(1)
   xic1 = xic1_gdark(1)
else
   teff_gdark = trad
   xic1_gdark = xic1
endif

!open(1, file='TRASH/gdark.dat', form='formatted')
!   do i=1, ntheta_gdark
!      write(1,'(10es20.8)') theta_gdark(i)*180.d0/pi, teff_gdark(i), xic1_gdark(i)
!   enddo
!close(1)
!
!xp=1.d0
!yp=1.d0
!do i=1, 50
!   zp = -3.d0 + (i-1)*6.d0/49.d0
!   rad = sqrt(xp**2+yp**2+zp**2)
!   theta_i = acos(abs(zp)/rad)
!   indx = log10(1.d0-1.8d0*theta_i/pi)*(1.d0-ntheta_gdark)+1.d0
!   write(*,'(4es20.8,2i5,6es20.8)') rad, zp, theta_i*180.d0/pi, indx, floor(indx), ceiling(indx), theta_gdark(floor(indx)), theta_i, theta_gdark(ceiling(indx)), &
!                                     xic1_gdark(floor(indx)), calc_icore_gdark(zp, rad), xic1_gdark(ceiling(indx))
!enddo
!stop 'go on in calc_bcondition'
!
end subroutine calc_bcondition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_bcondition_phot
!
!setting S_cont=B_nue(T) and S_line=B_nue(T) inside the photosphere
!    or  S_cont=B(T)     and S_line=B(T)     inside the photosphere (if grey atmosphere calculated)
!
use prog_type
use freq, only: xnue0, lfreqint
use dime3d, only: scont3d, sline3d, imask3d, ndxmax, ndymax, ndzmax
use params_input, only: trad
use bcondition, only: xic1
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: planck
!
! ... local functions
real(dp) :: bnue2
!
scont3d=0.d0
planck=bnue2(xnue0, 1.*trad, lfreqint)
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(imask3d(i,j,k).eq.4) then
            scont3d(i,j,k) = xic1
            sline3d(i,j,k) = xic1
         endif
      enddo
   enddo
enddo

!where statement is not working for large arrays
!where(imask3d.eq.4)
!   write(*,*) 'test'
!   scont3d = xic1  !planck
!   sline3d = xic1  !planck
!endwhere
!
!
!
end subroutine calc_bcondition_phot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_bintensity
!
!calculates the intensities of inner boundary points as a function
!of angle, and store them for each direction
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask3d, int3d
use angles, only: dim_omega, nodes_mu, n_x, n_y, n_z
use bcondition, only: xic1, n_inner, indx_xinner, indx_yinner, indx_zinner, int_inner
use options, only: opt_sol2d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, indx_omega, indx_ipoint, err
real(dp) :: nn_x, nn_y, nn_z, mueff
!
! ... local arrays
real(dp), dimension(:,:,:), allocatable :: int3d_test
!
! ... local characters
character(len=32) :: fname
!
write(*,*) '-------------calculating (angle dependent) inner boundary condition------------'
!
!setting arrays to describe the inner boundary condition
n_inner=0
do i=1, ndxmax 
   do j=1, ndymax
      do k=1, ndzmax
         if(imask3d(i,j,k).eq.4) n_inner=n_inner+1
      enddo
   enddo
enddo
!
write(*,*) 'number of inner points', n_inner
allocate(indx_xinner(n_inner), stat=err)
   if(err.ne.0) stop 'allocation error in calc_bintensity'
allocate(indx_yinner(n_inner), stat=err)
   if(err.ne.0) stop 'allocation error in calc_bintensity'
allocate(indx_zinner(n_inner), stat=err)
   if(err.ne.0) stop 'allocation error in calc_bintensity'
allocate(int_inner(n_inner), stat=err)
   if(err.ne.0) stop 'allocation error in calc_bintensity'
!
!setting boundary condition for inner points
do indx_omega=1, dim_omega
   nn_x=n_x(indx_omega)
   nn_y=n_y(indx_omega)
   nn_z=n_z(indx_omega)
!
   indx_ipoint=1
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            if(imask3d(i,j,k).eq.4) then
!inside the star, set intensities correctly (for spherical star)
               indx_xinner(indx_ipoint) = i
               indx_yinner(indx_ipoint) = j
               indx_zinner(indx_ipoint) = k
               mueff=(nn_x*x(i)+nn_y*y(j)+nn_z*z(k))/sqrt(x(i)**2+y(j)**2+z(k)**2)
               if(mueff.ge.0.d0) then
                  int_inner(indx_ipoint) = xic1
               else
                  int_inner(indx_ipoint) = 0.d0
               endif
               indx_ipoint=indx_ipoint+1
            endif
         enddo
      enddo
   enddo
!
   write(fname,'(a19,i4.4)') 'data_bcondition/iomega', indx_omega
   open(1, file=fname, form='unformatted')
      write(1) int_inner
      write(1) indx_xinner
      write(1) indx_yinner
      write(1) indx_zinner
   close(1)
!
enddo
!
!
end subroutine calc_bintensity
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine make_benchmark
!
use prog_type
use mod_benchmark, only: benchmark_mod
use options, only: opt_method
!
implicit none
!
! ... local scalars
integer(i4b) :: i
!
!calculate delx, dely, delz grids if not done yet
if(opt_method.ne.0) call grid_delxyz
!
select case(benchmark_mod)
   case(1)
      write(*,*) '------------------------performing benchmark model 1---------------------------'
      write(*,*)
!creating grids
      call benchmark01_grid
!calculating 1d solution
      call benchmark01_solution
!calculating 2d sc solution 9in formal_sc2db.f90, which needs to be included in makefile
!      call fsc_cont2d_debug(1)
!calculating 2d fvm solution
      call ffvm_cont2d_debug(1)
!interpolate 2d solution onto ray
      call benchmark01_ray
!output to file
      call output_benchmark01
!
   case(2)
      write(*,*) '-----------------performing benchmark model 2: searchlight 2d------------------'
      write(*,*)
!calculating solution
      call benchmark02_solution
!output to file
      call output_benchmark02
!
   case(3)
      write(*,*) '-----------------performing benchmark model 3: searchlight 2d------------------'
      write(*,*)
!calculating solution
      call benchmark03_solution
!output to file
      call output_benchmark03
!
   case(4)
      write(*,*) '----------performing benchmark model 4: 2d optically thin continuum------------'
      write(*,*)
!calculating solution
      call benchmark04_solution
!output to file
      call output_benchmark04
!
   case(5)
      write(*,*) '----performing benchmark model 5: 2d spherically symmetric continuum tests-----'
      write(*,*)
!calculating solution
      call benchmark05_solution
!output to file
      call output_benchmark05
!
   case(6)
      write(*,*) '------------performing benchmark model 6: 1d line transport test---------------'
      write(*,*)
!calculating solution
      call benchmark06_solution
!output to file
      call output_benchmark06
!
   case(7)
      write(*,*) '----performing benchmark model 7: spherically symmetric sobolev in 3d and 1d---'
      write(*,*)
!calculating solution
      call benchmark07_solution
!output to file
      call output_benchmark07
!
   case(8)
      write(*,*) '-------performing benchmark model 8: spherically symmetric line tests----------'
      write(*,*)
!calculating solution
      call benchmark08_solution
!output to file
      call output_benchmark08
!
   case(9)
      write(*,*) '---------performing benchmark model 9: timing for 2d line calculation----------'
      write(*,*)
!calculating solution
      call benchmark09
!output to file
!     call output_benchmark09
!
   case(10)
      write(*,*) '------------performing benchmark model 10: 3d searchlight beam test------------'
      write(*,*)
!calculating solution
      call benchmark10_solution
!output to file
      call output_benchmark10
!
   case(11)
      write(*,*) '----------performing benchmark model 11: 3d optically thin continuum-----------'
      write(*,*)
!calculating solution
      call benchmark11_solution
!output to file
      call output_benchmark11
!
   case(12)
      write(*,*) '----performing benchmark model 12: 3d spherically symmetric continuum tests----'
      write(*,*)
!calculating solution
!      call calc_startval_cont   !reasonable starting value required for nice convergence behaviour
      call benchmark12_solution
!output to file
      call output_benchmark12
!
   case(13)
      write(*,*) '---performing benchmark model 13: 3d spherically symmetric line transfer test--'
      write(*,*)
!calculating solution
      call benchmark13_solution
!output to file
      call output_benchmark13
!
   case(14)
      write(*,*) '---performing benchmark model 14: 3d sph. symm. line+continuum transfer test---'
      write(*,*)
!calculating solution
      call benchmark14_solution
!output to file
      call output_benchmark14

!
   case default
      write(*,*) '----------------------no benchmark is being performed--------------------------'
      write(*,*)
      return
end select
!
!
!
end subroutine make_benchmark

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine linetrans2d
!
!---------------------------line transport------------------------------
!
!---calculating intensities for different angles mu and frequencies-----
!--------------in x-z plane from short characteristics------------------
!-------performing angle and frequency integration on z-axis------------
!--back interpolation of newly calculated source function onto 3d grid--
!
use prog_type
use iter, only: itmaxl, devmaxl, epsmaxl_arr
use dimecr, only: n1d_cr, sline1d_cr, mintbar1d_cr, r1d_cr, aloline1d_diag_cr, &
                  norm1d_cr, epshistoryl1d_cr, epsl1d_cr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, sline3d, t3d, imask_innreg3d
use params_input, only: eps_line
use options, only: opt_ait_line, opt_ng_line, opt_method
use ng_extra, only: ng_const
use freq, only: xnue0
!use timing, only: ts_it, te_it, ttot_it, it_tot
use warnings, only: warn_itmaxl
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k, l, err
integer(i4b) :: s1, s2, s3, s4
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: dummy1, dummy2, rad
real(dp) :: eps_max
!
! ... local arrays
real(dp), dimension(:), allocatable :: bnue1d_cr
real(dp), dimension(:,:), allocatable :: sline1d_ng
!
!... local characters
!character(len=50) :: enter
!
! ... local functions
real(dp) :: bnue
!
! ... local logicals
!
!-----------------------------------------------------------------------
!
write(*,*) '-----------------------------line transport (cr)-------------------------------'
write(*,*)
!
allocate(bnue1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in linetrans2d'
allocate(sline1d_ng(4,n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in linetrans2d'
!
!-----------------------------------------------------------------------
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
!
!-----------------------------------------------------------------------
!
!calculating start values of line source function
call calc_startval_line
!
!calculating planck function along z direction
do i=1, n1d_cr
   bnue1d_cr(i) = bnue(xnue0, t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i))
   sline1d_cr(i) = sline3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
enddo
!
!setting start mean intensities and start deviations to zero
epsl1d_cr=0.d0
epshistoryl1d_cr=0.d0
mintbar1d_cr=0.d0
aloline1d_diag_cr=0.d0
!
!-----------------------------------------------------------------------
!
!ttot_it=0.d0
!
!************************start iteration scheme*************************
!
do i=1, itmaxl
!
!timing
!it_tot=i
!call cpu_time(ts_it)
!
!****************************output*************************************
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j_bar', 's_l', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mintbar1d_cr(j), sline1d_cr(j), epsl1d_cr(j), norm1d_cr(j), aloline1d_diag_cr(j)
   end do
   write(*,*)
!
!*************************output end************************************
!
   epsl1d_cr=mintbar1d_cr
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
   !
   select case(opt_method)
      case(-1)
         write(*,*) 'method: returning after initial condition'
         return
      case(0)
         call mintbar_fvm2d
      case(1)
         call mintbar_sc2d
      case(2)
         call mintbar_sc2d
      case default
         stop 'error in linetrans2d: opt_method not properly specified'
   end select      
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsl1d_cr, mintbar1d_cr, n1d_cr, eps_max)
   epshistoryl1d_cr(:,i)=epsl1d_cr
   epsmaxl_arr(i)=eps_max
!
   if(abs(eps_max).lt.devmaxl) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
!      call cpu_time(te_it)
!      ttot_it = ttot_it + te_it - ts_it
      exit
   else if(i.eq.itmaxl) then
      write(*,*) "no convergence after iteration no. ", i
      warn_itmaxl=.true.
   end if
!
!-------calculating alo-corrected source-functions on central ray-------
!
   write(*,*) '--------------calculating new iterate for source function (alo)----------------'
   write(*,*)
!
   dummy2=1.d0-eps_line
   do j=1, n1d_cr
      dummy1=1.d0-(1.d0-eps_line)*aloline1d_diag_cr(j)
      sline1d_cr(j)=(dummy2/dummy1)*mintbar1d_cr(j) - &
             (dummy2/dummy1)*aloline1d_diag_cr(j)*sline3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j) + &
              eps_line*bnue1d_cr(j)/dummy1
   enddo
!
!------------extrapolation of old subsequent source functions-----------
!
   if(opt_ng_line.or.opt_ait_line) then
!
!----------storing old source-functions for ng-extrapolation------------
!
      if(i.eq.s1) then
         write(*,*) '----------------------storing source fct. at step n-3--------------------------'
         write(*,*)
         sline1d_ng(1,:)=sline1d_cr
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         sline1d_ng(2,:)=sline1d_cr
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         sline1d_ng(3,:)=sline1d_cr
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         sline1d_ng(4,:)=sline1d_cr
         s4=s4+ng_const
         if(opt_ng_line) call ng_expol1d(sline1d_ng, n1d_cr)
         if(opt_ait_line) call ait_expol1d(sline1d_ng, n1d_cr)
         sline1d_cr=sline1d_ng(1,:)
      endif

   endif
!
!-----back-interpolation of central-ray-source-functions on 3d-grid-----
!
   do j=1, ndxmax
      do k=1, ndymax
         do l=1, ndzmax
            if(imask_innreg3d(j,k,l).eq.1) then
               sline3d(j,k,l) = 0.d0
            else
               rad=sqrt(x(j)**2 + y(k)**2 + z(l)**2)
               call find_index(rad, r1d_cr, n1d_cr, iim2, iim1, ii, iip1)
               sline3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), sline1d_cr(iim1), sline1d_cr(ii), rad)
            endif
         enddo
      enddo
   enddo
!
   call output_itstep_line(i)
!
!!timing
!   call cpu_time(te_it)
!   ttot_it=ttot_it + te_it - ts_it
!   write(*,*) 'time per iteration: ', te_it-ts_it
!
enddo
!
!
!
end subroutine linetrans2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine linetrans3d
!
!---------------------------line transport------------------------------
!
!-------calculating intensities for different angles theta, phi---------
!---------------------and different frequencies-------------------------
!-----------on a 3d cartesian grid using short characteristics----------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb, zero
use dimecr, only: n1d_cr!, r1d_cr
use dime3d, only: mintbar3d, sline3d, ssobo3d, normalization3d, aloline_nn3d, &
                  imask_totreg3d, imask3d, x, y, z, ndxmax, ndymax, ndzmax
use iter, only: itmaxl, devmaxl, epsmaxl_arr, it_start_line
use ng_extra, only: ng_const
use options, only: opt_ng_line, opt_ait_line, input_mod_dim, opt_opal, opt_method
use warnings, only: warn_itmaxl
use mod_interp2d, only: interpolation_threshold, wpa_interp2d, wpb_interp2d, wp_interp2d, &
                        wpa_interp1d, wpb_interp1d, wp_interp1d, wpa_integ1d, wpb_integ1d, wp_integ1d, &
                        lng_expol
use freq, only: xnue0
use bcondition, only: xic1
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, indx_threshold
integer(i4b) :: ix_epsmax, iy_epsmax, iz_epsmax
real(dp) :: s1, s2, s3, s4, s4b, eps_max, fdum, dummy1, dummy2, rad, trad
real(dp) :: ts, te
!
! ... local arrays
real(dp), dimension(:,:,:), allocatable :: eps3d
real(dp), dimension(:,:), allocatable :: sline3d_ng

!
! ... local logicals
logical :: check_fname
!
! ... local characters
character(len=3) :: fname_debug
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------line transport in 3d-------------------------------'
write(*,*)
!
call calc_startval_line
!
allocate(eps3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in linetrans3d'
allocate(sline3d_ng(4,ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in linetrans3d'
!
!
s1=it_start_line+1
s2=it_start_line+2
s3=it_start_line+3
s4=it_start_line+4
s4b=4
!
!-----------------------------------------------------------------------
!
!setting start mean intensities and start deviations to zero
eps3d=zero
normalization3d=zero
aloline_nn3d=zero
epsmaxl_arr=zero
!
!setting index for threshold
indx_threshold=it_start_line+5
!
!---------------------start iteration scheme----------------------------
!
do i=it_start_line+1, itmaxl
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 7(a20))') '#', 'z', 'jbar(rt)', 'sline(rt)', 'sline(sobo)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 7(e20.8))') j, z(j), mintbar3d(ndxmax/2+1,ndymax/2+1,j), sline3d(ndxmax/2+1,ndymax/2+1,j), &
                                    ssobo3d(ndxmax/2+1, ndymax/2+1, j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    aloline_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
   write(*,*)
   write(*,'(a25,es20.8)') 'alpha_min (interp1d)', wp_interp1d
   write(*,'(a25,es20.8)') 'alpha_min (interp2d)', wp_interp2d
   write(*,'(a25,es20.8)') 'alpha_min (integ1d)', wp_integ1d
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mintbar3d
!
   call output_itstep_line(i)
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '-------calculating frequency and angle integrated mean intensities in 3d-------'
   write(*,*) 'step', i
   write(*,*)
!   
   select case(opt_method)
      case(-1)
         write(*,*) 'method: returning after initial condition'
         return
      case(0)
         write(*,*) 'method: fvm'
         call mintbar_fvm3d
      case(1)
         write(*,*) 'method: sc-lin'
         call mintbar_sc3d_lin
      case(2)
         !use linear interpolations for the first s4b iteration steps
         !in order to get a first guess of a 'smooth' solution.
         !otherwise: convergence behaviour might be oscillatory
         if(i.le.s4b) then
            write(*,*) 'method: sc-lin (iteration < s4b)'
            call mintbar_sc3d_lin
         else
            write(*,*) 'method: sc-bez'
            call mintbar_sc3d
         endif
      case default
         stop 'error in linetrans3d: opt_method not properly specified'
   end select      
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)
!
   call calc_dev3d(eps3d, mintbar3d, imask_totreg3d, ndxmax, ndymax, ndzmax, eps_max, ix_epsmax, iy_epsmax, iz_epsmax)
   epsmaxl_arr(i)=eps_max
   write(*,'(a30, 3(i4), 3(f8.4), es18.8)') 'max (dev) at grid-point:', ix_epsmax, iy_epsmax, iz_epsmax, &
                                            x(ix_epsmax), y(iy_epsmax), z(iz_epsmax) , eps_max
   write(*,*)
!
   if(abs(eps_max).lt.devmaxl) then
      write(*,*) 'convergence after iteration no. ', i
      write(*,*) 'max (dev): ', eps_max
      write(*,*)
!for thin lines, ensure that higher order interpolation scheme
!   has been used, and not only linear approach
      if(i.gt.s4b) exit
   elseif(i.eq.itmaxl) then
      write(*,*) 'no convergence after iteration no. ', i
      write(*,*)
      warn_itmaxl=.true.
    endif
!
   if(i.gt.indx_threshold) then
      if(abs(epsmaxl_arr(i)).ge.abs(epsmaxl_arr(i-1)).and. &
         abs(epsmaxl_arr(i-1)).le.abs(epsmaxl_arr(i-2)).and. &
         abs(epsmaxl_arr(i-2)).ge.abs(epsmaxl_arr(i-3)).and. &
         abs(epsmaxl_arr(i-3)).le.abs(epsmaxl_arr(i-4)).and. &
         abs(epsmaxl_arr(i-4)).ge.abs(epsmaxl_arr(i-5))) then
!
         write(*,*) 'error in iteration scheme: oscillations!!!'
         write(*,*) 'max deviation at iteration i-5', epsmaxl_arr(i-5)
         write(*,*) 'max deviation at iteration i-4', epsmaxl_arr(i-4)
         write(*,*) 'max deviation at iteration i-3', epsmaxl_arr(i-3)
         write(*,*) 'max deviation at iteration i-2', epsmaxl_arr(i-2)
         write(*,*) 'max deviation at iteration i-1', epsmaxl_arr(i-1)
         write(*,*) 'max deviation at iteration i  ', epsmaxl_arr(i)
         write(*,*) 'possible solutions: '
         write(*,*) '   1. use linear interpolations for upwind/downwind source function'
         write(*,*) '      and for upwind intensities to have same solution procedure in'
         write(*,*) '      each iteration step (independent of the source function and intensities'
         write(*,*) '      themselves.'
         write(*,*) '   2. avoid monotonicity constraint in integration of source contribution'
         write(*,*) '      (try linear approximation of source function along ray)'
         write(*,*) '   3. increase the spatial grid resolution, in order to avoid'
         write(*,*) '      monotonicity constraints in quadratic interpolation procedures'
         write(*,*) '   4. increasing interpolation threshold in quadratic interpolation procedure'
         interpolation_threshold=min(interpolation_threshold+0.1d0,1.d0)
         write(*,*) 'setting interpolation threshold to', interpolation_threshold
         wpa_interp2d=wpa_interp2d+1.d0
         wpb_interp2d=wpb_interp2d+1.d0
         wp_interp2d=wpa_interp2d/wpb_interp2d
         wpa_interp1d=wpa_interp1d+1.d0
         wpb_interp1d=wpb_interp1d+1.d0
         wp_interp1d=wpa_interp1d/wpb_interp1d
         wpa_integ1d=wpa_integ1d+1.d0
         wpb_integ1d=wpb_integ1d+1.d0
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-1.d0)/(wpb_interp2d-1.d0), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-1.d0)/(wpb_interp1d-1.d0), wp_interp1d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-1.d0)/(wpb_integ1d-1.d0), wp_integ1d

!start ng-extrapolation from beginning
         s1=i+1
         s2=i+2
         s3=i+3
         s4=i+4
         indx_threshold=i+5
!
      endif
   endif
!
!--------------calculating alo-corrected source-functions---------------
!
   call sline_new3d
!
!------------extrapolation of old subsequent source functions-----------
!
   if(lng_expol) then
!start ng-extrapolation from beginning, when interpolation parameters have been adapted
         s1=i+1
         s2=i+2
         s3=i+3
         s4=i+4
         indx_threshold=i+5
   endif
!
   if(opt_ng_line.or.opt_ait_line) then
!
!----------storing old source-functions for ng-extrapolation------------
!
      if(i.eq.s1) then
         call store_ng3d(1,sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
         s1=s1+ng_const
      elseif(i.eq.s2) then
         call store_ng3d(2,sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
         s2=s2+ng_const
      elseif(i.eq.s3) then
         call store_ng3d(3,sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
         s3=s3+ng_const
      elseif(i.eq.s4) then
         call store_ng3d(4,sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
         s4=s4+ng_const
         if(opt_ng_line) call ng_expol3d(sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
         if(opt_ait_line) call ait_expol3d(sline3d_ng,ndxmax,ndymax,ndzmax,sline3d)
      endif
!
   endif
!
enddo
!
!
!
end subroutine linetrans3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conttrans2d
!
!-------------------------continuum transport---------------------------
!
!---------calculating intensities for different angles mu---------------
!--------------------------in xz-plane----------------------------------
!--------------performing angle integration on z-axis-------------------
!--back interpolation of newly calculated source function onto 3d grid--
!
use prog_type
use iter, only: itmaxc, devmaxc, epsmaxc_arr
use dimecr, only: n1d_cr, scont1d_cr, mint1d_cr, r1d_cr, alocont1d_diag_cr, &
                  norm1d_cr, epshistoryc1d_cr, epsc1d_cr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, scont3d, t3d, eps_cont3d, imask_innreg3d
use options, only: opt_ait_cont, opt_ng_cont, opt_method
use ng_extra, only: ng_const
use freq, only: xnue0, lfreqint
!use timing, only: ts_it, te_it, ttot_it, it_tot
use warnings, only: warn_itmaxc
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k, l, err
integer(i4b) :: s1, s2, s3, s4
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: dummy1, dummy2, rad
real(dp) :: eps_max
!
! ... local arrays
real(dp), dimension(:), allocatable :: bnue1d_cr!, r1d_rev, scont1d_rev
real(dp), dimension(:,:), allocatable :: scont1d_ng
!
!... local characters
!character(len=50) :: enter
!
! ... local functions
real(dp) :: bnue2
!
! ... local logicals
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------continuum transport (cr)------------------------------'
write(*,*)
!
allocate(bnue1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans2d'
allocate(scont1d_ng(4,n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans2d'
!
!-----------------------------------------------------------------------
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
!
!-----------------------------------------------------------------------
!
!calculating start values of continuum source function
call calc_startval_cont
!
!calculating log(r3d)
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!         log3d(i,j,k)=log10(sqrt(x(i)**2 + y(j)**2 + z(k)**2))
!      enddo
!   enddo
!enddo
!
!calculating planck function along z direction
do i=1, n1d_cr
   bnue1d_cr(i) = bnue2(xnue0, t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i), lfreqint)
   scont1d_cr(i) = scont3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
enddo
!
!setting start mean intensities and start deviations to zero
epsc1d_cr=0.d0
epshistoryc1d_cr=0.d0
mint1d_cr=0.d0
alocont1d_diag_cr=0.d0
!
!-----------------------------------------------------------------------
!
!ttot_it=0.d0
!
!************************start iteration scheme*************************
!
do i=1, itmaxc
!
!timing
!it_tot=i
!call cpu_time(ts_it)
!
!****************************output*************************************
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j', 's_c', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mint1d_cr(j), scont1d_cr(j), epsc1d_cr(j), norm1d_cr(j), alocont1d_diag_cr(j)
   end do
   write(*,*)
!
!*************************output end************************************
!
   epsc1d_cr=mint1d_cr
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
!
   select case(opt_method)
      case(-1)
         write(*,*) 'method: returning after initial condition'
         return
      case(0)
         call mint_fvm2d
      case(1)
         call mint_sc2d_lin
      case(2)
         call mint_sc2d
      case default
         stop 'opt_method not properly specified: 0 for FVM, 1 for SClin, 2 for SCbez'
   end select
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsc1d_cr, mint1d_cr, n1d_cr, eps_max)
   epshistoryc1d_cr(:,i)=epsc1d_cr
   epsmaxc_arr(i)=eps_max
!
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
!      call cpu_time(te_it)
!      ttot_it = ttot_it + te_it - ts_it
      exit
   else if(i.eq.itmaxc) then
      write(*,*) "no convergence after iteration no. ", i
      warn_itmaxc=.true.
   end if
!
!-------calculating alo-corrected source-functions on central ray-------
!
   write(*,*) '--------------calculating new iterate for source function (alo)----------------'
   write(*,*)
!
   do j=1, n1d_cr
      dummy2=1.d0-eps_cont3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)
      dummy1=1.d0-(1.d0-eps_cont3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j))*alocont1d_diag_cr(j)
      scont1d_cr(j)=(dummy2/dummy1)*mint1d_cr(j) - &
             (dummy2/dummy1)*alocont1d_diag_cr(j)*scont3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j) + &
              eps_cont3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*bnue1d_cr(j)/dummy1
   enddo
!
!------------extrapolation of old subsequent source functions-----------
!
   if(opt_ng_cont.or.opt_ait_cont) then
!
!----------storing old source-functions for ng-extrapolation------------
!
      if(i.eq.s1) then
         write(*,*) '----------------------storing source fct. at step n-3--------------------------'
         write(*,*)
         scont1d_ng(1,:)=scont1d_cr
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         scont1d_ng(2,:)=scont1d_cr
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         scont1d_ng(3,:)=scont1d_cr
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         scont1d_ng(4,:)=scont1d_cr
         s4=s4+ng_const
         if(opt_ng_cont) call ng_expol1d(scont1d_ng, n1d_cr)
         if(opt_ait_cont) call ait_expol1d(scont1d_ng, n1d_cr)
         scont1d_cr=scont1d_ng(1,:)
      endif

   endif
!
!-----back-interpolation of central-ray-source-functions on 3d-grid-----
!
!   do j=1, n1d_cr
!      scont1d_rev(j) = scont1d_cr(n1d_cr+1-j)
!      r1d_rev(j) = r1d_cr(n1d_cr+1-j)/sr
!   enddo
!!
!   scont1d_rev=log10(scont1d_rev*r1d_rev*r1d_rev)
!   r1d_rev=log10(r1d_rev)
!
   do j=1, ndxmax
      do k=1, ndymax
         do l=1, ndzmax
            if(imask_innreg3d(j,k,l).eq.1) then
               scont3d(j,k,l) = 0.d0
            else
               rad=sqrt(x(j)**2 + y(k)**2 + z(l)**2)
               call find_index(rad, r1d_cr, n1d_cr, iim2, iim1, ii, iip1)
               scont3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), scont1d_cr(iim1), scont1d_cr(ii), rad)
            endif
         enddo
      enddo
   enddo
!
   call output_itstep_cont(i)
!
!!timing
!   call cpu_time(te_it)
!   ttot_it=ttot_it + te_it - ts_it
!   write(*,*) 'time per iteration: ', te_it-ts_it
!
enddo
!
!
!
end subroutine conttrans2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conttrans3d
!
!---------------------------continuum transport-------------------------
!
!-------calculating intensities for different angles theta, phi---------
!-----------on a 3d cartesian grid using short characteristics----------
!------------------------or finite volume method------------------------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb, zero
use dimecr, only: n1d_cr!, r1d_cr
use dime3d, only: mint3d, scont3d, normalization3d, alocont_nn3d, &
                  imask_totreg3d, imask3d, x, y, z, ndxmax, ndymax, ndzmax, &
                  fcontx3d, fconty3d, fcontz3d, kcontxx3d, kcontyy3d, kcontzz3d, &
                  kcontxy3d, kcontxz3d, kcontyz3d
use iter, only: itmaxc, devmaxc, epsmaxc_arr, it_start_cont
use ng_extra, only: ng_const
use options, only: opt_ng_cont, opt_ait_cont, input_mod_dim, opt_method
use warnings, only: warn_itmaxc
use mod_interp2d, only: interpolation_threshold, wpa_interp2d, wpb_interp2d, wp_interp2d, &
                        wpa_interp1d, wpb_interp1d, wp_interp1d, wpa_integ1d, wpb_integ1d, wp_integ1d
use freq, only: xnue0
use bcondition, only: xic1
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, indx_threshold
integer(i4b) :: ix_epsmax, iy_epsmax, iz_epsmax
real(dp) :: s1, s2, s3, s4, s4b, eps_max, fdum, dummy1, dummy2, rad, trad
real(dp) :: ts, te
!
! ... local arrays
real(dp), dimension(:,:,:), allocatable :: eps3d
real(dp), dimension(:,:), allocatable :: scont3d_ng

!
! ... local logicals
logical :: check_fname
!
! ... local characters
character(len=3) :: fname_debug
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!
write(*,*) '--------------------------continuum transport in 3d----------------------------'
write(*,*)
!
call calc_startval_cont
!
allocate(eps3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans3d'
allocate(scont3d_ng(4,ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans3d'
!
!
s1=it_start_cont+1
s2=it_start_cont+2
s3=it_start_cont+3
s4=it_start_cont+4
s4b=4
!
!setting index for threshold
indx_threshold=it_start_cont+5
!
!-----------------------------------------------------------------------
!
!setting start mean intensities and start deviations to zero
fcontx3d=zero
fconty3d=zero
fcontz3d=zero
kcontxx3d=zero
kcontyy3d=zero
kcontzz3d=zero
kcontxy3d=zero
kcontxz3d=zero
kcontyz3d=zero
eps3d=zero
normalization3d=zero
alocont_nn3d=zero
epsmaxc_arr=zero
!
!setting index for threshold
indx_threshold=it_start_cont+5
!
!---------------------start iteration scheme----------------------------
!
do i=it_start_cont+1, itmaxc
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'z', 'j(sc)', 'scont(sc)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 6(e20.8))') j, z(j), mint3d(ndxmax/2+1,ndymax/2+1,j), scont3d(ndxmax/2+1,ndymax/2+1,j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    alocont_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
   write(*,*)
   write(*,*)
   write(*,'(a25,es20.8)') 'alpha_min (interp1d)', wp_interp1d
   write(*,'(a25,es20.8)') 'alpha_min (interp2d)', wp_interp2d
   write(*,'(a25,es20.8)') 'alpha_min (integ1d)', wp_integ1d
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mint3d
!
   call output_itstep_cont(i)
!
!----------------calculating mean intensities in 3d---------------------
!-------------------------------step i----------------------------------
!
   write(*,*) '----------------calculating angle integrated mean intensities in 3d------------'
   write(*,*) 'step', i
   write(*,*)
!
   select case(opt_method)
      case(-1)
         write(*,*) 'method: returning after initial condition'
         return
      case(0)
         write(*,*) 'method: fvm'
         call mint_fvm3d
      case(1)
         write(*,*) 'method: sc-lin'
         call mint_sc3d_lin
      case(2)
         if(i.le.s4b) then
         !use linear interpolations for the first s4b iteration steps
         !in order to get a first guess of a 'smooth' solution.
         !otherwise: convergence behaviour might be oscillatory
            write(*,*) 'method: sc-lin (iteration < s4b)'
            call mint_sc3d_lin
         else
            write(*,*) 'method: sc-bez'
            call mint_sc3d
         endif
   end select
   write(*,*)
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)
!
   call calc_dev3d(eps3d, mint3d, imask_totreg3d, ndxmax, ndymax, ndzmax, eps_max, ix_epsmax, iy_epsmax, iz_epsmax)
   epsmaxc_arr(i)=eps_max
   write(*,'(a30, 3(i4), 3(f8.4), es18.8)') 'max (dev) at grid-point:', ix_epsmax, iy_epsmax, iz_epsmax, &
                                            x(ix_epsmax), y(iy_epsmax), z(iz_epsmax) , eps_max
   write(*,*)
!
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) 'convergence after iteration no. ', i
      write(*,*) 'max (dev): ', eps_max
      write(*,*)
!for thin lines, ensure that higher order interpolation scheme
!   has been used, and not only linear approach
      if(i.gt.s4b) exit
!ensure that no ng-extrapolation is called for subsequent converged solutions (otherwise NaNs)
      s4=1
   elseif(i.eq.itmaxc) then
      write(*,*) 'no convergence after iteration no. ', i
      write(*,*)
      warn_itmaxc=.true.
    endif
!
   if(i.gt.indx_threshold) then
      if(abs(epsmaxc_arr(i)).ge.abs(epsmaxc_arr(i-1)).and. &
         abs(epsmaxc_arr(i-1)).le.abs(epsmaxc_arr(i-2)).and. &
         abs(epsmaxc_arr(i-2)).ge.abs(epsmaxc_arr(i-3)).and. &
         abs(epsmaxc_arr(i-3)).le.abs(epsmaxc_arr(i-4)).and. &
         abs(epsmaxc_arr(i-4)).ge.abs(epsmaxc_arr(i-5))) then
!
         write(*,*) 'error in iteration scheme: oscillations!!!'
         write(*,*) 'max deviation at iteration i-5', epsmaxc_arr(i-5)
         write(*,*) 'max deviation at iteration i-4', epsmaxc_arr(i-4)
         write(*,*) 'max deviation at iteration i-3', epsmaxc_arr(i-3)
         write(*,*) 'max deviation at iteration i-2', epsmaxc_arr(i-2)
         write(*,*) 'max deviation at iteration i-1', epsmaxc_arr(i-1)
         write(*,*) 'max deviation at iteration i  ', epsmaxc_arr(i)
         write(*,*) 'possible solutions: '
         write(*,*) '   1. use linear interpolations for upwind/downwind source function'
         write(*,*) '      and for upwind intensities to have same solution procedure in'
         write(*,*) '      each iteration step (independent of the source function and intensities'
         write(*,*) '      themselves.'
         write(*,*) '   2. avoid monotonicity constraint in integration of source contribution'
         write(*,*) '      (try linear approximation of source function along ray)'
         write(*,*) '   3. increase the spatial grid resolution, in order to avoid'
         write(*,*) '      monotonicity constraints in quadratic interpolation procedures'
         write(*,*) '   4. increasing interpolation threshold in quadratic interpolation procedure'
         interpolation_threshold=min(interpolation_threshold+0.1d0,1.d0)
         write(*,*) 'setting interpolation threshold to', interpolation_threshold
         wpa_interp2d=wpa_interp2d+1.d0
         wpb_interp2d=wpb_interp2d+1.d0
         wp_interp2d=wpa_interp2d/wpb_interp2d
         wpa_interp1d=wpa_interp1d+1.d0
         wpb_interp1d=wpb_interp1d+1.d0
         wp_interp1d=wpa_interp1d/wpb_interp1d
         wpa_integ1d=wpa_integ1d+1.d0
         wpb_integ1d=wpb_integ1d+1.d0
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-1.d0)/(wpb_interp2d-1.d0), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-1.d0)/(wpb_interp1d-1.d0), wp_interp1d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-1.d0)/(wpb_integ1d-1.d0), wp_integ1d

!start ng-extrapolation from beginning
         s1=i+1
         s2=i+2
         s3=i+3
         s4=i+4
         indx_threshold=i+5
!
      endif
   endif
!
!--------------calculating alo-corrected source-functions---------------
!
   call scont_new3d
!
!------------extrapolation of old subsequent source functions-----------
!

   if(opt_ng_cont.or.opt_ait_cont) then
!
!----------storing old source-functions for ng-extrapolation------------
!
      if(i.eq.s1) then
         call store_ng3d(1,scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         s1=s1+ng_const
      elseif(i.eq.s2) then
         call store_ng3d(2,scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         s2=s2+ng_const
      elseif(i.eq.s3) then
         call store_ng3d(3,scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         s3=s3+ng_const
      elseif(i.eq.s4) then
         call store_ng3d(4,scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         s4=s4+ng_const
         if(opt_ng_cont) call ng_expol3d(scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         if(opt_ait_cont) call ait_expol3d(scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
      endif
!
   endif
!
enddo
!
!
!
end subroutine conttrans3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine mint_sc2d
!
!-----------------------------------------------------------------------
!------------calculates mean intensities on z-axis using----------------
!------------------short characteristics method-------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, int3d
use dimecr, only: n1d_cr, alocont1d_diag_cr, mint1d_cr, norm1d_cr
use angles, only: dim_mu, weight_mu
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: muindx
!
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!------------calculating intensities for given nodes--------------------
!
mint1d_cr=0.d0
norm1d_cr=0.d0
!
write(*,*) '--------------------calculating intensities for all mu-------------------------'
write(*,*)
!
do muindx=1, dim_mu
!
   call fsc_cont2d(muindx)
!
   do j=1, n1d_cr
!on positive z-axis
!      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
!on negative z-axis
      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
   enddo
end do
!
!
!
end subroutine mint_sc2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine mint_sc2d_lin
!
!-----------------------------------------------------------------------
!------------calculates mean intensities on z-axis using----------------
!------------------short characteristics method-------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, int3d
use dimecr, only: n1d_cr, alocont1d_diag_cr, mint1d_cr, norm1d_cr
use angles, only: dim_mu, weight_mu
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: muindx
!
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!------------calculating intensities for given nodes--------------------
!
mint1d_cr=0.d0
norm1d_cr=0.d0
!
write(*,*) '--------------------calculating intensities for all mu-------------------------'
write(*,*)
!
do muindx=1, dim_mu
!
   call fsc_cont2d_lin(muindx)
!
   do j=1, n1d_cr
!on positive z-axis
!      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
!on negative z-axis
      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
   enddo
end do
!
!
!
end subroutine mint_sc2d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine mintbar_sc2d
!
!-----------------------------------------------------------------------
!---------calculates scattering integral on z-axis using----------------
!------------------short characteristics method-------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, int3d, velx3d, vely3d, velz3d, vth3d
use dimecr, only: n1d_cr, aloline1d_diag_cr, mintbar1d_cr, norm1d_cr
use angles, only: dim_mu, weight_mu, nodes_mu
use freq, only: nxobs, weight_xobs, nodes_xobs
use params_input, only: vth_fiducial
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: muindx, xobsindx
real(dp) :: velx_p, vely_p, velz_p, vel_p, vth_p
real(dp) :: xobs, n_x, n_y, n_z, phinorm
!
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!------------calculating intensities for given nodes--------------------
!
mintbar1d_cr=0.d0
norm1d_cr=0.d0
!
write(*,*) '----------calculating intensities for all mu and frequencies-------------------'
write(*,*)
!
do xobsindx=1, nxobs
   write(*,'(a55,i4, a1, i4)') 'calculating all angles for (freq-point/nxobs)', xobsindx, '/', nxobs
   xobs=nodes_xobs(xobsindx)
!
   do muindx=1, dim_mu
      n_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
      n_y=0.d0
      n_z=nodes_mu(muindx)
!
      call fsc_line2d(muindx,xobsindx)
!
      do j=1, n1d_cr
!on positive z-axis
!         mintbar1d_cr(j)=mintbar1d_cr(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)*weight_xobs(xobsindx)
!         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)*weight_xobs(xobsindx)
!on negative z-axis
         velx_p=velx3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)
         vely_p=vely3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)
         velz_p=velz3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)
         vth_p = vth3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)
         vel_p = n_x*velx_p + n_y*vely_p + n_z*velz_p
         call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
!
         mintbar1d_cr(j)=mintbar1d_cr(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)*weight_xobs(xobsindx)*phinorm
         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)*weight_xobs(xobsindx)*phinorm
      enddo
!
   enddo
enddo
write(*,*)
!
!
!
end subroutine mintbar_sc2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine mint_fvm2d
!
!-----------------------------------------------------------------------
!------------calculates mean intensities on z-axis using----------------
!----------------------finite volume method-----------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, int3d, alocont_nn3d
use dimecr, only: n1d_cr, mint1d_cr, norm1d_cr, alocont1d_diag_cr
use angles, only: dim_mu, weight_mu
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: muindx
!
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!------------calculating intensities for given nodes--------------------
!
mint1d_cr=0.d0
norm1d_cr=0.d0
alocont1d_diag_cr=0.d0
!
write(*,*) '---------------------calculating intensities for all mu------------------------'
write(*,*)
!
do muindx=1, dim_mu
!
   call ffvm_cont2d(muindx)
!
   do j=1, n1d_cr
!on positive z-axis
!      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
!      alocont1d_diag_cr(j)=alocont1d_diag_cr(j) + alocont_nn3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+j,4)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
!on negative z-axis
      mint1d_cr(j)=mint1d_cr(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
      alocont1d_diag_cr(j)=alocont1d_diag_cr(j) + alocont_nn3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j,4)*weight_mu(muindx)
      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
   enddo
end do
!
!
!
end subroutine mint_fvm2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine mintbar_fvm2d
!
!-----------------------------------------------------------------------
!---------calculates scattering integral on z-axis using----------------
!------------------short characteristics method-------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, int3d, velx3d, vely3d, velz3d, vth3d, mintbar3d, normalization3d
use dimecr, only: n1d_cr, aloline1d_diag_cr, mintbar1d_cr, norm1d_cr
use angles, only: dim_mu, weight_mu, nodes_mu
use freq, only: nxobs, weight_xobs, nodes_xobs
use params_input, only: vth_fiducial
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: muindx, xobsindx
real(dp) :: velx_p, vely_p, velz_p, vel_p, vth_p
real(dp) :: xobs, n_x, n_y, n_z, phinorm
!
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!------------calculating intensities for given nodes--------------------
!
mintbar1d_cr=0.d0
norm1d_cr=0.d0
mintbar3d = 0.d0
normalization3d = 0.d0
!
write(*,*) '----------calculating intensities for all mu and frequencies-------------------'
write(*,*)
!
do xobsindx=1, nxobs
   write(*,'(a55,i4, a1, i4)') 'calculating all angles for (freq-point/nxobs)', xobsindx, '/', nxobs
   do muindx=1, dim_mu
      call ffvm_line2d(muindx,xobsindx)
      do i=1, n1d_cr
         mintbar1d_cr(i) = mintbar3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)
         norm1d_cr(i) = normalization3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)
      enddo
   enddo
enddo
write(*,*)
!
!
!
end subroutine mintbar_fvm2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_startval_cont
!
!-----------------------------------------------------------------------
!---------calculates start value of continuum source function-----------
!-----------------------------------------------------------------------
!
use prog_type
use dimecr, only: n1d_cr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, scont3d, t3d, eps_cont3d, mint3d, imask_innreg3d, imask_totreg3d, opac3d, imask3d
use freq, only: xnue0, lfreqint
use options, only: opt_start_cont
use iter, only: it_start_cont
use params_stellar, only: sr
!
implicit none
!
! ... local scalars
integer(i4b) :: i,j,k, kk
integer(i4b) :: indx_therm3d, indx_therm1d
real(dp) :: opac1, opac2, r1, r2, dtau, tau, thdepth, rad
real(dp) :: r_therm, t_therm, i_core, dum_bnue, dilfac
real(dp) :: eps_cont
!
! ... local arrays
!
! ... local functions
real(dp) :: bnue2, bint
!
!-----------------------calculate thermalization depth------------------
!
!calculate a mean of eps_cont
call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)


if(eps_cont.eq.0.d0) then
   !set thermalization depth to arbitrary value
   thdepth = 10.d0
else
   thdepth = 1.d0/sqrt(eps_cont)
endif
!
!---------------------calculate spherical symmetric tau-----------------
!
!start at outer boundary with tau=0
tau=0.d0
indx_therm1d=1
!
do i=ndzmax-1, ndzmax-n1d_cr, -1
   opac1=opac3d(ndxmax/2+1,ndymax/2+1,i+1)
   opac2=opac3d(ndxmax/2+1,ndymax/2+1,i)
   r1=z(i+1)!*sr
   r2=z(i)!*sr
   dtau=(opac2+opac1)*(r1-r2)/2.d0
   tau=tau+dtau
   indx_therm3d=i
   indx_therm1d=indx_therm1d+1
!
!   write(*,'(7es20.8,3l5)') tau, opac1, opac2, r1, r2, dtau, thdepth, mask_innreg3d(ndxmax/2+1,ndymax/2+1,i), tau.gt.thdepth, .not.mask_innreg3d(ndxmax/2+1,ndymax/2+1,i)
   if(imask_innreg3d(ndxmax/2+1,ndymax/2+1,i).eq.1) exit
   if(tau.gt.thdepth) exit
enddo

r_therm=z(indx_therm3d)
t_therm=t3d(ndxmax/2+1,ndymax/2+1,indx_therm3d)
!
i_core=bnue2(xnue0, t_therm, lfreqint)
!
!-----------------------------------------------------------------------
!
if(opt_start_cont) then
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            if(imask_totreg3d(i,j,k).ne.0) then
               rad=sqrt(x(i)**2 + y(j)**2 + z(k)**2)
               if(rad.le.r_therm) then
                  !when thermalized, source function = planck function
                  scont3d(i,j,k) = bnue2(xnue0, t3d(i,j,k), lfreqint)
               else
                  !when not thermalized, calculate source function from
                  !mean intensity (i_core * dilution factor)
                  dilfac = r_therm**2 / rad**2
                  dilfac = 0.5d0*(1.d0 - sqrt(1.d0 - dilfac))
                  dum_bnue=bnue2(xnue0, t3d(i,j,k), lfreqint)
                  scont3d(i,j,k)=(1.d0-eps_cont)*dilfac*i_core + eps_cont*dum_bnue
               endif
            endif
         enddo
      enddo
   enddo
!   scont3d=0.d0
   mint3d=0.d0
   it_start_cont=0
else
   call input_itstep_cont(it_start_cont)
endif
!
!
!
end subroutine calc_startval_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_startval_line
!
!-----------------------------------------------------------------------
!-----------calculates start value of line source function--------------
!-----------------------------------------------------------------------
!
use prog_type
use dime3d, only: sline3d, ssobo3d, mintbar3d, imask3d, ndxmax, ndymax, ndzmax
use options, only: opt_start_line
use iter, only: it_start_line, epsmaxl_arr
!
implicit none
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
!
! ... local arrays
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(opt_start_line) then
!start value from sobolev approximation
   mintbar3d=0.d0
   epsmaxl_arr=0.d0
   it_start_line=0
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            select case(imask3d(i,j,k))
               case(1,2,3)
                  sline3d(i,j,k)=ssobo3d(i,j,k)
               case default
            endselect
         enddo
      enddo
   enddo
else
   call input_itstep_line(it_start_line)
endif
!
end subroutine calc_startval_line
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_cr
!
use prog_type
use dime3d, only: ndzmax, z
use dimecr, only: n1d_cr, r1d_cr, scont1d_cr, sline1d_cr, mint1d_cr, mintbar1d_cr, &
                  epsc1d_cr, epsl1d_cr, epshistoryc1d_cr, epshistoryl1d_cr, &
                  norm1d_cr, alocont1d_diag_cr, aloline1d_diag_cr

use iter, only: itmaxc
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
!--------------calculate number of data points on central ray-----------
!
n1d_cr=0
do i=ndzmax/2+1, ndzmax-1
   if(z(i).ge.1.d0) then
      n1d_cr=n1d_cr+1
   endif
enddo
!
if(n1d_cr.eq.0) stop 'error in grid_cr: n1d_cr = 0'
!
!-----------------------------------------------------------------------
!
allocate(r1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(norm1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(scont1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(sline1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(mint1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(mintbar1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(alocont1d_diag_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(aloline1d_diag_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
!
allocate(epsc1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(epshistoryc1d_cr(n1d_cr, itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(epsl1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
allocate(epshistoryl1d_cr(n1d_cr, itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in grid_cr'
!
!-----------------------------------------------------------------------
!
do i=1, n1d_cr
   r1d_cr(i) = z(ndzmax-n1d_cr-1+i)
enddo
!
norm1d_cr=0.d0
!
end subroutine grid_cr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine make_debug
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask3d, int3d, opalbar3d, velx3d, vely3d, velz3d, imask_totreg3d
use angles, only: dim_omega, n_x, n_y, n_z
use bcondition, only: xic1
use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, indx_omega, indx_ipoint, err
integer(i4b) :: n_inner
real(dp) :: nn_x, nn_y, nn_z, mueff
real(dp) :: ts, te
!
! ... local arrays
integer, dimension(:), allocatable :: indx_x, indx_y, indx_z
real(dp), dimension(:), allocatable :: iboundary
real(dp), dimension(:,:,:), allocatable :: int3d_test
!
! ... local characters
character(len=32) :: fname
!
write(*,*) '------------------------running dummy debug routine----------------------------'
!
!j=ndymax/2+1
!k=ndzmax/2+1
!i=ndxmax/2+1
!do i=38-2,38+2
!   do j=ndymax/2+1-2, ndymax/2+1+2
!      do k=ndzmax/2+1-2, ndzmax/2+1+2
!         write(*,'(4es20.8,i5)') x(i), y(j), z(k), (x(i)/smajorax_a)**2 + (y(j)/smajorax_b)**2 + (z(k)/smajorax_c)**2, imask3d(i,j,k)
!      enddo
!   enddo
!enddo
!stop
!do i=1, ndxmax
!   write(*,*) i, x(i), imask3d(i,j,k), mask_totreg3d(i,j,k)
!enddo
!stop 'go on in debug routine'
!
!setting arrays to describe the inner boundary condition
n_inner=0
do i=1, ndxmax 
   do j=1, ndymax
      do k=1, ndzmax
         if(imask3d(i,j,k).eq.4) n_inner=n_inner+1
      enddo
   enddo
enddo
write(*,*) 'number of inner points', n_inner
allocate(indx_x(n_inner), stat=err)
allocate(indx_y(n_inner), stat=err)
allocate(indx_z(n_inner), stat=err)
allocate(iboundary(n_inner), stat=err)
!
!setting boundary condition for inner points
do indx_omega=1, dim_omega
   nn_x = n_x(indx_omega)
   nn_y = n_y(indx_omega)
   nn_z = n_z(indx_omega)
!
   indx_ipoint=1
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            if(imask3d(i,j,k).eq.4) then
!inside the star, set intensities correctly
               indx_x(indx_ipoint) = i
               indx_y(indx_ipoint) = j
               indx_z(indx_ipoint) = k
               mueff=(nn_x*x(i)+nn_y*y(j)+nn_z*z(k))/sqrt(x(i)**2+y(j)**2+z(k)**2)
               if(mueff.ge.0.d0) then
                  iboundary(indx_ipoint) = xic1
               else
                  iboundary(indx_ipoint) = 0.d0
               endif
               indx_ipoint=indx_ipoint+1
            endif
         enddo
      enddo
   enddo
!
   write(fname,'(a19,i4.4)') 'data_bcondition/iomega', indx_omega
   open(1, file=fname, form='unformatted')
      write(1) iboundary
      write(1) indx_x
      write(1) indx_y
      write(1) indx_z
   close(1)
!
enddo
!
!setting boundary condition for inner points using 'standard' procedure
int3d=0.d0
call cpu_time(ts)
do indx_omega=1, dim_omega
   nn_x = n_x(indx_omega)
   nn_y = n_y(indx_omega)
   nn_z = n_z(indx_omega)

   do k=3, ndzmax-2
      do j=3, ndymax-2
         do i=3, ndxmax-2
            if(imask3d(i,j,k).eq.4) then
!inside the star, set intensities correctly
               mueff=(nn_x*x(i)+nn_y*y(j)+nn_z*z(k))/sqrt(x(i)**2+y(j)**2+z(k)**2)
               if(mueff.ge.0.d0) then
                  int3d(i,j,k) = xic1
               else
                  int3d(i,j,k) = 0.d0
               endif
            endif
         enddo
      enddo
   enddo
!
enddo

call cpu_time(te)

write(*,*) 'total time for standard', te-ts
!
!
!
!setting boundary condition for inner points by reading in everything
allocate(int3d_test(ndxmax,ndymax,ndzmax))
int3d_test=0.d0
call cpu_time(ts)
do indx_omega=1, dim_omega
!
   write(fname,'(a19,i4.4,a5,i4.4)') 'data_bcondition/iomega', indx_omega
!
   open(1, file=fname, form='unformatted')
      read(1) iboundary
      read(1) indx_x
      read(1) indx_y
      read(1) indx_z
   close(1)
    do i=1, n_inner
      int3d_test(indx_x(i),indx_y(i),indx_z(i))=iboundary(i)
   enddo
!
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            if(int3d_test(i,j,k).ne.int3d(i,j,k)) then
               write(*,'(3i5, i3, 3es20.8)') i, j, k, imask3d(i,j,k), int3d_test(i,j,k), int3d(i,j,k)
               mueff=(n_x(indx_omega)*x(i)+n_y(indx_omega)*y(j)+n_z(indx_omega)*z(k))/sqrt(x(i)**2+y(j)**2+z(k)**2)
               write(*,*) mueff
               stop 'error: test and standard different'
            endif
         enddo
      enddo
   enddo
!
enddo

call cpu_time(te)

write(*,*) 'total time for reading and storing', te-ts

stop
!
end subroutine make_debug
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine test_timing
!
!this routine calculates continuum and line transport
!with the FVM and SClin, SCbez methods to test the timing properties
!
use prog_type
use angles, only: dim_omega
use options, only: opt_method
use dime3d, only: ndxmax, ndymax, ndzmax, scont3d, int3d, aloline_on_nn3d, aloline_nn3d, mintbar3d, &
                  mintbar3d_tmp, aloline_nn3d_tmp, normalization3d_tmp, normalization3d, sline3d, &
                  x, y, z, ssobo3d, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, &
                  opalbar3d, t3d     
use mod_directories, only: output_dir_test
use bcondition, only: xic1
use freq, only: nxobs, xnue0
use params_input, only: vmin, vmax, beta, eps_line, vth_fiducial
!
implicit none
!
integer(i4b) :: oindx, xobsindx, i, k, err
integer(i4b) :: nticks_sec, nticks_max, nticks_initial, nticks_final
real(dp) :: ts_sclin, ts_scbez, ts_fvm, ts_sclin2, ts_scbez2, ts_fvm2, &
            te_sclin, te_scbez, te_fvm, te_sclin2, te_scbez2, te_fvm2, &
            tt_sclin, tt_scbez, tt_fvm, tt_sclin2, tt_scbez2, tt_fvm2
!
!calculate delx, dely, delz grids if not done yet (required for fvm)
if(opt_method.ne.0) call grid_delxyz
!
write(*,*) '------------------------testing timing properties------------------------------'
!
goto 20
!
!-------------------------timing of continuum---------------------------
!
10 continue
!
call calc_startval_cont
!
tt_sclin=0.d0
tt_scbez=0.d0
tt_fvm=0.d0
open(1, file=output_dir_test//'/test_timing_continuum.dat', form='formatted')
write(*,*) 
write(*,'(11a20)') 'oindx', 'dim_omega', 't_tot(fvm)', 't_tot(sclin)', 't_tot(scbez)', &
                   't_i(fvm)', 't_i(sclin)', 't_i(scbez)', '<t_i>(fvm)', '<t_i>(sclin)', '<t_i>(scbez)'
write(1,'(11a20)') 'oindx', 'dim_omega', 't_tot(fvm)', 't_tot(sclin)', 't_tot(scbez)', &
                   't_i(fvm)', 't_i(sclin)', 't_i(scbez)', '<t_i>(fvm)', '<t_i>(sclin)', '<t_i>(scbez)'

do oindx=1, dim_omega
   call cpu_time(ts_sclin)
   call fsc_cont3d_lin(oindx)
   call cpu_time(te_sclin)
   tt_sclin = tt_sclin + te_sclin-ts_sclin

   call cpu_time(ts_scbez)
   call fsc_cont3d(oindx)
   call cpu_time(te_scbez)
   tt_scbez = tt_scbez + te_scbez-ts_scbez
!
   call cpu_time(ts_fvm)
   call ffvm_cont3d(oindx)
   call cpu_time(te_fvm)
   tt_fvm = tt_fvm + te_fvm-ts_fvm
!
   write(*,'(2i20, 9es20.8)') oindx, dim_omega, &
                             tt_fvm, tt_sclin, tt_scbez, &
                             te_fvm-ts_fvm, te_sclin-ts_sclin, te_scbez-ts_scbez, &
                             tt_fvm/oindx, tt_sclin/oindx, tt_scbez/oindx
   write(1,'(2i20, 9es20.8)') oindx, dim_omega, &
                             tt_fvm, tt_sclin, tt_scbez, &
                             te_fvm-ts_fvm, te_sclin-ts_sclin, te_scbez-ts_scbez, &
                             tt_fvm/oindx, tt_sclin/oindx, tt_scbez/oindx
enddo
!
!
!
!running complete mint-routine directly
call cpu_time(ts_fvm)
call mint_fvm3d
call cpu_time(te_fvm)
tt_fvm=te_fvm-ts_fvm
!
call cpu_time(ts_sclin)
call mint_sc3d_lin
call cpu_time(te_sclin)
tt_sclin=te_sclin-ts_sclin
!
call cpu_time(ts_scbez)
call mint_sc3d
call cpu_time(te_scbez)
tt_scbez=te_scbez-ts_scbez
!
write(*,'(a60)') 'total time for calculating mean intensity'
write(*,'(3a20)' ) 'FVM', 'SClin', 'SCbez'
write(*,'(3es20.8)') tt_fvm, tt_sclin, tt_scbez
write(*,*)
!
write(1,'(a60)') 'total time for calculating mean intensity'
write(1,'(3a20)' ) 'FVM', 'SClin', 'SCbez'
write(1,'(3es20.8)') tt_fvm, tt_sclin, tt_scbez
!

close(1)

!
!-----------------------------timing of line----------------------------
!
20 continue
!
call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
call calc_startval_line
!
!deallocation of global (threadprivate) arrays
if(allocated(int3d)) deallocate(int3d)
if(allocated(aloline_on_nn3d)) deallocate(aloline_on_nn3d)
!
!set clock-ticks
call system_clock(count=nticks_initial, count_rate=nticks_sec, count_max=nticks_max)
!
open(1, file=output_dir_test//'/test_timing_line.dat', form='formatted')
write(*,'(13a20)') 'xobsindx', 'nxobs', 'oindx', 'dim_omega', 't_tot(fvm)', 't_tot(sclin)', 't_tot(scbez)', &
                   't_i(fvm)', 't_i(sclin)', 't_i(scbez)', '<t_i>(fvm)', '<t_i>(sclin)', '<t_i>(scbez)'
write(1,'(13a20)') 'xobsindx', 'nxobs', 'oindx', 'dim_omega', 't_tot(fvm)', 't_tot(sclin)', 't_tot(scbez)', &
                   't_i(fvm)', 't_i(sclin)', 't_i(scbez)', '<t_i>(fvm)', '<t_i>(sclin)', '<t_i>(scbez)'
!
!begin of parallel region
!$omp parallel &
!$omp private(err, oindx, xobsindx)
!
!allocation of global (threadprivate) arrays
allocate(aloline_on_nn3d(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error in ffvm_line3d: aloline_on_nn3d'
allocate(aloline_nn3d_tmp(ndxmax,ndymax,ndzmax,27), stat=err)
   if(err.ne.0) stop 'allocation error in ffvm_line3d: aloline_nn3d_tmp'
allocate(normalization3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error ffvm_line3d: normalization3d_tmp'
allocate(mintbar3d_tmp(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error ffvm_line3d: mintbar3d_tmp'
allocate(int3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error ffvm_line3d: int3d'
aloline_nn3d_tmp=0.d0
normalization3d_tmp=0.d0
mintbar3d_tmp=0.d0
!
tt_sclin=0.d0
tt_scbez=0.d0
tt_fvm=0.d0
k=1
!$omp do schedule(dynamic)
do xobsindx=1, nxobs, 20
!   write(*,'(a55,i4, a1, i4, a5)') 'calculating all angles for (freq-point/nxobs)', xobsindx, '/', nxobs
   do oindx=1, dim_omega, 40
      call cpu_time(ts_sclin)
      call fsc_line3d_lin(oindx,xobsindx)
      call cpu_time(te_sclin)
      tt_sclin = tt_sclin + te_sclin-ts_sclin

      call cpu_time(ts_scbez)
      call fsc_line3d(oindx,xobsindx)
      call cpu_time(te_scbez)
      tt_scbez = tt_scbez + te_scbez-ts_scbez
!
      call cpu_time(ts_fvm)
      call ffvm_line3d(oindx,xobsindx)
      call cpu_time(te_fvm)
      tt_fvm = tt_fvm + te_fvm-ts_fvm

      write(*,'(4i20, 9es20.8)') xobsindx, nxobs, oindx, dim_omega, &
                                 tt_fvm, tt_sclin, tt_scbez, &
                                 te_fvm-ts_fvm, te_sclin-ts_sclin, te_scbez-ts_scbez, &
                                 tt_fvm/k, tt_sclin/k, tt_scbez/k
      write(1,'(4i20, 9es20.8)') xobsindx, nxobs, oindx, dim_omega, &
                                 tt_fvm, tt_sclin, tt_scbez, &
                                 te_fvm-ts_fvm, te_sclin-ts_sclin, te_scbez-ts_scbez, &
                                 tt_fvm/k, tt_sclin/k, tt_scbez/k
      k=k+1
   enddo
enddo
!$omp enddo
!
!add up temporary arrays
!$omp critical
   aloline_nn3d = aloline_nn3d+aloline_nn3d_tmp
   normalization3d = normalization3d + normalization3d_tmp
   mintbar3d = mintbar3d + mintbar3d_tmp
!$omp end critical
!
!deallocation of global (threadprivate) arrays
deallocate(aloline_nn3d_tmp)
deallocate(normalization3d_tmp)
deallocate(mintbar3d_tmp)
deallocate(aloline_on_nn3d)
deallocate(int3d)
!
!$omp end parallel
!
!
call system_clock(count=nticks_final)
write(*,*)
write(*,*) 'time for all angles and frequencies', dble(nticks_final-nticks_initial)/nticks_sec
write(*,*) 'nticks_final', nticks_final
write(*,*) 'nticks_init', nticks_initial
write(*,*) 'nticks_sec', nticks_sec
!
!
!
!for complete mintbar-routines
!call system_clock(count=nticks_initial, count_rate=nticks_sec, count_max=nticks_max)
!call mintbar_fvm3d
!call system_clock(count=nticks_final)
!tt_fvm=dble(nticks_final-nticks_initial)/nticks_sec
!!
!call system_clock(count=nticks_initial, count_rate=nticks_sec, count_max=nticks_max)
!call mintbar_sc3d_lin
!call system_clock(count=nticks_final)
!tt_sclin=dble(nticks_final-nticks_initial)/nticks_sec
!!
!call system_clock(count=nticks_initial, count_rate=nticks_sec, count_max=nticks_max)
!call mintbar_sc3d
!call system_clock(count=nticks_final)
!tt_scbez=dble(nticks_final-nticks_initial)/nticks_sec
!!
!write(*,'(a60)') 'total time for calculating scattering integral'
!write(*,'(3a20)' ) 'FVM', 'SClin', 'SCbez'
!write(*,'(3es20.8)') tt_fvm, tt_sclin, tt_scbez
!write(*,*)
!write(1,'(a60)') 'total time for calculating scattering integral'
!write(1,'(3a20)' ) 'FVM', 'SClin', 'SCbez'
!write(1,'(3es20.8)') tt_fvm, tt_sclin, tt_scbez
!close(1)

end subroutine test_timing
