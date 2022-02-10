!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------------program----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!v1: formal solution of intensity along a given ray
!    several test opacities provided (see module test_formal_sol)
!
!v2: external module file: mod_forsol.f90
!      including subroutines to read in model atmosphere and source-functions
!         (read_dim3d, allocate_fs3d, read_mod3d, read_source3d)
!      test-routines decoupled from main program
!      setting up central ray ((p=0 => only radial grid)
!         -> including subroutines to transform coordinate systems
!         -> including subroutines to interpolate physical quantities onto ray
!
!v3: calculating grid of p-impact parameter (equidistant in log-space)
!    for now: z-rays are calculated with equidistant steps
!    including calculation of profile functions
!    including xobs-grid
!    including zeta-grid, with corresponding weights
!    including continuum description in ray-solution
!       -> problem: division by chi_tot (outside of 3-d grid with extrapolation)
!                   gives floating point error
!
!v4: calculating z-coordinates for each p-ray as z=sqrt(r_i^2 - p^2)
!       -> calculating reasonable radial grid (with more data points than input)
!    integration scheme
!    test procedure with analytic opacity and sourcefct - distributions
!    including error estimation of p- and zeta-integration
!
!v5: new subroutines to interpolate physical quantities onto ray
!        (old ones did not work properly)
!    including subroutines to read photspheric profile from external file
!        so far, external profile is only read in, not used
!
!v8: including photospheric profile in formal solver routines
!    complete decoupling from old module-files (for grid.eo)
!    debugging calculation of ovservers frequency grid:
!       del_xobs is given in v_thermal, whereas now, 
!       arbitrary vth_fiducial-values are allowed, without changing
!       grid (nor profile function)
!    parallelization of zeta-integration in formal-solver-routine
!
!v9: small debugs: make same units in complete code:
!     velocities in vthfiducial
!     integrated line opacities in 1/r_star
!     profile-fct as phi_x=1/sqrt(pi)/delta * exp(((x-v)/delta)^2)
!       with delta = vth/vthfiducial
!     => line opacity = integrated line opacity*profile
!     including adm model: interpolation from 2d-grid
!     including subroutine to calculate vthermal as a function
!        of position: vmicro increases linearly with abs(vel)
!     including options and read-in procedures to calculate only surface-brightness
!     including options and subroutines to calculate int2d and tau2d along a direction
!     note: bug found, when interpolating velocities and shear-velocities are present
!
!v10: including subroutine setup_ray_adm, which calculates the adm-velocity law
!         directly when setting up ray (instead of interpolation)
!      additional subroutines to check and compare differences of adm-interpolation
!         and direct calculation of adm-velocity law
!      consistency checks, such that input of 3d-solution and formal solver
!         are equal
!      including subroutine setup_ray_bvel, which calculates beta-velocits law
!         directly when setting up ray (instead of interpolation)
!      including subroutine setup_ray_spc, which creates the formal-solution
!         from a spherical grid (read in from mhd simulations); the source-function
!         is interpolated onto the spherical grid
!
!v11: new names for all f90-files:
!          formerly                   updated
!          formal_sol.f90             spec.f90
!          formal_ray.f90             formal_lc.f90
!          input_forsol.f90           spec_input.f90
!          output_forsol.f90          spec_output.f90
!          mod_formal.f90             mod_spec.f90
!          no model-calculations in this routine anymore, only read-in routines
!      including subroutine setup_ray3d_rot, which accounts for the
!          distortion of the stellar surface for large rotational velocities
!
!v12: including continuum along each ray
!      debugging subrutine formal_lc to combine continuum and line together
!      including vth_lim as minimum allowed thermal velocity (if T and vmicro low)
!
!vNEW: including the possibility to calculate several surface-brightnesses
!      including parallelization of surface-brightness calculations
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!todo: problem with weights for continuum integration:
!         if integration is only performed over core (i=ic for p<=1, i=0 else)
!         outermost weight is too large (step-function)
!
!note: normalized profile is restricted to +/- x_obs (max)
!         needs to be set to +/- 3 for puls diploma models!!!
!
!      input velocities in vthfiducial
!      input opacities in 1/sr
!
!-----------------------------------------------------------------------
!-----------------------start-------------------------------------------
!-----------------------------------------------------------------------
!
program spectrum3d
!
use prog_type
use fund_const, only: pi, cgs_clight
use options_spec, only: opt_photprof, input_mod, interp_photprof, opt_surface, opt_int2d, output_dir, opt_incl_gdark, opt_incl_sdist
use params_spec, only: vth_min, vth_lim, tmin, na, vth_fiducial, vmicro, xic1, vmax, teff, xlogg, yhe, xnue0, trad
use dime_spec, only: np, nxobs_fs, nzeta, nr, nalpha, ngamma
use mod_spectrum, only: del_vel, del_vel2, xobs, xnue, xic_nue, xicc_nue, &
                        acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic, rmax, rmin, &
                        flux_tot, flux_cont, normt, ftot_err, fcont_err, flux_emi, flux_abs, &
                        flux_tot_tmp, flux_cont_tmp, normt_tmp, ftot_err_tmp, fcont_err_tmp, flux_emi_tmp, flux_abs_tmp, &
                        ftot_p, fcont_p, ftot_errp, fcont_errp, normt_p, femi_p, fabs_p, &
                        relerr_tot, relerr_totp, relerr_cont, relerr_contp, &
                        zeta, zetaw, zetaw_err, p, pw, pw_err, &
                        lcore, phinorm, iin, iin_c, iem, iem_c, iemi, iabs, &
                        nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, temp_ray, profile_ray, vth_ray, velz_ray, &
                        alpha, gamma, alpha_arr, gamma_arr
use timing_spec, only: ts_tot, te_tot, t_tot, ts_obsdir, te_obsdir, ttot_obsdir, t_setup, t_setup_tmp, &
                       t_hunt, t_hunt_tmp, t_trilin, t_trilin_tmp, ticks_obsdir, &
                       ticks_initial_obsdir, ticks_final_obsdir, ttot_obsdir_sys
use mod_surfb, only: xobss_surface, alphas_surface, gammas_surface, nsurfb
use mod_gdark, only: xic1_factor
use omp_lib
!
implicit none
!
! ... local scalars
integer(i4b) :: indx_p, indx_zeta, indx_xobs, indx_zray, indx_alpha, indx_gamma, indx_num
integer(i4b) :: err
integer(i4b) :: ticks_sec, ticks_max
real(dp) :: err1, err2, sum, sum2
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal
!
!
!
call timing_init
call cpu_time(ts_tot)
!
call read_input
!
!
!------------------------read in model atmosphere-----------------------
!
select case(input_mod)
!
!---------------------------1d model------------------------------------
!
   case(0) 
      call read_model1d
      call print_model1d
      call calc_rotmat(0.d0, 0.d0)
      if(opt_incl_sdist.or.opt_incl_gdark) then
         write(*,*) 'error in main: 1d model including opt_incl_sdist/opt_incl_gdark'
         write(*,*) '   still to be implemented (see setup_ray3d_rot)'
         stop
      endif
!
!--------------------3d model in cartesian coordinates------------------
!
   case(1) 
      call read_model3d
      call print_model3d
      call calc_rotmat(0.d0, 0.d0)
!
!--------------------3d model in spherical coordinates------------------
!
   case(2)
      call read_model3d_spc
      call print_model3d_spc
      call calc_rotmat(0.d0, 0.d0)
      if(opt_incl_sdist.or.opt_incl_gdark) then
         write(*,*) 'error in main: 3d spc model including opt_incl_sdist/opt_incl_gdark'
         write(*,*) '   still to be implemented (see setup_ray3d_rot)'
         stop
      endif
!!
!!-----------------------halpha-model------------------------------------
!!
!   case(1)
!      call read_input_forsol
!!note: overwriting some input from read_input_forsol
!      call read_model1d_halpha
!!      call calc_rotmat(0.d0, 0.d0)
!!
!      rmax=100.d0
!      rmin=1.d0
!!
!      call print_physical
!      call print_model1d
!!      call check_interpolation1d
!!
!!--------------------halpha-model from petrenz & puls-------------------
!!
!   case(2)
!      call read_input_forsol
!!
!      rmax=100.d0
!      rmin=1.d0
!      call read_model1d_petrenz
!!
!      call calc_rotmat(0.d0, 0.d0)
!!
!      call print_physical
!      call print_model1d
!!      call check_interpolation1d
!!
!!--------------------resonance-line-model-------------------------------
!!
!   case(3,4)
!      call read_input_forsol
!      call read_model1d_res
!!
!      call calc_rotmat(0.d0, 0.d0)
!!
!      rmax=100.d0!12.d0!100.d0
!      rmin=1.d0
!!
!      call print_physical
!      call print_model1d
!!     call check_interpolation1d
!!
!!--------------------sobolev-theory-model-------------------------------
!!
!   case(5)
!      call read_input_forsol
!!
!      call calc_rotmat(0.d0, 0.d0)
!!
!      rmax=100.d0
!      rmin=1.d0
!!
!      call read_model1d_sobo
!!
!      call print_physical
!      call print_model1d
!!     call check_interpolation1d
!!
!!-------------------------adm-model (2d) for h-alpha--------------------
!!
!   case(6)
!      call read_input_forsol
!!
!      call calc_rotmat(0.d0, 0.d0)
!!
!      rmax=12.d0
!      rmin=1.d0
!      v_esc=6.d7
!!
!      call read_model2d_adm
!
!      call print_physical
!      call print_model2d
!!
!      call output_model2d
!!
!!---------------adm-model (3d) with s_line from main program------------
!!
!   case(7)   
!!
!      write(*,*) '------formal solution for 3d-(adm)-model--------'
!!at the very beginning: read all parameters from namelist (will be overwritten)
!      call read_input_forsol
!!
!      call read_dim3d
!      call allocate_fs3d
!      call read_input3d
!      call read_model3d
!      call read_source3d
!!
!      rmin=1.d0
!      rmax=20.d0
!!recalculate opacity model, and enlarge the grid
!      call recalc_adm_model3d
!!      call recalc_adm2_model3d
!!
!      call check_params
!!
!      call calc_rotmat(obliquity, obliquity_mod3d)
!!
!      call print_physical
!      call print_model3d
!
!!      stop 'check units: all velocities need to be in vthfiducial, opacity in 1/rstar, etc'!
!
!   case(10) 
!      write(*,*) '--------formal solution for 3d-mhd-model--------'
!!at the very beginning: read all parameters from namelist (will be overwritten)
!      call read_input_forsol
!!
!      call read_dim3d
!      call allocate_fs3d
!      call read_input3d
!      call read_model3d
!      call read_source3d
!!
!      rmin=1.d0
!      rmax=20.d0
!!
!      call check_params
!!
 !     call calc_rotmat(obliquity, obliquity_mod3d)
!!
!      call print_physical
!      call print_model3d
!!
!!now, read spherical system, and interpolate source-function onto that system
!      call read_model3d_spc
!!      call print_model3d_spc
!!
!!      stop 'check units: all velocities need to be in vthfiducial, opacity in 1/rstar, etc'
!!
   case default
      stop 'input model not specified'

!
end select
!
!--------------------test procedure for a single p-ray------------------
!
call test_pray
!
!------------------------setting up xobs grid---------------------------
!
!calculate thermal velocity (which is taken at tmin for given mass number)
vth_min = max(vth_lim,vthermal(vmicro, tmin, na))   !thermal velocity not allowed to be smaller than 5 km/s for computational reasons
!
!calculate maximum alowed velocity steps in fiducial vthermal
del_vel2=del_vel*vth_min/vth_fiducial
!
call grid_xobs
!
!--------------------setting up photospheric profile--------------------
!
select case(opt_photprof)
   case(0)
!calculate photospheric profile on xobs-grid (setting it to xic1)
      call calc_photprof(1, nxobs_fs, xic_nue, xicc_nue, xic1)
   case(1)
      call get_photprof_herrero(1, teff, xlogg, yhe, nxobs_fs, xnue, xic_nue, xicc_nue, xnue0, trad)
   case(2)
      call get_photprof(1, teff, xlogg, yhe, nxobs_fs, xnue, xic_nue, xicc_nue, xnue0, xic1)
   case default
      stop 'option opt_photprof not properly selected'
end select
!
!output photospheric profile
!write(*,*) output_dir
!write(*,*) nxobs_fs
!write(*,*) xobs*vth_fiducial/vmax
!write(*,*) xic_nue
!write(*,*) xicc_nue
!write(*,*) 't1'
call output_photprof(output_dir, 'photospheric_profile.dat', nxobs_fs, xobs*vth_fiducial/vmax, xic_nue, xicc_nue)
!
!
select case(interp_photprof)
   case(0)
      write(*,*) 'linear interpolation of photospheric profile'
   case(1)
      write(*,*) 'monotonic cubic spline interpolation of photospheric profile'
      allocate(acoeff_xic(nxobs_fs), stat=err)
         if(err.ne.0) stop 'error main: acoeff_xic'
      allocate(bcoeff_xic(nxobs_fs), stat=err)
         if(err.ne.0) stop 'error main: bcoeff_xic'
      allocate(ccoeff_xic(nxobs_fs), stat=err)
         if(err.ne.0) stop 'error main: ccoeff_xic'
      allocate(dcoeff_xic(nxobs_fs), stat=err)
         if(err.ne.0) stop 'error main: dcoeff_xic'
      call cube_mono_coeff(nxobs_fs, xobs, xic_nue, acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic)
   case default
      stop 'option interp_photprof not selected'
end select
!
!----------------------calculate gravity darkening----------------------
!
call calc_gdark
!
!------------------------setting up p-ray grid--------------------------
!
!for log-log increment
call grid_pray(.true.)
!
!for log increment
!call grid_pray(.false.)
!
!------------------------setting up zeta grid---------------------------
!
call grid_zeta
!
!------------------------setting radial grid----------------------------
!
!for log-log increment
!call grid_radial(.true.)
!
!for log increment
call grid_radial(.false.)
!
!---setting up alpha, gamma-arrays (specifying direction to observer)---
!
call grid_obsdir
!
!--------------check weights for analytical functions-------------------
!
call check_weights
!
!------calculate emergent intensity for a given direction and xobs------
!-------------(can not be included in parallel region!!!)---------------
!
if(opt_surface) then
   do indx_num = 1, nsurfb
      call calc_iem_surface(xobss_surface(indx_num), alphas_surface(indx_num), gammas_surface(indx_num))
      call output_surface(indx_num, xobss_surface(indx_num), alphas_surface(indx_num), gammas_surface(indx_num))
      write(*,*) indx_num, nsurfb
   enddo 
   stop
endif
!
!-----------calculate intensity and optical depths in 2d:---------------
!------------------slice along a given direction------------------------
!
if(opt_int2d) then
   do indx_num = 1, nsurfb
      call calc_int2d(xobss_surface(indx_num), alphas_surface(indx_num), gammas_surface(indx_num))
      call output_int2d(indx_num, xobss_surface(indx_num), alphas_surface(indx_num), gammas_surface(indx_num))
   enddo
   stop
endif
!
!--------------------allocate fluxes (not threadprivate)----------------
!
allocate(flux_tot(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_tot'
allocate(flux_cont(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_cont'
allocate(normt(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: normt'
allocate(ftot_err(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: ftot_err'
allocate(fcont_err(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: fcont_err'
allocate(flux_emi(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_emi'
allocate(flux_abs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_abs'
!   
!
indx_num=0
!
do indx_alpha=1, nalpha
   do indx_gamma=1, ngamma
!
      call timing_init
      call cpu_time(ts_obsdir)
      call system_clock(count_rate=ticks_sec, count_max=ticks_max)
      call system_clock(count=ticks_initial_obsdir)
!
!-------calculate direction to observer and transofrmation matrix-------
!-------------for given angles alpha, gamma, in order-------------------
!-----------to obtain p-rays in (x,y,z)-coordinate system---------------
!
      alpha=alpha_arr(indx_alpha)
      gamma=gamma_arr(indx_gamma)

!      alpha=27.d0*pi/180.d0!0.d0
!      gamma=13.d0*pi/180.d0
!
      write(*,*) '***************************calculating new direction***************************'
      write(*,*) 'alpha, gamma=', alpha, gamma
      write(*,*)

      call calc_transmat
!
!----------------set up a single p-ray from p-ray-grid------------------
!
      flux_tot=0.d0
      flux_cont=0.d0
      normt=0.d0
      ftot_err=0.d0
      fcont_err=0.d0
      flux_emi=0.d0
      flux_abs=0.d0
!
      !$omp parallel &
      !$omp private(indx_zeta, indx_p, indx_xobs, indx_zray, err1, err2), &
      !$omp copyin(xic1_factor)
!
!allocation of global (threadprivate) arrays
      call allocate_fluxes
!
      flux_tot_tmp=0.d0 
      flux_cont_tmp=0.d0
      flux_emi_tmp=0.d0
      flux_abs_tmp=0.d0
      normt_tmp=0.d0
      ftot_err_tmp=0.d0
      fcont_err_tmp=0.d0
      
      !$omp do schedule(dynamic)
      do indx_zeta=1, nzeta
!
         write(*,'(a30, i5, a2, i5)') 'calculating indx_zeta (nzeta)', indx_zeta, '/', nzeta
         ftot_p=0.d0
         fcont_p=0.d0
         ftot_errp=0.d0
         fcont_errp=0.d0
         normt_p=0.d0
         femi_p=0.d0
         fabs_p=0.d0
!
         do indx_p=1, np
!
            select case(input_mod)            
               case(0)
                  call setup_ray1d(zeta(indx_zeta), p(indx_p), lcore)
               case(1)
!                  call setup_ray3d_rot(zeta(indx_zeta), p(indx_p), lcore)
                  call setup_ray3d(zeta(indx_zeta), p(indx_p), lcore)
               case(2)
                  call setup_ray3d_spc(zeta(indx_zeta), p(indx_p), lcore)
               case default
                  stop 'error in spec main: input_mod not specified'
            end select
!
            do indx_xobs=1, nxobs_fs
!
!               call print_ray(xobs(indx_xobs), zeta(indx_zeta), p(indx_p), iin, iin_c, 'TRASH/testa.dat')
!
               !calculation of profile function for given xobs along ray
               do indx_zray=1, nz_ray
                  call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs), phinorm)
                  profile_ray(indx_zray)=phinorm
!                  write(*,'(10es20.8)') velz_ray(indx_zray), profile_ray(indx_zray), vth_ray(indx_zray), vth_fiducial
               enddo
!               stop 'go on bla'
!
               !setting intensity from core
               call get_iin(indx_xobs, lcore, iin, iin_c)
!               write(*,*) iin, iin_c, xic1
!               if(omp_get_thread_num().eq.1.and.indx_xobs.eq.1) write(*,'(i5,l5,5es20.8)') indx_xobs, lcore, iin, iin_c, p(indx_p)
!
               !formal solution along ray
               call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)

!               if(iem.lt.0.) then 
!                  write(*,*) p(indx_p), pw(indx_p), iem, iem_c, iabs, iemi
!                  write(*,*) indx_p, indx_zeta
!                  stop 'error: go on in spec'
!               endif
!
!p-integration
               ftot_p(indx_xobs) = ftot_p(indx_xobs) + iem * p(indx_p) * pw(indx_p)
               fcont_p(indx_xobs) = fcont_p(indx_xobs) + iem_c * p(indx_p) * pw(indx_p)
               normt_p(indx_xobs) = normt_p(indx_xobs) + 1.d0 * p(indx_p) * pw(indx_p)
               femi_p(indx_xobs) = femi_p(indx_xobs) + iemi * p(indx_p) * pw(indx_p)
               fabs_p(indx_xobs) = fabs_p(indx_xobs) + iabs * p(indx_p) * pw(indx_p)
!corresponding errors
               ftot_errp(indx_xobs) = ftot_errp(indx_xobs) + iem * p(indx_p) * pw_err(indx_p)
               fcont_errp(indx_xobs) = fcont_errp(indx_xobs) + iem_c * p(indx_p) * pw_err(indx_p)

               !debug open
!               sum=sum+iem_c * p(indx_p) * pw(indx_p)
!               write(*,*) ftot_p(indx_xobs), sum, p(indx_p)
!debug close
!
            !given p-ray is calculated for all xobs
            enddo
         !all p-rays are calculated for all xobs
         enddo
!         write(*,'(a6,3es20.8,2i5)') 'end', fcont_p(1), fcont_p(2), zeta(indx_zeta), indx_zeta, omp_get_thread_num()
!
!output error p-integration
         call calc_err(nxobs_fs, ftot_p, fcont_p, ftot_errp, fcont_errp, err1, err2)
         write(*,'(a60, 3es20.8)') 'max rel error in p-integration (tot, cont) for given zeta', zeta(indx_zeta), err1, err2
         relerr_totp=max(err1, relerr_totp)
         relerr_contp=max(err2, relerr_contp)
!
!zeta-integration
         flux_tot_tmp=flux_tot_tmp + zetaw(indx_zeta)*ftot_p
         flux_cont_tmp=flux_cont_tmp + zetaw(indx_zeta)*fcont_p
         normt_tmp=normt_tmp + zetaw(indx_zeta) * normt_p
         ftot_err_tmp=ftot_err_tmp + zetaw_err(indx_zeta) * ftot_errp
         fcont_err_tmp=fcont_err_tmp + zetaw_err(indx_zeta) * fcont_errp
         flux_emi_tmp=flux_emi_tmp + zetaw(indx_zeta)*femi_p
         flux_abs_tmp=flux_abs_tmp + zetaw(indx_zeta)*fabs_p
!debug open
!sum2=sum2+sum*zetaw(indx_zeta)
!write(*,*) flux_cont_tmp(1), sum2
!debug close
!
!      stop 'zeta0 done'
      !all zeta-angles are calculated for all p-rays and all xobs
      enddo
!      write(*,*) '1', omp_get_thread_num(), allocated(flux_tot_tmp), allocated(flux_cont_tmp)

!now: adding up all thread-private arrays
      !$omp critical
      flux_tot=flux_tot + flux_tot_tmp
      flux_cont=flux_cont + flux_cont_tmp
      normt=normt + normt_tmp
      flux_emi=flux_emi + flux_emi_tmp
      flux_abs=flux_abs + flux_abs_tmp
!corresponding errors
      ftot_err=ftot_err + ftot_err_tmp
      fcont_err=fcont_err + fcont_err_tmp
!timint
      t_setup=t_setup+t_setup_tmp
      t_hunt=t_hunt+t_hunt_tmp
      t_trilin=t_trilin+t_trilin_tmp
      !$omp end critical
!
!deallocate all arrays that have been allocated in parallel region
      call deallocate_fluxes
      call deallocate_fs1d
!
      !$omp end parallel
!
      !output error zeta-integration
      call calc_err(nxobs_fs, flux_tot, flux_cont, ftot_err, fcont_err, relerr_tot, relerr_cont)
      write(*,*)
      write(*,'(a50, 2es20.8)') 'max rel error in zeta-integration (tot, cont)', relerr_tot, relerr_cont
      write(*,*)
!
      call cpu_time(te_obsdir)
      ttot_obsdir = te_obsdir-ts_obsdir
      call system_clock(count=ticks_final_obsdir)
      ticks_obsdir=ticks_final_obsdir-ticks_initial_obsdir
      ttot_obsdir_sys = real(ticks_obsdir)/ticks_sec
!
      call print_flux
      call print_ferr
      indx_num=indx_num+1
      call output_fluxem(indx_num)
      call output_fluxem_debug(indx_num)
      call print_timing
!
   !all gamma angles are calculated
   enddo
!all alpha angles are calculated
enddo
!
call cpu_time(te_tot)
t_tot=te_tot-ts_tot
write(*,'(a40, e20.8)') 'total computation time (all angles)', t_tot
write(*,*) ''
!
!
!
end program spectrum3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_xobs
!
!-----------------calculates frequency integration nodes----------------
!
use prog_type
use fund_const, only: cgs_clight, cgs_kb, cgs_mp
use dime_spec, only: nxobs_fs
use mod_spectrum, only: del_xobs, xobs, xnue, xic_nue, xicc_nue
use params_spec, only: xnue0, vmax, vth_fiducial, vth_min, xlim
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: xmax
real(dp) :: deldop_fiducial, delta, dum_dxobs
!
!calculate ratio of actual (minimum) thermal velocity and scaling used in code
delta=vth_min/vth_fiducial
!
!calculate xmax (corresponding to vmax/vth_min + xlim, but in units of vth_fiducial)
xmax = vmax/vth_fiducial + xlim*delta
!
!calculate fiducial doppler width
deldop_fiducial = xnue0 * vth_fiducial/cgs_clight
!
!calculate number of needed frequency points 
!   note: (nue-nue0)/deldop_min = del_xobs
!         nxobs_fs needs to be odd in order to include xobs=0
dum_dxobs=del_xobs*delta
nxobs_fs=2*nint(xmax/dum_dxobs) + 1
!
!write(*,*) nxobs_fs
!write(*,*) vmax
!write(*,*) xmax
!write(*,*) xlim
!write(*,*) delta, vth_min, vth_fiducial
!write(*,*) 
!stop 'go on in spec xobs'
!
!allocate xobs-array, xnue-array, xic-nue-array
allocate(xobs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xobs'
allocate(xnue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xnue'
allocate(xic_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xic_nue'
allocate(xicc_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xicc_nue'
!
!calculate observers frame frequency grid
!
do i=1, nxobs_fs
   xobs(i) = -xmax + (i-1)*2*xmax/(nxobs_fs-1)
   xnue(i) = xobs(i)*deldop_fiducial + xnue0
!   write(*,*) i, xobs(i)
enddo
!
!-------------------------print out xobs-grid---------------------------
!
write(*,*) '-----------------------------calculating xobs-grid-----------------------------'
write(*,*)
!
write(*,*) 'xobs-grid:'
write(*,'(8es20.8)') xobs
write(*,*)
!
!
end subroutine grid_xobs
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_zeta
!
!-----------calculates polar angle integration nodes and weights--------
!
use prog_type
use fund_const
use dime_spec, only: nzeta, nzeta1, ni
use mod_spectrum, only: zeta, zetaw, zetaw1, zetaw_err
use mod_integ1d, only: precalc_weight_simps_err
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
real(dp) :: del, sum, sum1, sum_err, sumex
!
! ... local arrays
real(dp), dimension(:), allocatable :: dum_zeta
!
!------------------create grid with nzeta1 grid points------------------
!-------------------(do not have to be equidistant)---------------------
!
allocate(dum_zeta(nzeta1), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: dum_zeta'
!
if(nzeta.lt.2) stop 'error in grid_zeta: nzeta1 has to be ge 2'
!
do i=1, nzeta1
   dum_zeta(i) = 2.d0*pi*float(i-1)/(nzeta1-1)
!   write(*,*) dum_zeta(i)
enddo
!
!------------------insert ni equidistant subintervals-------------------
!
!allocate zeta-, zetaw-, zetaw1 and zetaw_err - arrays
allocate(zeta(nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: zeta'
allocate(zetaw(nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: zetaw'
allocate(zetaw1(nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: zetaw1'
allocate(zetaw_err(nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: zetaw_err'
!
k=1
!
do i=1, nzeta1-1
   del=(dum_zeta(i+1)-dum_zeta(i))/ni
   do j=1, ni
      zeta(k+j-1) = dum_zeta(i) + float(j-1)*del
   enddo
   k=k+ni
enddo
zeta(nzeta)=dum_zeta(nzeta1)
!
!-----------------------------------------------------------------------
!
!calculating weights for trapezoidal rule including error weights
!  (need to have at least 2 subsequent equidistant subintervals)
!call precalc_weight_trapez_err(zeta, nzeta, ni, zetaw, zetaw1, zetaw_err)
!
!calculating weights for simpson's rule including error weights
!  (need to have at least 4 subsequent equidistant subintervals)
call precalc_weight_simps_err(zeta, nzeta, ni, zetaw, zetaw1, zetaw_err)
!
!calculating weights (trapez or simpson - integration) without error weights
!call precalc_weight_trapez(zeta, nzeta, zetaw)
!call precalc_weight_simps(zeta, nzeta, zetaw)
!
!------------------------test integration-------------------------------
!
sumex=2.d0*pi
sum=0.d0
sum1=0.d0
sum_err=0.d0
do i=1, nzeta
   sum=sum+zetaw(i)
   sum1=sum1+zetaw1(i)
   sum_err=sum_err+zetaw_err(i)
enddo
!
if(abs(sumex-sum)/sumex.gt.1.d-14) then
   write(*,*) '------------error in zeta integration-----------'
   write(*,'(6a30)') 'integral(num)', 'integral(num), half step', 'integral(exact)', 'abs error_1(num)', 'abs error_2(num)', 'error(exact)'
   write(*,'(6(e30.8))') sum, sum1, sumex, sum_err, (sum-sum1)/6.d0, (sumex-sum)
!error_2(num) is (sum-sum1)/6.d0  for trapezoidal rule and
!                (sum-sum1)/15.d0 for simpson's rule
   stop
endif
!
!
!
end subroutine grid_zeta
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_pray(lloglog)
!
!-----------------sets up grid of impact parameter----------------------
!
use prog_type
use dime_spec, only: ni, np1_nc, np1_c, np
use mod_spectrum, only: p, pw, pw1, pw_err, rmax, hmax
use mod_integ1d, only: precalc_weight_simps_err
!
implicit none
!
! ... arguments
logical, intent(in) :: lloglog
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: np_dum
real(dp) :: del, sum, sum1, sum_err, sumex
!
! ... local arrays
real(dp), dimension(:), allocatable :: p_dum
!
!********************create grid with unequal steps*********************
!
np_dum = np1_nc + np1_c
allocate(p_dum(np_dum), stat=err)
   if(err.ne.0) stop 'allocation error grid_pray: p_dum'
p_dum=0.d0
!
!--------------------equidistant inside core----------------------------
!
if(np1_c.lt.2) stop 'error grid_pray: np1_c has to be ge 2'
!
do i=1, np1_c
   p_dum(i)=1.d0*float(i-1)/float(np1_c-1)
enddo
!
!------------------------outside core-----------------------------------
!
if(np1_nc.lt.2) stop 'error in grid_pray: np1_nc has to be ge 2'
!
!additional point right above photosphere
p_dum(np1_c+1)=p_dum(np1_c) + 1.d-4
!
if(.not.lloglog) then
!
!logarithmic increment
   del = log10(rmax/p_dum(np1_c+1))/float(np1_nc-1)
!
   do i=2, np1_nc
      p_dum(np1_c+i)=p_dum(np1_c+i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax
   p_dum(np1_nc+np1_c)=rmax
!
else
!
!log-log-increment
   del = log10(rmax)/log10(p_dum(np1_c+1))
   del = log10(del)
   del = del/(np1_nc-1)
!
   do i=2, np1_nc
      p_dum(np1_c+i) = 10.d0**(log10(p_dum(np1_c+i-1))*10.d0**del)
   enddo
   p_dum(np1_nc+np1_c)=rmax
endif
!
!-------------------include equidistant subintervals--------------------
!------------in order to calculate error estimation weights-------------
!
allocate(p(np), stat=err)
   if(err.ne.0) stop 'allocation error grid_pray: p'
!
allocate(pw(np), stat=err)
   if(err.ne.0) stop 'allocation error grid_pray: pw'
!
allocate(pw1(np), stat=err)
   if(err.ne.0) stop 'allocation error grid_pray: pw1'
!
allocate(pw_err(np), stat=err)
   if(err.ne.0) stop 'allocation error grid_pray: pw_err'
!
k=1
!
hmax=0.d0
!
do i=1, np1_c+np1_nc-1
!
   del=(p_dum(i+1)-p_dum(i))/ni
!
   if(del.gt.hmax) then
      hmax=del
   endif
!
   do j=1, ni
      p(k+j-1) = p_dum(i) + float(j-1)*del
   enddo
   k=k+ni
enddo
p(np)=p_dum(np1_c+np1_nc)
!
!-------------------------print out p-grid------------------------------
!
write(*,*) '-----------------------------calculating p-grid--------------------------------'
write(*,*)
!
write(*,*) 'p-grid:'
write(*,'(8es20.8)') p
write(*,*)
!
!-----------------------------------------------------------------------
!
!calculating weights for trapezoidal rule including error weights
!  (need to have at least 2 subsequent equidistant subintervals)
!call precalc_weight_trapez_err(p, np, ni, pw, pw1, pw_err)
!
!calculating weights for simpson's rule including error weights
!  (need to have at least 4 subsequent equidistant subintervals)
call precalc_weight_simps_err(p, np, ni, pw, pw1, pw_err)
!
!calculating weights for trapezoidal rule
!   (no need of equidistant subintervals since no error weights are calculated)
!call precalc_weight_trapez(p, np, pw)
!
!------------------------test integration-------------------------------
!
!note factor of 1/2
sumex=rmax*rmax/2.d0
sum=0.d0
sum1=0.d0
sum_err=0.d0
do i=1, np
   sum=sum+p(i)*pw(i)
   sum1=sum1+p(i)*pw1(i)
   sum_err=sum_err+p(i)*pw_err(i)
enddo
!
!
if(abs(sumex-sum)/sumex.gt.1.d-14) then
   write(*,*) '--------------error in p-integration-------------'
   write(*,'(6a30)') 'integral(num)', 'integral(num), half step', 'integral(exact)', 'abs error_1(num)', 'abs error_2(num)', 'error(exact)'
   write(*,'(6(e30.8))') sum, sum1, sumex, sum_err, (sum-sum1)/6.d0, (sumex-sum)
!error_2(num) is (sum-sum1)/6.d0  for trapezoidal rule and
!                (sum-sum1)/15.d0 for simpson's rule
   stop
endif
!
!
!
end subroutine grid_pray
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_radial(lloglog)
!
!----------set up radial grid for calculation of z-coordinates----------
!
use prog_type
use dime_spec, only: nr
use mod_spectrum, only: r, rmax
!
implicit none
!
! ... arguments
logical, intent(in) :: lloglog
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: del
!
!-----------------------------------------------------------------------
!
write(*,*) '-----------------------------calculating r-grid--------------------------------'
write(*,*)
!
if(allocated(r)) deallocate(r)
allocate(r(nr), stat=err)
   if(err.ne.0) stop 'allocation error grid_radial: r'
!
!-----------------------------------------------------------------------
!
if(.not.lloglog) then
!
!---------------equidistant in log-space outside core-------------------
!
!logarithmic increment
   del = log10(rmax/1.d0)/(nr)
!
   r(1)=1.d0
   do i=2, nr
      r(i)=r(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax
   r(nr)=rmax
!
else
!
!---------------------or equidsitant in log-logs------------------------
!
   r(1)=1.d0
   r(2)=1.d0+1.d-3
!
   del = log10(rmax)/log10(r(2))
   del = log10(del)
   del = del/(nr-2)
!
   do i=3, nr
      r(i) = 10.d0**(log10(r(i-1))*10.d0**del)
   enddo
   r(nr)=rmax
!
endif
!
write(*,*) 'r-grid:'
write(*,'(8es20.8)') r
!
!
end subroutine grid_radial
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_radial_rot(lloglog, rmin, r, nr)
!
!----------set up radial grid for calculation of z-coordinates----------
!-----from r_min to r_max, to account for the distortion of the---------
!----stellar surface when large rotational velocities are present-------
!
use prog_type
use mod_spectrum, only: rmax
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nr
real(dp), intent(in) :: rmin
logical, intent(in) :: lloglog
real(dp), dimension(nr) :: r
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: del
!
!-----------------------------------------------------------------------
!
!write(*,*) '-----------------------------calculating r-grid--------------------------------'
!write(*,*)
!
!-----------------------------------------------------------------------
!
if(.not.lloglog) then
!
!---------------equidistant in log-space outside core-------------------
!
!logarithmic increment
   del = log10(rmax/rmin)/(nr)
!
   r(1)=rmin
   do i=2, nr
      r(i)=r(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax
   r(nr)=rmax
!
else
!
!---------------------or equidsitant in log-logs------------------------
!
   r(1)=rmin
   r(2)=rmin+1.d-3
!
   del = log10(rmax)/log10(r(2))
   del = log10(del)
   del = del/(nr-2)
!
   do i=3, nr
      r(i) = 10.d0**(log10(r(i-1))*10.d0**del)
   enddo
   r(nr)=rmax
!
endif
!
!
end subroutine grid_radial_rot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_obsdir
!
!----------set up alpha, gamma arrays, for which formal solution--------
!--shall be performed (alpha, gamma specify direction towards observer)-
!
use prog_type
use fund_const, only: pi
use options_spec, only: opt_obsdir_read
use dime_spec, only: nalpha, ngamma
use mod_spectrum, only: alpha_arr, gamma_arr
!
implicit none
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: alpha_start, alpha_end
!
!
write(*,*) '----------------------------observers directions-------------------------------'
write(*,*)
!
!alpha_start=48.59d0*pi/180.d0
alpha_start=pi/2.d0
alpha_end=0.d0
!
!-----------------------------------------------------------------------
!
err=0
!
if(allocated(alpha_arr)) deallocate(alpha_arr, stat=err)
   if(err.ne.0) stop 'deallocation error grid_obsidr: alpha_arr'
if(allocated(gamma_arr)) deallocate(gamma_arr, stat=err)
   if(err.ne.0) stop 'deallocation error grid_obsidr: gamma_arr'
!
allocate(alpha_arr(nalpha), stat=err)
   if(err.ne.0) stop 'allocation error grid_obsdir: alpha_arr'
allocate(gamma_arr(ngamma), stat=err)
   if(err.ne.0) stop 'allocation error grid_obsdir: gamma_arr'
!
!---------------equidistant in log-space outside core-------------------
!
if(opt_obsdir_read) then
   open(1, file='in_alpha.dat', form='formatted')
      do i=1, nalpha
         read(1,*) alpha_arr(i)
      enddo
   close(1)
   open(1, file='in_gamma.dat', form='formatted')
      do i=1, ngamma
         read(1,*) gamma_arr(i)
      enddo
   close(1)
!transform to rad
   alpha_arr=alpha_arr*pi/180.d0
   gamma_arr=gamma_arr*pi/180.d0
else
   if(ngamma.eq.1) then
      !default: gamma=0.d0
      gamma_arr(1)=0.d0
   else
      do i=1, ngamma
         gamma_arr(i)=(i-1.d0)*2.d0*pi/(ngamma-1)
      enddo
   endif

   if(nalpha.eq.1) then
      !default: alpha_start
      alpha_arr(1)=alpha_start
   else

!for different angles
      do i=1, nalpha
!for constant alpha
         alpha_arr(i) = alpha_start + (i-1)*(alpha_end-alpha_start)/(nalpha-1)
!for constan sin(alpha)
!         alpha_arr(i) = asin(1.d0-float(i-1)/float(nalpha-1))
      enddo

   endif

endif
!
!
!
end subroutine grid_obsdir
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d(zeta, p, lcore)
!
!------sets up a p-ray for a given observers frame frequency xobs-------
!to calculate:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: zeta   angle-coordinate of cylindrical system
!       p      radial distance of point in cylindrical system
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
!attention: v_x, v_y, v_z need to be interpolated linearly
!              (otherwise: when rotational veloctities are included, 
!               extrapolation may overestimate correct velocities by a factor of 1000 or more!!!!)
!
use prog_type
use dime_spec, only:  nr
use dime_model3d, only: ndxmax, ndymax, ndzmax, &
                        x, y, z, r3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                        vth_ray, temp_ray,  profile_ray, velz_ray, r, transmat, rotmat_trans, nhat, &
                        del_vel2, rmin, vphot_proj, rmax
use params_spec, only: vth_fiducial, vth_min, sr, xic1, vmin, vmax, tmin, na, vmicro
use timing_spec, only: ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, t_hunt_tmp, &
                       ts_trilin, te_trilin, t_trilin_tmp
use mod_interp1d, only: interpol_yp
use mod_interp3d, only: get_xyz_indx, get_xyz_values1, get_xyz_values2, trilin_complete
!
implicit none
!
! ... arguments
real(dp), intent(in) :: zeta, p
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, indx
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_vmicro, dum_opac, dum_scont 
real(dp) :: x1, x2, y1, y2, z1, z2
real(dp) :: rada, radb, radc, radd, rade, radf, radg, radh, radp
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velx, vely, velz, velr, vel_abs, rad
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc, vec_vel
real(dp), dimension(:), allocatable :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(:), allocatable :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogx, llogy, llogz, llogf, lr2, ldum
!
! ... local functions
real(dp) :: calc_vmicro, vthermal
logical :: boundary
!
!
!
call cpu_time(ts_setup)
!
!might need to double check allocations
allocate(zdum_ray(nz_max), veldum_ray(nz_max), vthdum_ray(nz_max), opacdum_ray(nz_max), &
         opalbardum_ray(nz_max), scontdum_ray(nz_max), slinedum_ray(nz_max), tempdum_ray(nz_max))

allocate(zdum2_ray(nz_max), veldum2_ray(nz_max), vthdum2_ray(nz_max), opacdum2_ray(nz_max), &
         opalbardum2_ray(nz_max), scontdum2_ray(nz_max), slinedum2_ray(nz_max), tempdum2_ray(nz_max))
!
!
!
!
iz=0
!
do i=nr, 1, -1
   zdum=r(i)*r(i)-p*p
   if(zdum.gt.0.d0) then
      iz=iz+1
      zdum_ray(iz)=sqrt(zdum)
   else
      iz=iz+1
      zdum_ray(iz)=0.d0
      exit
   endif
enddo
!
!
!check if core ray 
lcore=.false.
!inner-most point in carthesian coordinates
vec_cyc(1)=p*cos(zeta)
vec_cyc(2)=p*sin(zeta)
vec_cyc(3)=zdum_ray(iz)
vec_cac=matmul(transmat, vec_cyc)
!
!note: function can in general calculate arbitrary boundary shapes, here: spherical
lcore=boundary(vec_cac)
!
!
if(.not.lcore) then
!set non core rays
   iz_dum=2*iz-1
   do i=1, iz-1
      zdum_ray(iz+i)=-zdum_ray(iz-i)
   enddo
   iz=iz_dum
   vphot_proj=0.d0
else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
   rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
   call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz)
   vphot_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
endif
!
!-----------------------------------------------------------------------
!
do i=1, iz
!
!calculate z_ray in carthesian coordinates
   vec_cyc(1)=p*cos(zeta)
   vec_cyc(2)=p*sin(zeta)
   vec_cyc(3)=zdum_ray(i)
   vec_cac=matmul(transmat, vec_cyc)
!
!check if point of ray lies within region where information is stored
   call info_region(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)
!
!interpolation only if point lies within region where information is stored
   if(linfo_phot.and.linfo_max) then
!
!timing for setting up interpolation-parameter
      call cpu_time(ts_hunt)
!
!search for indices of a cube surrounding the point of interest (or neighbouring cubes for extrapolation)
      call get_xyz_indx(vec_cac(1), vec_cac(2), vec_cac(3), x, y, z, ndxmax, ndymax, ndzmax, &
                        indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol, rmin, rmax)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
      call get_xyz_values1(vec_cac(1), vec_cac(2), vec_cac(3), x, y, z, ndxmax, ndymax, ndzmax, &
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                          x1, x2, y1, y2, z1, z2, &
                          rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                          llogx, llogy, llogz, expol)
!timing for setting up interpolation-parameter
      call cpu_time(te_hunt)
      t_hunt_tmp=t_hunt_tmp + te_hunt - ts_hunt
!
!timing for complete trilinear interpolation
      call cpu_time(ts_trilin)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, velx3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,*) 'at z', zdum_ray(i)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, velx)
!
!get vely on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, vely3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, vely)

!get velz on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, velz3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, velz)
!
!transform to formal-solver carthesian coordinates
      vec_vel=matmul(rotmat_trans, (/ velx, vely, velz /))
!
!calculation of velocity projected onto ray
      veldum_ray(i)=nhat(1)*vec_vel(1) + nhat(2)*vec_vel(2) + nhat(3)*vec_vel(3)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, opalbar3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!here, can apply only lin-lin interpolation by setting llogx, llogy, llogz, llogf=.false.
!      llogx=.false.
!      llogy=.false.
!      llogz=.false.
!      llogf=.false.
!here, can decide if values shall be interpolated by function*r^2
      lr2=.true.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
!line opacity in units of 1/rstar
      opalbardum_ray(i) = valp
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, sline3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogx=.false.
!      llogy=.false.
!      llogz=.false.
!      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
      slinedum_ray(i)=valp
!
!-------------------------thermal velocities----------------------------
!
!get thermal velocity on cube vertices
!      call get_xyz_values2(ndxmax, ndymax, ndzmax, vth3d, &
!                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
!                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      lr2=.false.
!!actual interpolation
!      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
!                           x1, x2, y1, y2, z1, z2, &
!                           vala, valb, valc, vald, vale, valf, valg, valh, &
!                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
!                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
!      vthdum_ray(i)=valp
!
!calculation of absolute value of velocity (needed to calculate vmicro)
      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial

!note: minimum microturbulent velocity set to input-vmicro
      dum_vmicro = calc_vmicro(vmicro, vmax, vel_abs)
      dum_vmicro = vmicro
!
!calculate corresponding thermal velocity (for a fixed tmin)
      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!
!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
      tempdum_ray(i)=0.d0
!
!-----------------------------------------------------------------------
!
!for now: without continuum
      scontdum_ray(i)=0.d0
      opacdum_ray(i)=0.d0
!
      call cpu_time(te_trilin)
      t_trilin_tmp=t_trilin_tmp + te_trilin - ts_trilin
!
!------set all quantities to zero, if outside information region--------
!
   else
      opalbardum_ray(i) = 0.d0
      slinedum_ray(i) = 0.d0
      veldum_ray(i) = 0.d0
      vthdum_ray(i) = 1.d-8
      tempdum_ray(i) = 0.d0
      scontdum_ray(i) = 0.d0
      opacdum_ray(i) = 0.d0
!
   endif
!
enddo
!
!-----------------------------------------------------------------------
!
t_trilin_tmp=t_trilin_tmp/(iz*1.d0)
t_hunt_tmp=t_hunt_tmp/(iz*1.d0)
!
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!
!velocity at beginning point
dum_vel1=veldum_ray(1)
!
veldum2_ray(1)=veldum_ray(1)
zdum2_ray(1)=zdum_ray(1)
opacdum2_ray(1)=opacdum_ray(1)
opalbardum2_ray(1)=opalbardum_ray(1)
scontdum2_ray(1)=scontdum_ray(1)
slinedum2_ray(1)=slinedum_ray(1)
tempdum2_ray(1)=tempdum_ray(1)
vthdum2_ray(1)=vthdum_ray(1)
!
k=0
iz_dum=iz
!
do i=2, iz
!
   dum_vel2=veldum_ray(i)
!
!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!   if(dum_vel1*dum_vel2.lt.0.d0) then
!      write(*,'(a92)') 'error in setup_ray: velocity jump along ray from positive velocities <-> negative velocities'
!      write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!      write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!      write(*,'(a82)') '   => 2. use an analytic description of the velocity law, instead of interpolation'
!      stop
!   endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up. note: ceiling is only f90
      nadd=ceiling(abs(dum_delv)/del_vel2)

!      write(*,*) dum_delv, del_vel2, nadd

      do j=1, nadd-1
!add nadd-1 additional points in velocity space
!         write(*,*) 'test', i-1, j, k, i+k, nadd
!         write(*,*) veldum_ray(i-1), veldum_ray(i)
!         write(*,*) veldum2_ray(i+k)
         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/nadd
!interpolation of all variables
         zdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
         opacdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
         opalbardum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
         scontdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
         slinedum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
         tempdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
         vthdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
      enddo
   endif

   veldum2_ray(i+k)=veldum_ray(i)
   zdum2_ray(i+k)=zdum_ray(i)
   opacdum2_ray(i+k)=opacdum_ray(i)
   opalbardum2_ray(i+k)=opalbardum_ray(i)
   scontdum2_ray(i+k)=scontdum_ray(i)
   slinedum2_ray(i+k)=slinedum_ray(i)
   tempdum2_ray(i+k)=tempdum_ray(i)
   vthdum2_ray(i+k)=vthdum_ray(i)

   dum_vel1=dum_vel2
! 
enddo
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum2_ray(nz_ray+1-i)
   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
   velz_ray(i) = veldum2_ray(nz_ray+1-i)

!  write(*,'(8es20.8)') z_ray(i), opalbar_ray(i), opac_ray(i), scont_ray(i), sline_ray(i), vth_ray(i), temp_ray(i), velz_ray(i)

enddo
!
!
!
call cpu_time(te_setup)
t_setup_tmp = t_setup_tmp + te_setup - ts_setup
!
!if(abs(p-0.25d0).lt.1.d-12) then
!   open(1, file='trash/vz_along.dat', form='formatted', position='append')
!      do i=1, nz_ray
!         write(1,'(4es20.8)') z_ray(i), opalbar_ray(i), sline_ray(i), velz_ray(i)
!      enddo
!   close(1)
!   read(*,*)
!endif
!
end subroutine setup_ray3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d_rot(zeta, p, lcore)
!
!------sets up a p-ray for a given observers frame frequency xobs-------
!to calculate:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: zeta   angle-coordinate of cylindrical system
!       p      radial distance of point in cylindrical system
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
!attention: v_x, v_y, v_z need to be interpolated linearly
!              (otherwise: when rotational veloctities are included, 
!               extrapolation may overestimate correct velocities by a factor of 1000 or more!!!!)
!
!this routine is the same routine as setup_ray3d, however accounts for the 
!rotational distortion of the stellar surface
!the stellar surface is approximated by an ellpsoid
!
use prog_type
use fund_const, only: cgs_grav, xmsu, rsu
use dime_spec, only:  nr
use dime_model3d, only: ndxmax, ndymax, ndzmax, &
                        x, y, z, r3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                        vth_ray, temp_ray,  profile_ray, velz_ray, r, transmat, rotmat_trans, nhat, &
                        del_vel2, rmin, vphot_proj, rmax
use params_spec, only: vth_fiducial, vth_min, sr, xic1, vmin, vmax, tmin, na, vmicro, vrot, rstar, xlogg, mstar
use timing_spec, only: ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, t_hunt_tmp, &
                       ts_trilin, te_trilin, t_trilin_tmp
use mod_gdark, only: smajorax_a, smajorax_b, smajorax_c, xic1_factor
use mod_interp1d, only: interpol_yp
use mod_interp3d, only: get_xyz_indx, get_xyz_values1, get_xyz_values2, trilin_complete
!
implicit none
!
! ... arguments
real(dp), intent(in) :: zeta, p
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, indx
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_vmicro, dum_opac, dum_scont 
real(dp) :: x1, x2, y1, y2, z1, z2
real(dp) :: rada, radb, radc, radd, rade, radf, radg, radh, radp
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velx, vely, velz, velr, vel_abs, rad
real(dp) :: m11, m12, m13, m21, m22, m23, m31, m32, m33, x_cyc, y_cyc, z_cyc, z_cyc1, z_cyc2, &
            fdum1, fdum2, fdum3, acoeff, bcoeff, ccoeff, disc, &
            r_boundary
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc, vec_vel
real(dp), dimension(:), allocatable :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(:), allocatable :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
real(dp), dimension(:), allocatable :: r_dum
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogx, llogy, llogz, llogf, lr2, ldum, linfo_phot2, linfo_max2
!
! ... local functions
real(dp) :: calc_vmicro, vthermal
logical :: boundary
!
!
!
call cpu_time(ts_setup)
!
!might need to double check allocations
allocate(zdum_ray(nz_max), veldum_ray(nz_max), vthdum_ray(nz_max), opacdum_ray(nz_max), &
         opalbardum_ray(nz_max), scontdum_ray(nz_max), slinedum_ray(nz_max), tempdum_ray(nz_max))

allocate(zdum2_ray(nz_max), veldum2_ray(nz_max), vthdum2_ray(nz_max), opacdum2_ray(nz_max), &
         opalbardum2_ray(nz_max), scontdum2_ray(nz_max), slinedum2_ray(nz_max), tempdum2_ray(nz_max))
!
!-----------calculate start-value of z-coordinate along ray-------------
!-----------------------for distorted surface---------------------------
!
allocate(r_dum(nr), stat=err)
   if(err.ne.0) stop 'error in setup_ray3d_rot: allocation'
!
!x,y-coordinates in cylindrical system
x_cyc=p*cos(zeta)
y_cyc=p*sin(zeta)
!
m11=transmat(1,1)
m12=transmat(1,2)
m13=transmat(1,3)
m21=transmat(2,1)
m22=transmat(2,2)
m23=transmat(2,3)
m31=transmat(3,1)
m32=transmat(3,2)
m33=transmat(3,3)
!
fdum1 = x_cyc**2*(m11**2/smajorax_a**2 + m21**2/smajorax_b**2 + m31**2/smajorax_c**2)
fdum2 = y_cyc**2*(m12**2/smajorax_a**2 + m22**2/smajorax_b**2 + m32**2/smajorax_c**2)
fdum3 = 2.d0*x_cyc*y_cyc*(m11*m12/smajorax_a**2 + m21*m22/smajorax_b**2 + m31*m32/smajorax_c**2)
ccoeff=fdum1+fdum2+fdum3-1.d0
!
fdum1 = x_cyc*(m11*m13/smajorax_a**2 + m21*m23/smajorax_b**2 + m31*m33/smajorax_c**2) + &
        y_cyc*(m12*m13/smajorax_a**2 + m22*m23/smajorax_b**2 + m32*m33/smajorax_c**2)
bcoeff=2.d0*fdum1
!
acoeff = m13**2/smajorax_a**2 + m23**2/smajorax_b**2 + m33**2/smajorax_c**2
!
!solve quadratic equation acoeff*z**2 + bcoeff*z + ccoeff = 0
disc = bcoeff**2-4.d0*acoeff*ccoeff
if(disc.lt.0.d0) then
   lcore=.false.
   r_dum = r
else
   z_cyc1 = (-bcoeff + sqrt(disc))/2.d0/acoeff
   z_cyc2 = (-bcoeff - sqrt(disc))/2.d0/acoeff
   z_cyc = max(z_cyc1,z_cyc2)
   lcore = .true.
   r_boundary = sqrt(p**2+z_cyc**2)
   call grid_radial_rot(.true., r_boundary, r_dum, nr)
!   write(*,*) 'p, zeta', p, zeta
!   write(*,*) 'lcore', lcore
!   write(*,*) 'z_cyc', z_cyc
!   write(*,*) smajorax_a, smajorax_b, smajorax_c, r_boundary, p, z_cyc
!   write(*,'(10f10.4)') r_dum
!   stop
endif
!
!
iz=0
!
do i=nr, 1, -1
   zdum=r_dum(i)*r_dum(i)-p*p
   if(zdum.gt.0.d0) then
      iz=iz+1
      zdum_ray(iz)=sqrt(zdum)
   else
      iz=iz+1
      zdum_ray(iz)=0.d0
      exit
   endif
enddo
!
!
vec_cyc(1)=p*cos(zeta)
vec_cyc(2)=p*sin(zeta)
vec_cyc(3)=zdum_ray(iz)
vec_cac=matmul(transmat, vec_cyc)
!
!
!
if(.not.lcore) then
!set non core rays
   iz_dum=2*iz-1
   do i=1, iz-1
      zdum_ray(iz+i)=-zdum_ray(iz-i)
   enddo
   iz=iz_dum
   vphot_proj=0.d0
else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
   rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
   call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz)
   vphot_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
!calculation of local effective temperature and scaling factor for global profile array
!   (assuming that the photospheric profile is constant over the stellar surface!!!)
   call calc_xic_factor(vec_cac(1),vec_cac(2),vec_cac(3))
endif

!write(*,'(15f10.5)') zdum_ray
!
!-----------------------------------------------------------------------
!
do i=1, iz
!
!calculate z_ray in carthesian coordinates
   vec_cyc(1)=p*cos(zeta)
   vec_cyc(2)=p*sin(zeta)
   vec_cyc(3)=zdum_ray(i)
   vec_cac=matmul(transmat, vec_cyc)
!
!check if point of ray lies within region where information is stored
!   call info_region(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)
   call info_region_rot(vec_cac(1), vec_cac(2), vec_cac(3), smajorax_a, smajorax_b, smajorax_c, rmax, linfo_phot, linfo_max, ldum)
!
!interpolation only if point lies within region where information is stored
   if(linfo_phot.and.linfo_max) then
!
!timing for setting up interpolation-parameter
      call cpu_time(ts_hunt)
!
!search for indices of a cube surrounding the point of interest (or neighbouring cubes for extrapolation)
      call get_xyz_indx(vec_cac(1), vec_cac(2), vec_cac(3), x, y, z, ndxmax, ndymax, ndzmax, &
                        indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol, rmin, rmax)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
      call get_xyz_values1(vec_cac(1), vec_cac(2), vec_cac(3), x, y, z, ndxmax, ndymax, ndzmax, &
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                          x1, x2, y1, y2, z1, z2, &
                          rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                          llogx, llogy, llogz, expol)
!timing for setting up interpolation-parameter
      call cpu_time(te_hunt)
      t_hunt_tmp=t_hunt_tmp + te_hunt - ts_hunt
!
!timing for complete trilinear interpolation
      call cpu_time(ts_trilin)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, velx3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,*) 'at z', zdum_ray(i)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, velx)
!
!get vely on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, vely3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, vely)

!get velz on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, velz3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      write(*,'(8es20.8)') vala, valb, valc, vald, vale, valf, valg, valh
      llogx=.false.
      llogy=.false.
      llogz=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, velz)
!
!transform to formal-solver carthesian coordinates
      vec_vel=matmul(rotmat_trans, (/ velx, vely, velz /))
!
!calculation of velocity projected onto ray
      veldum_ray(i)=nhat(1)*vec_vel(1) + nhat(2)*vec_vel(2) + nhat(3)*vec_vel(3)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, opalbar3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!here, can apply only lin-lin interpolation by setting llogx, llogy, llogz, llogf=.false.
!      llogx=.false.
!      llogy=.false.
!      llogz=.false.
!      llogf=.false.
!here, can decide if values shall be interpolated by function*r^2
      lr2=.true.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
!line opacity in units of 1/rstar
      opalbardum_ray(i) = valp
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
      call get_xyz_values2(ndxmax, ndymax, ndzmax, sline3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogx=.false.
!      llogy=.false.
!      llogz=.false.
!      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
                           x1, x2, y1, y2, z1, z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
      slinedum_ray(i)=valp
!
!-------------------------thermal velocities----------------------------
!
!get thermal velocity on cube vertices
!      call get_xyz_values2(ndxmax, ndymax, ndzmax, vth3d, &
!                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
!                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      lr2=.false.
!!actual interpolation
!      call trilin_complete(vec_cac(1), vec_cac(2), vec_cac(3), &
!                           x1, x2, y1, y2, z1, z2, &
!                           vala, valb, valc, vald, vale, valf, valg, valh, &
!                           rada, radb, radc, radd, rade, radf, radg, radh, radp, &
!                           expol, .false., llogx, llogy, llogz, llogf, lr2, valp)
!      vthdum_ray(i)=valp
!
!calculation of absolute value of velocity (needed to calculate vmicro)
      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial

!note: minimum microturbulent velocity set to input-vmicro
      dum_vmicro = calc_vmicro(vmicro, vmax, vel_abs)
      dum_vmicro = vmicro
!
!calculate corresponding thermal velocity (for a fixed tmin)
      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!
!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
      tempdum_ray(i)=0.d0
!
!-----------------------------------------------------------------------
!
!for now: without continuum
      scontdum_ray(i)=0.d0
      opacdum_ray(i)=0.d0
!
      call cpu_time(te_trilin)
      t_trilin_tmp=t_trilin_tmp + te_trilin - ts_trilin
!
!------set all quantities to zero, if outside information region--------
!
   else
      opalbardum_ray(i) = 0.d0
      slinedum_ray(i) = 0.d0
      veldum_ray(i) = 0.d0
      vthdum_ray(i) = 1.d-8
      tempdum_ray(i) = 0.d0
      scontdum_ray(i) = 0.d0
      opacdum_ray(i) = 0.d0
!
   endif
!
enddo
!
!-----------------------------------------------------------------------
!
t_trilin_tmp=t_trilin_tmp/(iz*1.d0)
t_hunt_tmp=t_hunt_tmp/(iz*1.d0)
!
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!
!velocity at beginning point
dum_vel1=veldum_ray(1)
!
veldum2_ray(1)=veldum_ray(1)
zdum2_ray(1)=zdum_ray(1)
opacdum2_ray(1)=opacdum_ray(1)
opalbardum2_ray(1)=opalbardum_ray(1)
scontdum2_ray(1)=scontdum_ray(1)
slinedum2_ray(1)=slinedum_ray(1)
tempdum2_ray(1)=tempdum_ray(1)
vthdum2_ray(1)=vthdum_ray(1)
!
k=0
iz_dum=iz
!
do i=2, iz
!
   dum_vel2=veldum_ray(i)
!
!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!   if(dum_vel1*dum_vel2.lt.0.d0) then
!      write(*,'(a92)') 'error in setup_ray: velocity jump along ray from positive velocities <-> negative velocities'
!      write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!      write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!      write(*,'(a82)') '   => 2. use an analytic description of the velocity law, instead of interpolation'
!      stop
!   endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up. note: ceiling is only f90
      nadd=ceiling(abs(dum_delv)/del_vel2)

!      write(*,*) dum_delv, del_vel2, nadd

      do j=1, nadd-1
!add nadd-1 additional points in velocity space
!         write(*,*) 'test', i-1, j, k, i+k, nadd
!         write(*,*) veldum_ray(i-1), veldum_ray(i)
!         write(*,*) veldum2_ray(i+k)
         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/nadd
!interpolation of all variables
         zdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
         opacdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
         opalbardum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
         scontdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
         slinedum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
         tempdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
         vthdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
      enddo
   endif

   veldum2_ray(i+k)=veldum_ray(i)
   zdum2_ray(i+k)=zdum_ray(i)
   opacdum2_ray(i+k)=opacdum_ray(i)
   opalbardum2_ray(i+k)=opalbardum_ray(i)
   scontdum2_ray(i+k)=scontdum_ray(i)
   slinedum2_ray(i+k)=slinedum_ray(i)
   tempdum2_ray(i+k)=tempdum_ray(i)
   vthdum2_ray(i+k)=vthdum_ray(i)

   dum_vel1=dum_vel2
! 
enddo
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum2_ray(nz_ray+1-i)
   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
   velz_ray(i) = veldum2_ray(nz_ray+1-i)

!  write(*,'(8es20.8)') z_ray(i), opalbar_ray(i), opac_ray(i), scont_ray(i), sline_ray(i), vth_ray(i), temp_ray(i), velz_ray(i)

enddo
!
!
!
call cpu_time(te_setup)
t_setup_tmp = t_setup_tmp + te_setup - ts_setup
!
!if(abs(p-0.25d0).lt.1.d-12) then
!   open(1, file='trash/vz_along.dat', form='formatted', position='append')
!      do i=1, nz_ray
!         write(1,'(4es20.8)') z_ray(i), opalbar_ray(i), sline_ray(i), velz_ray(i)
!      enddo
!   close(1)
!   read(*,*)
!endif
!
end subroutine setup_ray3d_rot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d_spc(zeta, p, lcore)
!
!------sets up a p-ray for a given observers frame frequency xobs-------
!-----------------from a spherical grid---------------------------------
!to calculate:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: zeta   angle-coordinate of cylindrical system
!       p      radial distance of point in cylindrical system
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
use prog_type
use dime_spec, only:  nr
use dime_model3d, only: nr_spc, ntheta_spc, nphi_spc, r_spc, theta_spc, phi_spc, &
                        opac3d, scont3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                        vth_ray, temp_ray,  profile_ray, velz_ray, r, transmat, rotmat_trans, nhat, &
                        del_vel2, rmin, vphot_proj, rmax
use params_spec, only: vth_fiducial, vth_min, sr, xic1, vmin, vmax, tmin, na, vmicro
use timing_spec, only: ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, t_hunt_tmp, &
                       ts_trilin, te_trilin, t_trilin_tmp
use mod_interp1d, only: interpol_yp, interpol_yp_spline
use mod_interp3d, only: get_rtp_indx, get_rtp_values1, get_rtp_values2, trilin_complete
!
implicit none
!
! ... arguments
real(dp), intent(in) :: zeta, p
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, indx
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 30000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_vmicro, dum_opac, dum_scont 
real(dp) :: r1, r2, theta1, theta2, phi1, phi2
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velx, vely, velz, velr, velphi, vel_abs
real(dp) :: radp, thetap, phip, fdum
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc, vec_vel
!real(dp), dimension(nz_max) :: zdum_ray, veldum_ray, vthdum_ray, &
!                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
!real(dp), dimension(nz_max) :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
!                               slinedum2_ray, tempdum2_ray

real(dp), dimension(:), allocatable :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(:), allocatable :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogr, llogt, llogp, llogf, lr2, ldum
!
! ... local functions
real(dp) :: calc_vmicro
real(dp) :: vthermal
logical :: boundary
!
!
!
call cpu_time(ts_setup)
!
!might need to double check allocations
!maybe  need to initialize these arrays with zeros
allocate(zdum_ray(nz_max), veldum_ray(nz_max), vthdum_ray(nz_max), opacdum_ray(nz_max), &
         opalbardum_ray(nz_max), scontdum_ray(nz_max), slinedum_ray(nz_max), tempdum_ray(nz_max))

allocate(zdum2_ray(nz_max), veldum2_ray(nz_max), vthdum2_ray(nz_max), opacdum2_ray(nz_max), &
     opalbardum2_ray(nz_max), scontdum2_ray(nz_max), slinedum2_ray(nz_max), tempdum2_ray(nz_max))
!

!
iz=0
!
do i=nr, 1, -1
   zdum=r(i)*r(i)-p*p
   if(zdum.gt.0.d0) then
      iz=iz+1
      zdum_ray(iz)=sqrt(zdum)
   else
      iz=iz+1
      zdum_ray(iz)=0.d0
      exit
   endif
enddo
!
!check if core ray 
lcore=.false.
!inner-most point in carthesian coordinates
vec_cyc(1)=p*cos(zeta)
vec_cyc(2)=p*sin(zeta)
vec_cyc(3)=zdum_ray(iz)
vec_cac=matmul(transmat, vec_cyc)
!
!note: function can in general calculate arbitrary boundary shapes, here: spherical
lcore=boundary(vec_cac)
!
!
if(.not.lcore) then
!set non core rays
   iz_dum=2*iz-1
   do i=1, iz-1
      zdum_ray(iz+i)=-zdum_ray(iz-i)
   enddo
   iz=iz_dum
   vphot_proj=0.d0
else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
   radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
   call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), radp, velx, vely, velz)
   vphot_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
endif
!
!-----------------------------------------------------------------------
!
do i=1, iz
!
!calculate z_ray in carthesian coordinates
   vec_cyc(1)=p*cos(zeta)
   vec_cyc(2)=p*sin(zeta)
   vec_cyc(3)=zdum_ray(i)
   vec_cac=matmul(transmat, vec_cyc)
!
!check if point of ray lies within region where information is stored
!note: info_region2 is used, because vec_cac on inner boundary is not exactly r=1
   call info_region2(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)
!   write(*,'(i5,2l5,10es40.20)') i, linfo_phot, linfo_max, vec_cac(1), vec_cac(2), vec_cac(3), sqrt(vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2), &
!                                 sqrt(vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2)-rmin, sqrt(vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2)-rmax
!
   if(i.eq.iz) then
      if(.not.linfo_phot) then
         write(*,*) vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2, linfo_phot, linfo_max
         stop 'error in setup_ray3d_spc: linfo_phot is false, i.e., point (x,y,z) is within rmin'
      endif
   endif
!interpolation only if point lies within region where information is stored
   if(linfo_phot.and.linfo_max) then
!
!calculate corresponding points in spherical coordinate system
      radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
      call get_angles_spc(vec_cac(1), vec_cac(2), vec_cac(3), thetap, phip)
!
!timing for setting up interpolation-parameter
      call cpu_time(ts_hunt)
!
!search for indices of a the surrounding grid-cell (or neighbouring grid-cells for extrapolation)
      call get_rtp_indx(radp, thetap, phip, r_spc, theta_spc, phi_spc, nr_spc, ntheta_spc, nphi_spc, &
                        indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, expol, rmin, rmax)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
      call get_rtp_values1(radp, thetap, phip, r_spc, theta_spc, phi_spc, nr_spc, ntheta_spc, nphi_spc, &
                          indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                          r1, r2, theta1, theta2, phi1, phi2, &
                          llogr, llogt, llogp)
!timing for setting up interpolation-parameter
      call cpu_time(te_hunt)
      t_hunt_tmp=t_hunt_tmp + te_hunt - ts_hunt
!
!timing for complete trilinear interpolation
      call cpu_time(ts_trilin)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cell vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, velx3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
      llogr=.false.
      llogt=.false.
      llogp=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, velx)
!
!get vely on cell vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, vely3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
      llogr=.false.
      llogt=.false.
      llogp=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, vely)
!
!get velz on cell vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, velz3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
      llogr=.false.
      llogt=.false.
      llogp=.false.
      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, velz)
!
!for debugging
!      fdum=0.d0
!      call vel_t3(radp*sin(thetap)*cos(phip), radp*sin(thetap)*sin(phip), radp*cos(thetap), radp, velx, vely, velz, fdum)

!transform to formal-solver carthesian coordinates
      vec_vel=matmul(rotmat_trans, (/ velx, vely, velz /))
!
!calculation of velocity projected onto ray
      veldum_ray(i)=nhat(1)*vec_vel(1) + nhat(2)*vec_vel(2) + nhat(3)*vec_vel(3)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, opalbar3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogr=.false.
!      llogt=.false.
!      llogp=.false.
!      llogf=.false.
      lr2=.true.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
!line opacity in units of 1/rstar
      opalbardum_ray(i) = valp
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, sline3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogr=.false.
!      llogt=.false.
!      llogp=.false.
!      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
      slinedum_ray(i)=valp
!
!----------------------interpolation of continuum opacity---------------
!
!get continuum-opacities on cube vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, opac3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogr=.false.
!      llogt=.false.
!      llogp=.false.
!      llogf=.false.
      lr2=.true.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
!line opacity in units of 1/rstar
      opacdum_ray(i) = valp
!
!---------------------continuum source function-------------------------
!
!get sline on cube vertices
      call get_rtp_values2(nr_spc, ntheta_spc, nphi_spc, scont3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!      llogr=.false.
!      llogt=.false.
!      llogp=.false.
!      llogf=.false.
      lr2=.false.
!actual interpolation
      call trilin_complete(radp, thetap, phip, &
                           r1, r2, theta1, theta2, phi1, phi2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, &
                           r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                           expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
      scontdum_ray(i)=valp
!
!----temporary solution: set v_micro, temperature, continuum to zero----
!
!calculation of absolute value of velocity (needed to calculate vmicro)
      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial
!
!note: minimum microturbulent velocity set to input-vmicro
      dum_vmicro = calc_vmicro(vmicro, vmax*1.d5, vel_abs)
      dum_vmicro = vmicro
!
!calculate corresponding thermal velocity (for a fixed tmin)
      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!      vthdum_ray(i) = 5.d5   !test low thermal velocities

!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
      tempdum_ray(i)=0.d0
!
!-----------------------------------------------------------------------
!
!for now: without continuum
!      scontdum_ray(i)=0.d0
!      opacdum_ray(i)=0.d0
!
      call cpu_time(te_trilin)
      t_trilin_tmp=t_trilin_tmp + te_trilin - ts_trilin
!
!------set all quantities to zero, if outside information region--------
!
   else
      opalbardum_ray(i) = 0.d0
      slinedum_ray(i) = 0.d0
      veldum_ray(i) = 0.d0
      vthdum_ray(i) = 1.d-8
      tempdum_ray(i) = 0.d0
      scontdum_ray(i) = 0.d0
      opacdum_ray(i) = 0.d0
!
   endif

!   write(*,'(5es20.8)') zdum_ray(i), slinedum_ray(i), opacdum_ray(i), opalbardum_ray(i), &
!                        veldum_ray(i)
                                                
!
enddo
!
!-----------------------------------------------------------------------
!
t_trilin_tmp=t_trilin_tmp/(iz*1.d0)
t_hunt_tmp=t_hunt_tmp/(iz*1.d0)
!
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!
!velocity at beginning point
dum_vel1=veldum_ray(1)
!
veldum2_ray(1)=veldum_ray(1)
zdum2_ray(1)=zdum_ray(1)
opacdum2_ray(1)=opacdum_ray(1)
opalbardum2_ray(1)=opalbardum_ray(1)
scontdum2_ray(1)=scontdum_ray(1)
slinedum2_ray(1)=slinedum_ray(1)
tempdum2_ray(1)=tempdum_ray(1)
vthdum2_ray(1)=vthdum_ray(1)
!
k=0
iz_dum=iz
!
!do i=2, iz
!   write(*,*) i, zdum_ray(i-1), zdum_ray(i), zdum_ray(i)-zdum_ray(i-1)
!enddo
!write(*,*) zdum_ray(635), zdum_ray(636)
!stop 'test'
!
!
do i=2, iz
!
   dum_vel2=veldum_ray(i)
!
!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!  if(dum_vel1*dum_vel2.lt.0.d0) then
!     write(*,'(a96)') 'error in setup_ray_scp: velocity jump along ray from positive velocities <-> negative velocities'
!     write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!     write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!     write(*,'(a82)') '   => 2. use an analytic description of the velocity law, instead of interpolation'
!     stop
!  endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up. note: ceiling is only f90
      nadd=ceiling(abs(dum_delv)/del_vel2)

      do j=1, nadd-1
!add nadd-1 additional points in velocity space
         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/float(nadd)
         !interpolation of all variables
         zdum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
         opacdum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
         opalbardum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
         scontdum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
         slinedum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
         tempdum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
         vthdum2_ray(i+k) = interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
!  
      enddo
!
   endif

   veldum2_ray(i+k)=veldum_ray(i)
   zdum2_ray(i+k) = zdum_ray(i)
   opacdum2_ray(i+k) = opacdum_ray(i)
   opalbardum2_ray(i+k) = opalbardum_ray(i)
   scontdum2_ray(i+k) = scontdum_ray(i)
   slinedum2_ray(i+k) = slinedum_ray(i)
   tempdum2_ray(i+k) = tempdum_ray(i)
   vthdum2_ray(i+k) = vthdum_ray(i)
!
   dum_vel1=dum_vel2
! 
enddo
!
if(iz_dum.gt.nz_max) stop 'error in setup_ray3d_spc: increase nz_max'
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum2_ray(nz_ray+1-i)
   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
   velz_ray(i) = veldum2_ray(nz_ray+1-i)
enddo
!
!
!stop 'go on in setup_ray_spc'
!
call cpu_time(te_setup)
t_setup_tmp = t_setup_tmp + te_setup - ts_setup
!
!
!
end subroutine setup_ray3d_spc
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine setup_ray_2d(zeta, p, lcore)
!!
!!------sets up a p-ray for a given observers frame frequency xobs-------
!!to calculate:
!!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!!   2. z-ray-coordinates in cartesian coordinate system
!!   3. opalbar_ray   bilinear interpolation from given 2d-grid
!!   4. sline_ray     bilinear interpolation from given 2d-grid
!!   5. velx(y,z)_ray bilinear interpolation from given 2d-grid
!!   5. velz_ray      projected velocity along the ray
!!
!!input: zeta   angle-coordinate of cylindrical system
!!       p      radial distance of point in cylindrical system
!!
!!output: lcore    logical to describe if core or non-core ray
!!        all physical quantities along ray
!!
!!note: interpolation is performed in 2d-coordinate system specified
!!      by (r, latitude), in which latitude=0 is aligned with
!!      magnetic pole axis
!!
!use prog_type
!use dime_forsol, only: nr
!use forsol, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
!                  vth_ray, temp_ray,  profile_ray, velz_ray, r, transmat, nhat, del_vel2, rmin, &
!                  vphot_proj, rmax
!use params_fs, only: obliquity, vmicro, vmax, tmin, na, vthfiducial
!use timing_forsol, only: ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, t_hunt_tmp, &
!                         ts_trilin, te_trilin, t_trilin_tmp
!use dime_model_1d, only: n1d, r1d
!use dime_model_2d, only: ntheta, latitude, vel_r2d, vel_theta2d, opalbar2d, sline2d
!!
!implicit none
!!
!! ... arguments
!real(dp), intent(in) :: zeta, p
!logical, intent(out) :: lcore
!!
!! ... local scalars
!integer(i4b) :: i,j,k
!integer(i4b) :: indx_theta1, indx_theta2, indx_rad1, indx_rad2
!integer(i4b) :: iz, iz_dum, nadd
!integer(i4b), parameter :: nz_max = 20000
!real(dp) :: zdum, dum_vel1, dum_vel2, dum_delv, dum_vmicro, vel_abs
!real(dp) :: rad1, rad2, theta1, theta2
!real(dp) :: vala, valb, valc, vald, valp
!real(dp) :: rad, theta, velx, vely, velz, xm, ym, zm, velx_m, vely_m, velz_m, vel_r, vel_theta
!real(dp) :: fdum, phi
!!
!! ... local arrays
!real(dp), dimension(3) :: vec_cac, vec_cyc
!real(dp), dimension(nz_max) :: zdum_ray, velxdum_ray, velydum_ray, velzdum_ray, veldum_ray, vthdum_ray, &
!                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
!real(dp), dimension(nz_max) :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
!                               slinedum2_ray, tempdum2_ray
!!
!! ... local logicals
!logical :: linfo_phot, linfo_max, llogr, llogt, llogf
!!
!! ... local functions
!real(dp) :: interpol_yp
!real(dp) :: calc_vmicro, vthermal
!logical :: boundary
!
!!
!call cpu_time(ts_setup)
!!
!iz=0
!!
!do i=nr, 1, -1
!   zdum=r(i)*r(i)-p*p
!   if(zdum.gt.0.d0) then
!      iz=iz+1
!      zdum_ray(iz)=sqrt(zdum)
!   else
!      iz=iz+1
!      zdum_ray(iz)=0.d0
!      exit
!   endif
!enddo
!!
!!check if core ray 
!lcore=.false.
!!inner-most point in carthesian coordinates
!vec_cyc(1)=p*cos(zeta)
!vec_cyc(2)=p*sin(zeta)
!vec_cyc(3)=zdum_ray(iz)
!vec_cac=matmul(transmat, vec_cyc)
!!
!!note: function can in general calculate arbitrary boundary shapes, here: spherical
!lcore=boundary(vec_cac)
!!
!!
!if(.not.lcore) then
!!set non core rays
!   iz_dum=2*iz-1
!   do i=1, iz-1
!      zdum_ray(iz+i)=-zdum_ray(iz-i)
!   enddo
!   iz=iz_dum
!   vphot_proj=0.d0
!else
!!calculation of photospheric velocity (only rotational velocity) projected onto ray
!   rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
!   call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz)
!   vphot_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
!endif
!!
!!-----------------------------------------------------------------------
!!
!!open(1, file='trash/model2d.dat', form='formatted', position='append')
!!
!do i=1, iz
!!
!!calculate z_ray in carthesian coordinates
!   vec_cyc(1)=p*cos(zeta)
!   vec_cyc(2)=p*sin(zeta)
!   vec_cyc(3)=zdum_ray(i)
!   vec_cac=matmul(transmat, vec_cyc)
!!
!!check if point of ray lies within region where information is stored
!   call info_region(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)
!!
!!interpolation only if point lies within region where information is stored
!   if(linfo_phot.and.linfo_max) then
!!
!!timing for setting up interpolation-parameter
!      call cpu_time(ts_hunt)
!!
!!calculate coordinates in magnetic coordinate system      
!      call transform_sigr_sigm(vec_cac(1), vec_cac(2), vec_cac(3), xm, ym, zm, obliquity)
!      rad=sqrt(xm**2 + ym**2 + zm**2)
!      theta=acos(zm/rad)
!!
!!search for theta index surrounding point of interest
!      call find_index(theta, latitude, ntheta, , iim2, iim1, ii, iip1)
!      indx_theta1=indx_theta2-1
!!search for radial index surrounding point of interest
!      call find_index(rad, r1d, n1d, iim2, iim1, ii, iip1)
!      indx_rad1=indx_rad2-1
!!
!      call get_xy_values1(rad, theta, r1d, latitude, n1d, ntheta, &
!                         indx_rad1, indx_rad2, indx_theta1, indx_theta2, &
!                         rad1, rad2, theta1, theta2, llogr, llogt)
!!
!!timing for setting up interpolation-parameter
!      call cpu_time(te_hunt)
!      t_hunt_tmp=t_hunt_tmp + te_hunt - ts_hunt
!!
!!timing for complete bilinear interpolation
!      call cpu_time(ts_trilin)
!!
!!----------------------interpolation of velocity components-------------
!!
!!r-component
!      call get_xy_values2(n1d, ntheta, r1d, latitude, vel_r2d, &
!                          indx_rad1, indx_rad2, indx_theta1, indx_theta2, &
!                          vala, valb, valc, vald, llogf)
!!
!!      write(*,'(4es20.8)') vala, valb, valc, vald
!!      llogx=.false.
!!      llogy=.false.
!!      llogz=.false.
!!      llogf=.false.
!!actual interpolation
!      call bilin(rad, theta, rad1, rad2, theta1, theta2, &
!                 vala, valb, valc, vald, llogr, llogt, llogf, vel_r)
!!
!!theta-component
!      call get_xy_values2(n1d, ntheta, r1d, latitude, vel_theta2d, &
!                          indx_rad1, indx_rad2, indx_theta1, indx_theta2, &
!                          vala, valb, valc, vald, llogf)
!!
!!      write(*,'(4es20.8)') vala, valb, valc, vald
!!      llogx=.false.
!!      llogy=.false.
!!      llogz=.false.
!!      llogf=.false.
!!actual interpolation
!      call bilin(rad, theta, rad1, rad2, theta1, theta2, &
!                 vala, valb, valc, vald, llogr, llogt, llogf, vel_theta)
!!
!!transform velocities (given 2d-plane) to carthesian coordinates with z_m=z      
!      call get_angles_spc(xm, ym, zm, fdum, phi)
!      velx_m = vel_r * sin(theta)*cos(phi) + vel_theta*cos(theta)*cos(phi)
!      vely_m = vel_r * sin(theta)*sin(phi) + vel_theta*cos(theta)*sin(phi)
!      velz_m = vel_r * cos(theta) - vel_theta*sin(theta)
!!
!!transform velocities from carthesian coordinates with z_m=z to
!!   own carthesian system with z_r=z
!      call transform_sigm_sigr(velx_m, vely_m, velz_m, velx, vely, velz, obliquity)
!      velxdum_ray(i)=velx
!      velydum_ray(i)=vely
!      velzdum_ray(i)=velz
!!
!!calculation of velocity projected onto ray
!      veldum_ray(i)=nhat(1)*velxdum_ray(i) + nhat(2)*velydum_ray(i) + nhat(3)*velzdum_ray(i)
!!
!!----------------------interpolation of line opacity--------------------
!!
!      call get_xy_values2(n1d, ntheta, r1d, latitude, opalbar2d, &
!                          indx_rad1, indx_rad2, indx_theta1, indx_theta2, &
!                          vala, valb, valc, vald, llogf)
!!
!!      write(*,'(4es20.8)') vala, valb, valc, vald
!!      llogx=.false.
!!      llogy=.false.
!!      llogz=.false.
!!      llogf=.false.
!!actual interpolation
!      call bilin(rad, theta, rad1, rad2, theta1, theta2, &
!                 vala, valb, valc, vald, llogr, llogt, llogf, valp)
!      opalbardum_ray(i) = valp
!!
!!------------------------line source function---------------------------
!!
!      call get_xy_values2(n1d, ntheta, r1d, latitude, sline2d, &
!                          indx_rad1, indx_rad2, indx_theta1, indx_theta2, &
!                          vala, valb, valc, vald, llogf)
!!
!!      write(*,'(4es20.8)') vala, valb, valc, vald
!!      llogx=.false.
!!      llogy=.false.
!!      llogz=.false.
!!      llogf=.false.
!!actual interpolation
!      call bilin(rad, theta, rad1, rad2, theta1, theta2, &
!                 vala, valb, valc, vald, llogr, llogt, llogf, valp)
!      slinedum_ray(i) = valp
!!
!!----temporary solution: set temperature, continuum to zero-------------
!!
!!calculation of absolute value of velocity (needed to calculate vmicro)
!      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vthfiducial
!!
!!note: minimum microturbulent velocity set to input-vmicro
!!      dum_vmicro = calc_vmicro(vmicro, vmax*1.d5, vel_abs)
!      dum_vmicro=vmicro
!!
!!calculate corresponding thermal velocity (for a fixed tmin)
!      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!!
!!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
!      tempdum_ray(i)=0.d0
!!
!!for now: without continuum
!      scontdum_ray(i)=0.d0
!      opacdum_ray(i)=0.d0
!!
!      call cpu_time(te_trilin)
!      t_trilin_tmp=t_trilin_tmp + te_trilin - ts_trilin
!
!
!!     write(*,'(17es13.4)') rad, theta, vec_cac(1), vec_cac(2), vec_cac(3), xm, ym, zm, velx_m, vely_m, velz_m, velxdum_ray(i), velydum_ray(i), velzdum_ray(i), &
!!                           vel_r, vel_theta, opalbardum_ray(i)
!!
!!    write(1,'(7es20.8)') rad, theta, slinedum_ray(i), opalbardum_ray(i), velx_m, vely_m, velz_m
!!
!!*******************************debug begin*****************************
!!
!!!      if(expol) then
!!!test if direct interpolation makes things better
!!!         rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
!!!         call find_index(rad, r1d, n1d, iim2, iim1, ii, iip1)
!!!         slinedum_ray(i) = interpol_yp_spline(acoeff_sline1d(indx), bcoeff_sline1d(indx), &
!!!                                              ccoeff_sline1d(indx), dcoeff_sline1d(indx), &
!!!                                              r1d(indx), rad)
!!!         opalbardum_ray(i) = interpol_yp_spline(acoeff_opalbar1d(indx), bcoeff_opalbar1d(indx), &
!!!                                              ccoeff_opalbar1d(indx), dcoeff_opalbar1d(indx), &
!!!                                              r1d(indx), rad)
!!!
!!!         call bvel3d(vmin*1.d5, vmax*1.d5, beta, vec_cac(1), vec_cac(2), vec_cac(3), dum_velx, dum_vely, dum_velz, dum_gradv)
!!!         velxdum_ray(i)=dum_velx
!!!         velydum_ray(i)=dum_vely
!!!         velzdum_ray(i)=dum_velz
!!!         dum_vel=sqrt(dum_velx*dum_velx + dum_vely*dum_vely + dum_velz*dum_velz)
!!!calculation of velocity projected onto ray
!!!         veldum_ray(i)=nhat(1)*velxdum_ray(i) + nhat(2)*velydum_ray(i) + nhat(3)*velzdum_ray(i)
!!!!calculation of line opacity
!          stop 'include stuff for opalbar_model_hamann'
!          opalbardum_ray(i) = opalbar_model_hamann(sr, vinf, mdot, kappa0, alpha, vth_fid, rad, rho)
!!!!line source function
!         slinedum_ray(i) = sobo1d(r1d(i)*sr, velr/vth, opalbar1d(i)/sr, tmin, xic1, xnue0, eps_line)
!         stop 'adapt sobolev input'
!!!      endif
!!!
!!!*******************************debug end*******************************
!!!
!!!
!!!------set all quantities to zero, if outside information region--------
!!!
!   else
!      opalbardum_ray(i) = 0.d0
!      slinedum_ray(i) = 0.d0
!      velxdum_ray(i) = 0.d0
!      velydum_ray(i) = 0.d0
!      velzdum_ray(i) = 0.d0
!      veldum_ray(i) = 0.d0
!      vthdum_ray(i) = 0.d0
!      tempdum_ray(i) = 0.d0
!      scontdum_ray(i) = 0.d0
!      opacdum_ray(i) = 0.d0
!   endif
!!
!!   write(*,'(2l4, 3es20.8, 3es20.8, 2es60.40)') linfo_phot, linfo_max, vec_cac, slinedum_ray(i), opacdum_ray(i), opalbardum_ray(i), &
!!                                                 sqrt(vec_cac(1)*vec_cac(1)+vec_cac(2)*vec_cac(2)+vec_cac(3)*vec_cac(3)), &
!!                                                 sqrt(vec_cac(1)*vec_cac(1)+vec_cac(2)*vec_cac(2)+vec_cac(3)*vec_cac(3))-1.d0
!!
!
!!stop
!!
!enddo
!
!!close(1)
!!
!!stop
!!
!!-----------------------------------------------------------------------
!!
!t_trilin_tmp=t_trilin_tmp/(iz*1.d0)
!t_hunt_tmp=t_hunt_tmp/(iz*1.d0)
!!
!!----------------check if resonance zones are resolved------------------
!!-----------if not resolved: add grid points and interpolate------------
!!
!!velocity at beginning point
!dum_vel1=veldum_ray(1)
!!
!veldum2_ray(1)=veldum_ray(1)
!zdum2_ray(1)=zdum_ray(1)
!opacdum2_ray(1)=opacdum_ray(1)
!opalbardum2_ray(1)=opalbardum_ray(1)
!scontdum2_ray(1)=scontdum_ray(1)
!slinedum2_ray(1)=slinedum_ray(1)
!tempdum2_ray(1)=tempdum_ray(1)
!vthdum2_ray(1)=vthdum_ray(1)
!!
!k=0
!iz_dum=iz
!!
!do i=2, iz
!!
!   dum_vel2=veldum_ray(i)
!!
!!calculation of velocity steps
!   dum_delv=dum_vel2-dum_vel1
!!
!!   if(dum_vel1*dum_vel2.lt.0.d0) then
!!      write(*,'(a92)') 'error in setup_ray: velocity jump along ray from positive velocities <-> negative velocities'
!!      write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!!      write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!!      write(*,'(a79)') '   => 2. use an analytic description of the velocity law, instead interpolation'
!!      stop
!!   endif
!!
!   if(abs(dum_delv).gt.del_vel2) then
!!always round up. note: ceiling is only f90
!      nadd=ceiling(abs(dum_delv)/del_vel2)
!      do j=1, nadd-1
!!add nadd-1 additional points in velocity space
!         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/nadd
!!interpolation of all variables
!         zdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
!         opacdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
!         opalbardum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
!         scontdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
!         slinedum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
!         tempdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
!         vthdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!!
!         k=k+1
!         iz_dum=iz_dum+1
!      enddo
!   endif
!!
!   veldum2_ray(i+k)=veldum_ray(i)
!   zdum2_ray(i+k)=zdum_ray(i)
!   opacdum2_ray(i+k)=opacdum_ray(i)
!   opalbardum2_ray(i+k)=opalbardum_ray(i)
!   scontdum2_ray(i+k)=scontdum_ray(i)
!   slinedum2_ray(i+k)=slinedum_ray(i)
!   tempdum2_ray(i+k)=tempdum_ray(i)
!   vthdum2_ray(i+k)=vthdum_ray(i)
!
!   dum_vel1=dum_vel2
!! 
!enddo
!
!!------------------store everything in global arrays--------------------
!!
!nz_ray=iz_dum
!!
!call allocate_fs1d
!!
!do i=1, nz_ray
!   z_ray(i) = zdum2_ray(nz_ray+1-i)
!   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
!   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
!   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
!   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
!   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
!   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
!   velz_ray(i) = veldum2_ray(nz_ray+1-i)
!enddo
!!
!!
!!
!call cpu_time(te_setup)
!t_setup_tmp = t_setup_tmp + te_setup - ts_setup
!!
!!
!end subroutine setup_ray_2d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine setup_ray1d(zeta, p, lcore)

!------sets up a p-ray for a given observers frame frequency xobs-------
!to calculate:
!   1. z-ray          from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray    monotonic spline interpolation from given 1d-grid
!   4. sline_ray      monotonic spline interpolation from given 1d-grid
!   5. velr_ray       monotonic spline interpolation from given 1d-grid
!   6. vel-components from velr_ray and coordinates as well as rotation-law
!
!input: zeta   angle-coordinate of cylindrical system
!       p      radial distance of point in cylindrical system
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
use prog_type
use fund_const, only: pi
use dime_spec, only: nr
use dime_model1d, only: n1d, r1d, sline1d, opalbar1d, velr1d, & 
                        acoeff_sline1d, bcoeff_sline1d, ccoeff_sline1d, dcoeff_sline1d, &
                        acoeff_scont1d, bcoeff_scont1d, ccoeff_scont1d, dcoeff_scont1d, &
                        acoeff_velr1d, bcoeff_velr1d, ccoeff_velr1d, dcoeff_velr1d, &
                        acoeff_opalbar1d, bcoeff_opalbar1d, dcoeff_opalbar1d, ccoeff_opalbar1d, &
                        acoeff_opac1d, bcoeff_opac1d, ccoeff_opac1d, dcoeff_opac1d, &
                        acoeff_t1d, bcoeff_t1d, dcoeff_t1d, ccoeff_t1d
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                        vth_ray, temp_ray,  profile_ray, velz_ray, r, transmat, nhat, del_vel2, rmin, rmax, &
                        vphot_proj
use params_spec, only: vth_fiducial, vth_min, sr, xic1, vmicro, vmax, tmin, na
use timing_spec, only: ts_setup, te_setup, t_setup_tmp, ts_hunt, te_hunt, t_hunt_tmp, &
                       ts_trilin, te_trilin, t_trilin_tmp
use mod_interp1d, only: find_index, interpol_yp, interpol_yp_spline
!
implicit none
!
! ... arguments
real(dp), intent(in) :: zeta, p
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv, dum_vmicro
real(dp) :: velx, vely, velz, velr, vel_abs
real(dp) :: rad
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc
real(dp), dimension(:), allocatable :: zdum_ray, velxdum_ray, velydum_ray, velzdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(:), allocatable :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogx, llogy, llogz, llogf, lr2, ldum
!
! ... local functions
logical :: boundary
real(dp) :: calc_vmicro, vthermal
!
call cpu_time(ts_setup)
!
!might need to double check allocations
!maybe  need to initialize these arrays with zeros
allocate(zdum_ray(nz_max), veldum_ray(nz_max), vthdum_ray(nz_max), opacdum_ray(nz_max), &
         opalbardum_ray(nz_max), scontdum_ray(nz_max), slinedum_ray(nz_max), tempdum_ray(nz_max))

allocate(zdum2_ray(nz_max), veldum2_ray(nz_max), vthdum2_ray(nz_max), opacdum2_ray(nz_max), &
     opalbardum2_ray(nz_max), scontdum2_ray(nz_max), slinedum2_ray(nz_max), tempdum2_ray(nz_max))
!
!
!
iz=0
!
do i=nr, 1, -1
   zdum=r(i)**2-p**2
   if(zdum.gt.0.d0) then
      iz=iz+1
      zdum_ray(iz)=sqrt(zdum)
   else
      iz=iz+1
      zdum_ray(iz)=0.d0
      exit
   endif
enddo
!
!check if core ray 
lcore=.false.
!inner-most point in carthesian coordinates
vec_cyc(1)=p*cos(zeta)
vec_cyc(2)=p*sin(zeta)
vec_cyc(3)=zdum_ray(iz)
vec_cac=matmul(transmat, vec_cyc)
!
!note: function can in general calculate arbitrary boundary shapes, here: spherical
lcore=boundary(vec_cac)
!
if(.not.lcore) then
!set non core rays
   iz_dum=2*iz-1
   do i=1, iz-1
      zdum_ray(iz+i)=-zdum_ray(iz-i)
   enddo
   iz=iz_dum
   vphot_proj=0.d0
else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
   rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
   call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz)
   vphot_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
endif
!
!write(*,'(13es20.8)') p, zeta*180.d0/pi, vec_cac, nhat, velx, vely, velz, vphot_proj, rad
!
!-----------------------------------------------------------------------
!
do i=1, iz
!
!calculate z_ray in carthesian coordinates
   vec_cyc(1)=p*cos(zeta)
   vec_cyc(2)=p*sin(zeta)
   vec_cyc(3)=zdum_ray(i)
   vec_cac=matmul(transmat, vec_cyc)

!   write(*,'(6es20.8)') vec_cac(1), vec_cac(2), vec_cac(3), zdum_ray(i), p*cos(zeta), p*sin(zeta)

   rad=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
!
!check if point of ray lies within region where information is stored
   call info_region2(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)

!   write(*,'( 5e30.20 , 2l4)') zdum_ray(i), vec_cac(1), vec_cac(2), vec_cac(3), sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2), linfo_phot, linfo_max
!
!interpolation only if point lies within region where information is stored
   if(linfo_phot.and.linfo_max) then
!
!timing for setting up interpolation-parameter
      call cpu_time(ts_hunt)

      call find_index(rad, r1d, n1d, iim2, iim1, ii, iip1)
!
!timing for setting up interpolation-parameter
      call cpu_time(te_hunt)
      t_hunt_tmp=t_hunt_tmp + te_hunt - ts_hunt
!
!timing for complete interpolation
      call cpu_time(ts_trilin)
!
!----------------------interpolation of velocity components-------------
!
      velr = interpol_yp_spline(acoeff_velr1d(ii), bcoeff_velr1d(ii), &
                                ccoeff_velr1d(ii), dcoeff_velr1d(ii), &
                                r1d(ii), rad)
!
      call vel_t3(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz, velr)
      velxdum_ray(i)=velx
      velydum_ray(i)=vely
      velzdum_ray(i)=velz
!
!calculation of velocity projected onto ray
      veldum_ray(i)=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
!
!----------------------interpolation of continuum opacity---------------

      opacdum_ray(i) = interpol_yp_spline(acoeff_opac1d(ii), bcoeff_opac1d(ii), &
                                          ccoeff_opac1d(ii), dcoeff_opac1d(ii), &
                                          r1d(ii), rad)
!
!----------------------interpolation of line opacity--------------------

      opalbardum_ray(i) = interpol_yp_spline(acoeff_opalbar1d(ii), bcoeff_opalbar1d(ii), &
                                             ccoeff_opalbar1d(ii), dcoeff_opalbar1d(ii), &
                                             r1d(ii), rad)
!
!------------------------continnuum source function---------------------
!
      scontdum_ray(i) = interpol_yp_spline(acoeff_scont1d(ii), bcoeff_scont1d(ii), &
                                           ccoeff_scont1d(ii), dcoeff_scont1d(ii), &
                                           r1d(ii), rad)
!
!------------------------line source function---------------------------
!
      slinedum_ray(i) = interpol_yp_spline(acoeff_sline1d(ii), bcoeff_sline1d(ii), &
                                           ccoeff_sline1d(ii), dcoeff_sline1d(ii), &
                                           r1d(ii), rad)
!
!------------------------temperature------------------------------------
!
      tempdum_ray(i) = interpol_yp_spline(acoeff_t1d(ii), bcoeff_t1d(ii), &
                                          ccoeff_t1d(ii), dcoeff_t1d(ii), &
                                          r1d(ii), rad)
!
!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
      tempdum_ray(i)=0.d0
!
!--------------------------microturbulent velocity----------------------
!
!calculation of absolute value of velocity (needed to calculate vmicro)
      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial
!
!note: minimum microturbulent velocity set to input-vmicro
      dum_vmicro = calc_vmicro(vmicro, vmax, vel_abs)
      dum_vmicro = vmicro
!
!----------------------------thermal velocity---------------------------
!calculate corresponding thermal velocity (for a fixed tmin)
!      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
      vthdum_ray(i) = vthermal(dum_vmicro, tempdum_ray(i), na)
!
!for now: without continuum
!      scontdum_ray(i)=0.d0
!      opacdum_ray(i)=0.d0
!
      call cpu_time(te_trilin)
      t_trilin_tmp=t_trilin_tmp + te_trilin - ts_trilin

!
!------set all quantities to zero, if outside information region--------
!------(extrapolation to zero when additional points are added)---------
!
   else
      opalbardum_ray(i) = 0.d0
      slinedum_ray(i) = 0.d0
!      velxdum_ray(i) = 0.d0
!      velydum_ray(i) = 0.d0
!      velzdum_ray(i) = 0.d0
!      veldum_ray(i) = 0.d0
!      vthdum_ray(i) = 0.d0
      tempdum_ray(i) = 0.d0
      scontdum_ray(i) = 0.d0
      opacdum_ray(i) = 0.d0
!
!velocity is extrapolated linearly in log-space (spline-extrapolation is very dangerous)
      call find_index(rad, r1d, n1d, iim2, iim1, ii, iip1)
      velr = interpol_yp(log10(r1d(ii)), log10(r1d(iim1)), &
                         log10(velr1d(ii)), log10(velr1d(iim1)), log10(rad))
      velr = 10.d0**velr

      call vel_t3(vec_cac(1), vec_cac(2), vec_cac(3), rad, velx, vely, velz, velr)
      velxdum_ray(i)=velx
      velydum_ray(i)=vely
      velzdum_ray(i)=velz
      veldum_ray(i)=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz

      vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial
      dum_vmicro = calc_vmicro(vmicro, vmax*1.d5, vel_abs)
      dum_vmicro = vmicro
      vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
      write(*,*) tmin
stop 'go on in setup_ray1d'
!
   endif

!   write(*,'(i5, 2l4, 3es20.8, 3es20.8, 2es60.40)') i, linfo_phot, linfo_max, vec_cac, slinedum_ray(i), veldum_ray(i), opalbardum_ray(i), &
!                                                    sqrt(vec_cac(1)*vec_cac(1)+vec_cac(2)*vec_cac(2)+vec_cac(3)*vec_cac(3)), &
!                                                    sqrt(vec_cac(1)*vec_cac(1)+vec_cac(2)*vec_cac(2)+vec_cac(3)*vec_cac(3))-1.d0
!
enddo
!
!-----------------------------------------------------------------------
!
t_trilin_tmp=t_trilin_tmp/(iz*1.d0)
t_hunt_tmp=t_hunt_tmp/(iz*1.d0)
!
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!
!velocity at beginning point
dum_vel1=veldum_ray(1)
!
veldum2_ray(1)=veldum_ray(1)
zdum2_ray(1)=zdum_ray(1)
opacdum2_ray(1)=opacdum_ray(1)
opalbardum2_ray(1)=opalbardum_ray(1)
scontdum2_ray(1)=scontdum_ray(1)
slinedum2_ray(1)=slinedum_ray(1)
tempdum2_ray(1)=tempdum_ray(1)
vthdum2_ray(1)=vthdum_ray(1)
!
k=0
iz_dum=iz
!
do i=2, iz
!
   dum_vel2=veldum_ray(i)
!
!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!   if(dum_vel1*dum_vel2.lt.0.d0) then
!      write(*,'(a92)') 'error in setup_ray: velocity jump along ray from positive velocities <-> negative velocities'
!      write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!      write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!      write(*,'(a79)') '   => 2. use an analytic description of the velocity law, instead interpolation'
!      stop
!   endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up. note: ceiling is only f90
      nadd=ceiling(abs(dum_delv)/del_vel2)
      do j=1, nadd-1
!add nadd-1 additional points in velocity space
         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/nadd
!interpolation of all variables
         zdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
         opacdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
         opalbardum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
         scontdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
         slinedum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
         tempdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
         vthdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
      enddo
   endif

   veldum2_ray(i+k)=veldum_ray(i)
   zdum2_ray(i+k)=zdum_ray(i)
   opacdum2_ray(i+k)=opacdum_ray(i)
   opalbardum2_ray(i+k)=opalbardum_ray(i)
   scontdum2_ray(i+k)=scontdum_ray(i)
   slinedum2_ray(i+k)=slinedum_ray(i)
   tempdum2_ray(i+k)=tempdum_ray(i)
   vthdum2_ray(i+k)=vthdum_ray(i)

   dum_vel1=dum_vel2
! 
enddo
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum2_ray(nz_ray+1-i)
   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
   velz_ray(i) = veldum2_ray(nz_ray+1-i)
enddo
!
!
!do i=1, nz_ray
!   write(*,'(10es20.8)') z_ray(i), opalbar_ray(i), opac_ray(i), scont_ray(i), sline_ray(i), vth_ray(i), temp_ray(i), velz_ray(i)
!enddo
!stop
!stop
!
call cpu_time(te_setup)
t_setup_tmp = t_setup_tmp + te_setup - ts_setup
!
end subroutine setup_ray1d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine calc_transmat
!
!---------------calculation of transformation matrix--------------------
!--------for system (ex, ey, ez) to system (eex, eey, eez)--------------
!
use prog_type
use mod_spectrum, only: alpha, gamma, nhat, transmat, rotmat
!
implicit none
!
! ... local scalars
real(dp) :: check1, check2, check3
!
! ... local arrays
real(dp), dimension(3) :: ex, ey, ez, eex, eey, eez
!
write(*,*) '----------calculate transformation matrix--------------'
write(*,*) '(alpha, gamma): ', alpha, gamma
write(*,*)
!
!calculate nhat (aligned with z-axis of cylindrical coordinate system)
nhat(1)=sin(alpha)*cos(gamma)
nhat(2)=sin(alpha)*sin(gamma)
nhat(3)=cos(alpha)
if(abs(nhat(1)).lt.1.d-14) nhat(1)=0.d0
if(abs(nhat(2)).lt.1.d-14) nhat(2)=0.d0
if(abs(nhat(3)).lt.1.d-14) nhat(3)=0.d0
!
!calculate carthesian coordinate system that is aligned with 
!cylindrical coordinates
!
!ez predefined: direction to observer
ez=nhat
!
!calculate orthogonal ex from dot_product(ex,ez)=0
if(ez(1).eq.0.d0) then
   if(ez(2).eq.0.d0) then
      ex = (/ 1.d0, 0.d0, 0.d0 /)
   else if(ez(3).eq.0.d0) then
      ex = (/ 1.d0, 0.d0, 0.d0 /)
   else
      ex(1) = 1.d0
      ex(2) = 1.d0
      ex(3) = -ex(2)*ez(2)/ez(3)
   endif
else if (ez(2).eq.0.d0) then
   if(ez(1).eq.0.d0) then
      ex = (/ 1.d0, 0.d0, 0.d0 /)
   else if(ez(3).eq.0.d0) then
      ex = (/ 0.d0, 1.d0, 0.d0 /)
   else
      ex(3) = 1.d0
      ex(2) = 1.d0
      ex(1) = -ex(3)*ez(3)/ez(1)
   endif
else if (ez(3).eq.0.d0) then
   if(ez(1).eq.0.d0) then
      ex = (/ 1.d0, 0.d0, 0.d0 /)
   else if(ez(2).eq.0.d0) then
      ex = (/ 0.d0, 1.d0, 0.d0 /)
   else
      ex(3) = 1.d0
      ex(2) = 1.d0
      ex(1) = -ex(2)*ez(2)/ez(1)
   endif
else
   ex(1) = 1.d0
   ex(2) = 1.d0
   ex(3) = (-ex(1)*ez(1)-ex(2)*ez(2))/ez(3)
endif
!
!calculate orthogonal ey from cross-product
call cross_product(ez, ex, ey)
!
!normalize unit vectors
ex = ex/sqrt(dot_product(ex,ex))
ey = ey/sqrt(dot_product(ey,ey))
ez = ez/sqrt(dot_product(ez,ez))
!
!check for orthogonality
check1=dot_product(ex,ey)
check2=dot_product(ex,ez)
check3=dot_product(ey,ez)
!
if(abs(check1).gt.1.d-14) stop 'error in calc_transmat: ex, ey not orthogonal'
if(abs(check2).gt.1.d-14) stop 'error in calc_transmat: ex, ez not orthogonal'
if(abs(check3).gt.1.d-14) stop 'error in calc_transmat: ey, ez not orthogonal'
!
!----transformation matrix from system (ex,ey,ez) to (eex, eey, eez)----
!
eex = (/ 1.d0, 0.d0, 0.d0 /)
eey = (/ 0.d0, 1.d0, 0.d0 /)
eez = (/ 0.d0, 0.d0, 1.d0 /)
!
transmat = reshape((/ dot_product(eex,ex), dot_product(eex,ey), dot_product(eex,ez), &
                      dot_product(eey,ex), dot_product(eey,ey), dot_product(eey,ez), &
                      dot_product(eez,ex), dot_product(eez,ey), dot_product(eez,ez) /), shape(transmat))
transmat = transpose(transmat)
!
transmat = matmul(rotmat, transmat)
!
write(*,'(a5, 3es20.8)') 'e_x', ex
write(*,'(a5, 3es20.8)') 'e_y', ey
write(*,'(a5, 3es20.8)') 'e_z', ez
write(*,'(3es20.8)') transmat
write(*,*)
!
!
end subroutine calc_transmat
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine calc_rotmat(obliquity_fs, obliquity_3d)
!
!---------------calculation of transformation matrix--------------------
!---to match two carthesian coordinate systems, that shall be aligned---
!--------------------along a given vector-------------------------------
!
!   e.g.: for obliquity beta1 defined for 3d-input system
!         and different obliquity beta2 defined here, in formal solver:
!         need to rotate 3d-input-system to obtain formal solution
!         for beta2
!
use prog_type
use fund_const, only: pi
use mod_spectrum, only: rotmat, rotmat_trans
!
implicit none
!
! ... arguments
real(dp), intent(in) :: obliquity_fs, obliquity_3d
!
! ... local scalars
real(dp) :: del_beta, cbeta, sbeta
!
! ... local arrays
real(dp), dimension(3) :: vec
!
!
write(*,*) '-----------------------calculate rotation matrix-------------------------------'
write(*,*) '(obliquity (3d-system): ', obliquity_3d*180.d0/pi
write(*,*) '(obliquity (fs-system): ', obliquity_fs*180.d0/pi
write(*,*)
!
del_beta=obliquity_fs-obliquity_3d
!
!calculate cos and sine of del_beta
cbeta=cos(del_beta)
sbeta=sin(del_beta)
!
if(abs(cbeta).lt.1.d-14) cbeta=0.d0
if(abs(sbeta).lt.1.d-14) sbeta=0.d0
!
!to rotate from 3d coordinates to formal solver system
rotmat_trans = reshape((/ cbeta, 0.d0, -sbeta, &
                           0.d0, 1.d0,   0.d0, &
                          sbeta, 0.d0, cbeta /), shape(rotmat_trans))
!
!to rotate from formal solver system to 3d coordinates
rotmat = transpose(rotmat_trans)
!
write(*,*) rotmat
write(*,*)
!
!vec=(/ sin(obliquity_fs), 0.d0, cos(obliquity_fs) /)
!write(*,*) vec
!write(*,*) matmul(rotmat, vec)
!
!
end subroutine calc_rotmat
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_gdark
!
!----------------------calculate gravity darkening law------------------
!-------------following petrenz/puls 1996, cranmer/owocki 1995----------
!---------------------(neglecting eddington factor!!!)------------------
!
use prog_type
use fund_const, only: cgs_grav, xmsu, xlsu, pi, cgs_sb
use mod_gdark
use params_spec, only: xlogg, sr, rstar, mstar, lstar, vrot, trad, xic1, xnue0, teff
use options_spec, only: opt_incl_gdark, opt_incl_sdist
use mod_spectrum, only: xic_nue, xicc_nue
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: ntheta_fine=101
integer(i4b) :: i, iim2, iim1, ii, iip1
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, l_star_cgs, w0, omega, sigma, sigma_fit, c1, c2, sum
!
! ... local arrays
real(dp), dimension(ntheta_fine) :: theta_fine, rsurface_fine, gr_fine, gtheta_fine, gperp_fine, integrand_fine, teff_fine
!
! ... local functions
real(dp) :: bnue, calc_req
!
!create theta-grid
do i=1, ntheta_fine
   theta_fine(i) = 1.d-5 + (i-1)*(pi/2.d0-1.d-5)/(ntheta_fine-1)
enddo
!
!
r_pole_cgs = sr
m_star_cgs = sr**2 * 10.d0**xlogg/cgs_grav
mstar = m_star_cgs/xmsu
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
!
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
   xic1_gdark(i) = bnue(xnue0, teff_gdark(i))
enddo
!
if(.not.opt_incl_gdark.and..not.opt_incl_sdist) then
!if neither surface distortion nor gravity darkening is included, overwrite arrays
   teff_gdark = trad
   xic1_gdark = xic1
elseif(.not.opt_incl_gdark) then
!if surface distortion is included and gravity darkening is excluded, calculate average effective temperatures
!
!------------------------------old version------------------------------
!
!   do i=1, ntheta_fine
!      integrand_fine(i) = 4.d0*pi*cgs_sb*(rsurface_fine(i)*r_pole_cgs)**2*sin(theta_fine(i))
!!      write(*,'(i5,10es20.8)') i, rsurface_fine(i), teff_fine(i), theta_fine(i), integrand_fine(i)
!   enddo
!!calculate average effective temperature
!   teff=0.d0
!   do i=2, ntheta_fine
!      teff = teff + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
!   enddo
!   write(*,*) 'teff (old)', (l_star_cgs/teff)**0.25d0
!
!------------------------------old version------------------------------
!
!account for distorted surface when calculating the flux integral
   do i=1, ntheta_fine
      integrand_fine(i) = 4.d0*pi*cgs_sb*(rsurface_fine(i)*r_pole_cgs)**2*sin(theta_fine(i))*gperp_fine(i)/(-gr_fine(i))
!      write(*,'(i5,10es20.8)') i, rsurface_fine(i), teff_fine(i), theta_fine(i), integrand_fine(i)
   enddo
!calculate average effective temperature
   teff=0.d0
   do i=2, ntheta_fine
      teff = teff + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
   enddo
   write(*,*) 'teff (new)', (l_star_cgs/teff)**0.25d0
!
!
!
   teff=(l_star_cgs/teff)**0.25d0
   teff_fine=teff
   trad=teff
   xic_nue=xic_nue/xic1
   xicc_nue=xicc_nue/xic1
   xic1=bnue(xnue0,teff)
   xic_nue=xic_nue*xic1
   xicc_nue=xicc_nue*xic1
   xic1_gdark=xic1
   teff_gdark=teff
endif
!
!-------------------calculate stellar surface distortion----------------
!
!calculate radius at equator
if(opt_incl_sdist) then 
   smajorax_a = calc_req(rstar, mstar, vrot)
   smajorax_b = smajorax_a
   smajorax_c = 1.d0
else
   smajorax_a = 1.d0
   smajorax_b = 1.d0
   smajorax_c = 1.d0
endif
!
!set scaling factor for surface intensity to one
!(required if non-rotating routines are used, and is adapted anyways when rotating routines are used)
xic1_factor=1.d0
!
!
!calculate luminosity (as test)
!do i=1, ntheta_fine
!   integrand_fine(i) = 4.d0*pi*cgs_sb*(rsurface_fine(i)*r_pole_cgs)**2.d0*sin(theta_fine(i))*teff_fine(i)**4
!enddo
!sum=0.d0
!do i=2, ntheta_fine
!   sum=sum + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
!enddo
!write(*,*) sum/l_star_cgs, sum, l_star_cgs
!stop
!
!
end subroutine calc_gdark
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_xic_factor(x, y, z)
!
!calculates a scaling factor for core intensity:
!   I_c(theta) = xic_factor * I_c(pole)
!   => xic_factor = Bnue(Teff(theta))/Bnue(Teff(pole))
!
!on input:   x, y, z describes the coordinates for a rotating star
!                    on the stellar surface
!
use prog_type
use mod_gdark, only: xic1_factor, ntheta_gdark, theta_gdark, xic1_gdark
use params_spec, only: xic1, xnue0
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: x, y, z
!
! ... local scalars
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: rad, theta, xic1_local
!
! ... local functions
real(dp) :: bnue
!
!calculate co-latitude
rad = sqrt(x**2+y**2+z**2)
theta = acos(abs(z/rad))
!
!interpolate to obtain local xic1
call find_index(theta, theta_gdark, ntheta_gdark, iim2, iim1, ii, iip1)
xic1_local = interpol_yp(theta_gdark(iim1), theta_gdark(ii), xic1_gdark(iim1), xic1_gdark(ii), theta)
xic1_factor = xic1_local/xic1
!
!
!
end subroutine calc_xic_factor
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_iin(indx_xobs, lcore, iin, iin_c)
!
use prog_type
use fund_const, only: cgs_clight
use options_spec, only: interp_photprof
use dime_spec, only: nxobs_fs
use params_spec, only: xnue0, xic1
use mod_spectrum, only: xic_nue, xicc_nue, xobs, xnue, vphot_proj, acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic
use mod_gdark, only: xic1_factor
use omp_lib
use mod_interp1d, only: find_index, interpol_yp, interpol_yp_spline
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_xobs
logical, intent(in) :: lcore
real(dp), intent(out) :: iin, iin_c
!
! ... local scalars
integer(i4b) ::  iim2, iim1, ii, iip1
real(dp) :: nue_cmf, xcmf, delta
!
! ... local functions


if(lcore) then
!note: need to interpolate from photospheric profile, since this is only given in cmf, 
!      whereas i am calculating in observers frame (shift of profile due to rotational velocities)
!calculate comoving frame frequency
   xcmf = xobs(indx_xobs) - vphot_proj
!
!constant continuum
   iin_c=xic1*xic1_factor
!
!find index within xobs-grid, for which xcmf matches and interpolate
   call find_index(xcmf, xobs, nxobs_fs, iim2, iim1, ii, iip1)
!   if(omp_get_thread_num().eq.1.and.indx_xobs.eq.1) write(*,'(i5,l5,5es20.8)') indx_xobs, lcore, xic1, xic1_factor, iin, iin_c
!
   select case(interp_photprof)
      case(0) 
         iin=interpol_yp(xobs(iim1), xobs(ii), xic_nue(iim1), xic_nue(ii), xcmf) * xic1_factor
         iin_c=interpol_yp(xobs(iim1), xobs(ii), xicc_nue(iim1), xicc_nue(ii), xcmf)  * xic1_factor !for frequency dependent continuum
      case(1)
         iin=interpol_yp_spline(acoeff_xic(ii), bcoeff_xic(ii), &
                                ccoeff_xic(ii), dcoeff_xic(ii), &
                                xobs(ii), xcmf) * xic1_factor
      case default
         stop 'error get_iin: interp_photprof not defined'
   end select
!
!
else
!set non-core ray
   iin=0.d0
   iin_c=0.d0
endif

!
end subroutine get_iin
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine print_physical
!!
!use prog_type
!use fund_const
!use params_fs, only: vthfiducial, sr, xnue0, xic1, xlogg, xmloss, yhe, &
!                     vmax, vmin, vmicro, teff, trad, tmin, eps_line, vrot, &
!                     v_esc, rhoc_star, t_inf, obliquity, ralfven, chi_inf
!!
!implicit none
!!
!write(*,*) '-------physical parameter used for calculations--------'
!!
!write(*,'(a20, e20.8)') 'r_star [cm]', sr
!write(*,'(a20, e20.8)') 'v_max [km/s]', vmax
!write(*,'(a20, e20.8)') 'v_min [km/s]', vmin
!write(*,'(a20, e20.8)') 'vth_fiducial [cm/s]', vthfiducial
!write(*,'(a20, e20.8)') 'v_micro [km/s]', vmicro
!write(*,'(a20, e20.8)') 'v_rot [km/s]', vrot
!write(*,'(a20, e20.8)') 't_eff [k]', teff
!write(*,'(a20, e20.8)') 't_rad [k]', trad
!write(*,'(a20, e20.8)') 't_min [k]', tmin
!write(*,'(a20, e20.8)') 'nue_0 [1/s]', xnue0
!write(*,'(a20, e20.8)') 'lambda [a]', cgs_clight/xnue0*1.d8
!write(*,'(a20, e20.8)') 'epsilon_l', eps_line
!write(*,'(a20, e20.8)') 'xic1', xic1
!write(*,'(a20, e20.8)') 'xlogg', xlogg
!write(*,'(a20, e20.8)') 'mdot', xmloss
!write(*,'(a20, e20.8)') 'v_esc [cm/s]', v_esc
!write(*,'(a20, e20.8)') 'rhoc_star [g/cm^3]', rhoc_star
!write(*,'(a20, e20.8)') 't_inf [k]', t_inf
!write(*,'(a20, e20.8)') 'obliquity', obliquity
!write(*,'(a20, e20.8)') 'ralfven [r_star]', ralfven
!write(*,'(a20, e20.8)') 'chi_inf', chi_inf
!write(*,*)
!
!end subroutine print_physical
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine print_model2d
!!
!use prog_type
!use fund_const, only: pi
!use dime_model_1d, only: n1d, r1d
!use dime_model_2d, only: ntheta, latitude, opalbar2d, vel_r2d, vel_theta2d, sline2d
!!
!implicit none
!!
!! ... arguments
!integer(i4b) :: i, j
!!
!write(*,*) '----------atmospheric structure from 2d-model----------'
!!
!j=1
!write(*,*) '      for theta=', latitude(j)*180./pi
!write(*,'(a4, 7(a20))') '#', 'r [r_star]', 'opath [1/r_star]', &
!                        'opalbar [1/r_star]', 'vel_r [vth_fiducial]', &
!                        'vel_theta [vth_fid]', 's_cont', 's_line [xic1]'
!do i=1, n1d
!   write(*,'(i4, e20.8, a20, 3e20.8, a20, e20.8)')  i, &
!                         r1d(i), '-----', opalbar2d(i,j), &
!                         vel_r2d(i,j), vel_theta2d(i,j), &
!                         '-----', sline2d(i,j)
!enddo
!write(*,*)
!!
!j=ntheta/2
!write(*,*) '      for theta=', latitude(j)*180./pi
!write(*,'(a4, 7(a20))') '#', 'r [r_star]', 'opath [1/r_star]', &
!                        'opalbar [1/r_star]', 'vel_r [vth_fiducial]', &
!                        'vel_theta [vth_fid]', 's_cont', 's_line [xic1]'
!do i=1, n1d
!   write(*,'(i4, e20.8, a20, 3e20.8, a20, e20.8)')  i, &
!                         r1d(i), '-----', opalbar2d(i,j), &
!                         vel_r2d(i,j), vel_theta2d(i,j), &
!                         '-----', sline2d(i,j)
!enddo
!write(*,*)
!!
!!
!!
!end subroutine print_model2d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine print_ray(xobs, zeta, p, iin, iin_c, fname)
!
use prog_type
use mod_spectrum, only: nhat, nz_ray, z_ray, opac_ray, opalbar_ray, velz_ray, &
                        scont_ray, sline_ray, temp_ray, profile_ray
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, p, zeta, iin, iin_c
character(len=*), intent(in) :: fname
!
! ... local scalars
integer(i4b) :: i
!
write(*,*) '-----------atmospheric structure along ray-------------'
!
open(1,file=trim(fname))
   write(1,'(a15, 3e20.8)') 'for nhat:    ', nhat
   write(1,'(a15, 2e20.8)') 'for p, zeta: ', p ,zeta
   write(1,'(a15, e20.8)')  'for xobs:    ', xobs
   write(1,'(a15, e20.8)')  'boundary iin ', iin
   write(1,'(a15, e20.8)')  'boundary iin_c ', iin_c
   write(1,*)
   write(1,'(a4, 8(a20))') '#', 'z [r_star]', 'opac [1/cm]', 'opalbar [1/cm]', 'vel_z [cm/s]', 's_cont', 's_line', 'temp', 'profile'
   do i=1, nz_ray
      write(1,'(i4, 8(e20.8))')  i, z_ray(i), opac_ray(i), opalbar_ray(i), velz_ray(i), scont_ray(i), sline_ray(i), temp_ray(i), profile_ray(i)
   enddo
   write(1,*)
close(1)
!
write(*,'(a15, 3e20.8)') 'for nhat:    ', nhat
write(*,'(a15, 2e20.8)') 'for p, zeta: ', p ,zeta
write(*,'(a15, e20.8)')  'for xobs:    ', xobs
write(*,'(a15, e20.8)')  'boundary iin ', iin
write(*,'(a15, e20.8)')  'boundary iin_c ', iin_c
write(*,*)
write(*,'(a4, 8(a20))') '#', 'z [r_star]', 'opac [1/cm]', 'opalbar [1/cm]', 'vel_z [cm/s]', 's_cont', 's_line', 'temp', 'profile'
do i=1, nz_ray
   write(*,'(i4, 8(e20.8))')  i, z_ray(i), opac_ray(i), opalbar_ray(i), velz_ray(i), scont_ray(i), sline_ray(i), temp_ray(i), profile_ray(i)
enddo
write(*,*)

end subroutine print_ray
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine print_flux
!
use prog_type
use fund_const, only: pi
use dime_spec, only: nxobs_fs
use mod_spectrum, only: flux_tot, flux_cont, xobs, normt
use params_spec, only: vmax, vth_fiducial
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
!
!xobs in units of vmax, insted of vthfiducial
!
write(*,*) '-------------------------emergent flux profile---------------------------------'
!
write(*,'(a4, 5(a20))') '#', 'xobs', 'flux_tot', 'flux_cont', 'f_tot/f_cont', 'normt'
do i=1, nxobs_fs
   write(*,'(i4, 5(e20.8))')  i, xobs(i)*vth_fiducial/vmax, flux_tot(i), flux_cont(i), &
                              flux_tot(i)/flux_cont(i), normt(i)
enddo
write(*,*)
!
end subroutine print_flux
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine print_ferr
!
use prog_type
use dime_spec, only: np, nzeta
use mod_spectrum, only: relerr_cont, relerr_tot, relerr_contp, relerr_totp, hmax
!
implicit none
!
! ... local scalars
!

write(*,*) '---------------------error estimation of flux integral-------------------------'
write(*,'(a50, 2es20.8)') 'max rel error p-integ f_tot, f_cont', relerr_totp, relerr_contp
write(*,'(a50, 2es20.8)') 'max rel error zeta-integ f_tot, f_cont', relerr_tot, relerr_cont
write(*,'(a50, i10, i10, es20.8)') '# p/zeta-grid-points, maximum step size', np, nzeta, hmax
write(*,*)
!
!
end subroutine print_ferr
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine calc_err(nd, sol1, sol2, abserr1, abserr2, relerr1, relerr2)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp), dimension(nd), intent(in) :: sol1, sol2, abserr1, abserr2
real(dp) :: relerr1, relerr2
!
! ... local scalars
integer(i4b) :: i
real(dp) :: dum_err1, dum_err2
!
relerr1=0.d0
relerr2=0.d0
!
do i=1, nd
!
   if(sol1(i).eq.0.d0) then
      dum_err1=0.d0
   else
      dum_err1=abs(abserr1(i)/sol1(i))
   endif
!
   if(dum_err1.gt.relerr1) then
      relerr1=dum_err1
   endif
!
   if(sol2(i).eq.0.d0) then
      dum_err2=0.d0
   else
      dum_err2=abs(abserr2(i)/sol2(i))
   endif
!
   if(dum_err2.gt.relerr2) then
      relerr2=dum_err2
   endif
!
enddo
!
!
!
end subroutine calc_err
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine allocate_fluxes
!
use prog_type
use dime_spec, only: nxobs_fs, nzeta
use mod_spectrum, only: flux_tot_tmp, flux_cont_tmp, ftot_p, fcont_p, ftot_errp, fcont_errp, &
                        normt_p, normt_tmp, ftot_err_tmp, fcont_err_tmp, flux_emi_tmp, flux_abs_tmp, &
                        femi_p, fabs_p
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
if(allocated(flux_tot_tmp)) deallocate(flux_tot_tmp)
if(allocated(flux_cont_tmp)) deallocate(flux_cont_tmp)
if(allocated(normt_tmp)) deallocate(normt_tmp)
if(allocated(ftot_err_tmp)) deallocate(ftot_err_tmp)
if(allocated(fcont_err_tmp)) deallocate(ftot_err_tmp)
if(allocated(ftot_p)) deallocate(ftot_p)
if(allocated(fcont_p)) deallocate(fcont_p)
if(allocated(ftot_errp)) deallocate(ftot_errp)
if(allocated(fcont_errp)) deallocate(fcont_errp)
if(allocated(normt_p)) deallocate(normt_p)
if(allocated(flux_emi_tmp)) deallocate(flux_emi_tmp)
if(allocated(flux_abs_tmp)) deallocate(flux_abs_tmp)
if(allocated(femi_p)) deallocate(femi_p)
if(allocated(fabs_p)) deallocate(fabs_p)
!
!
allocate(flux_tot_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_tot_tmp'
!
allocate(flux_cont_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_cont_p'
!
allocate(normt_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: normt_tmp'
!
allocate(ftot_err_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_err_tmp'
!
allocate(fcont_err_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_err_tmp'
!
allocate(ftot_p(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_p'
!
allocate(fcont_p(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_p'
!
allocate(ftot_errp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_errp'
!
allocate(fcont_errp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_errp'
!
allocate(normt_p(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: normt_p'
!
allocate(flux_emi_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_emi_tmp'
!
allocate(flux_abs_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_abs_tmp'
!
allocate(femi_p(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: femi_p'
!
allocate(fabs_p(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fabs_p'
!
!
!
end subroutine allocate_fluxes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine deallocate_fluxes
!
use prog_type
use dime_spec, only: nxobs_fs, nzeta
use mod_spectrum, only: flux_tot_tmp, flux_cont_tmp, ftot_p, fcont_p, ftot_errp, fcont_errp, &
                        normt_p, normt_tmp, ftot_err_tmp, fcont_err_tmp, flux_emi_tmp, flux_abs_tmp, &
                        femi_p, fabs_p
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
deallocate(flux_tot_tmp)
deallocate(flux_cont_tmp)
deallocate(normt_tmp)
deallocate(ftot_err_tmp)
deallocate(fcont_err_tmp)
deallocate(ftot_p)
deallocate(fcont_p)
deallocate(ftot_errp)
deallocate(fcont_errp)
deallocate(normt_p)
deallocate(flux_emi_tmp)
deallocate(flux_abs_tmp)
deallocate(femi_p)
deallocate(fabs_p)
!
!
!
end subroutine deallocate_fluxes
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine allocate_fs3d
!!
!use prog_type
!use dime_model_3d, only: ndxmax_fs, ndymax_fs, ndzmax_fs, x, y, z, &
!                         sline3d, scont3d, t3d, opath3d, opalbar3d, &
!                         vel3d, velx3d, vely3d, velz3d, mask_totreg3d, &
!                         bmask3d, r3d
!!
!implicit none
!!
!! ... arguments
!!
!! ... local scalars
!integer(i4b) :: err
!!
!allocate(x(ndxmax_fs))
!allocate(y(ndymax_fs))
!allocate(z(ndzmax_fs))
!!
!allocate(r3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: r3d'
!allocate(sline3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: sline3d'
!allocate(scont3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: scont3d'
!allocate(t3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: t3d'
!allocate(opath3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: opath3d'
!allocate(opalbar3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: opalbar3d'
!allocate(vel3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: vel3d'
!allocate(velx3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: velx3d'
!allocate(vely3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: vely3d'
!allocate(velz3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: velz3d'
!allocate(mask_totreg3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: mask_totreg3d'
!allocate(bmask3d(ndxmax_fs, ndymax_fs, ndzmax_fs), stat=err)
!   if(err.ne.0) stop 'allocation error allocate_fs3d: bmask3d'
!!
!!
!!
!end subroutine allocate_fs3d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine timing_init
!
use prog_type
use timing_spec, only: t_tot, t_setup, t_setup_tmp, t_hunt, t_hunt_tmp, &
                       t_trilin, t_trilin_tmp, ttot_obsdir
!
implicit none
!
t_tot=0.d0
ttot_obsdir=0.d0
t_setup=0.d0
t_setup_tmp=0.d0
t_hunt=0.d0
t_hunt_tmp=0.d0
t_trilin=0.d0
t_trilin_tmp=0.d0
!
end subroutine timing_init
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine print_timing
!
use prog_type
use dime_spec, only: nxobs_fs, nzeta, np
use timing_spec, only: ttot_obsdir, t_setup, t_hunt, t_trilin, ttot_obsdir_sys
!
implicit none
!
!
write(*,*) '----------------------------------timing---------------------------------------'
write(*,*) 'nxobs_fs: ', nxobs_fs
write(*,*) 'nzeta:    ', nzeta
write(*,*) 'np:       ', np
write(*,*)
write(*,'(a40, 2e20.8)') 'total computation time (cpu/sys)', ttot_obsdir, ttot_obsdir_sys
write(*,'(a40, e20.8)') 'total time for hunting', t_hunt
write(*,'(a40, e20.8)') 'total time for interpolation', t_trilin
write(*,'(a40, e20.8)') 'total time for setting up all rays', t_setup
write(*,'(a40, e20.8)') 'average time for one hunt', t_hunt/(nxobs_fs*nzeta*np)
write(*,'(a40, e20.8)') 'average time for interpolation', t_trilin/(nxobs_fs*nzeta*np)
write(*,'(a40, e20.8)') 'average time for setting up each ray', t_setup/(nxobs_fs*nzeta*np)
write(*,*)

end subroutine print_timing
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine check_weights
!
!-------------------testing weights and error weights-------------------
!----------------------------for p-grid---------------------------------
!
use prog_type
use fund_const, only: pi
use dime_spec, only: np
use mod_spectrum, only: p, pw, pw1, pw_err
!
implicit none
!
! ... local scalars
real(dp) :: sum1, sum2, sumerr, sumex
!
! ... local arrays
real(dp), dimension(np) :: yvalue
!
!
write(*,*) '--------------------checking weights and error weights-------------------------'
write(*,*)
!
write(*,'(6a20)') 'function', 'half step', 'full step', 'exact', 'error(num)', 'error(exact)'
!
!linear function
yvalue=p
sumex=0.5d0*(p(np)**2.d0 - p(1)**2.d0)
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'linear', sum1, sum2, sumex, sumerr, sumex-sum1
!
!quadratic function
yvalue=p*p
sumex=1.d0/3.d0 * (p(np)**3.d0 - p(1)**3.d0)
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'quadratic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cubic function
yvalue=p*p*p
sumex=1.d0/4.d0 * (p(np)**4.d0 - p(1)**4.d0)
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cubic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!to the four
yvalue=p*p*p*p / p(np)
sumex=1.d0/5.0/p(np) * (p(np)**5.d0 - p(1)**5.d0)
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'to the four', sum1, sum2, sumex, sumerr, sumex-sum1
!
!exponential
yvalue=exp(p/p(np)*5.d0)
sumex=(exp(5.d0)-exp(p(1)/p(np)*5.d0))*p(np)/5.d0
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'exponential', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cosine function
yvalue=cos(p/p(np) * 2.d0*pi)
sumex=(sin(2.d0*pi)-sin(p(1)/p(np) * 2.d0*pi))*p(np)/(2.d0*pi)
sum1=sum(pw*yvalue)
sum2=sum(pw1*yvalue)
sumerr=sum(pw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cos', sum1, sum2, sumex, sumerr, sumex-sum1
write(*,*)
!
!
!
end subroutine check_weights
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine check_t1
!!
!!------------------testing analytic soblev calculations, etc------------
!!
!use prog_type
!use options_forsol, only: outputdir_fs
!use dime_forsol, only: nzeta, np
!use dime_model_3d, only: ndxmax_fs, ndymax_fs, ndzmax_fs, x, y, z, &
!                         r3d, opalbar3d, sline3d, velx3d, vely3d, velz3d
!use forsol, only: nz_ray, z_ray, opalbar_ray, sline_ray, vth_ray, &
!                  temp_ray,  profile_ray, velz_ray, r, transmat, nhat, p, zeta, rmin, rmax
!use params_fs, only: vthfiducial, sr
!!
!implicit none
!!
!! ... arguments
!!
!! ... local scalars
!integer(i4b) :: i, indx_zeta, indx_p
!integer(i4b) :: err
!integer(i4b) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2
!integer(i4b) :: iz, indx_max
!integer(i4b), parameter :: nz_max=20000
!real(dp), parameter :: delz=1.d-1, zmax=20.d0
!real(dp) :: zdum, val
!real(dp) :: x1, x2, y1, y2, z1, z2
!real(dp) :: rad, dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_sline, dum_opalbar
!!
!! ... local arrays
!real(dp), dimension(3) :: vec_cac, vec_cyc
!real(dp), dimension(:), allocatable :: vvelx_ray, vvely_ray, vvelz_ray, zdum_ray
!real(dp), dimension(:,:,:),  allocatable :: rad3d_ray, velx3d_ray, vely3d_ray, velz3d_ray, sline3d_ray, opalbar3d_ray
!!
!! ... local logicals
!logical :: expol, core, linfo_phot, linfo_max, llogx, llogy, llogz, llogf, lr2
!!
!! ... local functions
!real(dp) :: opalbar_model_hamann
!!
!!---------calculate maximum z-indx to allocate arrays properly----------
!!
!indx_max=0
!!
!do indx_zeta=1, nzeta
!   do indx_p=1, np
!!
!      core=.false.
!      iz=1
!      zdum=zmax
!      do i=1, nz_max
!         zdum=zdum-delz
!         vec_cyc(1)=p(indx_p)*cos(zeta(indx_zeta))
!         vec_cyc(2)=p(indx_p)*sin(zeta(indx_zeta))
!         vec_cyc(3)=zdum
!         vec_cac=matmul(transmat, vec_cyc)
!         if(sqrt(vec_cac(1)*vec_cac(1) + vec_cac(2)*vec_cac(2) + vec_cac(3)*vec_cac(3)).lt.1.d0) then
!            core=.true.
!            exit
!         endif
!         if(zdum.lt.-zmax) exit
!         iz=iz+1
!      enddo
!      if(indx_max.lt.iz) indx_max=iz
!   enddo
!enddo
!!
!!---------------------allocate arrays for output------------------------
!!
!allocate(rad3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: rad3d_ray'
!rad3d_ray=0.d0
!!
!allocate(sline3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: sline3d_ray'
!sline3d_ray=0.d0
!!
!allocate(velx3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: velx3d_ray'
!velx3d_ray=0.d0
!!
!allocate(vely3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: vely3d_ray'
!vely3d_ray=0.d0
!!
!allocate(velz3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: velz3d_ray'
!velz3d_ray=0.d0
!!
!allocate(opalbar3d_ray(np, nzeta, indx_max), stat=err)
!   if(err.ne.0) stop 'allocation error in check_interpolation: opalbar3d_ray'
!opalbar3d_ray=0.d0
!!
!!-----------------------------------------------------------------------
!!
!do indx_zeta=1, nzeta
!   do indx_p=1, np
!!
!!--------------------calculate z-grid along ray-------------------------
!!
!      core=.false.
!      iz=1
!      zdum=zmax
!      do i=1, nz_max
!         zdum=zdum-delz
!         vec_cyc(1)=p(indx_p)*cos(zeta(indx_zeta))
!         vec_cyc(2)=p(indx_p)*sin(zeta(indx_zeta))
!         vec_cyc(3)=zdum
!         vec_cac=matmul(transmat, vec_cyc)
!         if(sqrt(vec_cac(1)*vec_cac(1) + vec_cac(2)*vec_cac(2) + vec_cac(3)*vec_cac(3)).lt.1.d0) then
!            core=.true.
!            exit
!         endif
!         if(zdum.lt.-zmax) exit
!         iz=iz+1
!      enddo
!      if(indx_max.lt.iz) indx_max=iz
!!
!      nz_ray=iz
!!
!!------------------------allocate arrays--------------------------------
!!
!      call allocate_fs1d
!!
!      if(allocated(zdum_ray)) deallocate(zdum_ray)
!      allocate(zdum_ray(nz_ray))
!
!      z_ray(1)=zmax
!      do i=2, nz_ray
!         z_ray(i)=z_ray(i-1)-delz
!      enddo
!      !reverse array
!      do i=1, nz_ray
!         zdum_ray(i)=z_ray(nz_ray+1-i)
!      enddo
!      z_ray=zdum_ray
!!
!      if(allocated(vvelx_ray)) deallocate(vvelx_ray)
!      if(allocated(vvely_ray)) deallocate(vvely_ray)
!      if(allocated(vvelz_ray)) deallocate(vvelz_ray)
!      allocate(vvelx_ray(nz_ray), stat=err)
!         if(err.ne.0) stop 'allocation error check_interpolation: vvelx_ray'
!      allocate(vvely_ray(nz_ray), stat=err)
!         if(err.ne.0) stop 'allocation error check_interpolation: vvvely_ray'
!      allocate(vvelz_ray(nz_ray), stat=err)
!         if(err.ne.0) stop 'allocation error check_interpolation: vvelz_ray'
!!
!!-----------------------------------------------------------------------
!!
!      do i=1, nz_ray
!!
!      !calculate z_ray in carthesian coordinates
!         vec_cyc(1)=p(indx_p)*cos(zeta(indx_zeta))
!         vec_cyc(2)=p(indx_p)*sin(zeta(indx_zeta))
!         vec_cyc(3)=z_ray(i)
!         vec_cac=matmul(transmat, vec_cyc)
!         rad=sqrt(vec_cac(1)*vec_cac(1)+vec_cac(2)*vec_cac(2)+vec_cac(3)*vec_cac(3))
!         rad3d_ray(indx_p, indx_zeta, i)=rad
!!
!      !check if point of ray lies within region where information is stored
!         call info_region(vec_cac(1), vec_cac(2), vec_cac(3), rmin, rmax, linfo_phot, linfo_max, ldum)
!!
!      !analytic calculation only if point lies within region where information is stored
!         if(linfo_phot.and.linfo_max) then
!!
!!velocity components
!            call bvel3d(vmin*1.d5, vmax*1.d5, beta, vec_cac(1), vec_cac(2), vec_cac(3), dum_velx, dum_vely, dum_velz, dum_gradv)
!            dum_vel=sqrt(dum_velx*dum_velx + dum_vely*dum_vely + dum_velz*dum_velz)
!!
!            vvelx_ray(i)=dum_velx
!            velx3d_ray(indx_p, indx_zeta, i)=dum_velx
!            vvely_ray(i)=dum_vely
!            vely3d_ray(indx_p, indx_zeta, i)=dum_vely
!            vvelz_ray(i)=dum_velz
!            velz3d_ray(indx_p, indx_zeta, i)=dum_velz
!!
!!line opacity in units of 1/rstar
!            dum_opalbar=opalbar_model_hamann(sr, vinf, mdot, kappa0, alpha, vth_fid, rad, rho)
!          stop 'include stuff for opalbar_model_hamann'
!            opalbar_ray(i) = dum_opalbar
!            opalbar3d_ray(indx_p,indx_zeta,i) = dum_opalbar
!!
!!line source function
!         sline_ray(i) = sobo1d(r1d(i)*sr, velr/vth, opalbar1d(i)/sr, tmin, xic1, xnue0, eps_line)
!         stop 'adapt sobolev input'
!            sline3d_ray(indx_p, indx_zeta, i)=dum_sline
!!
!!------set all quantities to zero, if outside information region--------
!!
!         else
!            sline_ray(i)=0.d0
!            sline3d_ray(indx_p, indx_zeta,i)=0.d0
!            opalbar_ray(i)=0.d0
!            opalbar3d_ray(indx_p, indx_zeta,i)=0.d0
!            velz_ray(i)=0.d0
!            velx3d_ray(indx_p, indx_zeta, i)=0.d0
!            vely3d_ray(indx_p, indx_zeta, i)=0.d0
!            velz3d_ray(indx_p, indx_zeta, i)=0.d0
!         endif
!!
!      enddo
!!
!   enddo
!enddo
!!
!open(1, file=trim(outputdir_fs)//'/forsol_dime.dat', form='formatted')
!   write(1,'(3a15)') 'np', 'nzeta', 'indx_max'
!   write(1,'(3i15)') np, nzeta, indx_max
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_pgrid.dat', form='unformatted')
!   write(1) p
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_zetagrid.dat', form='unformatted')
!   write(1) zeta
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_sline.dat', form='unformatted')
!   write(1) sline3d_ray
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_opalbar.dat', form='unformatted')
!   write(1) opalbar3d_ray
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_velx.dat', form='unformatted')
!   write(1) velx3d_ray
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_vely.dat', form='unformatted')
!   write(1) vely3d_ray
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_velz.dat', form='unformatted')
!   write(1) velz3d_ray
!close(1)
!!
!open(1, file=trim(outputdir_fs)//'/forsol_rad.dat', form='unformatted')
!   write(1) rad3d_ray
!close(1)
!!
!
!
!!
!end subroutine check_t1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function boundary(vec)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), dimension(3), intent(in) :: vec
logical :: boundary
!
! ... local scalars
real(dp) :: rad
!
rad=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
!
!note: possible rounding errors
if(rad.lt.1.d0+1.d-14) then
   boundary=.true.
else
   boundary=.false.
endif
!
end function boundary
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine vel_t3(xcoord, ycoord, zcoord, rad, velx, vely, velz, velr)
!
!-----------calculate velocity components for velocity-law--------------
!----------------v_phi(r,theta) = vrot * sin(theta)/r-------------------
!----------------v_r(r) = beta-velocity-law-----------------------------
!
!note: xcoord = r*sin(theta)*cos(phi)
!      vel_x = velr*sin(theta)*cos(phi) - velphi*sin(phi)
!      vel_y = velr*sin(theta)*sin(phi) + velphi*cos(phi)
!      vel_z = velr*cos(theta)
!   can be reduced to velx = velr/r * xcoord - vrot/r^2 * ycoord
!
!
use prog_type
use params_spec, only: vrot, vth_fiducial, vmin, vmax
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xcoord, ycoord, zcoord, rad, velr
real(dp) :: velx, vely, velz
!
! ... local scalars
real(dp) :: dum_r, dum_phi, dum_velr
!
! ... local functions
!
!in units of vthfiducial
dum_r = velr / rad
dum_phi = vrot*1.d5/rad/rad/vth_fiducial
!
!for rigid body rotation:
!dum_phi = vrot*1.d5/vthfiducial
!
velx = dum_r*xcoord - dum_phi*ycoord
vely = dum_r*ycoord + dum_phi*xcoord
velz = dum_r*zcoord
!
!write(*,*) rad, sqrt(velx**2 + vely**2 + velz**2)
!
end subroutine vel_t3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine vel_phot(xcoord, ycoord, zcoord, rad, velx, vely, velz)
!
!-----------calculate velocity components for velocity-law--------------
!----------------vel_r = vmin
!----------------v_phi(r,theta) = vrot * sin(theta)/r-------------------
!
use prog_type
use params_spec, only: vmin, vrot, vth_fiducial
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xcoord, ycoord, zcoord, rad
real(dp) :: velx, vely, velz
!
! ... local scalars
real(dp) :: dum_phi, dum_r
!
! ... local functions
!
!in units of vthfiducial
!dum_r = vmin/rad/vth_fiducial   !if radial velocity is included for photosphere
dum_r = 0.d0                    !if radial velocity is neglected for photosphere
dum_phi = vrot*1.d5/rad/rad/vth_fiducial
!
velx = dum_r*xcoord - dum_phi*ycoord
vely = dum_r*ycoord + dum_phi*xcoord
velz = dum_r*zcoord
!
!
end subroutine vel_phot
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine test_sobolev
!!
!use prog_type
!use forsol, only: rmax, rmin
!use params_fs, only: vthfiducial, xnue0, xic1, sr, beta
!!
!!
!implicit none
!!
!! ... local scalars
!integer(i4b), parameter :: nd=100
!integer(i4b) :: i
!real(dp) :: xcoord, ycoord, zcoord, velx, vely, velz, velr, gradv, opalbar
!real(dp) :: delr
!!
!! ... local arrays
!real(dp), dimension(nd) :: r, vel1d, gradv1d, opalbar1d, ssobo1d
!!
!! ... local functions
!real(dp) :: opalbar_model_hamann
!!
!!--------make arbitrary radial grid (here: equidistant in log(r))-------
!!
!delr=(log10(rmax)-log10(rmin))/(nd-1)
!r(1)=rmin
!do i=2, nd
!   r(i) = r(i-1) * 10.d0**delr
!enddo
!!
!!--------------------calculate physical values--------------------------
!!
!!calculate opacities, velocities and sline along radial grid
!do i=1, nd
!   xcoord=0.d0
!   ycoord=0.d0
!   zcoord=r(i)
!!
!!velocities in cm/s, gradv in 1/s
!   call bvel3d(vmin*1.d5, vmax*1.d5, beta, xcoord, ycoord, zcoord, velx, vely, velz, gradv)
!   velr=sqrt(velx*velx + vely*vely + velz*velz)
!   vel1d(i)=velr
!   gradv1d(i)=gradv
!!
!!line opacity in units of 1/rstar
!   opalbar = opalbar_model_hamann(sr, vinf, mdot, kappa0, alpha, vth_fid, rad, rho)
!   stop 'include vinf etc in model_hamann'
!   opalbar1d(i) = opalbar/sr
!!
!!line source function
!         ssobo1d(i) = sobo1d(r1d(i)*sr, velr/vth, opalbar1d(i)/sr, tmin, xic1, xnue0, eps_line)
!         stop 'adapt sobolev input'
!!
!enddo
!
!open(1, file='trash/sobo_forsol.dat')
!do i=1, nd
!   write(1,fmt='(5(e20.8))') r(i), vel1d(i)/vthfiducial, gradv1d(i)/vthfiducial, opalbar1d(i), &
!                             ssobo1d(i)*r(i)*r(i)/xic1
!enddo
!close(1)
!!
!end subroutine test_sobolev
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine subroutine_debug(indx_xobs)
!!
!use prog_type
!use fund_const, only: pi
!use dime_forsol, only:  nxobs_fs, np, nzeta, nr, nalpha, ngamma
!use dime_model_3d, only: ndxmax_fs, ndymax_fs, ndzmax_fs, x, y, z, opalbar3d
!use forsol, only: opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, &
!                  z_ray, nz_ray, p, pw, pw1, pw_err, zeta, zetaw, zetaw_err, xobs, xnue, xic_nue, &
!                  vth_ray, velz_ray, temp_ray, iin, iin_c, iem, iem_c, flux_tot, &
!                  flux_cont, ftot_errp, ftot_err, fcont_errp, fcont_err, ftot_p, &
!                  fcont_p, normt, normt_p, relerr_totp, relerr_contp, relerr_cont, &
!                  relerr_tot, r, alpha, gamma, alpha_arr, gamma_arr, lcore, del_vel, del_vel2, vphot_proj
!use timing_forsol, only: ts_tot, te_tot, t_tot, ts_obsdir, te_obsdir, ttot_obsdir
!use params_fs, only: vthfiducial, vmicro, vth_min, sr, beta, xnue0, xic1, vmax, tmin, teff, xlogg, &
!                     yhe
!!
!implicit none
!!
!! ... arguments
!integer(i4b), intent(in) :: indx_xobs
!!
!! ... local scalars
!integer(i4b) :: indx_zeta, indx_p, indx_zray
!real(dp) :: phinorm
!real(dp) :: iemi, iabs
!!
!! ... local arrays
!real(dp), dimension(:,:), allocatable :: iem_debug, iemc_debug, iin_debug, vz_proj
!real(dp), dimension(:), allocatable :: prof_debug
!!
!
!!calculate transformation matrix for given alpha, gamma
!call calc_transmat
!!
!!call test_sobolev
! 
!call allocate_fluxes
!!
!flux_tot=0.d0
!flux_cont=0.d0
!normt=0.d0
!ftot_err=0.d0
!!
!!
!allocate(iem_debug(np, nzeta))
!allocate(iemc_debug(np, nzeta))
!allocate(iin_debug(np, nzeta))
!allocate(vz_proj(np, nzeta))
!allocate(prof_debug(np))
!!
!!write(*,*) 'calculating xobs=', xobs(indx_xobs)
!!
!open(10, file='trash/debug_model01.dat', form='formatted')
!open(11, file='trash/debug_profile.dat', form='formatted')
!open(12, file='trash/zeta_integ.dat', form='formatted')
!!
!do indx_zeta=1, nzeta !39
!!
!   write(*,'(a30, i5, a2, i5)') 'calculating indx_zeta (nzeta)', indx_zeta, '/', nzeta
!   write(*,*) zeta(indx_zeta)
!   ftot_p=0.d0
!   fcont_p=0.d0
!   ftot_errp=0.d0
!   fcont_errp=0.d0
!   normt_p=0.d0
!!
!   do indx_p=1, np!146, 146 !np
!!      call setup_ray_1d(zeta(indx_zeta), p(indx_p), lcore)
!      call setup_ray(zeta(indx_zeta), p(indx_p), lcore)
!!      call setup_ray_2d(zeta(indx_zeta), p(indx_p), lcore)
!!
!      do indx_zray=1, nz_ray
!          call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs), phinorm)
!         profile_ray(indx_zray)=phinorm
!!         write(10,'(6es20.8)') z_ray(indx_zray), sqrt(z_ray(indx_zray)**2 + p(indx_p)**2), sline_ray(indx_zray), velz_ray(indx_zray), opalbar_ray(indx_zray), profile_ray(indx_zray)
!      enddo
!!
!      call get_iin(indx_xobs, lcore, iin, iin_c)
!!
!      call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)
!!      if(lcore) then 
!!         write(*,'(3es20.8)') iin, iem, iin-iem
!!      endif
!!p-integration
!      ftot_p(indx_xobs) = ftot_p(indx_xobs) + iem * p(indx_p) * pw(indx_p)
!      fcont_p(indx_xobs) = fcont_p(indx_xobs) + iem_c * p(indx_p) * pw(indx_p)
!      normt_p(indx_xobs) = normt_p(indx_xobs) + 1.d0 * pw(indx_p)
!!corresponding errors
!      ftot_errp(indx_xobs) = ftot_errp(indx_xobs) + iem * p(indx_p) * pw_err(indx_p)
!      fcont_errp(indx_xobs) = fcont_errp(indx_xobs) + iem_c * p(indx_p) * pw_err(indx_p)
!!
!      iem_debug(indx_p,indx_zeta) = iem
!      iemc_debug(indx_p, indx_zeta) = iem_c
!      iin_debug(indx_p, indx_zeta) = iin
!      vz_proj(indx_p, indx_zeta) = vphot_proj/2.7d8
!
!!      if(indx_p.gt.35.and.indx_p.lt.39) then
!!         write(*,'(i5, 4es20.8)') indx_p, p(indx_p), pw(indx_p), iin, ftot_p(indx_xobs)
!!      endif
!      write(11,*) p(indx_p), ftot_p(indx_xobs)
!!      write(*,'(4es20.8)') p(indx_p), iem_c, iem_c*p(indx_p)*pw(indx_p), fcont_p(indx_xobs)
!!      write(*,'(6es20.8)') zeta(indx_zeta)*180.d0/pi, p(indx_p), iem, iem_c, iem*p(indx_p)*pw(indx_p), iin
!!      write(*,'(9es20.8, l3)') p(indx_p), iem, iem_c, iem*p(indx_p), xic_nue(indx_xobs), xic1, ftot_p(indx_xobs), fcont_p(indx_xo!bs), iin, lcore
!   enddo
!!
!   write(12,*) zeta(indx_zeta), ftot_p(indx_xobs)
!   write(*,'(4es20.8)') ftot_p(indx_xobs), fcont_p(indx_xobs), ftot_p(indx_xobs)/fcont_p(indx_xobs), normt_p(indx_xobs)
!!
!!zeta-integration
!   flux_tot=flux_tot + zetaw(indx_zeta) * ftot_p
!   flux_cont=flux_cont + zetaw(indx_zeta) * fcont_p
!   normt=normt + zetaw(indx_zeta) * normt_p
!!corresponding errors
!   ftot_err=ftot_err + zetaw_err(indx_zeta) * ftot_errp
!   fcont_err=fcont_err + zetaw_err(indx_zeta) * fcont_errp
!!
!enddo
!!
!close(10)
!close(11)
!close(12)
!!
!write(*,*) 'at xobs', xobs(indx_xobs)
!write(*,'(3es20.8)') flux_tot(indx_xobs)/pi, flux_cont(indx_xobs)/pi, flux_tot(indx_xobs)/flux_cont(indx_xobs)
!write(*,*) np, nzeta
!!
!open(1, file='trash/pgrid.dat', form='unformatted')
!   write(1), p
!close(1)
!!
!open(1, file='trash/zgrid.dat', form='unformatted')
!   write(1), zeta
!close(1)
!!
!open(1, file='trash/iem_debug.dat', form='unformatted')
!   write(1), iem_debug
!close(1)
!open(1, file='trash/iemc_debug.dat', form='unformatted')
!   write(1), iemc_debug
!close(1)
!open(1, file='trash/iin_debug.dat', form='unformatted')
!   write(1), iin_debug
!close(1)
!open(1, file='trash/vz_proj.dat', form='unformatted')
!   write(1), vz_proj
!close(1)
!open(1, file='trash/dime_forsol.dat', form='formatted')
!   write(1,*) np, nzeta
!   write(1,*) xobs(indx_xobs)
!   write(1,*) alpha, gamma   
!close(1)
!!
!deallocate(iem_debug)
!deallocate(iemc_debug)
!deallocate(iin_debug)
!deallocate(vz_proj)
!!
!!
!end subroutine subroutine_debug
!
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine subroutine_debug_vproj(indx_xobs, indx_zeta, indx_p)
!!
!use prog_type
!use forsol, only: z_ray, nz_ray, velz_ray, p, zeta
!!
!implicit none
!!
!! ... arguments
!integer(i4b), intent(in) :: indx_xobs, indx_zeta, indx_p
!!
!! ... local scalars
!integer(i4b) :: i
!!
!! ... local arrays
!!
!! ... local logicals
!logical :: lcore
!
!!calculate transformation matrix for given alpha, gamma
!call calc_transmat
!!
!!call setup_ray(zeta(indx_zeta), p(indx_p), lcore)
!call setup_ray_2d_adm(zeta(indx_zeta), p(indx_p), lcore)
!!
!open(1, file='trash/vn.dat', form='formatted')
!   do i=1, nz_ray
!      write(1,*), z_ray(i), velz_ray(i)
!   enddo
!close(1)
!!
!end subroutine subroutine_debug_vproj
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine calc_iem_surface(xobs_surface, alpha_surface, gamma_surface)
!
use prog_type
use options_spec, only: input_mod
use dime_spec, only: np, nzeta, nxobs_fs
use mod_spectrum, only: alpha, gamma, xobs, zeta, p, vth_ray, velz_ray, z_ray, temp_ray, &
                  profile_ray, nz_ray, opalbar_ray, opac_ray, sline_ray, scont_ray
use params_spec, only: vth_fiducial, xic1
use mod_surfb, only: iem_surface, iemi_surface, iabs_surface, icont_surface
use mod_interp1d, only: find_index, interpol_yp
use mod_gdark, only: xic1_factor
!
implicit none
!
! ... arguments
real(dp), intent(in) :: alpha_surface, gamma_surface, xobs_surface
!
! ... local scalars
integer(i4b) :: indx_zeta, indx_p, indx_zray, indx_xobs1, indx_xobs2, iim2, iim1, ii, iip1, err
real(dp) :: iin, iin_c, iem1, iem2, iem_c1, iem_c2, phinorm, iemi, iemi1, iemi2, iabs, iabs1, iabs2
real(dp) :: xobs_surface1, xobs_surface2
!
! ... local arrays
!
! ... local logicals
logical :: lcore
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!
alpha=alpha_surface
gamma=gamma_surface
!
call find_index(xobs_surface, xobs, nxobs_fs, iim2, iim1, ii, iip1)
!
indx_xobs1=ii
indx_xobs2=iim1
!
xobs_surface1=xobs(indx_xobs1)
xobs_surface2=xobs(indx_xobs2)
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------calculating emergent intensity on surface-------------------'
write(*,*) 'for xobs', xobs_surface1, xobs_surface2
write(*,*) '   alpha', alpha
write(*,*) '   gamma', gamma
write(*,*) 'interpolation onto', xobs_surface
write(*,*)
!
!-----------------------------------------------------------------------
!
if(allocated(iem_surface)) deallocate(iem_surface)
allocate(iem_surface(np, nzeta), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iem_surface'
!
if(allocated(iemi_surface)) deallocate(iemi_surface)   
allocate(iemi_surface(np, nzeta), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iemi_surface'
!
if(allocated(iabs_surface)) deallocate(iabs_surface)   
allocate(iabs_surface(np, nzeta), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iabs_surface'
!
if(allocated(icont_surface)) deallocate(icont_surface)
allocate(icont_surface(np, nzeta), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: icont_surface'
!
!-----------------------------------------------------------------------
!
call calc_transmat
!
!-----------------------------------------------------------------------
!
!$omp parallel &
!$omp private(indx_zeta, indx_p, indx_zray), &
!$omp private(lcore, phinorm, iin, iin_c, iem1, iem_c1, iemi1, iabs1), &
!$omp private(iem2, iem_c2, iemi2, iabs2), &
!$omp copyin(xic1_factor)
!
!
!$omp do schedule(static, 4)
do indx_zeta=1, nzeta
!
   write(*,'(a30, i5, a2, i5)') 'calculating indx_zeta (nzeta)', indx_zeta, '/', nzeta
!
   do indx_p=1, np
!
!select case statement is allowed here, because not much do-loops are present, 
!   which might lower execution time considerable
      select case(input_mod)
         case(0)
            call setup_ray1d(zeta(indx_zeta), p(indx_p), lcore)
         case(1)
            call setup_ray3d(zeta(indx_zeta), p(indx_p), lcore)
!            call setup_ray3d_rot(zeta(indx_zeta), p(indx_p), lcore)
         case(2)
            call setup_ray3d_spc(zeta(indx_zeta), p(indx_p), lcore)            
         case default
            stop 'error in calc_iem_surface: input_mod not specified'
      end select
!
!------------------------calculation of indx_xobs1----------------------
!
!calculation of profile function for given xobs along ray
      do indx_zray=1, nz_ray
         call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs1), phinorm)
         profile_ray(indx_zray)=phinorm
      enddo
!
!setting intensity from core
      call get_iin(indx_xobs1, lcore, iin, iin_c)
!
!formal solution along ray
      call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem1, iem_c1, iemi1, iabs1)
!      write(*,*) z_ray
!      write(*,'(es20.8, l5, 3es20.8)') p(indx_p), lcore, iin, iin_c, iem1
!
!------------------------calculation of indx_xobs2----------------------
!
!calculation of profile function for given xobs along ray
      do indx_zray=1, nz_ray
         call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs2), phinorm)
         profile_ray(indx_zray)=phinorm
      enddo
!
!setting intensity from core
      call get_iin(indx_xobs2, lcore, iin, iin_c)
!
!formal solution along ray
      call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem2, iem_c2, iemi2, iabs2)
!
!------------------interpolation of both xobs-solutions-----------------
!
      iem_surface(indx_p,indx_zeta)=interpol_yp(xobs_surface1, xobs_surface2, iem1, iem2, xobs_surface)
      iemi_surface(indx_p,indx_zeta)=interpol_yp(xobs_surface1, xobs_surface2, iemi1, iemi2, xobs_surface)
      iabs_surface(indx_p,indx_zeta)=interpol_yp(xobs_surface1, xobs_surface2, iabs1, iabs2, xobs_surface)
      icont_surface(indx_p,indx_zeta)=interpol_yp(xobs_surface1, xobs_surface2, iem_c1, iem_c2, xobs_surface)

   enddo
!
enddo
!$omp end do
!$omp end parallel
!
!
!
end subroutine calc_iem_surface
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_int2d(xobs_2d, alpha_2d, gamma_2d)
!
use prog_type
use fund_const, only: pi
use options_spec, only: input_mod
use dime_spec, only: np, nzeta, nxobs_fs
use params_spec, only: vth_fiducial
use mod_spectrum, only: alpha, gamma, xobs, zeta, p, vth_ray, velz_ray, z_ray, temp_ray, &
                        profile_ray, nz_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, nz_ray
use mod_int2d, only: zcoord_2d, xcoord_2d, int_2d, tau_2d, nz_ray_max, iemi_2d, iabs_2d, vn_2d
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: alpha_2d, gamma_2d, xobs_2d
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: indx_zeta1, indx_zeta2, indx_p, indx_zray, indx_xobs1, indx_xobs2, nz_rest, &
                iim2, iim1, ii, iip1, err
real(dp) :: xobs_2d1, xobs_2d2, zeta_2d1, zeta_2d2, zeta_2d, del
real(dp) :: iin, iin_c, iem1, iem2, iem_c, phinorm, iemi, iemi1, iemi2, iabs, iabs1, iabs2
!
! ... local arrays
real(dp), dimension(:), allocatable :: int1d_1, iemi1d_1, iabs1d_1, tau1d_1, &
                                       int1d_2, iemi1d_2, iabs1d_2, tau1d_2
!
! ... local logicals
logical :: lcore
!
! ... local functions
!
!-----------------------------------------------------------------------
!
alpha=alpha_2d
gamma=gamma_2d
!
!----------------------------find xobs----------------------------------
!
call find_index(xobs_2d, xobs, nxobs_fs, iim2, iim1, ii, iip1)
!
indx_xobs1=ii
indx_xobs2=iim1
!
xobs_2d1=xobs(indx_xobs1)
xobs_2d2=xobs(indx_xobs2)
!
!----------------------------find zeta----------------------------------
!
zeta_2d=pi/4.d0
!
call find_index(zeta_2d, zeta, nzeta, iim2, iim1, ii, iip1)
indx_zeta1=ii
if(zeta_2d.eq.0.d0) then
   indx_zeta1=iim1
endif
zeta_2d1=zeta(indx_zeta1)
!
if(zeta_2d1.gt.pi) then
   stop 'error calc_int2d: zeta_2d1 is not allowed to be greater then 180'
endif
!
!zeta_2d=zeta_2d+pi
!
!call find_index(zeta_2d, zeta, nzeta, iim2, iim1, ii, iip1)
!indx_zeta2=ii
!zeta_2d2=zeta(indx_zeta2)
!!
!if(abs(zeta_2d2-zeta_2d1-pi).gt.1.d-8) then
!   stop 'error calc_int2d: zeta_2d2-zeta_2d1 neq 180'
!endif
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------calculating int2d and tau2d along ray-----------------------'
write(*,*) 'for xobs', xobs_2d1, xobs_2d2
write(*,*) '   alpha', alpha
write(*,*) '   gamma', gamma
write(*,*) 'interpolation onto', xobs_2d
write(*,*) '   zeta1', zeta_2d1*180.d0/pi
write(*,*) '   zeta2', zeta_2d2*180.d0/pi
write(*,*)
!
!-----------------------------------------------------------------------
!
call calc_transmat
!
!-----------------------------------------------------------------------
!
nz_ray_max=0
!
do indx_p=1, np
!
   select case(input_mod)
      case(0)
         call setup_ray1d(zeta(indx_zeta1), p(indx_p), lcore)
      case(1)
         call setup_ray3d(zeta(indx_zeta1), p(indx_p), lcore)
      case(2)
         call setup_ray3d_spc(zeta(indx_zeta1), p(indx_p), lcore)         
      case default
         stop 'error in calc_int2d: input_mod not specified'
   end select

   if(nz_ray.gt.nz_ray_max) then
      nz_ray_max=nz_ray
   endif
!
enddo
!
!-----------------------------------------------------------------------
!
if(allocated(zcoord_2d)) deallocate(zcoord_2d)
allocate(zcoord_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: zcoord_2d'
if(allocated(xcoord_2d)) deallocate(xcoord_2d)
allocate(xcoord_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: xcoord_2d'
if(allocated(int_2d)) deallocate(int_2d)
allocate(int_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: int_2d'
if(allocated(iemi_2d)) deallocate(iemi_2d)
allocate(iemi_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: iemi_2d'
if(allocated(iabs_2d)) deallocate(iabs_2d)
allocate(iabs_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: iabs_2d'
if(allocated(tau_2d)) deallocate(tau_2d)
allocate(tau_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: tau_2d'
if(allocated(vn_2d)) deallocate(vn_2d)
allocate(vn_2d(np,nz_ray_max), stat=err)
   if(err.ne.0) stop 'allocation error calc_int2d: vn_2d'
!
zcoord_2d=0.d0
int_2d=0.d0
iemi_2d=0.d0
iabs_2d=0.d0
tau_2d=0.d0
vn_2d=0.d0
!
!-----------------------------------------------------------------------
!
do indx_p=1, np
!
   select case(input_mod)
      case(0)
         call setup_ray1d(zeta(indx_zeta1), p(indx_p), lcore)
      case(1)
         call setup_ray3d(zeta(indx_zeta1), p(indx_p), lcore)
      case(2)
         call setup_ray3d_spc(zeta(indx_zeta1), p(indx_p), lcore)         
      case default
         stop 'error in calc_int2d: input_mod not specified'
   end select
!
!-----------------calculate everything for indx_xobs1-------------------
!
   do indx_zray=1, nz_ray
       call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs1), phinorm)
       profile_ray(indx_zray)=phinorm
   enddo
!
   call get_iin(indx_xobs1, lcore, iin, iin_c)
   !
   allocate(int1d_1(nz_ray), stat=err)
   allocate(iemi1d_1(nz_ray), stat=err)
   allocate(iabs1d_1(nz_ray), stat=err)
   allocate(tau1d_1(nz_ray), stat=err)
!
   call formal_ray2(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, int1d_1, iemi1d_1, iabs1d_1, tau1d_1)
!
!-----------------calculate everything for indx_xobs2-------------------
!
   do indx_zray=1, nz_ray
       call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs2), phinorm)
       profile_ray(indx_zray)=phinorm
   enddo
!
   call get_iin(indx_xobs2, lcore, iin, iin_c)
!
   allocate(int1d_2(nz_ray), stat=err)
   allocate(iemi1d_2(nz_ray), stat=err)
   allocate(iabs1d_2(nz_ray), stat=err)
   allocate(tau1d_2(nz_ray), stat=err)
!
   call formal_ray2(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, int1d_2, iemi1d_2, iabs1d_2, tau1d_2)
!
!--------------------interpolate onto xobs2d----------------------------
!
   do j=1, nz_ray
      int_2d(indx_p, nz_ray_max-nz_ray+j) = interpol_yp(xobs_2d1, xobs_2d2, int1d_1(j), int1d_2(j), xobs_2d)
      iemi_2d(indx_p, nz_ray_max-nz_ray+j) = interpol_yp(xobs_2d1, xobs_2d2, iemi1d_1(j), iemi1d_2(j), xobs_2d)
      iabs_2d(indx_p, nz_ray_max-nz_ray+j) = interpol_yp(xobs_2d1, xobs_2d2, iabs1d_1(j), iabs1d_2(j), xobs_2d)
      tau_2d(indx_p, nz_ray_max-nz_ray+j) = interpol_yp(xobs_2d1, xobs_2d2, tau1d_1(j), tau1d_2(j), xobs_2d)
   enddo
!
   deallocate(int1d_1)
   deallocate(iemi1d_1)
   deallocate(iabs1d_1)
   deallocate(tau1d_1)

   deallocate(int1d_2)
   deallocate(iemi1d_2)
   deallocate(iabs1d_2)
   deallocate(tau1d_2)

   xcoord_2d(indx_p, :) = p(indx_p)
   zcoord_2d(indx_p, nz_ray_max-nz_ray+1:nz_ray_max) = z_ray(:)
   vn_2d(indx_p, nz_ray_max-nz_ray+1:nz_ray_max) = velz_ray(:)
!
!note: need to set zcoord_2d to some dummy values, in order that
!      triangulation in idl-routine works
   zcoord_2d(indx_p, 1:nz_ray_max-nz_ray) = z_ray(1)-0.0001d0*abs(z_ray(1))

   if(p(indx_p).le.1.d0) then
      nz_rest=nz_ray_max-nz_ray
      del=(z_ray(nz_ray)+z_ray(1)-0.0001d0*z_ray(1))/float(nz_rest-2)
      zcoord_2d(indx_p, nz_rest) = z_ray(1)-0.0001d0*z_ray(1)
!      write(*,*) nz_rest, zcoord_2d(indx_p, nz_rest), zcoord_2d(indx_p, nz_rest+1)
      do j=1, nz_rest-2
!         write(*,*) j, z_ray(nz_ray)
         zcoord_2d(indx_p, j) = -z_ray(nz_ray) + (j-1)*del
      enddo
!      write(*,'(4es20.8)') zcoord_2d(indx_p, nz_rest-1), zcoord_2d(indx_p, nz_rest), zcoord_2d(indx_p, nz_rest+1), z_ray(1)
   endif
!
enddo
!
!
!
end subroutine calc_int2d
!
!
