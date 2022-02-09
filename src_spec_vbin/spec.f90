!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------------program----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!v02: line profiles for binary star systems
!
!vnew: including photospheric profiles for halpha, hbeta from coelho05
!      including global tilts of individual coordinate systems (calc_rotmat)
!      (including rotation axis of individual stars)
!      including options for setting up p-grids and r-grids
!      debugging subrutine formal_lc to combine continuum and line together
!
!-----------------------------------------------------------------------
!-----------------------start-------------------------------------------
!-----------------------------------------------------------------------
!
program spectrum3d
! 
use prog_type
use fund_const
use timing_spec
use options_spec, only: opt_photprof1, opt_photprof2, interp_photprof, opt_surface, opt_int2d, &
                        output_dir, input_mod, opt_pgrid01, opt_pgrid02, opt_rgrid01, opt_rgrid02
use dime_spec, only: nxobs_fs, nalpha, ngamma, ndirs
use mod_surfb, only: xobs_surface, alpha_surface, gamma_surface
use mod_spectrum, only: del_vel, del_vel2, lcore1, lcore2, xobs, xnue, &
                        acoeff_xic1, bcoeff_xic1, ccoeff_xic1, dcoeff_xic1, &
                        acoeff_xic2, bcoeff_xic2, ccoeff_xic2, dcoeff_xic2, &
                        xic1_nue, xicc1_nue, xic2_nue, xicc2_nue, &
                        flux_tot, flux_emi, flux_abs, flux_cont, normt, &
                        phinorm, iin, iin_c, iem, iem_c, iemi, iabs, &
                        nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, temp_ray, profile_ray, vth_ray, velz_ray, &
                        alpha, gamma, alpha_arr, gamma_arr
use mod_triangles, only: npoints, points_xcoord, points_ycoord, points_weight
use params_spec, only: logg1, yhe1, fehe1, aenh1, teff1, trad1, xic1, vmicro1, &
                       logg2, yhe2, fehe2, aenh2, teff2, trad2, xic2, vmicro2, &
                       xnue0, tmin, iline, na, vth_fiducial, vth_min
use mod_gdark, only: xic1_factor, xic2_factor
use omp_lib
!
implicit none
!
! ... local scalars
integer(i4b) :: indx_point, indx_xobs, indx_zray, indx_alpha, indx_gamma, indx_dir
integer(i4b) :: i, err
integer(i4b) :: ticks_sec, ticks_max
real(dp) :: err1, err2, sum, sum2
real(dp) :: cs1_x, cs1_y, cs2_x, cs2_y
real(dp) :: ts_tot, te_tot, ts_obsdir, te_obsdir
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
ts_tot = omp_get_wtime()
!
!
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
!
!--------------------3d model in cartesian coordinates------------------
!
   case(1) 
!
!--------------------3d model in spherical coordinates------------------
!
   case(2)
      call read_model3d_spc
      call print_model3d_spc
!
   case default
      stop 'input model not specified'

!
end select
!
!-------------------calculating rotation matrices-----------------------
!
call calc_rotmat
!
!--------------------test procedure for a single p-ray------------------
!
!call test_pray
!
!------------------------setting up xobs grid---------------------------
!
!calculate thermal velocity (which is taken at tmin for given mass number)
vth_min=vthermal(min(vmicro1,vmicro2), tmin, na)
!
!calculate maximum alowed velocity steps in fiducial vthermal
del_vel2=del_vel*vth_min/vth_fiducial
!
call grid_xobs
!
!--------------------setting up photospheric profile--------------------
!
select case(opt_photprof1)
   case(0)
!calculate photospheric profile on xobs-grid (setting it to xic1)
      call calc_photprof(iline, nxobs_fs, xic1_nue, xicc1_nue, xic1)
   case(1)
      call get_photprof_herrero(iline, teff1, logg1, yhe1, nxobs_fs, xnue, xic1_nue, xicc1_nue, xnue0, trad1)
!   case(2)
!      call get_photprof(teff1, logg1, yhe1, nxobs_fs, xnue, xic1_nue, xicc1_nue, xnue0, xic1)
   case(3)
      call get_photprof_fastwind(iline, teff1, logg1, yhe1, nxobs_fs, xnue, xic1_nue, xicc1_nue, xnue0, trad1)
   case(4)
      call get_photprof_coelho05(iline, teff1, logg1, fehe1, aenh1, nxobs_fs, xnue, xic1_nue, xicc1_nue, xnue0, trad1)
   case(5)
      call get_photprof_coelho14(iline, teff1, logg1, fehe1, aenh1, nxobs_fs, xnue, xic1_nue, xicc1_nue, xnue0, trad1)                  
   case default
      stop 'option opt_photprof1 not properly selected'
end select
!
select case(opt_photprof2)
   case(0)
!calculate photospheric profile on xobs-grid (setting it to xic1)
      call calc_photprof(iline, nxobs_fs, xic2_nue, xicc2_nue, xic2)
   case(1)
      call get_photprof_herrero(iline, teff2, logg2, yhe2, nxobs_fs, xnue, xic2_nue, xicc2_nue, xnue0, trad2)
!   case(2)
!      call get_photprof(teff2, logg2, yhe2, nxobs_fs, xnue, xic2_nue, xicc2_nue, xnue0, xic2)
   case(3)
      call get_photprof_fastwind(iline, teff2, logg2, yhe2, nxobs_fs, xnue, xic2_nue, xicc2_nue, xnue0, trad2)
   case(4)
      call get_photprof_coelho05(iline, teff2, logg2, fehe2, aenh2, nxobs_fs, xnue, xic2_nue, xicc2_nue, xnue0, trad2)
   case(5)
      call get_photprof_coelho14(iline, teff2, logg2, fehe2, aenh2, nxobs_fs, xnue, xic2_nue, xicc2_nue, xnue0, trad2)               
   case default
      stop 'option opt_photprof2 not properly selected'
end select
!
!output photospheric profile to be done
call output_photprof(output_dir, 'photprof_star01.dat', nxobs_fs, xobs, xic1_nue, xicc1_nue)
call output_photprof(output_dir, 'photprof_star02.dat', nxobs_fs, xobs, xic2_nue, xicc2_nue)
!
!
select case(interp_photprof)
   case(0)
      write(*,*) 'linear interpolation of photospheric profile'
   case(1)
      write(*,*) 'monotonic cubic spline interpolation of photospheric profile'
      allocate(acoeff_xic1(nxobs_fs), stat=err)
      allocate(acoeff_xic2(nxobs_fs), stat=err)
      allocate(bcoeff_xic1(nxobs_fs), stat=err)
      allocate(bcoeff_xic2(nxobs_fs), stat=err)
      allocate(ccoeff_xic1(nxobs_fs), stat=err)
      allocate(ccoeff_xic2(nxobs_fs), stat=err)      
      allocate(dcoeff_xic1(nxobs_fs), stat=err)
      allocate(dcoeff_xic2(nxobs_fs), stat=err)      
      call cube_mono_coeff(nxobs_fs, xobs, xic1_nue, acoeff_xic1, bcoeff_xic1, ccoeff_xic1, dcoeff_xic1)
      call cube_mono_coeff(nxobs_fs, xobs, xic2_nue, acoeff_xic2, bcoeff_xic2, ccoeff_xic2, dcoeff_xic2)      
   case default
      stop 'option interp_photprof not selected'
end select
!
!----------------------test if xobs grid large enough-------------------
!
call check_grid_xobs
!
!----------------------calculate gravity darkening----------------------
!
call calc_gdark
!
!------------------------setting up p-ray grid--------------------------
!-------------------(individually for each star)------------------------
!
call grid_pray(opt_pgrid01, opt_pgrid02)
!
!------------------------setting up zeta grid---------------------------
!
call grid_zeta
!
!------------------------setting radial grid----------------------------
!
call grid_radial(opt_rgrid01, opt_rgrid02)
!
!---setting up alpha, gamma-arrays (specifying direction to observer)---
!
call grid_obsdir
!
!
!------calculate emergent intensity for a given direction and xobs------
!-------------(can not be included in parallel region!!!)---------------
!
if(opt_surface) then
   call calc_iem_surface(xobs_surface, alpha_surface, gamma_surface)      
   call output_surface(xobs_surface, alpha_surface, gamma_surface)
   call output_triangles(1)   
   stop
endif
!
!-----------calculate intensity and optical depths in 2d:---------------
!------------------slice along a given direction------------------------
!
!if(opt_int2d) then
!   call calc_int2d(xobs_surface, alpha_surface, gamma_surface)
!   call output_int2d(xobs_surface, alpha_surface, gamma_surface)
!   stop
!endif
!-----------------------------------------------------------------------
!
!$omp parallel &
!$omp private(indx_dir, indx_alpha, indx_gamma, indx_point, indx_xobs, indx_zray, &
!$omp         cs1_x, cs1_y, cs2_x, cs2_y, ts_obsdir, te_obsdir), &
!$omp copyin(xic1_factor,xic2_factor)
!$omp do schedule(dynamic)
do indx_dir=1,ndirs
!
!--------------------------set up timing--------------------------------
!
   ts_obsdir=omp_get_wtime()
   call timing_init   
!
!----------calculate the transformation matrix and triangulation--------
!-----------------------for each direction------------------------------
!
   call conv_indx_1d_to_2d (indx_dir, ngamma, indx_gamma, indx_alpha)      
!
   alpha=alpha_arr(indx_alpha)
   gamma=gamma_arr(indx_gamma)
!
   write(*,*) '***************************calculating new direction***************************'
   write(*,*) 'alpha, gamma=', alpha, gamma
   write(*,*)
      
!define transformation matrix
   call calc_transmat
!
!get the triangulation
   call grid_triangles
!
!test the integration weights      
   call check_weights
!
!----------------set up a fluxes for a given direction------------------
!
!allocation of global (threadprivate) arrays
   call allocate_fluxes
!
!   do indx_point=15584, 15584!npoints
   do indx_point=1, npoints
!
!      write(*,'(a35, i8, a2, i8)') 'calculating indx_point (npoints)', indx_point, '/', npoints
!
      select case(input_mod)            
!         case(0)
!         case(1)
         case(2)
!set up the z-ray in system of star 1, star 2 and in global system
            call setup_ray3d_cs1spc(indx_point, cs1_x, cs1_y, lcore1)
            call setup_ray3d_cs2spc(indx_point, cs2_x, cs2_y, lcore2)   
            call setup_ray3d_cs0spc(indx_point, cs1_x, cs1_y, cs2_x, cs2_y, lcore1, lcore2)
         case default
            stop 'error in spec main: input_mod not specified'
      end select
!
!debug start
!      do indx_zray=1, nz_ray
!         write(*,*) indx_zray, nz_ray, opalbar_ray(indx_zray)
!         if(opalbar_ray(indx_zray).ne.opalbar_ray(indx_zray)) stop 'error opalbar'
!         if(opac_ray(indx_zray).ne.opac_ray(indx_zray)) stop 'error opac'      
!         if(vth_ray(indx_zray).ne.vth_ray(indx_zray)) stop 'error vth'
!         if(velz_ray(indx_zray).ne.velz_ray(indx_zray)) stop 'error velz'
!         if(sline_ray(indx_zray).ne.sline_ray(indx_zray)) stop 'error sline'
!         if(scont_ray(indx_zray).ne.scont_ray(indx_zray)) stop 'error scont'
!         if(z_ray(indx_zray).ne.z_ray(indx_zray)) stop 'error z'
!      enddo
!debug end
      
      do indx_xobs=1, nxobs_fs
!
         !calculation of profile function for given xobs along ray
         do indx_zray=1, nz_ray
            call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs), phinorm)
            profile_ray(indx_zray)=phinorm
!            write(*,'(10es20.8)') velz_ray(indx_zray), profile_ray(indx_zray), vth_ray(indx_zray), vth_fiducial
         enddo
!
        !setting intensity from core
         call get_iin(indx_xobs, lcore1, lcore2, iin, iin_c)
!         write(*,*) indx_point, lcore1, lcore2, iin, iin_c
!
         !formal solution along ray
         call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)

!debug
!         if(opalbar_ray(1).gt.1.d0) then
!            do i=1, nz_ray
!               write(*,*) i, z_ray(i), opalbar_ray(i), opac_ray(i), sline_ray(i), scont_ray(i)
!            enddo
!            stop
!         endif
!
!p-integration
         flux_tot(indx_xobs) = flux_tot(indx_xobs) + iem * points_weight(indx_point)
         flux_cont(indx_xobs) = flux_cont(indx_xobs) + iem_c * points_weight(indx_point)
         normt(indx_xobs) = normt(indx_xobs) + one * points_weight(indx_point)
         flux_emi(indx_xobs) = flux_emi(indx_xobs) + iemi * points_weight(indx_point)
         flux_abs(indx_xobs) = flux_abs(indx_xobs) + iabs * points_weight(indx_point)
!
      !given p-ray is calculated for all xobs
      enddo
!
   !all p-rays are calculated for all xobs
   enddo
!
   call print_flux
   call output_triangles(indx_dir)
   call output_fluxem(indx_dir)
   call output_fluxem_debug(indx_dir)
   
!deallocate all arrays that have been allocated in parallel region
   call deallocate_fluxes
   call deallocate_fs1d
!
!deallocate triangles that have been allocated in parallel region
   call deallocate_triangles      

   te_obsdir=omp_get_wtime()
   ttot_obsdir=te_obsdir-ts_obsdir
   call print_timing

   
!all angles calculated
enddo
!$omp end do
!$omp end parallel   

!
te_tot=omp_get_wtime()
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
use mod_spectrum, only: del_xobs, xobs, xnue, xic1_nue, xicc1_nue, xic2_nue, xicc2_nue
use params_spec, only: xnue0, vmax, vth_fiducial, vth_min, xlim, vmicro1, vmicro2
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
!
write(*,*) '---------------------calculating frequency grid--------------------------------'
write(*,*)
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
!write(*,*) 
!stop 'go on in spec'
!
!allocate xobs-array, xnue-array, xic-nue-array
allocate(xobs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xobs'
allocate(xnue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xnue'
allocate(xic1_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xic1_nue'
allocate(xicc1_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xicc1_nue'
allocate(xic2_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xic2_nue'
allocate(xicc2_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xicc2_nue'
!
!calculate observers frame frequency grid
!
do i=1, nxobs_fs
   xobs(i) = -xmax + (i-1)*2*xmax/(nxobs_fs-1)
   xnue(i) = xobs(i)*deldop_fiducial + xnue0
enddo
!
write(*,'(a20,es20.8)') 'vth_min', vth_min
write(*,'(a20,es20.8)') 'vth_fiducial', vth_fiducial
write(*,'(a20,es20.8)') 'delta', delta
write(*,'(a20,es20.8)') 'xmax', xmax
write(*,'(a20,i20)') 'nxobs_fs', nxobs_fs
write(*,*)
write(*,*) 'xobs-grid:'
write(*,'(8es20.8)') xobs
write(*,*)
!
!write(*,*) nxobs_fs, xobs(1)
!stop 'go on in grid_xobs'
!
end subroutine grid_xobs

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_grid_xobs
!
!------------check if photospheric profile is completely resolved-------
!
use prog_type
use fund_const, only: cgs_clight, cgs_kb, cgs_mp
use dime_spec, only: nxobs_fs
use mod_spectrum, only: del_xobs, xobs, xnue, xic1_nue, xicc1_nue, xic2_nue, xicc2_nue
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
!
!
!can uncomment these lines if the left and right edges are contaminated by another line
if(xic1_nue(1)/xicc1_nue(1).lt.0.95d0) stop 'error in check_grid_xobs: xlim too small'
if(xic1_nue(nxobs_fs)/xicc1_nue(nxobs_fs).lt.0.95d0) stop 'error in check_grid_xobs: xlim too small'
if(xic2_nue(1)/xicc2_nue(1).lt.0.95d0) stop 'error in check_grid_xobs: xlim too small'
if(xic2_nue(nxobs_fs)/xicc2_nue(nxobs_fs).lt.0.95d0) stop 'error in check_grid_xobs: xlim too small'


!write(*,*) xobs
!stop
!
end subroutine check_grid_xobs
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
use dime_spec, only: cs1_nzeta, cs2_nzeta
use mod_spectrum, only: cs1_zeta, cs2_zeta
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
!
!
allocate(cs1_zeta(cs1_nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: cs1_zeta'
allocate(cs2_zeta(cs2_nzeta), stat=err)
   if(err.ne.0) stop 'allocation error grid_zeta: cs2_zeta'
!
!equidistant
do i=1, cs1_nzeta
   cs1_zeta(i) = float(i-1)*two*pi/float(cs1_nzeta-1)
enddo
!
do i=1, cs2_nzeta
   cs2_zeta(i) = float(i-1)*two*pi/float(cs2_nzeta-1)
enddo
!
write(*,*) '--------------------------calculating zeta grids-------------------------------'
write(*,*)
write(*,*) 'cs1_zeta-grid:'
write(*,'(8es20.8)') cs1_zeta*180./pi
write(*,*)
write(*,*) 'cs2_zeta-grid:'
write(*,'(8es20.8)') cs2_zeta*180./pi
write(*,*)
!
!
!
end subroutine grid_zeta
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_pray(opt_pgrid01, opt_pgrid02)
!
!-----------------sets up grid of impact parameter----------------------
!
use prog_type
use fund_const
use dime_spec, only: cs1_np, cs1_np_nc, cs1_np_c, &
                     cs2_np, cs2_np_nc, cs2_np_c
use mod_spectrum, only: cs1_p, cs2_p
use params_spec, only: rmax1, rmax2
!
implicit none
!
! ... arguments
character(len=5), intent(in) :: opt_pgrid01, opt_pgrid02
!
! ... local scalars
integer(i4b) :: i,j,k,ip
integer(i4b) :: err
integer(i4b) :: np_dum
real(dp) :: del, sum, sum1, sum_err, sumex
!
! ... local arrays
real(dp), dimension(:), allocatable :: p_dum
!
!
write(*,*) '--------------------------calculating p-grid: star1----------------------------'
write(*,*)
!
!---------------------equidistant inside core---------------------------
!
if(cs1_np_c.lt.2) stop 'error grid_pray: cs1_np_c has to be ge 2'
!
allocate(cs1_p(cs1_np))
!
np_dum = cs1_np_c
allocate(p_dum(np_dum), stat=err)
!
do i=1, cs1_np_c
   p_dum(i)=one*float(i-1)/float(cs1_np_c-1)
enddo
!
ip=1
do i=1, cs1_np_c
   cs1_p(ip) = p_dum(i)
   ip=ip+1
enddo
!
!---------------------outside the core----------------------------------
!
if(cs1_np_nc.lt.2) stop 'error in grid_pray: np1_nc has to be ge 2'
!
deallocate(p_dum)
np_dum=cs1_np_nc
allocate(p_dum(np_dum), stat=err)
!
!additional point right above photosphere
p_dum(1) = one
p_dum(2) = one + 1.d-4
!
!linear grid spacing outside the core
if(trim(opt_pgrid01).eq.'lin') then
   del = (rmax1-p_dum(2))/float(cs1_np_nc-1)
   do i=3, cs1_np_nc
      p_dum(i) = p_dum(i-1)+del
   enddo
   p_dum(np_dum)=rmax1

!logarithmic increment outside core
elseif(trim(opt_pgrid01).eq.'log') then
   del = log10(rmax1/p_dum(2))/float(cs1_np_nc)
!
   do i=3, cs1_np_nc
      p_dum(i)=p_dum(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax1
   p_dum(np_dum)=rmax1
   !
!log-log-increment outside core
elseif(trim(opt_pgrid01).eq.'llog') then
   del = log10(rmax1)/log10(p_dum(2))
   del = log10(del)
   del = del/(cs1_np_nc)
!
   do i=3, cs1_np_nc
      p_dum(i) = 10.d0**(log10(p_dum(i-1))*10.d0**del)
   enddo
   p_dum(np_dum)=rmax1
else
   stop 'error in grid_pray: opt_pgrid01 not properly specified'
endif
!
do i=2, np_dum
   cs1_p(ip) = p_dum(i)
   ip=ip+1;
enddo
!
write(*,*) 'cs1_p-grid:'
write(*,'(8es20.8)') cs1_p
write(*,*)
!
!
write(*,*) '--------------------------calculating p-grid: star2----------------------------'
write(*,*)
!
!---------------------equidistant inside core---------------------------
!
if(cs2_np_c.lt.2) stop 'error grid_pray: cs2_np_c has to be ge 2'
!
allocate(cs2_p(cs2_np))
!
np_dum = cs2_np_c
deallocate(p_dum)
allocate(p_dum(np_dum), stat=err)
!
do i=1, cs2_np_c
   p_dum(i)=one*float(i-1)/float(cs2_np_c-1)
enddo
!
ip=1
do i=1, cs2_np_c
   cs2_p(ip) = p_dum(i)
   ip=ip+1
enddo
!
!---------------------logarithmic outside core--------------------------
!
if(cs2_np_nc.lt.2) stop 'error in grid_pray: np1_nc has to be ge 2'
!
deallocate(p_dum)
np_dum=cs2_np_nc
allocate(p_dum(np_dum), stat=err)
!
!additional point right above photosphere
p_dum(1) = one
p_dum(2) = one + 1.d-4
!
if(trim(opt_pgrid02).eq.'lin') then
   del = (rmax2-p_dum(2))/float(cs2_np_nc-1)
   do i=3, cs2_np_nc
      p_dum(i) = p_dum(i-1)+del
   enddo
   p_dum(np_dum)=rmax2

!logarithmic increment outside core
elseif(trim(opt_pgrid02).eq.'log') then
!
!logarithmic increment
   del = log10(rmax2/p_dum(2))/float(cs2_np_nc)
!
   do i=3, cs2_np_nc
      p_dum(i)=p_dum(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax2
   p_dum(np_dum)=rmax2
!
elseif(trim(opt_pgrid02).eq.'llog') then
!
!log-log-increment
   del = log10(rmax2)/log10(p_dum(2))
   del = log10(del)
   del = del/(cs2_np_nc)
!
   do i=3, cs2_np_nc
      p_dum(i) = 10.d0**(log10(p_dum(i-1))*10.d0**del)
   enddo
   p_dum(np_dum)=rmax2
!
else
   stop 'error in grid_pray: opt_pgrid02 not properly specified'
endif
!
do i=2, np_dum
   cs2_p(ip) = p_dum(i)
   ip=ip+1
enddo
!
write(*,*) 'cs2_p-grid:'
write(*,'(8es20.8)') cs2_p
write(*,*)
!
!-----------------------------------------------------------------------
!
!
!
end subroutine grid_pray
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_radial(opt_rgrid01, opt_rgrid02)
!
!----------set up radial grid for calculation of z-coordinates----------
!
use prog_type
use fund_const  
use dime_spec, only: cs1_nr, cs2_nr
use mod_spectrum, only: cs1_r, cs2_r
use params_spec, only: rmax1, rmax2
!
implicit none
!
! ... arguments
character(len=5), intent(in) :: opt_rgrid01, opt_rgrid02
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: del
!
!-----------------------------------------------------------------------
!
write(*,*) '--------------------------calculating r-grid: star 1---------------------------'
write(*,*)
!
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error grid_radial: cs1_r'
!
!-----------------------------------------------------------------------
!
!linear increment
if(trim(opt_rgrid01).eq.'lin') then
   del = (rmax1-one)/float(cs1_nr-1)
   do i=1, cs1_nr
      cs1_r(i) = one + (i-1)*del
   enddo
   cs1_r(cs1_nr)=rmax1

!logarithmic increment
elseif(trim(opt_rgrid01).eq.'log') then
!
!---------------equidistant in log-space outside core-------------------
!
!logarithmic increment
   del = log10(rmax1/one)/(cs1_nr)
!
   cs1_r(1)=one
   do i=2, cs1_nr
      cs1_r(i)=cs1_r(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax1
   cs1_r(cs1_nr)=rmax1
!
elseif(trim(opt_rgrid01).eq.'llog') then
!
!---------------------or equidsitant in log-logs------------------------
!
   cs1_r(1)=one
   cs1_r(2)=one+1.d-3
!
   del = log10(rmax1)/log10(cs1_r(2))
   del = log10(del)
   del = del/(cs1_nr-2)
!
   do i=3, cs1_nr
      cs1_r(i) = 10.d0**(log10(cs1_r(i-1))*10.d0**del)
   enddo
   cs1_r(cs1_nr)=rmax1
!
else
   stop 'error in grid_radial: opt_rgrid01 not properly specified'
endif
!
write(*,*) 'cs1_r-grid:'
write(*,'(8es20.8)') cs1_r
write(*,*)
!
!
!
write(*,*) '--------------------------calculating r-grid: star 2---------------------------'
write(*,*)
!
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error grid_radial: cs2_r'
!
!-----------------------------------------------------------------------
!
!linear increment
if(trim(opt_rgrid02).eq.'lin') then
   del = (rmax2-one)/float(cs2_nr-1)
   do i=1, cs2_nr
      cs2_r(i) = one + (i-1)*del
   enddo
   cs2_r(cs2_nr)=rmax2

!logarithmic increment
elseif(trim(opt_rgrid02).eq.'log') then
!
!---------------equidistant in log-space outside core-------------------
!
!logarithmic increment
   del = log10(rmax2/one)/(cs2_nr)
!
   cs2_r(1)=one
   do i=2, cs2_nr
      cs2_r(i)=cs2_r(i-1)*10.d0**del
   enddo
!
!set outermost point to exactly rmax2
   cs2_r(cs2_nr)=rmax2
!
elseif(trim(opt_rgrid02).eq.'llog') then
!
!---------------------or equidsitant in log-logs------------------------
!
   cs2_r(1)=one
   cs2_r(2)=one+1.d-3
!
   del = log10(rmax2)/log10(cs2_r(2))
   del = log10(del)
   del = del/(cs2_nr-2)
!
   do i=3, cs2_nr
      cs2_r(i) = 10.d0**(log10(cs2_r(i-1))*10.d0**del)
   enddo
   cs2_r(cs2_nr)=rmax2
!
else
   stop 'error in grid_radial: opt_rgrid02 not properly specified'
endif
!
write(*,*) 'cs2_r-grid:'
write(*,'(8es20.8)') cs2_r
write(*,*)
!
!
end subroutine grid_radial
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_radial_rot(lloglog, rmin, rmax, r, nr)
!
!----------set up radial grid for calculation of z-coordinates----------
!-----from r_min to r_max, to account for the distortion of the---------
!----stellar surface when large rotational velocities are present-------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nr
real(dp), intent(in) :: rmin, rmax
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
use dime_spec, only: nalpha, ngamma, ndirs
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
!total number of directions
ndirs=nalpha*ngamma
!
end subroutine grid_obsdir
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
use fund_const
use mod_gdark
use params_spec, only: logg1, sr1, rstar1, lstar1, vrot1, xnue0, teff1, trad1, xic1
use options_spec, only: opt_incl_gdark1, opt_incl_sdist1
use mod_spectrum, only: xic1_nue, xicc1_nue
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: ntheta_fine=101
integer(i4b) :: i, iim2, iim1, ii, iip1
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, l_star_cgs, w0, omega, sigma, sigma_fit, c1, c2, sum, mstar
!
! ... local arrays
real(dp), dimension(ntheta_fine) :: theta_fine, rsurface_fine, gr_fine, gtheta_fine, gperp_fine, integrand_fine, teff_fine
!
! ... local functions
real(dp) :: interpol_yp, bnue, calc_req
!
!
!
xic1_factor=one
xic2_factor=one
!
return

stop 'calc_gdark needs to be adapted for binary system'

!
!-----------------------------for star 1--------------------------------
!
!create theta-grid
do i=1, ntheta_fine
   theta_fine(i) = 1.d-5 + (i-1)*(pi/2.d0-1.d-5)/(ntheta_fine-1)
enddo
!
!
r_pole_cgs = sr1
m_star_cgs = sr1**2 * 10.d0**logg1/cgs_grav
mstar = m_star_cgs/xmsu
v_rot_cgs = vrot1
l_star_cgs = lstar1*xlsu
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
do i=1, ntheta1_gdark
   theta1_gdark(i) = c1*10.d0**(-float(i-1)/float(ntheta1_gdark-1)) + c2
   call find_index(theta1_gdark(i), theta_fine, ntheta_fine, iim2, iim1, ii, iip1)
   teff1_gdark(i) = interpol_yp(theta_fine(iim1), theta_fine(ii), teff_fine(iim1), teff_fine(ii), theta1_gdark(i))
   xic1_gdark(i) = bnue(xnue0, teff1_gdark(i))
enddo
!
if(.not.opt_incl_gdark1.and..not.opt_incl_sdist1) then
!if neither surface distortion nor gravity darkening is included, overwrite arrays
   teff1_gdark = trad1
   xic1_gdark = xic1
elseif(.not.opt_incl_gdark1) then
!if surface distortion is included and gravity darkening is excluded, calculate average effective temperatures
!
!------------------------------old version------------------------------
!
!   do i=1, ntheta_fine
!      integrand_fine(i) = 4.d0*pi*cgs_sb*(rsurface_fine(i)*r_pole_cgs)**2*sin(theta_fine(i))
!!      write(*,'(i5,10es20.8)') i, rsurface_fine(i), teff_fine(i), theta_fine(i), integrand_fine(i)
!   enddo
!!calculate average effective temperature
!   teff1=0.d0
!   do i=2, ntheta_fine
!      teff1 = teff1 + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
!   enddo
!   write(*,*) 'teff1 (old)', (l_star_cgs/teff1)**0.25d0
!
!------------------------------old version------------------------------
!
!account for distorted surface when calculating the flux integral
   do i=1, ntheta_fine
      integrand_fine(i) = 4.d0*pi*cgs_sb*(rsurface_fine(i)*r_pole_cgs)**2*sin(theta_fine(i))*gperp_fine(i)/(-gr_fine(i))
!      write(*,'(i5,10es20.8)') i, rsurface_fine(i), teff_fine(i), theta_fine(i), integrand_fine(i)
   enddo
!calculate average effective temperature
   teff1=zero
   do i=2, ntheta_fine
      teff1 = teff1 + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
   enddo
   write(*,*) 'teff1 (new)', (l_star_cgs/teff1)**0.25d0
!
!
!
   teff1=(l_star_cgs/teff1)**0.25d0
   teff_fine=teff1
   trad1=teff1
   xic1_nue=xic1_nue/xic1
   xicc1_nue=xicc1_nue/xic1
   xic1=bnue(xnue0,teff1)
   xic1_nue=xic1_nue*xic1
   xicc1_nue=xicc1_nue*xic1
   xic1_gdark=xic1
   teff1_gdark=teff1
endif
!
!-------------------calculate stellar surface distortion----------------
!
!calculate radius at equator
if(opt_incl_sdist1) then 
   smajorax1_a = calc_req(rstar1, mstar, vrot1)
   smajorax1_b = smajorax1_a
   smajorax1_c = one
else
   smajorax1_a = one
   smajorax1_b = one
   smajorax1_c = one
endif
!
!set scaling factor for surface intensity to one
!(required if non-rotating routines are used, and is adapted anyways when rotating routines are used)
xic1_factor=one
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
use mod_gdark, only: xic1_factor, ntheta1_gdark, theta1_gdark, xic1_gdark, &
                     xic2_factor, ntheta2_gdark, theta2_gdark, xic2_gdark
use params_spec, only: xic1, xic2, xnue0
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: x, y, z
!
! ... local scalars
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: rad, theta, xic1_local, xic2_local
!
! ... local functions
real(dp) :: bnue
!
!calculate co-latitude
rad = sqrt(x**2+y**2+z**2)
theta = acos(abs(z/rad))
!
!interpolate to obtain local xic1
call find_index(theta, theta1_gdark, ntheta1_gdark, iim2, iim1, ii, iip1)
xic1_local = interpol_yp(theta1_gdark(iim1), theta1_gdark(ii), xic1_gdark(iim1), xic1_gdark(ii), theta)
xic1_factor = xic1_local/xic1

!interpolate to obtain local xic2
call find_index(theta, theta2_gdark, ntheta2_gdark, iim2, iim1, ii, iip1)
xic2_local = interpol_yp(theta2_gdark(iim1), theta2_gdark(ii), xic2_gdark(iim1), xic2_gdark(ii), theta)
xic2_factor = xic2_local/xic2
!
!
!
end subroutine calc_xic_factor
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_iin(indx_xobs, lcore1, lcore2, iin, iin_c)
!
use prog_type
use fund_const, only: cgs_clight, zero
use options_spec, only: interp_photprof
use dime_spec, only: nxobs_fs
use params_spec, only: xnue0, xic1, xic2
use mod_spectrum, only: xobs, xnue, vphot_proj, &
                        xic1_nue, xicc1_nue, acoeff_xic1, bcoeff_xic1, ccoeff_xic1, dcoeff_xic1, &
                        xic2_nue, xicc2_nue, acoeff_xic2, bcoeff_xic2, ccoeff_xic2, dcoeff_xic2
use mod_gdark, only: xic1_factor, xic2_factor
use omp_lib
use mod_interp1d, only: find_index, interpol_yp, interpol_yp_spline
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_xobs
logical, intent(in) :: lcore1, lcore2
real(dp), intent(out) :: iin, iin_c
!
! ... local scalars
integer(i4b) ::  iim2, iim1, ii, iip1
real(dp) :: nue_cmf, xcmf, delta
!
! ... local functions


if(lcore1) then
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
!
   select case(interp_photprof)
      case(0) 
         iin=interpol_yp(xobs(iim1), xobs(ii), xic1_nue(iim1), xic1_nue(ii), xcmf) * xic1_factor
         iin_c=interpol_yp(xobs(iim1), xobs(ii), xicc1_nue(iim1), xicc1_nue(ii), xcmf)  * xic1_factor !for frequency dependent continuum
      case(1)
         iin=interpol_yp_spline(acoeff_xic1(ii), bcoeff_xic1(ii), &
                                ccoeff_xic1(ii), dcoeff_xic1(ii), &
                                xobs(ii), xcmf) * xic1_factor
      case default
         stop 'error get_iin: interp_photprof not defined'
   end select
!
!
elseif(lcore2) then
!note: need to interpolate from photospheric profile, since this is only given in cmf, 
!      whereas i am calculating in observers frame (shift of profile due to rotational velocities)
!calculate comoving frame frequency
   xcmf = xobs(indx_xobs) - vphot_proj
!
!constant continuum
   iin_c=xic2*xic2_factor
!
!find index within xobs-grid, for which xcmf matches and interpolate
   call find_index(xcmf, xobs, nxobs_fs, iim2, iim1, ii, iip1)
!
   select case(interp_photprof)
      case(0) 
         iin=interpol_yp(xobs(iim1), xobs(ii), xic2_nue(iim1), xic2_nue(ii), xcmf) * xic2_factor
         iin_c=interpol_yp(xobs(iim1), xobs(ii), xicc2_nue(iim1), xicc2_nue(ii), xcmf)  * xic2_factor !for frequency dependent continuum
      case(1)
         iin=interpol_yp_spline(acoeff_xic2(ii), bcoeff_xic2(ii), &
                                ccoeff_xic2(ii), dcoeff_xic2(ii), &
                                xobs(ii), xcmf) * xic2_factor
      case default
         stop 'error get_iin: interp_photprof not defined'
   end select

else
!set non-core ray
   iin=zero
   iin_c=zero
endif

!
end subroutine get_iin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_ray(xobs, cs0_x, cs0_y, iin, iin_c, fname)
!
use prog_type
use mod_spectrum, only: nhat, nz_ray, z_ray, opac_ray, opalbar_ray, velz_ray, &
                        scont_ray, sline_ray, temp_ray, profile_ray
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, cs0_x, cs0_y, iin, iin_c
character(len=*), intent(in) :: fname
!
! ... local scalars
integer(i4b) :: i
!
write(*,*) '-----------atmospheric structure along ray-------------'
!
open(1,file=trim(fname))
   write(1,'(a15, 3e20.8)') 'for nhat:    ', nhat
   write(1,'(a15, 2e20.8)') 'for x, y:    ', cs0_x, cs0_y
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
write(*,'(a15, 2e20.8)') 'for x, y:    ', cs0_x, cs0_y
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
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_flux
!
use prog_type
use fund_const, only: pi, half, two
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
real(dp) :: eqwidth
!
!xobs in units of vmax, insted of vthfiducial
!
write(*,*) '-------------------------emergent flux profile---------------------------------'
!
write(*,'(a4, 5(a20))') '#', 'xobs[vth*]', 'flux_tot', 'flux_cont', 'f_tot/f_cont', 'normt'
write(*,'(i4, 5(e20.8))')  i, xobs(1), flux_tot(1), flux_cont(1), &
                           flux_tot(1)/flux_cont(1), normt(1)
eqwidth = 0.d0
do i=2, nxobs_fs
   write(*,'(i4, 5(e20.8))')  i, xobs(i), flux_tot(i), flux_cont(i), &
                              flux_tot(i)/flux_cont(i), normt(i)
   eqwidth = eqwidth + half*(flux_tot(i-1)/flux_cont(i-1) + flux_tot(i)/flux_cont(i) - two)*(xobs(i)-xobs(i-1))
enddo
write(*,*) 'equivalent width', eqwidth
write(*,*)
!
end subroutine print_flux
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
use fund_const  
use dime_spec, only: nxobs_fs
use mod_spectrum, only: flux_tot, flux_cont, flux_emi, flux_abs, &
                        normt
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
if(allocated(flux_tot)) deallocate(flux_tot)
if(allocated(flux_cont)) deallocate(flux_cont)
if(allocated(normt)) deallocate(normt)
if(allocated(flux_emi)) deallocate(flux_emi)
if(allocated(flux_abs)) deallocate(flux_abs)
!
!
allocate(flux_tot(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_tot'
!
allocate(flux_cont(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_cont'
!
allocate(normt(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: normt'
!
allocate(flux_emi(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_emi'
!
allocate(flux_abs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_abs'
!
flux_tot=zero
flux_cont=zero
flux_emi=zero
flux_abs=zero      
normt=zero
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
use dime_spec, only: nxobs_fs
use mod_spectrum, only: flux_tot, flux_cont, flux_emi, flux_abs, normt
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
deallocate(flux_tot)
deallocate(flux_cont)
deallocate(normt)
deallocate(flux_emi)
deallocate(flux_abs)
!
!
!
end subroutine deallocate_fluxes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine timing_init
!
use prog_type
use fund_const
use timing_spec, only: ttot_obsdir, t_setup1, t_setup2, t_setup3, t_triangles
!
implicit none
!
ttot_obsdir=zero
t_setup1=zero
t_setup2=zero
t_setup3=zero
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
use dime_spec, only: nxobs_fs
use mod_triangles, only: npoints
use timing_spec, only: ttot_obsdir, t_setup1, t_setup2, t_setup3
!
implicit none
!
!
write(*,*) '----------------------------------timing---------------------------------------'
write(*,*) 'nxobs_fs: ', nxobs_fs
write(*,*) 'npoints:  ', npoints
write(*,*)
write(*,'(a40, 2e20.8)') 'total computation time', ttot_obsdir
write(*,'(a40, e20.8)') 'total time for setting ray of star1', t_setup1
write(*,'(a40, e20.8)') 'total time for setting ray of star2', t_setup2
write(*,'(a40, e20.8)') 'total time for setting ray globally', t_setup3
write(*,'(a40, e20.8)') 'total time for setting up ray', t_setup1+t_setup2+t_setup3
write(*,'(a40, e20.8)') 'average time for setting ray of star1', t_setup1/npoints
write(*,'(a40, e20.8)') 'average time for setting ray of star2', t_setup2/npoints
write(*,'(a40, e20.8)') 'average time for setting ray globally', t_setup3/npoints
write(*,'(a40, e20.8)') 'average time for setting up ray', (t_setup1+t_setup2+t_setup3)/npoints
write(*,*)

end subroutine print_timing
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine check_weights
!
!-------------------testing weights for triangular integration----------
!
use prog_type
use fund_const, only: pi
use mod_triangles, only: npoints, points_weight
use params_spec, only: rmax1, rstar1, rmax2, rstar2
use mod_spectrum, only: unit_length
!
implicit none
!
! ... local scalars
real(dp) :: sum0, sum1, sum2
!
! ... local arrays
!
!
write(*,*) '--------------------checking weights and error weights-------------------------'
write(*,*)
!
!integrating over the area of both stars
sum0=sum(points_weight)
!
!area of star 1
sum1=rmax1*rmax1*pi*(rstar1/unit_length)**2
!
!area of star 2
sum2=rmax2*rmax2*pi*(rstar2/unit_length)**2
!
write(*,'(a20, 5es20.8)') 'total area', sum0
write(*,'(a20, 5es20.8)') 'area star 1', sum1
write(*,'(a20, 5es20.8)') 'area star 2', sum2
write(*,*)
!
!
!
end subroutine check_weights
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_iem_surface(xobs_surface, alpha_surface, gamma_surface)
!
use prog_type
use options_spec, only: input_mod
use dime_spec, only: nxobs_fs
use mod_spectrum, only: alpha, gamma, xobs, vth_ray, velz_ray, z_ray, temp_ray, &
                  profile_ray, nz_ray, opalbar_ray, opac_ray, sline_ray, scont_ray
use params_spec, only: vth_fiducial, xic1, xic2
use mod_surfb, only: iem_surface, iemi_surface, iabs_surface
use mod_triangles, only: npoints, points_xcoord, points_ycoord
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: alpha_surface, gamma_surface
real(dp) :: xobs_surface
!
! ... local scalars
integer(i4b) :: i, indx_point, indx_zray, indx_xobs1, indx_xobs2, iim2, iim1, ii, iip1, err
real(dp) :: iin, iin_c, iem1, iem2, iem_c, phinorm, iemi, iemi1, iemi2, iabs, iabs1, iabs2
real(dp) :: xobs_surface1, xobs_surface2
real(dp) :: cs1_x, cs1_y, cs2_x, cs2_y
!
! ... local arrays
!
! ... local logicals
logical :: lcore1, lcore2
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
!
!define transformation matrix
call calc_transmat
!
!-----------------------------------------------------------------------
!
!get the triangulation for the given alpha and gamma
call grid_triangles
!
!
!-----------------------------------------------------------------------
!
allocate(iem_surface(npoints), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iem_surface'
!
allocate(iemi_surface(npoints), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iemi_surface'
!
allocate(iabs_surface(npoints), stat=err)
   if(err.ne.0) stop 'allocation-error in calc_iem_surface: iabs_surface'
!
!-----------------------------------------------------------------------
   !
!   do i=1, npoints
!      if(points_xcoord(i).lt.4.2.and.points_xcoord(i).gt.3.8) then
!         if(points_ycoord(i).lt.0.5.and.points_ycoord(i).gt.-0.5) then
!            write(*,*) i, points_xcoord(i), points_ycoord(i)
!         endif
!      endif
!   enddo
!   stop

!   do indx_point=3027, 3027!npoints
   do indx_point=1, npoints   
!
!   write(*,'(a35, i8, a2, i8)') 'calculating indx_point (npoints)', indx_point, '/', npoints
!
!set up the z-ray in system if star 1
   call setup_ray3d_cs1spc(indx_point, cs1_x, cs1_y, lcore1)
!set up the z-ray in system if star 2   
   call setup_ray3d_cs2spc(indx_point, cs2_x, cs2_y, lcore2)   
!set up the z-ray in global
   call setup_ray3d_cs0spc(indx_point, cs1_x, cs1_y, cs2_x, cs2_y, lcore1, lcore2)
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
   call get_iin(indx_xobs1, lcore1, lcore2, iin, iin_c)
!
!formal solution along ray
   call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem1, iem_c, iemi1, iabs1)
!   if(abs(points_xcoord(indx_point)+10.d0).lt.1.d-1.and.&
!          sqrt(points_xcoord(indx_point)**2+points_ycoord(indx_point)**2).gt.10.d0 ) then 
!      write(*,*) indx_point, lcore1, lcore2, iin, iin_c, iem1/xic1, points_xcoord(indx_point), points_ycoord(indx_point), sqrt(points_xcoord(indx_point)**2+points_ycoord(indx_point)**2), nz_ray
!   endif
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
   call get_iin(indx_xobs2, lcore1, lcore2, iin, iin_c)
!
!formal solution along ray
   call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem2, iem_c, iemi2, iabs2)
!
!------------------interpolation of both xobs-solutions-----------------
!
   iem_surface(indx_point)=interpol_yp(xobs_surface1, xobs_surface2, iem1, iem2, xobs_surface)
   iemi_surface(indx_point)=interpol_yp(xobs_surface1, xobs_surface2, iemi1, iemi2, xobs_surface)
   iabs_surface(indx_point)=interpol_yp(xobs_surface1, xobs_surface2, iabs1, iabs2, xobs_surface)

!
enddo
!
!do i=1, npoints
!   if(abs(points_xcoord(i)+18.477591d0).lt.1.d-3) then
!      write(*,*) i, points_xcoord(i), points_ycoord(i)
!   endif
!enddo
!-18.477591    -0.023894237

!
!stop 'go on in calc_iem_surface'
!
end subroutine calc_iem_surface
!
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine calc_int2d(xobs_2d, alpha_2d, gamma_2d)
!
use prog_type
use fund_const, only: pi
use options_spec, only: input_mod
use dime_spec, only: nxobs_fs
use params_spec, only: vth_fiducial
use mod_spectrum, only: alpha, gamma, xobs, vth_ray, velz_ray, z_ray, temp_ray, &
     profile_ray, nz_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, nz_ray
use mod_interp1d, only: find_index, interpol_yp
use mod_int2d, only: zcoord_2d, xcoord_2d, int_2d, tau_2d, nz_ray_max, iemi_2d, iabs_2d, vn_2d
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
!call find_index(zeta_2d, zeta, nzeta, iim2, iim1, ii, iip1)
indx_zeta1=ii
if(zeta_2d.eq.0.d0) then
   indx_zeta1=iim1
endif
!zeta_2d1=zeta(indx_zeta1)
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
!do indx_p=1, np
!
!   select case(input_mod)
!      case(0)
!         call setup_ray1d(zeta(indx_zeta1), p(indx_p), lcore)
!      case(1)
!         call setup_ray3d(zeta(indx_zeta1), p(indx_p), lcore)
!      case default
!         stop 'error in calc_int2d: input_mod not specified'
!   end select

   if(nz_ray.gt.nz_ray_max) then
      nz_ray_max=nz_ray
   endif
!
!enddo
!
!-----------------------------------------------------------------------
!
!allocate(zcoord_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: zcoord_2d'
!allocate(xcoord_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: xcoord_2d'
!allocate(int_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: int_2d'
!allocate(iemi_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: iemi_2d'
!allocate(iabs_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: iabs_2d'
!allocate(tau_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: tau_2d'
!allocate(vn_2d(np,nz_ray_max), stat=err)
!   if(err.ne.0) stop 'allocation error calc_int2d: vn_2d'
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
!do indx_p=1, np
!
!   select case(input_mod)
!      case(0)
!         call setup_ray1d(zeta(indx_zeta1), p(indx_p), lcore)
!      case(1)
!         call setup_ray3d(zeta(indx_zeta1), p(indx_p), lcore)
!      case default
!         stop 'error in calc_int2d: input_mod not specified'
!   end select
!
!-----------------calculate everything for indx_xobs1-------------------
!
   do indx_zray=1, nz_ray
       call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs1), phinorm)
       profile_ray(indx_zray)=phinorm
   enddo
!
   call get_iin(indx_xobs1, lcore, lcore, iin, iin_c)
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
   call get_iin(indx_xobs2, lcore, lcore, iin, iin_c)
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

!   xcoord_2d(indx_p, :) = p(indx_p)
   zcoord_2d(indx_p, nz_ray_max-nz_ray+1:nz_ray_max) = z_ray(:)
   vn_2d(indx_p, nz_ray_max-nz_ray+1:nz_ray_max) = velz_ray(:)
!
!note: need to set zcoord_2d to some dummy values, in order that
!      triangulation in idl-routine works
   zcoord_2d(indx_p, 1:nz_ray_max-nz_ray) = z_ray(1)-0.0001d0*abs(z_ray(1))

!   if(p(indx_p).le.1.d0) then
!      nz_rest=nz_ray_max-nz_ray
!      del=(z_ray(nz_ray)+z_ray(1)-0.0001d0*z_ray(1))/float(nz_rest-2)
!      zcoord_2d(indx_p, nz_rest) = z_ray(1)-0.0001d0*z_ray(1)
!!      write(*,*) nz_rest, zcoord_2d(indx_p, nz_rest), zcoord_2d(indx_p, nz_rest+1)
!      do j=1, nz_rest-2
!!         write(*,*) j, z_ray(nz_ray)
!         zcoord_2d(indx_p, j) = -z_ray(nz_ray) + (j-1)*del
!      enddo
!!      write(*,'(4es20.8)') zcoord_2d(indx_p, nz_rest-1), zcoord_2d(indx_p, nz_rest), zcoord_2d(indx_p, nz_rest+1), z_ray(1)
!   endif
!
!enddo
!
!
!
end subroutine calc_int2d
!
!
