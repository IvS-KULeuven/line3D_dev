subroutine calc_mod2d_abla
!
! calculates 2d model atmosphere from Dylan Kee's snapshots
! converts the hydro-model to own units, etc.
!
use prog_type
use fund_const, only: pi
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, t_modext2d, rho_modext2d, velr_modext2d, &
                   velth_modext2d, velphi_modext2d, vth_modext2d, eps_cont_modext2d
use params_input, only: teff, vmicro, na, eps_cont
use params_stellar, only: sr
use hdf5
use mod_grid, only: grid_loglog, grid_equi
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2da=380, ntheta_2da=72
integer(i4b) :: i, j
integer(i4b) :: err
!
! ... local arrays
real(dp), dimension(:), allocatable :: theta_dum
real(dp), dimension(:,:), allocatable :: rho_dum, t_dum
!
! ... local characters
!
! ... local functions
real(dp) :: vthermal
!
! ... for hdf5-file
integer(hid_t) :: file_id, density_id, pressure_id, radius_id, lat_id
integer(hsize_t), dimension(1) :: dims_rad
integer(hsize_t), dimension(1) :: dims_lat
integer(hsize_t), dimension(2) :: dims
character(len=7), parameter :: density_name='Density'
character(len=8), parameter :: pressure_name='Pressure'
character(len=6), parameter :: radius_name='Radius'
character(len=8), parameter :: lat_name='Latitude'
!
! ... namelist
character(len=500) :: fname_model
namelist / input_usr / fname_model

!default values
fname_model='tuebingen/test_tau1000.h5'
!fname_model='tuebingen/test_tau10.h5'
!
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_2da
ntheta_modext=2*ntheta_2da

!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
!
allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
!
allocate(theta_dum(ntheta_2da), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(rho_dum(ntheta_2da, nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
allocate(t_dum(ntheta_2da, nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_abla'
!
!---------------read everything from dylans files-----------------------
!
write(*,*) '-------------------------read model from dylan kee-----------------------------'
write(*,*) 'file name: ', trim(fname_model)
write(*,*)
!
dims_rad= (/ nr_modext /)
dims_lat= (/ ntheta_2da /)
dims=(/ ntheta_2da, nr_modext /)
!
call h5open_f(err)
call h5fopen_f(trim(fname_model), h5f_acc_rdonly_f, file_id, err)
!
!open all datasets
call h5dopen_f(file_id, radius_name, radius_id, err)
call h5dopen_f(file_id, lat_name, lat_id, err)
call h5dopen_f(file_id, density_name, density_id, err)
call h5dopen_f(file_id, pressure_name, pressure_id, err)
!
!read all datasets
call h5dread_f(radius_id, h5t_native_double, r_modext2d, dims_rad, err)
call h5dread_f(lat_id, h5t_native_double, theta_dum, dims_lat, err)
call h5dread_f(density_id, h5t_native_double, rho_dum, dims, err)
call h5dread_f(pressure_id, h5t_native_double, t_dum, dims, err)
!
!close everything
call h5dclose_f(radius_id, err)
call h5dclose_f(lat_id, err)
call h5dclose_f(density_id, err)
call h5dclose_f(pressure_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------mirror theta to obtain complete range [0,180]---------------
!
!convert pressure to temperature
t_dum=t_dum*9.95d-9/1.38d0/rho_dum
!
!now mirror grids and store everything in own arrays
do i=1, ntheta_2da
   theta_modext2d(i) = theta_dum(i)
   theta_modext2d(ntheta_2da+i) = pi - theta_dum(ntheta_2da+1-i)
enddo
!
do i=1, nr_modext
   do j=1, ntheta_2da
      rho_modext2d(i,j) = rho_dum(j,i)
      rho_modext2d(i,ntheta_modext+1-j) = rho_dum(j,i)
      t_modext2d(i,j) = t_dum(j,i)
      t_modext2d(i, ntheta_modext+1-j) = t_dum(j,i)
      vth_modext2d(i,j) = vthermal(vmicro*1.d5, teff, na)
      eps_cont_modext2d(i,j) = eps_cont
   enddo
enddo
!
!no velocity information
velr_modext2d=0.d0
velth_modext2d=0.d0
velphi_modext2d=0.d0
!
r_modext2d=r_modext2d/sr
!
!
end subroutine calc_mod2d_abla
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_ablb
!
!----------sets up 2d atmospheric structure from wind ablation model----
!-------------------(by dylan kee, see Kee et al 2016)------------------
!
use prog_type
use fund_const, only: cgs_grav, cgs_mp, cgs_kb, xmsu, pi, sigmae
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, xmloss, vmax, hei, yhe, xlogg, vmin, &
                        teff, tmin, na, vmicro, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_loglog, grid_equi
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2d=1000, ntheta_2d=201
integer(i4b) :: i, j, err
real(dp) :: vkep, vsound, kappae, rhod, rhow, xmloss_cgs, &
            scale_height, b
real(dp) :: rad, theta, rho, dum1, errf, rsin, vel_r, vel_phi, vth
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal
!
! ... namelist
real(dp) :: theta_d, dtheta_abl, beta_accr, tau_d
namelist / input_usr / theta_d, dtheta_abl, beta_accr, tau_d

!default values
theta_d = 12.d0    !in degree
dtheta_abl = 4.d0  !in degree
beta_accr = 1.5d0
tau_d = 1.d3
!
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)

!in rad
theta_d=theta_d*pi/180.d0
dtheta_abl=dtheta_abl*pi/180.d0

write(*,'(a20,es20.8)') 'theta_d', theta_d
write(*,'(a20,es20.8)') 'dtheta_abl', dtheta_abl
write(*,'(a20,es20.8)') 'beta_accr', beta_accr
write(*,'(a20,es20.8)') 'tau_d', tau_d
write(*,*)
!
!---------------allocate grids, and set up spatial coordinates----------
!set up radial grid (equidistant in log-space)
!
nr_modext=nr_2d
ntheta_modext=ntheta_2d
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_ablb'
!
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(rmin, rmin+1.d-4, rlim, nr_modext, r_modext2d)
r_modext2d(nr_modext)=rlim
!
!theta-grid equidistant
call grid_equi(0.d0, pi, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameter from input file--------------
!
!wind density
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*3600.d0)
rhow=xmloss_cgs/4.d0/pi/sr/sr/vmax/1.d5
!
!kappae: mass absorption coefficient
kappae=sigmae*(1.d0+hei*yhe)/(1.d0+4.d0*yhe)/cgs_mp
!
!disc-density
rhod=tau_d/SR/kappae
!
!kepler velocity (in km/s)
vkep=sqrt(10.d0**xlogg * SR)/1.d5
!
!scale height
scale_height=vmin/vkep
!
!velocity law
b=1.d0-vmin/vmax
!
!(constant) thermal velocity
vth = vthermal(vmicro*1.d5, teff, na)
!
!--------------------print all parameters for ablation model------------
!
write(*,*) '---------------------------ablation parameters---------------------------------'
write(*,*)
!
write(*,'(a20, es20.8)') 'v_kepler [km/s]', vkep
write(*,'(a20, es20.8)') 'kappa_e [cm^2/g]', kappae
write(*,'(a20, es20.8)') 'rho_d [g/cm^3]', rhod
write(*,'(a20, es20.8)') 'v_thermal [km/s]', vth*1.d-5
write(*,'(a20, es20.8)') 'v_min [km/s]', vmin
write(*,'(a20, es20.8)') 'rho_w [g/cm^3]', rhow
write(*,'(a20, es20.8)') 'scale h [sqrt(cm)]', scale_height
write(*,'(a20, 2es20.8)') 'rhow/rhod', rhow/rhod
write(*,*)
!
!-------------------calculation of model atmosphere---------------------
!
do i=1, nr_modext
   do j=1, ntheta_modext
      rad=r_modext2d(i)
      theta=theta_modext2d(j)
!avoid theta=0, 180
      if(abs(theta).lt.1.d-8) theta=1.d-8
      if(abs(theta-pi).lt.1.d-8) theta=pi-1.d-8
!
      dum1=1.d0-b/rad
      errf=erfc((theta_d - abs(theta-pi/2.d0))/dtheta_abl)
      rsin=rad*sin(theta)
      vel_r = vmax*1.d5*dum1*errf/2.d0
      vel_phi = vkep*1.d5*(1.d0-errf/2.d0) / sqrt(rsin)
      rho = rhod * rsin**(-1.d0*beta_accr) * &
            exp(scale_height**2 * (1.d0/rad - 1.d0/rsin)) * (1.d0-errf/2.d0) + &
            rhow * vmax*1.d5 / (rad**2 * vel_r) * errf/2.d0
!
      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = tmin
      vth_modext2d(i,j) = vth      
      velr_modext2d(i,j)=vel_r
      velth_modext2d(i,j)=0.d0
      velphi_modext2d(i,j)=vel_phi
      eps_cont_modext2d(i,j) = eps_cont
!
   enddo
enddo
!
!
!r_modext2d=r_modext2d*sr
!
end subroutine calc_mod2d_ablb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_be
!
!----------sets up 2d atmospheric structure for a Be-star disc----------
!
use prog_type
use fund_const
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, hei, yhe, xlogg, &
                        teff, tmin, na, vmicro, vrot, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_log
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, ii, err
integer(i4b) :: idum1, idum2
real(dp) :: rad, theta, sint, cost, vth, velphi, rho, hi, temp
real(dp) :: mstar_cgs, theta_discl, theta_discu, csound, rho0_disc, mmw
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal, vsound, mean_molecular_weight

! ... namelist
real(dp) :: mdisc, rdisc, tdisc, dtheta, slope
namelist / input_usr / mdisc, rdisc, tdisc, dtheta, slope

!default values
mdisc = 4.5d-8    !in msun
rdisc = 13.5d0  !in rstar
tdisc = 30.d3
dtheta = 45.d0 !in degree
slope=1.5d0
!
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)

write(*,'(a20, es20.8)') 'mdisc [m_sun]', mdisc
write(*,'(a20, es20.8)') 'rdisc [rstar]', rdisc
write(*,'(a20, es20.8)') 'tdisc [K]', tdisc
write(*,'(a20, es20.8)') 'dtheta [degree]', dtheta
write(*,'(a20, es20.8)') 'slope', slope
write(*,*)
!
!---------------allocate grids, and set up spatial coordinates----------
!
nr_modext=101
idum1=12
idum2=20
ntheta_modext=2*idum1+2*idum2+1
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_be'
!
!
!radial grid equidistant in log-log-space
call grid_log(rmin, rlim, nr_modext,  r_modext2d)
r_modext2d(nr_modext)=rlim

if(r_modext2d(nr_modext-1).le.rdisc) stop 'error in calc_mod2d_be: r(nr-1) < rdisc'
!
!theta-grid equidistant
dtheta=dtheta*pi/180.d0
theta_discl = pi/two - dtheta
theta_discu = pi/two + dtheta
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, theta_discu-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   theta_modext2d(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   theta_modext2d(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
theta_modext2d(ii)=pi/two
ii=ii+1
do i=1, idum1
   theta_modext2d(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   theta_modext2d(ii)=fdum2_arr(i)
   ii=ii+1
enddo
theta_modext2d(1)=zero
theta_modext2d(ntheta_modext)=pi

!call grid_equi(0.d0, pi, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameter from input file--------------
!
hi = 1.d0   !number free electrons for each hydrogen atom
mmw = mean_molecular_weight(hi,hei,yhe)  !mean molecular weight
mstar_cgs = sr**2*ten**xlogg/cgs_grav
vrot = sqrt(cgs_grav*mstar_cgs/sr)    !breakup velocity (consistent with v_phi of disc model)
!
mdisc=mdisc*xmsu
!
csound = vsound(tdisc,mmw)            
rho0_disc = mdisc*sqrt(cgs_grav*mstar_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr**3.5d0
!
!--------------------print all parameters for Be-disc model-------------
!
write(*,*) '---------------------------Be-disc parameters----------------------------------'
write(*,*)
!
write(*,'(a20, es20.8)') 'vrot [km/s]', vrot/1.d5
write(*,'(a20, es20.8)') 'mstar [msun]', mstar_cgs/xmsu
write(*,'(a20, es20.8)') 'rho0 [g/cm^3]', rho0_disc
write(*,*)
!
!-------------------calculation of model atmosphere---------------------
!
do i=1, nr_modext
   do j=1, ntheta_modext
!
      rad=r_modext2d(i)
      theta=theta_modext2d(j)
  
      sint=sin(theta)
      cost=cos(theta)      
!
      if(rad.le.rdisc.and. &
         theta.gt.theta_discl.and.&
         theta.lt.theta_discu) then
!
         velphi = sqrt(cgs_grav*mstar_cgs/rad/sr/sint)
         temp = tdisc
         vth = vthermal(vmicro*1.d5, temp, na)

         rho = rho0_disc*(rad*sint)**(-slope) * exp(cgs_grav*mstar_cgs/csound**2 * (one/sr/rad-one/sr/rad/sint))
      else
         velphi=zero
         temp=tmin
         vth = vthermal(vmicro*1.d5, temp, na)
         rho=zero
      endif

      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = temp
      vth_modext2d(i,j) = vth      
      velr_modext2d(i,j) = zero
      velth_modext2d(i,j) = zero
      velphi_modext2d(i,j) = velphi
      eps_cont_modext2d(i,j) = eps_cont
!
!      if(j.eq.(ntheta_modext/2+1)) then
!          write(*,'(10es20.8)') theta_modext2d(j), rad, rho, rho0_disc, temp
!      endif

   enddo
!   j=30
!   write(*,*) theta_modext2d(j), rho_modext2d(i,j)
enddo
!
!
r_modext2d=r_modext2d*sr
!
end subroutine calc_mod2d_be
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_adma
!
!-----------sets up 2d atmospheric structure from ADM-model-------------
!---------------------(see Owocki et al 2016)---------------------------
!
use prog_type
use fund_const
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, xmloss, vmax, hei, yhe, xlogg, vmin, &
                        teff, tmin, na, vmicro, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_loglog, grid_equi
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2d=1000, ntheta_2d=201
integer(i4b) :: i, j, err
real(dp) :: mmw, xmloss_cgs, v_esc, v_inf, t_inf, vth, rhoc_star, rhow_star
real(dp) :: rad, mu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
real(dp) :: rho, temp, vel, vel_r, vel_theta
!
! ... namelist
real(dp) :: ralfven, delta, chi_inf, obliquity
namelist / input_usr / ralfven, delta, chi_inf, obliquity
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal
!
!default values
ralfven = 2.7d0   !in rstar
delta = 0.5d0
chi_inf = 1.d-1
obliquity = 0.d0  !in degree
!
write(*,*) '---------------------calculate 2D ADM model------------------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
write(*,'(a20,es20.8)') 'ralfven', ralfven
write(*,'(a20,es20.8)') 'delta', delta
write(*,'(a20,es20.8)') 'chi_inf', chi_inf
write(*,*)

!---------------allocate grids, and set up spatial coordinates----------
!set up radial grid (equidistant in log-space)
!
nr_modext=nr_2d
ntheta_modext=ntheta_2d
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_adma'
!
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(rmin, rmin+1.d-4, rlim, nr_modext, r_modext2d)
!
!theta-grid equidistant
call grid_equi(0.d0, pi, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameters from input file-------------
!
v_esc=sqrt(2.d0*sr*10.d0**xlogg)
v_inf=vmax*1.d5
xmloss_cgs=xmloss*xmsu/(3600.d0*24.d0*365.25d0)
rhoc_star=xmloss_cgs/(4.d0*PI*SR*SR*V_ESC)
rhow_star=xmloss_cgs/(4.d0*PI*SR*SR*V_INF)
!
!mean  molecular weight (for completely ionized hydrogen and helium abundance yhe)
mmw=2.d0*(1.d0-yhe)+ (1.d0+hei)*yhe/4.d0
mmw=1.d0/mmw
!
t_inf=3.d0*mmw*cgs_mp*v_inf*v_inf / cgs_kb / 16.d0
!
!(constant) thermal velocity
vth = vthermal(vmicro*1.d5, teff, na)
!
!--------------------print all parameters for ablation model------------
!
write(*,*) '-----------------------------ADM parameters------------------------------------'
write(*,*)
!
write(*,'(a20, es20.8)') 'v_escape [km/s]', v_esc/1.d5
write(*,'(a20, es20.8)') 'rhoc_star[g/cm^3]', rhoc_star
write(*,'(a20, es20.8)') 'rhow_star[g/cm^3]', rhow_star
write(*,'(a20, es20.8)') 'v_thermal [km/s]', vth*1.d-5
write(*,'(a20, es20.8)') 't_inf [K]', t_inf
write(*,*)
!
!-------------------calculation of model atmosphere---------------------
!
do i=1, nr_modext
      rad=r_modext2d(i)
   do j=1, ntheta_modext
!
      if(rad.gt.1.d0) then
!can calculate all mu_star, and mu_shock only if radius larger than 1 for numerical reasons
!
!set position coordinates (rad, mu), with mu=cos(theta)
         mu=cos(theta_modext2d(j))
!for current position, calculate intersection point of the attached closed loop
!   with the photosphere, as well as the apex radius of this field line
         mu_star=sqrt(1.d0-(1.d0-mu*mu)/rad)
!set maximum of mu_star to be equal to 1.-1.d-9 (otherwise: r_apex is not defined on z-axis)
         mu_star=min(mu_star, 0.999999999d0)               
         r_apex=1.d0/(1.d0-mu_star*mu_star)
!for current position, calculate the shock-radius r_shock, mu_shock (with initial guess for 
!   interval surrounding mu_shock: [mu_lower, mu_upper]
         mu_lower=0.d0
         mu_upper=mu_star
         call get_mu_shock(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
         r_shock=r_apex*(1.d0-mu_shock*mu_shock)
!
!calculate wind upflow component everywhere
         call component_w(rad, mu, rho, temp, vel, vel_r, vel_theta, &
                          teff, v_esc, v_inf, rhoc_star)
!
         if(r_apex.lt.ralfven) then
!in that regime: calculate cool downflow component and overwrite variables
            call component_c(rad, mu, r_apex, rho, temp, vel, vel_r, vel_theta, &
                             teff, v_esc, rhoc_star, delta)
!
            if(rad.gt.r_shock) then
!in that regime: calculate post-shock component and overwrite variables
               CALL COMPONENT_S(rad, mu, mu_star, mu_shock, rho, temp, vel, vel_r, vel_theta, &
                                teff, t_inf, v_esc, v_inf, rhoc_star, delta)
            endif
         endif
!
      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = temp
      vth_modext2d(i,j) = vth      
      velr_modext2d(i,j)=vel_r
      velth_modext2d(i,j)=vel_theta
      velphi_modext2d(i,j)=0.d0
      eps_cont_modext2d(i,j) = eps_cont
!
      endif

   enddo
enddo
!
!replace all values at r=1 with neighbouring ones
do i=1, ntheta_modext
   rho_modext2d(1,i)=rho_modext2d(2,i)
   t_modext2d(1,i)=t_modext2d(2,i)
   vth_modext2d(1,i)=vth_modext2d(2,i)
   velr_modext2d(1,i)=velr_modext2d(2,i)
   velth_modext2d(1,i)=velth_modext2d(2,i)
   velphi_modext2d(1,i)=velphi_modext2d(2,i)
   eps_cont_modext2d(1,i)=eps_cont_modext2d(2,i)
enddo
!
!
r_modext2d=r_modext2d*sr
!
!j=ntheta_modext/2+1
!do i=1, nr_modext
!   write(*,'(7es20.8)') r_modext2d(i), rho_modext2d(i,j), t_modext2d(i,j), vth_modext2d(i,j), velr_modext2d(i,j), velth_modext2d(i,j), velphi_modext2d(i,j)
!enddo
!
call calc_mdot_adm(ralfven, chi_inf, delta)
!
end subroutine calc_mod2d_adma
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_admb
!
!-----------sets up 2d atmospheric structure from ADM-model-------------
!---------------------(see Owocki et al 2016)---------------------------
!-----------spherical beta-velocity law is used in wind-region----------
!
use prog_type
use fund_const
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, xmloss, vmax, hei, yhe, xlogg, vmin, &
                        teff, tmin, na, vmicro, beta, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_loglog, grid_equi
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2d=1000, ntheta_2d=201
integer(i4b) :: i, j, err
real(dp) :: mmw, xmloss_cgs, v_esc, v_inf, t_inf, vth, rhoc_star, rhow_star
real(dp) :: rad, mu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
real(dp) :: rho, temp, vel, vel_r, vel_theta, b
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal, bvel
!
! ... namelist
real(dp) :: ralfven, delta, chi_inf, obliquity
namelist / input_usr / ralfven, delta, chi_inf, obliquity

!default values
ralfven = 2.7d0   !in rstar
delta = 0.5d0
chi_inf = 1.d-1
obliquity = 0.d0  !in degree
!
write(*,*) '---------------------calculate 2D ADM model (with spherical beta law)----------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
write(*,'(a20,es20.8)') 'ralfven', ralfven
write(*,'(a20,es20.8)') 'delta', delta
write(*,'(a20,es20.8)') 'chi_inf', chi_inf
write(*,*)
!
!---------------allocate grids, and set up spatial coordinates----------
!set up radial grid (equidistant in log-space)
!
nr_modext=nr_2d
ntheta_modext=ntheta_2d
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_admb'
!
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(rmin, rmin+1.d-4, rlim, nr_modext, r_modext2d)
!
!theta-grid equidistant
call grid_equi(0.d0, pi, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameters from input file-------------
!
v_esc=sqrt(2.d0*sr*10.d0**xlogg)
v_inf=vmax*1.d5
xmloss_cgs=xmloss*xmsu/(3600.d0*24.d0*365.25d0)
rhoc_star=xmloss_cgs/(4.d0*PI*SR*SR*V_ESC)
rhow_star=xmloss_cgs/(4.d0*PI*SR*SR*V_INF)
!
!mean  molecular weight (for completely ionized hydrogen and helium abundance yhe)
mmw=2.d0*(1.d0-yhe)+ (1.d0+hei)*yhe/4.d0
mmw=1.d0/mmw
!
t_inf=3.d0*mmw*cgs_mp*v_inf*v_inf / cgs_kb / 16.d0
!
!(constant) thermal velocity
vth = vthermal(vmicro*1.d5, teff, na)
!
!b-factor for beta-velocity law
b=1.d0-(vmin/vmax)**(1.d0/beta)
!
!--------------------print all parameters for ablation model------------
!
write(*,*) '------------------------------ADM parameters-----------------------------------'
write(*,*)
!
write(*,'(a20, es20.8)') 'v_escape [km/s]', v_esc/1.d5
write(*,'(a20, es20.8)') 'rhoc_star[g/cm^3]', rhoc_star
write(*,'(a20, es20.8)') 'rhow_star[g/cm^3]', rhow_star
write(*,'(a20, es20.8)') 'v_thermal [km/s]', vth*1.d-5
write(*,'(a20, es20.8)') 't_inf [K]', t_inf
write(*,*)
!
!-------------------calculation of model atmosphere---------------------
!
do i=1, nr_modext
      rad=r_modext2d(i)
   do j=1, ntheta_modext
!
      if(rad.gt.1.d0) then
!can calculate all mu_star, and mu_shock only if radius larger than 1 for numerical reasons
!
!set position coordinates (rad, mu), with mu=cos(theta)
         mu=cos(theta_modext2d(j))
!for current position, calculate intersection point of the attached closed loop
!   with the photosphere, as well as the apex radius of this field line
         mu_star=sqrt(1.d0-(1.d0-mu*mu)/rad)
!set maximum of mu_star to be equal to 1.-1.d-9 (otherwise: r_apex is not defined on z-axis)
         mu_star=min(mu_star, 0.999999999d0)               
         r_apex=1.d0/(1.d0-mu_star*mu_star)
!for current position, calculate the shock-radius r_shock, mu_shock (with initial guess for 
!   interval surrounding mu_shock: [mu_lower, mu_upper]
         mu_lower=0.d0
         mu_upper=mu_star
         call get_mu_shock(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
         r_shock=r_apex*(1.d0-mu_shock*mu_shock)
!
!calculate wind upflow component everywhere
         call component_w(rad, mu, rho, temp, vel, vel_r, vel_theta, &
                          teff, v_esc, v_inf, rhoc_star)
         vel_theta=0.d0
         vel_r=bvel(rad, v_inf, b, beta)
!
         if(r_apex.lt.ralfven) then
!in that regime: calculate cool downflow component and overwrite variables
            call component_c(rad, mu, r_apex, rho, temp, vel, vel_r, vel_theta, &
                             teff, v_esc, rhoc_star, delta)
!
            if(rad.gt.r_shock) then
!in that regime: calculate post-shock component and overwrite variables
               CALL COMPONENT_S(rad, mu, mu_star, mu_shock, rho, temp, vel, vel_r, vel_theta, &
                                teff, t_inf, v_esc, v_inf, rhoc_star, delta)
            endif
         endif
!
      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = temp
      vth_modext2d(i,j) = vth      
      velr_modext2d(i,j)=vel_r
      velth_modext2d(i,j)=vel_theta
      velphi_modext2d(i,j)=0.d0
      eps_cont_modext2d(i,j) = eps_cont
!
      endif

   enddo
enddo
!
!replace all values at r=1 with neighbouring ones
do i=1, ntheta_modext
   rho_modext2d(1,i)=rho_modext2d(2,i)
   t_modext2d(1,i)=t_modext2d(2,i)
   vth_modext2d(1,i)=vth_modext2d(2,i)
   velr_modext2d(1,i)=velr_modext2d(2,i)
   velth_modext2d(1,i)=velth_modext2d(2,i)
   velphi_modext2d(1,i)=velphi_modext2d(2,i)
   eps_cont_modext2d(1,i)=eps_cont_modext2d(2,i)
enddo
!
r_modext2d=r_modext2d*sr
!
!j=ntheta_modext/2+1
!do i=1, nr_modext
!   write(*,'(7es20.8)') r_modext2d(i), rho_modext2d(i,j), t_modext2d(i,j), vth_modext2d(i,j), velr_modext2d(i,j), velth_modext2d(i,j), velphi_modext2d(i,j)
!enddo
!
call calc_mdot_adm(ralfven, chi_inf, delta)
!
end subroutine calc_mod2d_admb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mdot_adm(ralfven, chi_inf, delta)
!
!estimating mdot by:
!   mdot/mdto_b0 = 2 * integral from mu_c to 1 of feeding rate
! mu_c is footpoint of closed loop for apex-radius=r_alfven
!
use prog_type
use fund_const, ONLY: pi, xmsu, cgs_mp, cgs_kb
use params_input, only: xmloss, teff, yhe, hei, vmax, xlogg
use params_stellar, only: sr
use mod_integ1d, only: precalc_weight_trapez
!
implicit none
!
! ... arguments
real(dp), intent(in) :: ralfven, chi_inf, delta
!
! ... local scalars
integer(i4b), parameter :: nmu=101
integer(i4b) :: i, j, nr_test
real(dp) :: mdot_true1, mdot_true2, sqrt3, mu_crit, cc
real(dp) :: r1, r2, ralfven_test
real(dp) :: sum1, sum2
real(dp) :: rho, vel, velr, velth, temp
real(dp) :: rad, mu, mu_star, r_apex, mu_lower, mu_upper, mu_shock, r_shock, del
real(dp) :: flux_w, flux_c, flux_s, flux_tot, fw, fc, fs, ft
real(dp) :: v_esc, v_inf, xmloss_cgs, rhoc_star, rhow_star, mmw, t_inf
!
! ... local arrays
real(dp), dimension(nmu) :: mu_grid1, mu_grid2, mu_weight1, mu_weight2
real(dp), dimension(:), allocatable :: r_grid, r_weight
!
! ... local characters
!
! ... local functions
!
!
write(*,*) '-----------------calculating true mass-loss rate for ADM model-----------------'
write(*,*)
!
!-------------calculate physical parameters from input file-------------
!
v_esc=sqrt(2.d0*sr*10.d0**xlogg)
v_inf=vmax*1.d5
xmloss_cgs=xmloss*xmsu/(3600.d0*24.d0*365.25d0)
rhoc_star=xmloss_cgs/(4.d0*PI*SR*SR*V_ESC)
rhow_star=xmloss_cgs/(4.d0*PI*SR*SR*V_INF)
!
!mean  molecular weight (for completely ionized hydrogen and helium abundance yhe)
mmw=2.d0*(1.d0-yhe)+ (1.d0+hei)*yhe/4.d0
mmw=1.d0/mmw
!
t_inf=3.d0*mmw*cgs_mp*v_inf*v_inf / cgs_kb / 16.d0
!
!--------------------------from theoretical formulae--------------------
!
sqrt3=sqrt(3.d0)
!
!calculating mdot with feeding rate defined along magnetic loops ('correct ansatz')
mu_crit=sqrt(1.d0-1.d0/ralfven)
cc=1.d0-mu_crit+atan(sqrt3*mu_crit)/sqrt3 - atan(sqrt3)/sqrt3
mdot_true1=4.d0*cc/3.d0 * xmloss
!
mdot_true2=(1.d0-mu_crit) * xmloss
!
!---------------check if densities and velocities are correct-----------
!
r1=1.d0+1.d-8
r2=ralfven

mu_crit=sqrt(1.d0-1.d0/ralfven)
!
!calculate mu-grid from mu_crit to 1 => integration at stellar surface
!calculate mu_grid from 0 to 1 => integration at alfven radius
do i=1, nmu
   mu_grid1(i)=mu_crit+(i-1)*(1d0-mu_crit)/float(nmu-1)
   mu_grid2(i)=(i-1)/float(nmu-1)
enddo
!
call precalc_weight_trapez(mu_grid1, nmu, mu_weight1)
call precalc_weight_trapez(mu_grid2, nmu, mu_weight2)
!
sum1=0.d0
sum2=0.d0
!
do i=1, nmu
!
   call component_w(r1, mu_grid1(i), rho, temp, vel, velr, velth, &
                    teff, v_esc, v_inf, rhoc_star)
   sum1 = sum1 + rho*velr*r1*r1*mu_weight1(i)
!
   call component_w(r2, mu_grid2(i), rho, temp, vel, velr, velth, &
                    teff, v_esc, v_inf, rhoc_star)
   sum2 = sum2 + rho*velr*r2*r2*mu_weight2(i)
!
enddo
!
!in cgs
sum1=sum1*2.d0*pi*SR*SR*2.d0
sum2=sum2*2.d0*pi*SR*SR*2.d0
!
!in mdot/yr
sum1=sum1*365.25d0*24.d0*3600.d0/xmsu
sum2=sum2*365.25d0*24.d0*3600.d0/xmsu
!
!-----------------------------print out everything----------------------
!
write(*,'(a55, 2es20.8)') 'mdot (using feeding rate along b-field line)', mdot_true1, mdot_true1/xmloss
write(*,'(a55, 2es20.8)') 'mdot (using scaling relation (e.g. Petit et al 2017)', mdot_true2, mdot_true2/xmloss
write(*,'(a55, 2es20.8)') 'mdot (integrating over mu numerically at r=1)', sum1, sum1/xmloss
write(*,'(a55, 2es20.8)') 'mdot (integrating over mu numerically AT r=r_A)', sum2, sum2/xmloss
write(*,*)
!
if(abs(sum1-sum2)/sum1.gt.1.d-3) stop 'error1 in calc_mdot_adm: error to be found in rho_w'
if(abs(sum1-mdot_true1)/sum1.gt.1.d-3) stop 'error2 in calc_mdot_adm: error to be found in rho_w'
!
!
!
end subroutine calc_mdot_adm
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_bc
!
!-----------sets up 2d atmospheric structure from WCD-model-------------
!--------------------(Bjorkman&Cassinelli 1993)-------------------------
!
use prog_type
use fund_const
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, xmloss, vmax, hei, yhe, xlogg, vmin, &
                              teff, tmin, na, vmicro, beta, vrot, eps_cont, mstar
use params_stellar, only: sr
use mod_grid, only: grid_loglog, grid_equi
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2d=1000, ntheta_2d=201
integer(i4b), parameter :: itmax=1000
integer(i4b) :: i, j, k, indx_osc, err
real(dp) :: vinf_pole, vinf, vcrit, vesc, vth, gedd, lstar, mstar_eff
real(dp) :: theta, theta0, theta0_ii, theta0_iim1, phi0, phi0_new, &
            st, st0, ct0, sp0, db_dst0, dc0_db, dc0_dst0, dphi0_dst0, dmu_dmu0
real(dp) :: rad, rho, rho1d, temp, vel, vel_r, vel_phi, vel_theta, b, bbeta
real(dp) :: xmloss_cgs, xmloss_theta0, xmloss1d
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal, bvel
!
! ... namelist
real(dp) :: zeta, gamma, xi, alpha_cak, delta_cak
namelist / input_usr / mstar, zeta, gamma, xi, alpha_cak, delta_cak

!default values
mstar = 52.5 !in msun
zeta = 4.18d0
gamma = 0.35d0
xi = -0.43d0
alpha_cak = 0.66d0
delta_cak = 0.07d0
!
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
write(*,'(a20,es20.8)') 'mstar', mstar
write(*,'(a20,es20.8)') 'zeta', zeta
write(*,'(a20,es20.8)') 'gamma', gamma
write(*,'(a20,es20.8)') 'xi', xi
write(*,'(a20,es20.8)') 'alpha_cak', alpha_cak
write(*,'(a20,es20.8)') 'delta_cak', delta_cak
write(*,*)
!
!
!---------------allocate grids, and set up spatial coordinates----------
!set up radial grid (equidistant in log-space)
!
nr_modext=nr_2d
ntheta_modext=ntheta_2d
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_bc'
!
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(rmin, rmin+1.d-4, rlim, nr_modext, r_modext2d)
!
!theta-grid equidistant (theta=0 and theta=pi not allowed for numerical reasons later on)
call grid_equi(1.d-5, pi-1.d-5, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameters from input file-------------
!
write(*,*) '------------------------------BC parameters------------------------------------'
write(*,*)
!
!stellar luminosity
lstar = 4.d0*pi*sr**2*teff**4*cgs_sb
write(*,'(a20,es20.8)') 'L_star [L_sun]', lstar/xlsu
!
!eddington factor
gedd = sigmae * lstar/4.d0/pi/cgs_clight/cgs_grav/mstar/xmsu/cgs_mp * (1.d0+hei*yhe)/(1.d0+4.d0*yhe)
write(*,'(a20,es20.8)') 'eddington factor', gedd
!
!effective mass
mstar_eff = mstar*(1.d0-gedd)
write(*,'(a20,2es20.8)') 'effective mass', mstar_eff, mstar_eff/mstar
!
!escape velocity
vesc=sqrt(2.d0*cgs_grav*mstar_eff*xmsu/sr)
write(*,'(a20,es20.8)') 'vesc [km/s]', vesc/1.d5
!
!critical (break up) velocity
vcrit = vesc/sqrt(2.d0)
write(*,'(a20,es20.8)') 'vcrit [km/s]', vcrit/1.d5
!
!terminal velocity at the pole
vinf_pole=vmax*1.d5
vrot = vrot*1.d5
vmin = vmin*1.d5
write(*,'(a20,es20.8)') 'vinf [km/s] at pole', vinf_pole/1.d5
!
write(*,'(a20,es20.8)') 'mdot [msun/yr]', xmloss
xmloss_cgs=xmloss*xmsu/(3600.d0*24.d0*365.25d0)
!
write(*,'(a20,es20.8)') 'zeta (input)', zeta
zeta=vinf_pole/vesc
write(*,'(a20,es20.8)') 'zeta (vinf_p/vesc)', zeta
write(*,'(a20,es20.8)') 'gamma', gamma
write(*,'(a20,es20.8)') 'xi', xi
write(*,'(a20,es20.8)') 'omega', vrot/vcrit
write(*,*)
!
!ab-factor for beta-velocity law
!b=1.d0-(vmin/vmax)**(1.d0/beta)
!
!-------------------calculation of model atmosphere---------------------
!
bbeta=1.d0-beta
!
do i=1, nr_modext
!
   rad=r_modext2d(i)
!
   do j=1, ntheta_modext
!
      theta=theta_modext2d(j)
!
!---------------------------first step----------------------------------
!
!first guess of phi0 (=phi' in BC paper)
      phi0 = min(theta,1.d-5)
!calculate theta0
      theta0 = acos(cos(theta)/cos(phi0))
!calculate terminal velocity and b-factor for given theta0
      call vinf_bc(vmin, vrot, vesc, vcrit, beta, zeta, gamma, theta0, vinf, b)
      theta0_iim1 = theta0
!
!start iteration
      do k=1, itmax
         if(beta.eq.1.d0) then
            phi0_new = (log(1.d0-b/rad)-log(1.d0-b))/b * sin(theta0_iim1)*vrot/vinf
         else
            phi0_new = ((1.d0-b/rad)**bbeta - (1.d0-b)**bbeta)/b/bbeta * sin(theta0_iim1)*vrot/vinf
         endif
         theta0_ii = acos(cos(theta)/cos(phi0_new))
         call vinf_bc(vmin, vrot, vesc, vcrit, beta, zeta, gamma, theta0_ii, vinf, b)
         if(abs(theta0_ii-theta0_iim1).lt.1.d-10) exit
         theta0_iim1 = theta0_ii
      enddo
!
      theta0 = theta0_ii
      phi0 = phi0_new
!
      if(k.eq.itmax+1) then
         write(*,*) 'no converged solution found in calc_model2d_bc'
         write(*,'(a20,es20.8)') 'at rad', rad
         write(*,'(a20,2es20.8)') 'at theta', theta*180.d0/pi, theta
         write(*,'(a20,es20.8)') 'vinf', vinf
         write(*,'(a20,es20.8)') 'b', b
         write(*,'(a20,es20.8)') 'theta0', theta0
         write(*,'(a20,es20.8)') 'phi0', phi0
         stop
      endif
!
      st = sin(theta)
      sp0 = sin(phi0)
      st0 = sin(theta0)
      ct0 = cos(theta0)
!calculate velocity components
      vel_r = vinf*(1.d0-b/rad)**beta
      vel_theta = vrot*st0*ct0*sp0/st/rad
      vel_phi = vrot*st0**2/st/rad
!calculate density
      db_dst0 = -gamma*vrot/beta/(1.d0-st0*vrot/vcrit)/vcrit * (vmin/vinf)**(1.d0/beta)
      if(beta.eq.1.d0) then
         dc0_db = (b/(1.d0-b) - b/(rad-b) - log(1.d0-b/rad) + log(1.d0-b)) / b**2
      else
         dc0_db = (b/rad*(beta-1.d0)*(1.d0-b/rad)**(-beta) - &
                  (1.d0-b/rad)**(1.d0-beta) + &
                  b*(1.d0-beta)*(1.d0-b)**(-beta) + &
                  (1.d0-b)**(1.d0-beta))/(1.d0-beta)/b**2
!         stop 'to include'
      endif
      dc0_dst0 = dc0_db*db_dst0
      dphi0_dst0 = (1.d0+gamma*vrot*st0/(vcrit-vrot*st0))*phi0/st0 + vrot*st0/vinf * dc0_dst0
      dmu_dmu0 = cos(phi0) + ct0**2*sp0/st0 * dphi0_dst0
!      xmloss_theta0 = xmloss_cgs*(1.d0-st0*vrot/vcrit)**xi   !accroding to petrenz&puls1996
      xmloss_theta0 = xmloss_cgs*(1.d0-(st0*vrot/vcrit)**2)**(1.d0-1.d0/(alpha_cak-delta_cak))   !according to puls review 2008

      rho = xmloss_theta0/4.d0/pi/(rad*sr)**2/dmu_dmu0/vel_r
      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = 0.75d0*teff
      vth_modext2d(i,j) = vthermal(vmicro*1.d5, 0.75d0*teff, na)
      velr_modext2d(i,j) = vel_r
      velth_modext2d(i,j) = vel_theta !0.d0 !vel_theta
      velphi_modext2d(i,j) = vel_phi
      eps_cont_modext2d(i,j) = eps_cont
!
   enddo
enddo
!
!-----------------------------------------------------------------------
!
!calculate 1d mass loss rates at all radii (integrated mass flux)
do i=1, nr_modext
   xmloss1d = 0.d0
   do j=2, ntheta_modext
      xmloss1d = xmloss1d + 2.d0*pi*(r_modext2d(i)*sr)**2*(rho_modext2d(i,j)*velr_modext2d(i,j)*sin(theta_modext2d(j)) + rho_modext2d(i,j-1)*velr_modext2d(i,j-1)*sin(theta_modext2d(j)))*(theta_modext2d(j)-theta_modext2d(j-1))/2.d0
   enddo
!   write(*,*) r_modext2d(i), xmloss1d*(r_modext2d(i)*sr)**2 *(3600.d0*24.d0*365.25d0)/xmsu
enddo
!
write(*,*)
write(*,'(a20,es20.8)') 'mdot(1d) [msun/yr]', xmloss1d*(3600.d0*24.d0*365.25d0)/xmsu
!
!calculate 1d densities (only for compariston in plot routine, scaling...)
do i=1, nr_modext
   rho1d = xmloss1d/4.d0/pi/velr_modext2d(i,1)/(r_modext2d(i)*sr)**2
   do j=1, ntheta_modext
!      rho_modext2d(i,j) = rho_modext2d(i,j)/rho1d
   enddo
enddo
!
!
end subroutine calc_mod2d_bc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_rot
!
!-----------sets up 2d atmospheric structure for rotation---------------
!------radial beta velocity law
!------azimuthal velocities from vel_phi=vrot*sin(phi)/r
!------1d density stratification
!(e.g. petrenz 1995, paragraph 3.3)
!
use prog_type
use fund_const
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, velr_modext2d, velth_modext2d, &
                   velphi_modext2d, t_modext2d, vth_modext2d, rho_modext2d, eps_cont_modext2d
use params_input, only: rmin, rlim, xmloss, vmax, hei, yhe, xlogg, vmin, &
                        teff, tmin, na, vmicro, beta, mstar, vrot, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_loglog, grid_equi
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_2d=1000, ntheta_2d=201
integer(i4b) :: i, j, k, indx_osc, err
real(dp) :: vinf, vth
real(dp) :: rad, rho, temp, vel_r, vel_phi, vel_theta, b
real(dp) :: xmloss_cgs
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal, bvel
!
!---------------allocate grids, and set up spatial coordinates----------
!set up radial grid (equidistant in log-space)
!
nr_modext=nr_2d
ntheta_modext=ntheta_2d
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'

allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_rot'
!
!
!radial grid equidistant in log-log-space
!call grid_log(rmin, rlim, nr_modext,  r_modext2d)
call grid_loglog(rmin, rmin+1.d-4, rlim, nr_modext, r_modext2d)
!
!theta-grid equidistant (theta=0 and theta=pi not allowed for numerical reasons later on)
call grid_equi(1.d-5, pi-1.d-5, ntheta_modext, theta_modext2d)
!
!-------------calculate physical parameters from input file-------------
!
write(*,*) '--------------------------model parameters-------------------------------------'
write(*,*)
!
vrot = vrot*1.d5
vmin = vmin*1.d5
vinf = vmax*1.d5
write(*,'(a20,es20.8)') 'vmin [km/s] ', vmin/1.d5
write(*,'(a20,es20.8)') 'vinf [km/s] ', vinf/1.d5
write(*,'(a20,es20.8)') 'vrot [km/s] ', vrot/1.d5
write(*,'(a20,es20.8)') 'beta        ', beta
!
write(*,'(a20,es20.8)') 'mdot [msun/yr]', xmloss
xmloss_cgs=xmloss*xmsu/(3600.d0*24.d0*365.25d0)
!
!ab-factor for beta-velocity law
b=1.d0-(vmin/vinf)**(1.d0/beta)
!
!-------------------calculation of model atmosphere---------------------
!
!
do i=1, nr_modext
!
   rad=r_modext2d(i)

   vel_r = vinf*(1.d0-b/rad)**beta
   vel_theta = 0.d0
   rho = xmloss_cgs/4.d0/pi/(rad*sr)**2/vel_r

   do j=1, ntheta_modext
!
      vel_phi = vrot*sin(theta_modext2d(j))/rad

      rho_modext2d(i,j) = rho
      t_modext2d(i,j) = tmin
      vth_modext2d(i,j) = vthermal(vmicro*1.d5, tmin, na)
      velr_modext2d(i,j) = vel_r
      velth_modext2d(i,j) = vel_theta
      velphi_modext2d(i,j) = vel_phi
      eps_cont_modext2d(i,j) = eps_cont
!
   enddo
enddo
!
!
!
end subroutine calc_mod2d_rot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_vh1
!
! calculates 2d model atmosphere from VH-1 snapshots
!   (aloready preprocessed by idl routine)
! basically, only resets the thermal velocity...
!
use prog_type
use fund_const, only: pi
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, t_modext2d, rho_modext2d, velr_modext2d, &
                   velth_modext2d, velphi_modext2d, vth_modext2d, eps_cont_modext2d
use params_input, only: teff, vmicro, na, eps_cont
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: err
!
! ... local arrays
!
! ... local characters
character(len=24), parameter :: model2d_vh1='models/rotvh1/model2d.h5'
!model2d_vh1: input-file for vh1 input file (already preprocessed in idl)
!
! ... local functions
real(dp) :: vthermal
!
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, dset_id, attr_id
integer(hsize_t), dimension(1) :: dims_scalars
integer(hsize_t), dimension(1) :: dims_rad
integer(hsize_t), dimension(1) :: dims_lat
integer(hsize_t), dimension(2) :: dims
!
!---------------read everything from dylans files-----------------------
!
write(*,*) '-------------------------read model from dylan kee-----------------------------'
write(*,*) 'file name: ', model2d_vh1
write(*,*)
!
!
call h5open_f(err)
call h5fopen_f(model2d_vh1, h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!-----------------------allocate arrays---------------------------------
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
!
allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1 '
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1 '
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1 '
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_vh1'
!
dims_rad= (/ nr_modext /)
dims_lat= (/ ntheta_modext /)
dims=(/ nr_modext, ntheta_modext /)
!
!-------------------------------------------------------------
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext2d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext2d, dims_rad, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims, err)
   call h5dclose_f(dset_id, err)

call h5gclose_f(group_id, err)

!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!--------------------------calculate new thermal velocity---------------
!
do i=1, nr_modext
   do j=1, ntheta_modext
      vth_modext2d(i,j) = vthermal(vmicro*1.d5, teff, na)
      eps_cont_modext2d(i,j) = eps_cont
   enddo
enddo
!
!write(*,*) ntheta_modext, theta_modext2d(ntheta_modext/2)*180/pi
!write(*,*) velr_modext2d(:,1)
!write(*,*)
!write(*,*) velr_modext2d(:,ntheta_modext/2)
!
!
end subroutine calc_mod2d_vh1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod2d_florian
!
! reads in Florian's MHD simulations
!
use prog_type
use fund_const, only: pi
use dime_modext, only: nr_modext, ntheta_modext
use model2d, only: r_modext2d, theta_modext2d, t_modext2d, rho_modext2d, velr_modext2d, &
                   velth_modext2d, velphi_modext2d, vth_modext2d, eps_cont_modext2d
use params_input, only: teff, vmicro, na, eps_cont
use mod_directories, only: indat_file
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
integer(i4b) :: err
!
! ... local arrays
!
! ... local characters
character(len=500)  :: fname_model
!
! ... local logicals
logical :: lcheck
!
! ... local functions
real(dp) :: vthermal
!
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, dset_id, attr_id
integer(hsize_t), dimension(1) :: dims_scalars
integer(hsize_t), dimension(1) :: dims_rad
integer(hsize_t), dimension(1) :: dims_lat
integer(hsize_t), dimension(2) :: dims
!
! ... namelists
integer(i4b) :: opt_bvel
integer(i4b), dimension(:), allocatable :: nd
namelist / input_usr / fname_model
!
!
write(*,*) '-----------------------read 2d model from florian------------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
write(*,*) 'file name: ', trim(fname_model)
write(*,*)

inquire(file=trim(fname_model), exist=lcheck)
if(.not.lcheck) then
   write(*,*) 'error in calc_mod2d_florian: file "', trim(fname_model), '" does not exist'
   stop
endif
!
!
call h5open_f(err)
call h5fopen_f(trim(fname_model), h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!-----------------------allocate arrays---------------------------------
!
allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
!
allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod2d_florian'
!
dims_rad= (/ nr_modext /)
dims_lat= (/ ntheta_modext /)
dims=(/ nr_modext, ntheta_modext /)
!
!-------------------------------------------------------------
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext2d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext2d, dims_rad, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims, err)
   call h5dclose_f(dset_id, err)

call h5gclose_f(group_id, err)

!
call h5fclose_f(file_id, err)
call h5close_f(err)

!write(*,*) r_modext2d(11)
!
!--------------------------calculate new thermal velocity---------------
!
do i=1, nr_modext
   do j=1, ntheta_modext
      vth_modext2d(i,j) = vthermal(vmicro*1.d5, teff, na)
      eps_cont_modext2d(i,j) = eps_cont
   enddo
enddo
!
!do i=1, nr_modext
!   write(*,*) sqrt(velr_modext2d(i,1)**2+velth_modext2d(i,1)**2)/1.d5
!enddo
!stop
!
end subroutine calc_mod2d_florian
