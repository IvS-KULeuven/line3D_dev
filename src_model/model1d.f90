subroutine calc_mod1d
!
!-------calculated model atmosphere from beta velocity law--------------
!
use prog_type
use fund_const, only: cgs_mp, pi, xmsu, sigmae
use dime_modext, only: nr_modext
use model1d, only: r_modext1d, t_modext1d, rho_modext1d, velr_modext1d, vth_modext1d, eps_cont_modext1d
use params_input, only: vmax, vmin, beta, xmloss, yhe, hei, rlim, rmin, teff, vmicro, na, eps_cont, kcont
use params_stellar, only: sr
use mod_opacities, only: opac_thomson
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_1d=1000
integer(i4b) :: i, j
integer(i4b) :: err
real(dp) :: bconst, cconst1, cconst2, xmloss_cgs
real(dp) :: taumin, taumax, deltau, tau
real(dp) :: sigem, delr, velr
!
! ... local characters
character(len=9), parameter :: model1d_jo='models/jo'
!
! ... local arrays
real(dp), dimension(:), allocatable :: opac_modext1d
! ... local functions
real(dp) :: bvel, vthermal
!
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_1d
!
allocate(r_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(velr_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(rho_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(opac_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(t_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(vth_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(eps_cont_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
!
!------------------calculate/define required constants------------------
!
!b-factor for beta-velocity-law
bconst=1.d0-(vmin/vmax)**(1.D0/beta)
!
cconst1=(1.d0+4.d0*yhe)*cgs_mp
cconst2=(1.d0+hei*yhe)/cconst1
!
!mass-loss rate in cgs
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!calculate tau_min, set tau_max
taumin=(xmloss_cgs*sigmae*cconst2/(4.d0*pi)) / (vmax*1.d5 * rlim*sr)
taumax=20.d0
!
!calculate equidistant tau-steps
deltau=(taumax-taumin)/(nr_modext)
!
rmin=1.d0
!
!-----------------------------------------------------------------------
!
!calculate nr_modext radial grid points that are equidistant in tau
!
r_modext1d=0.d0
r_modext1d(1)=rmin
!
!write(*,*) nr_modext
!
outer: do j=1, 100000
!outer cyvle adapts the optical depth-steps deltau (controlled by taumax)
!            until maximum radial grid point is less than rlim
   deltau=(taumax-taumin)/(nr_modext-1)
!   write(*,*) deltau
!
   do i=1, nr_modext-1
!calculate velocity at given radius
      velr = bvel(r_modext1d(i), vmax*1.d5, bconst, beta)
!calculate opacity at given radius (in cgs)
      sigem = sigmae * cconst2 * xmloss_cgs / (4.d0*pi*velr*sr**2 * r_modext1d(i)**2)
!calculate radial grid step at given radius (in r_star)
      delr = -deltau/sigem / sr
!calculate radius for next grid point
      r_modext1d(i+1)=r_modext1d(i)-delr
!go to outer cycle if the calculated radius exceeds rlim
      if(r_modext1d(i+1).gt.rlim) then
         taumax=taumax*999.d0/1000.d0
         exit
      endif
   enddo

!   write(*,'(i5,20es20.8)') j, r_modext1d(1), r_modext1d(nr_modext), maxval(r_modext1d), rlim
!
!stop outer cycle if the maximum calculated raidus is less than rlim
   if(maxval(r_modext1d).lt.rlim) then
      write(*,*) 'found most reasonable radial grid, with maximum values'
      write(*,*) j, 'max', maxval(r_modext1d)
      exit outer
   endif
!
enddo outer
!
!
!
if(maxval(r_modext1d).gt.rlim) then
   stop 'error in calc_mod1d: adapt cycle'
endif
!
!-----------------------------------------------------------------------
!
!overwrite innermost and outermost grid points
!
r_modext1d(nr_modext)=rlim
r_modext1d(1)=rmin
!
!---------------------setting up all other grids------------------------
!
do i=1, nr_modext
!velocitiy
   velr_modext1d(i) = bvel(r_modext1d(i), vmax*1.d5, bconst, beta)
!
!thermal velocity
   vth_modext1d(i) =  vthermal(vmicro*1.d5, teff, na)
!
!density
   rho_modext1d(i) = xmloss_cgs/(4.d0*pi*r_modext1d(i)**2 * sr**2 * velr_modext1d(i))
!
!thermalization parameter
   eps_cont_modext1d(i) = eps_cont
!
!opacity (only locally for lucy temperature structure
   opac_modext1d(i) = opac_thomson(yhe, hei, rho_modext1d(i), kcont)*sr
!   write(*,*) yhe, hei, rho_modext1d(i)*r_modext1d(i)**2, kcont, opac_modext1d(i), opac_modext1d(i)/rho_modext1d(i)/sr
!
enddo
!
!tau=0.d0
!do i=nr_modext-1, 1, -1
!   tau = tau + 0.5d0*(opac_modext1d(i)+opac_modext1d(i+1))*(r_modext1d(i+1)-r_modext1d(i))
!   write(*,*) tau
!enddo
!stop
!write(*,*) rho_modext1d(1), velr_modext1d(1), xmloss_cgs/4.d0/pi/r_modext1d(1)**2/sr**2/velr_modext1d(1) * 0.34*sr/3.
!stop
!read temperature stratification from jo's model and interpolate
!call read_t1d(model1d_jo)
!t_modext1d=23.d3
call calc_tlucy1d(sr, teff, nr_modext, opac_modext1d, r_modext1d, t_modext1d)
!write(*,*) t_modext1d
!stop 'go on in model1d'
!
!
!set radial grid to cgs
r_modext1d=r_modext1d*sr
!
end subroutine calc_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_t1d(input_dir)
!
!--------read in model atmosphere from file 'continuum_jo.dat'----------
!----------to obtain temperature stratification from JOs grid-----------
!
use prog_type
use fund_const, only: rsu
use dime_modext, only: nr_modext
use model1d, only: r_modext1d, t_modext1d
use params_input, only: vmax, vmin, beta, xmloss, yhe, hei, rstar, teff, opt_incl_line
use mod_interp1d
!
implicit none
!
! ... arguments
character(len=*), intent(in) :: input_dir
!
! ... local scalars
integer(i4b) :: nr_jo
integer(i4b) :: i, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: idum
real(dp) :: fdum
real(dp) :: vmax_jo, vmin_jo, beta_jo, xmloss_jo, yhe_jo, hei_jo, sr_jo, teff_jo, &
            eps_cont_jo, eps_line_jo, kcont_jo, kline_jo, kappa0_jo, alpha_jo, &
            corrfc_jo, xic1_jo, lambda_jo, vth_jo
!
! ... local arrays
real(dp), dimension(:), allocatable :: r_jo, t_jo, rho_jo, opac_jo, tau_jo, mintray_jo, mintmom_jo, scontray_jo, scontmom_jo, &
                                       sline_jo, ssoboc_jo, ssobo_jo, scont_jo, velr_jo, opalbar_jo
!
! ... local logicals
logical :: check_fname
!
! ... local functions
!
! ... local characters
character(len=300) :: dum_line
!
check_fname=.false.
!
!----------reading temperature and solution from JOs output-------------
!
if(opt_incl_line) then
   inquire(file=trim(input_dir)//'/line_jo.dat', exist=check_fname)
   if(.not.check_fname) stop 'error in read_t1d: input file does not exist'
   write(*,*) 'reading temperature structure from '//trim(input_dir)//'/line_jo.dat'
   write(*,*)
!
   open(1, file=trim(input_dir)//'/line_jo.dat', form='formatted')
      read(1,*)
      read(1,'(i5, 6es20.8, i5, 12es20.8)') nr_jo, eps_cont_jo, eps_line_jo, kcont_jo, kline_jo, kappa0_jo, alpha_jo, &
                                            idum, corrfc_jo, teff_jo, xic1_jo, lambda_jo, vth_jo, xmloss_jo, vmin_jo, &
                                            vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
      read(1,*)
      read(1,*)
!
      allocate(r_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(t_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(rho_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(tau_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(opac_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(opalbar_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(velr_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(sline_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(ssoboc_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(ssobo_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(scont_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
!
      do i=1, nr_jo
      read(1,'(11es20.8)') r_jo(i), t_jo(i), rho_jo(i), tau_jo(i), opac_jo(i), opalbar_jo(i), velr_jo(i), &
                           sline_jo(i), ssoboc_jo(i), ssobo_jo(i), scont_jo(i)
      enddo
   close(1)


else

!reading continuum solution
   inquire(file=trim(input_dir)//'/continuum_jo.dat', exist=check_fname)
   if(.not.check_fname) stop 'error in read_t1d: input file does not exist'
   write(*,*) 'reading temperature structure from '//trim(input_dir)//'/continuum_jo.dat'
   write(*,*)
!
   open(1, file=trim(input_dir)//'/continuum_jo.dat', form='formatted')
      read(1,*)
      read(1,'(i5, 13es20.8)') nr_jo, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xmloss_jo, &
                               vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
      read(1,*)
      read(1,*)
      allocate(r_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(t_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(rho_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(opac_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(tau_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(mintray_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(mintmom_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(scontray_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
      allocate(scontmom_jo(nr_jo), stat=err)
         if(err.ne.0) stop 'allocation error in read_t1d'
!
      do i=1, nr_jo
         read(1,'(11es20.8)') r_jo(i), t_jo(i), rho_jo(i), opac_jo(i), tau_jo(i), mintray_jo(i), mintmom_jo(i), &
                              scontray_jo(i), scontmom_jo(i), fdum, fdum
      enddo
!
   close(1)
!
endif
!
!----check if input data is consistent (otherwise, wrong comparison)----
!
   write(*,'(a10, 2a30)') 'quantity', 'input from JOs program', 'input from own program'
   write(*,'(a10, 2es30.8)') 'teff', teff_jo, teff
   if(abs(teff-teff_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,'(a10, 2es30.8)') 'mdot', xmloss_jo, xmloss
   if(abs(xmloss-xmloss_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent' 
   write(*,'(a10, 2es30.8)') 'vmin', vmin_jo, vmin*1.d5
   if(abs(vmin*1.d5-vmin_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,'(a10, 2es30.8)') 'vmax', vmax_jo, vmax*1.d5
   if(abs(vmax*1.d5-vmax_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent' 
   write(*,'(a10, 2es30.8)') 'beta', beta_jo, beta
   if(abs(beta-beta_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,'(a10, 2es30.8)') 'yhe', yhe_jo, yhe
   if(abs(yhe-yhe_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,'(a10, 2es30.8)') 'hei', hei_jo, hei
   if(abs(hei-hei_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,'(a10, 2es30.8)') 'sr', sr_jo, rstar*rsu
   if(abs(rstar*rsu-sr_jo).gt.1.d-10) stop 'error in read_t1d: input data not consistent'
   write(*,*)
!
!
!-----interpolating temperature from JOs grid onto my own grid----------
!
!interpolation
do i=1, nr_modext
   call find_index(r_modext1d(i), r_jo, nr_jo, iim2, iim1, ii, iip1)
   t_modext1d(i)=interpol_yp(r_jo(ii), r_jo(iim1), t_jo(ii), t_jo(iim1), r_modext1d(i))
enddo
!
!
!
end subroutine read_t1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod1d_cis
!
!---calculate model atmosphere from average of cis's ldi simulations----
!
use prog_type
use fund_const, only: cgs_mp, cgs_kb, xmsu, rsu, yr, pi
use dime_modext, only: nr_modext
use model1d, only: r_modext1d, t_modext1d, rho_modext1d, velr_modext1d, vth_modext1d, eps_cont_modext1d
use params_input, only: rlim, rmin, teff, vmicro, na, tmin, eps_cont, xmloss, vmin, vmax, beta
use params_stellar, only: sr
!
implicit none
!
! ... local scalars
!for 1d-models from cis (including minimum and maximum number of timesteps)
integer(i4b), parameter :: nr_1d_cis=2000
!
integer(i4b) :: i, j
integer(i4b) :: err, idum
real(dp) :: cconst1
real(dp) :: fdum, vmin_cgs, vmax_cgs, bconst, mdot_cgs

integer(i4b), dimension(:), allocatable :: imask_smooth1d, ngauss_smooth1d
real(dp), dimension(:), allocatable :: sigma_smooth1d
!
! ... local characters
!for b dwarf
!character(len=41), parameter :: model1d_cis='models/cis/output_bstar_cooling_Feldmeier'
!for o dwarf (no cooling correction)
!character(len=35), parameter :: model1d_cis='models/cis/output_ostar_100xQmax_rc'
!for o dwarf (with cooling correction)
character(len=45), parameter :: model1d_cis='models/cis/output_ostar_100xQmax_rc_feldmeier'
!for o star
!character(len=34), parameter :: model1d_cis='models/cis/output_Ostar_rc_15xTeff'
character(len=500) :: fname
character(len=4) :: chdum
!
! ... local logicals
logical :: check1

! ... local functions
real(dp) :: vthermal
!
! ... local arrays
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_1d_cis
!
allocate(r_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(velr_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(rho_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(t_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(vth_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(eps_cont_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
!
!-----------------------------------------------------------------------
!
r_modext1d=0.d0
!
velr_modext1d=0.d0
rho_modext1d=0.d0
t_modext1d=0.d0
vth_modext1d=0.d0
!
write(fname,*) model1d_cis//'/LDI_statvars.vac'
fname=adjustl(fname)
!write(*,*) trim(fname)
inquire(file=trim(fname), exist=check1)
if(.not.check1) then
   write(*,*) 'error in calc_mod1d_cis: file "', trim(fname), '" does not exist'
   stop
endif
!
open(1, file=trim(fname), form='formatted')
   read(1,*) chdum
   do i=1, nr_modext
      read(1,'(i5, 9e22.14)') idum, r_modext1d(i), rho_modext1d(i), fdum, velr_modext1d(i), fdum, fdum, fdum, t_modext1d(i), fdum
   enddo
close(1)
!   
!normalize to lowermost point
r_modext1d = r_modext1d/r_modext1d(1)
!
!write(*,*) r_modext1d(nr_modext), rlim, maxval(r_modext1d)
if(maxval(r_modext1d).gt.rlim) then
   stop 'error in calc_mod1d_cis: increase rlim'
endif
!
!-----------------------------------------------------------------------
!
!overwrite innermost and outermost grid points
!
rmin=1.d0
r_modext1d(nr_modext)=rlim
r_modext1d(1)=rmin
!
!-------smooth radial velocity to avoid weird behaviour outside--------
!
!smooth over delta_r = 0.5d0 with maximum 5 grid points
!allocate(sigma_smooth1d(nr_modext))
!allocate(ngauss_smooth1d(nr_modext))
!allocate(imask_smooth1d(nr_modext))
!!
!do i=1, nr_modext
!   if(r_modext1d(i).gt.2.) then
!      sigma_smooth1d(i) = 0.5d0
!      ngauss_smooth1d(i) = 91
!      imask_smooth1d(i) = 1
!   else
!      sigma_smooth1d(i) = 1.d-5
!      ngauss_smooth1d(i) = 3
!      imask_smooth1d(i) = 0
!   endif
!enddo
!
!call smooth_gauss(nr_modext, r_modext1d, velr_modext1d, sigma_smooth1d, ngauss_smooth1d, imask_smooth1d)
!
!do i=1, nr_modext
!   write(*,*) i, r_modext1d(i), velr_modext1d(i)
!enddo
!stop
!
!test beta velocity law
vmin_cgs=vmin*1.d5
vmax_cgs=vmax*1.d5
mdot_cgs=xmloss*xmsu/yr
bconst=1.-(vmin_cgs/vmax_cgs)**(1./beta)
velr_modext1d = vmax_cgs*(1.-bconst/r_modext1d)**beta
rho_modext1d = mdot_cgs/4.d0/pi/(r_modext1d*sr)**2/velr_modext1d
!
!write(*,*) r_modext1d
!write(*,*) velr_modext1d
!
!---------------------setting up all other grids------------------------
!
!set temperature to (constant) minimum wind temperature
!write(*,*) tmin
t_modext1d = tmin
!
do i=1, nr_modext
!thermal velocity
   vth_modext1d(i) =  vthermal(vmicro*1.d5, teff, na)
!
!thermalization parameter
   eps_cont_modext1d(i) = eps_cont   
   
!
enddo
!
!
!set radial grid to cgs
r_modext1d=r_modext1d*sr
!
end subroutine calc_mod1d_cis
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod1d_luka
!
!---calculate model atmosphere from Lukas input file (r_v_dvdr)---------
!
use prog_type
use fund_const, only: cgs_mp, cgs_kb, xmsu, pi
use dime_modext, only: nr_modext
use model1d, only: r_modext1d, t_modext1d, rho_modext1d, velr_modext1d, vth_modext1d, eps_cont_modext1d
use params_input, only: rlim, rmin, teff, vmicro, na, tmin, vmax, xmloss, eps_cont
use params_stellar, only: sr
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_1d_luka=2000
integer(i4b) :: i, j
integer(i4b) :: err, idum
real(dp) :: cconst1
real(dp) :: fdum, xmloss_cgs
!
! ... local characters
character(len=11), parameter :: model1d_luka='models/luka'
character(len=500) :: fname
character(len=4) :: chdum
!
! ... local logicals
logical :: check1

! ... local functions
real(dp) :: vthermal
!
! ... local arrays
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_1d_luka
!
allocate(r_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
allocate(velr_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
allocate(rho_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
allocate(t_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
allocate(vth_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
allocate(eps_cont_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_luka'
!
!-----------------------------------------------------------------------
!
r_modext1d=0.d0
!
velr_modext1d=0.d0
rho_modext1d=0.d0
t_modext1d=0.d0
vth_modext1d=0.d0
!
write(fname,*) model1d_luka//'/input_r_v_dvdr_Luka.dat'
!write(*,*) trim(fname)
inquire(file=trim(fname), exist=check1)
if(.not.check1) then
   write(*,*) 'error in calc_mod1d_luka: file "', trim(fname), '" does not exist'
   stop
endif
!
open(1, file=trim(fname), form='formatted')
   do i=1, nr_modext
      read(1,*) r_modext1d(i), velr_modext1d(i), fdum
   enddo
close(1)
!   
!normalize to lowermost point
r_modext1d = r_modext1d/r_modext1d(1)
!
!-----------------------------------------------------------------------
!
!overwrite innermost grid points
!
rmin=1.d0
r_modext1d(1)=rmin
!
!---------------------setting up all other grids------------------------
!
!set temperature to (constant) minimum wind temperature
!write(*,*) tmin
t_modext1d = tmin
!
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
do i=1, nr_modext
!
!velocity in cgs
   velr_modext1d(i)=velr_modext1d(i)*vmax*1.d5
!
!radius in cgs
   r_modext1d(i)=r_modext1d(i)*sr
!
!density in cgs
   rho_modext1d(i) = xmloss_cgs/4.d0/pi/r_modext1d(i)**2/velr_modext1d(i)
!  
!thermal velocity
   vth_modext1d(i) =  vthermal(vmicro*1.d5, teff, na)  
!
!thermalization parameter
   eps_cont_modext1d(i) = eps_cont
!
enddo
!
!
end subroutine calc_mod1d_luka
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod1d_test
!
!---------------calculate model atmosphere for 1d test models-----------
!
use prog_type
use fund_const, only: cgs_mp, cgs_kb, pi, sigmae, xmsu
use dime_modext, only: nr_modext
use model1d, only: r_modext1d, t_modext1d, rho_modext1d, velr_modext1d, vth_modext1d, eps_cont_modext1d
use params_input, only: vmax, vmin, beta, xmloss, yhe, hei, rlim, rmin, teff, tmin, vmicro, na, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_log
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nr_1d=1000
integer(i4b) :: i, j
integer(i4b) :: err
real(dp) :: bconst, cconst1, cconst2, xmloss_cgs
real(dp) :: sigem, delr, velr
!
! ... local characters
character(len=9), parameter :: model1d_jo='models/jo'
!
! ... local functions
real(dp) :: bvel, vthermal
!
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_1d
!
allocate(r_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(velr_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(rho_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(t_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(vth_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
allocate(eps_cont_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d'
!
!------------------calculate/define required constants------------------
!
!b-factor for beta-velocity-law
bconst=1.d0-(vmin/vmax)**(1.d0/beta)
!
cconst1=(1.d0+4.d0*yhe)*cgs_mp
cconst2=(1.d0+hei*yhe)/cconst1
!
!mass-loss rate in cgs
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!
rmin=1.d0
!
r_modext1d=0.d0
call grid_log(rmin, rlim, nr_modext, r_modext1d)
r_modext1d(1)=rmin
r_modext1d(nr_modext)=rlim
!
!-----------------------------------------------------------------------
!
do i=1, nr_modext
   velr_modext1d(i) = bvel(r_modext1d(i), vmax*1.d5, bconst, beta)
   vth_modext1d(i) =  vthermal(vmicro*1.d5, teff, na)
   rho_modext1d(i) = xmloss_cgs/(4.d0*pi*r_modext1d(i)**2 * sr**2 * velr_modext1d(i))
   eps_cont_modext1d(i) = eps_cont!1.d0/r_modext1d(i) !eps_cont
   t_modext1d(i) = teff !tmin
enddo
!
!-----------------------------------------------------------------------
!
!set radial grid to cgs
r_modext1d=r_modext1d*sr
!
end subroutine calc_mod1d_test
