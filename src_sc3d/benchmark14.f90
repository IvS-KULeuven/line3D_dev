subroutine benchmark14_solution
!
!---------calculating 3d line + continuum radiative transfer------------
!-----------------for spherically symmetric problems--------------------
!------with the theoretical solution given from JOs 1d programs---------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb
use dimecr, only: n1d_cr, r1d_cr
use dime3d, only: t3d, eps_cont3d, mint3d, mintbar3d, int3d, opac3d, opalbar3d, &
                  scont3d, sline3d, ssobo3d, velx3d, vely3d, velz3d, vth3d, &
                  ndxmax, ndymax, ndzmax, x, y, z, imask_innreg3d, &
                  alocont_nn3d, aloline_nn3d, normalization3d, imask3d, imask_totreg3d
use angles, only: dim_mu, dim_phi, dim_omega
use bcondition, only: xic1, xic2
use mod_benchmark, only: epsmaxc_sc, epsmaxc_fvm, epsmaxl_sc, epsmaxl_fvm, opalbar1d_cr, &
                         opalbar1d_jo, ssobo1d_jo, sline1d_jo, scont1d_jo, &
                         t1d_cr, t1d_jo, opac1d_cr, opac1d_jo, sline1d_sc, sline1d_fvm, ssobo1d_cr, &
                         mintbar1d_sc, mintbar1d_fvm, velr1d_cr, &
                         mint3d_sc, mint3d_fvm, scont3d_sc, scont3d_fvm, &
                         mintbar3d_sc, mintbar3d_fvm, sline3d_sc, sline3d_fvm
use freq, only: xnue0
use iter, only: itmaxl, itmaxc, devmaxc, devmaxl, epsmaxl_arr, epsmaxc_arr, it_start_line
use ng_extra, only: ng_const
use options, only: opt_ng_line, opt_ait_line, opt_ng_cont, opt_ait_cont, input_mod_dim, opt_opal, opt_alo_line
use params_input, only: teff, rstar, beta, vmin, vmax, hei, yhe, xmloss, kcont, eps_line, kline, &
                        kappa0, alpha, vth_fiducial
use params_stellar, only: sr
use mod_debug, only: indx1, indx2, indxx, indxy, indxz
use warnings, only: warn_itmaxl, warn_itmaxc
use mod_interp1d, only: find_index, interpol_yp
use mod_interp2d, only: interpolation_threshold, wpa_interp2d, wpb_interp2d, wp_interp2d, &
     wpa_interp1d, wpb_interp1d, wp_interp1d, wpa_integ1d, wpb_integ1d, wp_integ1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, muindx, err
integer(i4b) :: iim2, iim1, ii, iip1, indx_threshold
integer(i4b) :: ix_epsmax, iy_epsmax, iz_epsmax
integer(i4b) :: nr_jo, opt_opal_jo
real(dp) :: eps_cont, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xnue0_jo, xmloss_jo, vth_jo, &
            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo, eps_line_jo, kline_jo, kappa0_jo, alpha_jo
real(dp) :: s1, s2, s3, s4, s4b, eps_max, fdum, dummy1, dummy2, rad, trad
real(dp) :: ts, te
!
! ... local arrays
real(dp), dimension(:,:,:), allocatable :: eps3d
real(dp), dimension(:,:), allocatable :: sline3d_ng, scont3d_ng
real(dp), dimension(:), allocatable :: r_jo, t_jo, rho_jo, opac_jo, tau_jo, opalbar_jo, &
                                       sline_jo, ssoboc_jo, ssobo_jo, scont_jo, velr_jo
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
!---------------------check if input model is consistent----------------
!
call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
inquire(file='./models/jo/line_jo.dat', exist=check_fname)
if(.not.check_fname) stop 'error in benchmark14: 1d solution file does not exist'
!
open(1, file='./models/jo/line_jo.dat', form='formatted')
   read(1,*)
   read(1,'(i5, 6es20.8, i5, 12es20.8)') nr_jo, eps_cont_jo, eps_line_jo, kcont_jo, kline_jo, kappa0_jo, alpha_jo, &
                                         opt_opal_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, vth_jo, xmloss_jo, vmin_jo, &
                                         vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
   read(1,*)
   read(1,*)

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

   do i=1, nr_jo
   read(1,'(11es20.8)') r_jo(i), t_jo(i), rho_jo(i), tau_jo(i), opac_jo(i), opalbar_jo(i), velr_jo(i), &
                        sline_jo(i), ssoboc_jo(i), ssobo_jo(i), scont_jo(i)
   enddo
close(1)
!
write(*,'(3a20)') 'quantity', 'own input', 'JOs input'
write(*,'(a20, 2es20.8)') 'teff', teff, teff_jo
if(abs(teff-teff_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'mdot', xmloss, xmloss_jo
if(abs(xmloss-xmloss_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmin', vmin, vmin_jo
if(abs(vmin*1.d5-vmin_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmax', vmax, vmax_jo
if(abs(vmax*1.d5-vmax_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'beta', beta, beta_jo
if(abs(beta-beta_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'yhe', yhe, yhe_jo
if(abs(yhe-yhe_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'hei', hei, hei_jo
if(abs(hei-hei_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'rstar', sr, sr_jo
if(abs(sr-sr_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'eps_cont', eps_cont, eps_cont_jo
if(abs(eps_cont-eps_cont_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'kcont', kcont, kcont_jo
if(abs(kcont-kcont_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
xnue0_jo=cgs_clight/lambda_jo/1.d-8
write(*,'(a20, 2es20.8)') 'xnue', xnue0, xnue0_jo
if(abs(xnue0-xnue0_jo)/xnue0.gt.1.d-6) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2es20.8)') 'xic1', xic1, xic1_jo
if(abs(xic1-xic1_jo).gt.1.d-5) then
   write(*,*) 'error in benchmark14: input data not consistent'
   trad = log(2.d0*cgs_planck*xnue0**3/cgs_clight**2/xic1_jo + 1.d0)
   trad = cgs_planck*xnue0/trad/cgs_kb
   write(*,*) 'use as radiation temperature on boundary:', trad
   stop
endif
write(*,'(a20, 2es20.8)') 'eps_line', eps_line, eps_line_jo
if(abs(eps_line-eps_line_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
write(*,'(a20, 2i20)') 'opt_opal', opt_opal, opt_opal_jo
if(opt_opal.ne.opt_opal_jo) stop 'error in benchmark14: input data not consistent'
if(opt_opal.eq.0) then
   write(*,'(a20, 2es20.8)') 'kline', kline, kline_jo
   if(abs(kline-kline_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
elseif(opt_opal.eq.1) then
   write(*,'(a20, 2es20.8)') 'kappa0', kappa0, kappa0_jo
   if(abs(kappa0-kappa0_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
   write(*,'(a20, 2es20.8)') 'alpha', alpha, alpha_jo
   if(abs(alpha-alpha_jo).gt.1.d-10) stop 'error in benchmark14: input data not consistent'
endif

write(*,*)
!
!interpolate everything onto own grid
allocate(sline1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(scont1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(ssobo1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(ssobo1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(opac1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(opalbar1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(velr1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(opac1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(opalbar1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(t1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(t1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
!
do i=1, n1d_cr
   t1d_cr(i) = t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   opac1d_cr(i) = opac3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   opalbar1d_cr(i) = opalbar3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   ssobo1d_cr(i) = ssobo3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   velr1d_cr(i) = vth_fiducial*sqrt(velx3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)**2 + &
                       vely3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)**2 + &
                       velz3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)**2)
   rad = z(ndzmax-n1d_cr-1+i)
   call find_index(rad, r_jo, nr_jo, iim2, iim1, ii, iip1)
   opac1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), opac_jo(iim1), opac_jo(ii), rad)
   opalbar1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), opalbar_jo(iim1), opalbar_jo(ii), rad)
   t1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), t_jo(iim1), t_jo(ii), rad)
   sline1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), sline_jo(iim1), sline_jo(ii), rad)
   ssobo1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), ssobo_jo(iim1), ssobo_jo(ii), rad)
   scont1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), scont_jo(iim1), scont_jo(ii), rad)
enddo
!
!in cgs
opac1d_cr=opac1d_cr/sr
!
!-----------------------------------------------------------------------
!
allocate(mint3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(mint3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(scont3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(scont3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(mintbar3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(mintbar3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(sline3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(sline3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
!
allocate(epsmaxc_sc(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(epsmaxc_fvm(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(epsmaxl_sc(itmaxl), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(epsmaxl_fvm(itmaxl), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
!
allocate(eps3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(sline3d_ng(4,ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
allocate(scont3d_ng(4,ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark14'
!
epsmaxl_sc=0.d0
epsmaxl_fvm=0.d0
epsmaxc_sc=0.d0
epsmaxc_fvm=0.d0
!
s1=1
s2=2
s3=3
s4=4
s4b=4
!
!***********************************************************************
!***********************************************************************
!
!                    SC CONTINUUM TRANSPORT
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!
!setting start mean intensities and start deviations to zero
mint3d=0.d0
eps3d=0.d0
normalization3d=0.d0
alocont_nn3d=0.d0
epsmaxc_arr=0.d0
!
call calc_startval_cont
!
!setting index for threshold
indx_threshold=5
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxc
!
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'z', 'j(sc)', 'scont(sc)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 6(e20.8))') j, z(j), mint3d(ndxmax/2+1,ndymax/2+1,j), scont3d(ndxmax/2+1,ndymax/2+1,j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    alocont_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mint3d
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities in 3d--------------'
   write(*,*) 'step', i
   write(*,*)
!
   if(i.le.s4b) then
!use linear interpolations for the first s4b iteration steps
!in order to get a first guess of a 'smooth' solution.
!otherwise: convergence behaviour might be oscillatory
      call mint_sc3d_lin
   else
      call mint_sc3d_lin
   endif
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
   call output_itstep_cont(i)
!
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) 'convergence after iteration no. ', i
      write(*,*) 'max (dev): ', eps_max
      write(*,*)
!for optically thin models, ensure that higher order interpolation scheme
!   has been used, and not only linear
      if(i.gt.s4b) exit  
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
         wpa_integ1d=wpa_integ1d+1.d0
         wpb_integ1d=wpb_integ1d+1.d0
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-1.d0)/(wpb_interp2d-1.d0), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-1.d0)/(wpb_integ1d-1.d0), wp_integ1d
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
!-----------------------------------------------------------------------
!
!store everything in benchmark12-arrays
epsmaxc_sc = epsmaxc_arr
mint3d_sc = mint3d
scont3d_sc = scont3d
!
!***********************************************************************
!***********************************************************************
!
!                    SC LINE + CONTINUUM TRANSPORT
!
!***********************************************************************
!***********************************************************************
!
!---------------------calculate sobolev solution------------------------
!
call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
!
call calc_startval_line
!
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!         select case(imask3d(i,j,k))
!            case(1,2,3)
!               sline3d(i,j,k) = ssobo3d(i,j,k)
!            case default
!         end select
!      enddo
!   enddo
!enddo
!
!***debug start
!opac3d=1.d-20
!***debug end
!
!-----------------------------------------------------------------------
!
!setting start mean intensities and start deviations to zero
mintbar3d=0.d0
eps3d=0.d0
normalization3d=0.d0
aloline_nn3d=0.d0
epsmaxl_arr=0.d0
!
!setting index for threshold
indx_threshold=it_start_line+5
!
!---------------------start iteration scheme----------------------------
!
do i=it_start_line+1, itmaxl
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 8(a20))') '#', 'z', 'jbar(sc)', 'sline(sc)', 'sline(sobo)', 'scont(sc)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 8(e20.8))') j, z(j), mintbar3d(ndxmax/2+1,ndymax/2+1,j), sline3d(ndxmax/2+1,ndymax/2+1,j), &
                                    ssobo3d(ndxmax/2+1, ndymax/2+1, j), scont3d(ndxmax/2+1, ndymax/2+1, j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    aloline_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mintbar3d
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '-------calculating frequency and angle integrated mean intensities in 3d-------'
   write(*,*) 'step', i
   write(*,*)
!
   if(i.le.s4b) then
!use linear interpolations for the first s4b iteration steps
!in order to get a first guess of a 'smooth' solution.
!otherwise: convergence behaviour might be oscillatory,
!           or worse, first correction might give negative values
      call mintbar_sc3d_lin
   else
      call mintbar_sc3d_lin
   endif
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
   call output_itstep_line(i)
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
!
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
!-----------------------------------------------------------------------
!
!store everything in benchmark14-arrays
epsmaxl_sc = epsmaxl_arr
mintbar3d_sc = mintbar3d
sline3d_sc = sline3d
!
!***********************************************************************
!***********************************************************************
!
!                    FVM CONTINUUM TRANSPORT
!
!***********************************************************************
!***********************************************************************
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
!setting start mean intensities and start deviations to zero
mint3d=0.d0
eps3d=0.d0
normalization3d=0.d0
alocont_nn3d=0.d0
epsmaxc_arr=0.d0
!
scont3d=0.d0
call calc_startval_cont
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxc
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'z', 'j(fvm)', 'scont(fvm)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 6(e20.8))') j, z(j), mint3d(ndxmax/2+1,ndymax/2+1,j), scont3d(ndxmax/2+1,ndymax/2+1,j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    alocont_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mint3d
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities in 3d--------------'
   write(*,*) 'step', i
   write(*,*)
!
   call mint_fvm3d
!
!   write(*,*) '----------------------calculating alo in sparse formats------------------------'
!   call alo_sparse
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
   call output_itstep_cont(i)
!
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      write(*,*)
      exit
   elseif(i.eq.itmaxc) then
      write(*,*) "no convergence after iteration no. ", i
      write(*,*)
      warn_itmaxc=.true.
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
!-----------------------------------------------------------------------
!
!store everything in benchmark14-arrays
epsmaxc_fvm = epsmaxc_arr
mint3d_fvm = mint3d
scont3d_fvm = scont3d
!
!***********************************************************************
!***********************************************************************
!
!                    FVM LINE + CONTINUUM TRANSPORT
!
!***********************************************************************
!***********************************************************************
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
call calc_startval_line
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               sline3d(i,j,k) = ssobo3d(i,j,k)
            case default
         end select
      enddo
   enddo
enddo
!
!setting start mean intensities and start deviations to zero
mintbar3d=0.d0
eps3d=0.d0
normalization3d=0.d0
aloline_nn3d=0.d0
epsmaxl_arr=0.d0
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxl
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 7(a20))') '#', 'z', 'jbar(fvm)', 'sline(fvm)', 'sline(sobo)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=ndzmax-n1d_cr-3, ndzmax
      write(*,fmt='(i5, 7(e20.8))') j, z(j), mintbar3d(ndxmax/2+1,ndymax/2+1,j), sline3d(ndxmax/2+1,ndymax/2+1,j), &
                                    ssobo3d(ndxmax/2+1,ndymax/2+1,j), &
                                    eps3d(ndxmax/2+1,ndymax/2+1,j), normalization3d(ndxmax/2+1,ndymax/2+1,j), &
                                    aloline_nn3d(ndxmax/2+1,ndymax/2+1,j,14)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   eps3d=mintbar3d
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '-------calculating frequency and angle integrated mean intensities in 3d-------'
   write(*,*) 'step', i
   write(*,*)
!
   call mintbar_fvm3d
!
!   write(*,*) '----------------------calculating alo in sparse formats------------------------'
!   call alo_sparse
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
   call output_itstep_line(i)
!
   if(abs(eps_max).lt.devmaxl) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      write(*,*)
      exit
   elseif(i.eq.itmaxl) then
      write(*,*) "no convergence after iteration no. ", i
      write(*,*)
      warn_itmaxl=.true.
    endif
!
!--------------calculating alo-corrected source-functions---------------
!
   call sline_new3d
!
!------------extrapolation of old subsequent source functions-----------
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
!-----------------------------------------------------------------------
!
!store everything in benchmark14-arrays
epsmaxl_fvm = epsmaxl_arr
mintbar3d_fvm = mintbar3d
sline3d_fvm = sline3d
!
end subroutine benchmark14_solution
