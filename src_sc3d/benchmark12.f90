subroutine benchmark12_solution
!
!------calculating 3d continuum radiative transfer for spherically------
!----symmetric problems, with the theoretical solution given from-------
!---------------------------JOs 1d programs-----------------------------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb, zero
use dimecr, only: n1d_cr, r1d_cr
use dime3d, only: t3d, eps_cont3d, mint3d, int3d, opac3d, scont3d, ndxmax, ndymax, ndzmax, x, y, z, &
                  alocont_nn3d, alocont_o_nn3d, normalization3d, imask3d, imask_totreg3d, &
                  fcontx3d, fconty3d, fcontz3d, &
                  kcontxx3d, kcontyy3d, kcontzz3d, kcontxy3d, kcontxz3d, kcontyz3d
use angles, only: dim_mu, dim_phi, dim_omega
use bcondition, only: xic1, xic2
use mod_benchmark, only: mint3d_sc, mint3d_fvm, mint1d_joray, mint1d_jomom, fcont1d_joray, fcont1d_jomom, epsmaxc_sc, epsmaxc_fvm, &
                         fcontr3d_sc, fcontr3d_fvm, fcontth3d_sc, fcontth3d_fvm, fcontphi3d_sc, fcontphi3d_fvm, &
                         kcontrr3d_sc, kcontthth3d_sc, kcontphiphi3d_sc, kcontrth3d_sc, kcontrphi3d_sc, kcontthphi3d_sc, &
                         kcontrr3d_fvm, kcontthth3d_fvm, kcontphiphi3d_fvm, kcontrth3d_fvm, kcontrphi3d_fvm, kcontthphi3d_fvm, &
                         t1d_jo, opac1d_jo
use freq, only: xnue0
use iter, only: itmaxc, devmaxc, epsmaxc_arr
use ng_extra, only: ng_const
use options, only: opt_ng_cont, opt_ait_cont, input_mod_dim
use params_input, only: teff, rstar, beta, vmin, vmax, hei, yhe, xmloss, kcont
use params_stellar, only: sr
use mod_debug, only: indx1, indx2, indxx, indxy, indxz
use warnings, only: warn_itmaxc
use mod_interp1d, only: find_index, interpol_yp
use mod_interp2d, only: interpolation_threshold, wpa_interp2d, wpb_interp2d, wp_interp2d, &
                        wpa_integ1d, wpb_integ1d, wp_integ1d
use timing, only: ttot_it_sc, ttot_it_fvm
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, muindx, err
integer(i4b) :: iim2, iim1, ii, iip1, indx_threshold
integer(i4b) :: nr_jo
integer(i4b) :: ix_epsmax, iy_epsmax, iz_epsmax
real(dp) :: eps_cont, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xnue0_jo, xmloss_jo, &
            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
real(dp) :: s1, s2, s3, s4, s4b, eps_max, fdum, dummy1, dummy2, rad, trad, theta, phi
real(dp) :: ts, te
real(dp) :: sint, cost, sinp, cosp, a11, a12, a13, a21, a22, a23, a31, a32, a33
!
! ... local arrays
real(dp), dimension(:,:,:), allocatable :: eps3d
real(dp), dimension(:,:), allocatable :: scont3d_ng
real(dp), dimension(:), allocatable :: r_jo, t_jo, rho_jo, opac_jo, tau_jo, mintray_jo, mintmom_jo, scontray_jo, scontmom_jo, &
                                       fcontray_jo, fcontmom_jo
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
ttot_it_sc=0.d0
ttot_it_fvm=0.d0
!
!---------------------check if input model is consistent----------------
!
call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
if(input_mod_dim.ne.1) stop 'error in benchmark12_solution: input dimension ne 1'
!
check_fname=.false.
inquire(file='./models/jo/continuum_jo.dat', exist=check_fname)
if(.not.check_fname) stop 'error in benchmark12: 1d solution file does not exist'
!
!read 1d solution
open(1, file='./models/jo/continuum_jo.dat', form='formatted')
   read(1,*)
   read(1,'(i5, 13es20.8)') nr_jo, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xmloss_jo, &
                            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
   read(1,*)
   read(1,*)
   allocate(r_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(t_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(rho_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(opac_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(tau_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(mintray_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(mintmom_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(fcontray_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(fcontmom_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(scontray_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
   allocate(scontmom_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark12'
!
   do i=1, nr_jo
      read(1,'(11es20.8)') r_jo(i), t_jo(i), rho_jo(i), opac_jo(i), tau_jo(i), mintray_jo(i), mintmom_jo(i), &
                           scontray_jo(i), scontmom_jo(i), fcontray_jo(i), fcontmom_jo(i)
   enddo
!
close(1)
!
write(*,'(3a20)') 'quantity', 'own input', 'JOs input'
write(*,'(a20, 2es20.8)') 'teff', teff, teff_jo
if(abs(teff-teff_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'mdot', xmloss, xmloss_jo
if(abs(xmloss-xmloss_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmin', vmin, vmin_jo
if(abs(vmin*1.d5-vmin_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmax', vmax, vmax_jo
if(abs(vmax*1.d5-vmax_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'beta', beta, beta_jo
if(abs(beta-beta_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'yhe', yhe, yhe_jo
if(abs(yhe-yhe_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'hei', hei, hei_jo
if(abs(hei-hei_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'rstar', sr, sr_jo
if(abs(sr-sr_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'eps_cont', eps_cont, eps_cont_jo
if(abs(eps_cont-eps_cont_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'kcont', kcont, kcont_jo
if(abs(kcont-kcont_jo).gt.1.d-10) stop 'error in benchmark12: input data not consistent'
xnue0_jo=cgs_clight/lambda_jo/1.d-8
write(*,'(a20, 2es20.8)') 'xnue', xnue0, xnue0_jo
if(abs(xnue0-xnue0_jo)/xnue0.gt.1.d-8) stop 'error in benchmark12: input data not consistent'
write(*,'(a20, 2es20.8)') 'xic1', xic1, xic1_jo
if(abs(xic1-xic1_jo).gt.1.d-5) then
   write(*,*) 'error in benchmark12: input data not consistent'
   trad = log(2.d0*cgs_planck*xnue0**3/cgs_clight**2/xic1_jo + 1.d0)
   trad = cgs_planck*xnue0/trad/cgs_kb
   write(*,*) 'use as radiation temperature on boundary:', trad
   stop
endif
if(abs(xic1-xic1_jo).gt.1.d-5) stop 'error in benchmark12: input data not consistent'
write(*,*)
!
!interpolate everything onto own grid
allocate(mint1d_joray(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(mint1d_jomom(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcont1d_joray(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcont1d_jomom(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(opac1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(t1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
!
do i=1, n1d_cr
   rad = z(ndzmax-n1d_cr-1+i)
   call find_index(rad, r_jo, nr_jo, iim2, iim1, ii, iip1)
   opac1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), opac_jo(iim1), opac_jo(ii), rad)
   t1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), t_jo(iim1), t_jo(ii), rad)
   mint1d_joray(i) = interpol_yp(r_jo(iim1), r_jo(ii), mintray_jo(iim1), mintray_jo(ii), rad)
   mint1d_jomom(i) = interpol_yp(r_jo(iim1), r_jo(ii), mintmom_jo(iim1), mintmom_jo(ii), rad)
   fcont1d_joray(i) = interpol_yp(r_jo(iim1), r_jo(ii), fcontray_jo(iim1), fcontray_jo(ii), rad)
   fcont1d_jomom(i) = interpol_yp(r_jo(iim1), r_jo(ii), fcontmom_jo(iim1), fcontmom_jo(ii), rad)
enddo
!
!-----------------------------------------------------------------------
!
allocate(mint3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(mint3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
!
allocate(epsmaxc_sc(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(epsmaxc_fvm(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
!
allocate(eps3d(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(scont3d_ng(4,ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
!
allocate(fcontr3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcontr3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcontth3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcontth3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcontphi3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(fcontphi3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'

allocate(kcontrr3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontthth3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontphiphi3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
 allocate(kcontrth3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontrphi3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontthphi3d_fvm(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontrr3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontthth3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontphiphi3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontrth3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontrphi3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
allocate(kcontthphi3d_sc(ndxmax,ndymax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark12'
!
!
!
epsmaxc_sc=0.d0
epsmaxc_fvm=0.d0
!
!*************start short characteristics method************************
!
!--------------------------for debugging the alo------------------------
!
!calculate alo at certain indices
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!         if(imask3d(i,j,k).eq.2) then
!            write(*,*) imask3d(i,j,k), i, j, k
!         endif
!      enddo
!   enddo
!enddo
!stop
!
!indxx=ndxmax/2+1
!indxy=ndymax/2+1
!indxz=37
!!
!i=indxx
!j=indxy
!k=indxz
!!
!do l=1, dim_omega
!   call fsc_cont3d(l)
!!  
!   write(*,*) '----------------angle indx', l, '-----------------'
!
!      if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,11)- int3d(i,j-1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.1.d-14.or. &
!         abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)).gt.1.d-14) then
!
!    write(*,*) 'level z(k-1)', z(k-1)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), int3d(i-1,j-1,k-1), int3d(i,j-1,k-1), int3d(i+1,j-1,k-1)
!    write(*,'(f8.4, 3es25.14)') y(j), int3d(i-1,j,k-1), int3d(i,j,k-1), int3d(i+1,j,k-1)
!    write(*,'(f8.4, 3es25.14)') y(j+1), int3d(i-1,j+1,k-1), int3d(i,j+1,k-1), int3d(i+1,j+1,k-1)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,3)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,6)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,9)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)
!    write(*,*)
!
!    write(*,*) 'level z(k)', z(k)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), int3d(i-1,j-1,k), int3d(i,j-1,k), int3d(i+1,j-1,k)
!    write(*,'(f8.4, 3es25.14)') y(j), int3d(i-1,j,k), int3d(i,j,k), int3d(i+1,j,k)
!    write(*,'(f8.4, 3es25.14)') y(j+1), int3d(i-1,j+1,k), int3d(i,j+1,k), int3d(i+1,j+1,k)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,12)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,15)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,18)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k), alocont_o_nn3d(i,j,k,11)- int3d(i,j-1,k), alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k), alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)
!    write(*,*)
!
!    write(*,*) 'level z(k+1)', z(k+1)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), int3d(i-1,j-1,k+1), int3d(i,j-1,k+1), int3d(i+1,j-1,k+1)
!    write(*,'(f8.4, 3es25.14)') y(j), int3d(i-1,j,k+1), int3d(i,j,k+1), int3d(i+1,j,k+1)
!    write(*,'(f8.4, 3es25.14)') y(j+1), int3d(i-1,j+1,k+1), int3d(i,j+1,k+1), int3d(i+1,j+1,k+1)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,21)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,24)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,27)
!    write(*,*)
!    write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!    write(*,'(f8.4, 3es25.14)') y(j-1), alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)
!    write(*,'(f8.4, 3es25.14)') y(j), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
!    write(*,'(f8.4, 3es25.14)') y(j+1), alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)
!    write(*,*)
!   endif
!enddo
!write(*,*) alocont_nn3d(i,j,k,14), mint3d(i,j,k)
!stop
!!write(*,*) 'level z(k-1)', z(k-1)
!!write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!!write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k-1), mint3d(i,j-1,k-1), mint3d(i+1,j-1,k-1)
!!write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k-1), mint3d(i,j,k-1), mint3d(i+1,j,k-1)
!!write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k-1), mint3d(i,j+1,k-1), mint3d(i+1,j+1,k-1)
!!write(*,*)
!!write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!!write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1), alocont_nn3d(i,j,k,2), alocont_nn3d(i,j,k,3)
!!write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4), alocont_nn3d(i,j,k,5), alocont_nn3d(i,j,k,6)
!!write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7), alocont_nn3d(i,j,k,8), alocont_nn3d(i,j,k,9)
!write(*,*)
!write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1)-mint3d(i-1,j-1,k-1), alocont_nn3d(i,j,k,2)-mint3d(i,j-1,k-1), alocont_nn3d(i,j,k,3)-mint3d(i+1,j-1,k-1)
!write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4)-mint3d(i-1,j,k-1), alocont_nn3d(i,j,k,5)-mint3d(i,j,k-1), alocont_nn3d(i,j,k,6)-mint3d(i+1,j,k-1)
!write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7)-mint3d(i-1,j+1,k-1), alocont_nn3d(i,j,k,8)-mint3d(i,j+1,k-1), alocont_nn3d(i,j,k,9)-mint3d(i+1,j+1,k-1)
!write(*,*)
!stop


!do i=indxx, indxx
!   do j=indxy, indxy
!      do k=indxz, indxz
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!!
!         if(mask_totreg3d(i,j,k)) then
!            indxx=i
!            indxy=j
!            indxz=k
!            write(*,*) 'calculating scont=1 at', i, j, k
!            call mint_sc3d
!
!            write(*,*) imask3d(i,j,k)
!            write(*,*) 'level z(k-1)', z(k-1)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k-1), mint3d(i,j-1,k-1), mint3d(i+1,j-1,k-1)
!            write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k-1), mint3d(i,j,k-1), mint3d(i+1,j,k-1)
!            write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k-1), mint3d(i,j+1,k-1), mint3d(i+1,j+1,k-1)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1), alocont_nn3d(i,j,k,2), alocont_nn3d(i,j,k,3)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4), alocont_nn3d(i,j,k,5), alocont_nn3d(i,j,k,6)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7), alocont_nn3d(i,j,k,8), alocont_nn3d(i,j,k,9)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1)-mint3d(i-1,j-1,k-1), alocont_nn3d(i,j,k,2)-mint3d(i,j-1,k-1), alocont_nn3d(i,j,k,3)-mint3d(i+1,j-1,k-1)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4)-mint3d(i-1,j,k-1), alocont_nn3d(i,j,k,5)-mint3d(i,j,k-1), alocont_nn3d(i,j,k,6)-mint3d(i+1,j,k-1)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7)-mint3d(i-1,j+1,k-1), alocont_nn3d(i,j,k,8)-mint3d(i,j+1,k-1), alocont_nn3d(i,j,k,9)-mint3d(i+1,j+1,k-1)
!            write(*,*)
!
!            write(*,*) 'level z(k)', z(k)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k), mint3d(i,j-1,k), mint3d(i+1,j-1,k)
!            write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k), mint3d(i,j,k), mint3d(i+1,j,k)
!            write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k), mint3d(i,j+1,k), mint3d(i+1,j+1,k)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,10), alocont_nn3d(i,j,k,11), alocont_nn3d(i,j,k,12)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,13), alocont_nn3d(i,j,k,14), alocont_nn3d(i,j,k,15)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,16), alocont_nn3d(i,j,k,17), alocont_nn3d(i,j,k,18)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,10)-mint3d(i-1,j-1,k), alocont_nn3d(i,j,k,11)- mint3d(i,j-1,k), alocont_nn3d(i,j,k,12)-mint3d(i+1,j-1,k)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,13)-mint3d(i-1,j,k), alocont_nn3d(i,j,k,14)-mint3d(i,j,k), alocont_nn3d(i,j,k,15)-mint3d(i+1,j,k)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,16)-mint3d(i-1,j+1,k), alocont_nn3d(i,j,k,17)-mint3d(i,j+1,k), alocont_nn3d(i,j,k,18)-mint3d(i+1,j+1,k)
!            write(*,*)
!
!            write(*,*) 'level z(k+1)', z(k+1)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k+1), mint3d(i,j-1,k+1), mint3d(i+1,j-1,k+1)
!            write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k+1), mint3d(i,j,k+1), mint3d(i+1,j,k+1)
!            write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k+1), mint3d(i,j+1,k+1), mint3d(i+1,j+1,k+1)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,19), alocont_nn3d(i,j,k,20), alocont_nn3d(i,j,k,21)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,22), alocont_nn3d(i,j,k,23), alocont_nn3d(i,j,k,24)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,25), alocont_nn3d(i,j,k,26), alocont_nn3d(i,j,k,27)
!            write(*,*)
!            write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!            write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,19)-mint3d(i-1,j-1,k+1), alocont_nn3d(i,j,k,20)-mint3d(i,j-1,k+1), alocont_nn3d(i,j,k,21)-mint3d(i+1,j-1,k+1)
!            write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,22)-mint3d(i-1,j,k+1), alocont_nn3d(i,j,k,23)-mint3d(i,j,k+1), alocont_nn3d(i,j,k,24)-mint3d(i+1,j,k+1)
!            write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,25)-mint3d(i-1,j+1,k+1), alocont_nn3d(i,j,k,26)-mint3d(i,j+1,k+1), alocont_nn3d(i,j,k,27)-mint3d(i+1,j+1,k+1)
!            write(*,*)

!            if(abs(alocont_nn3d(i,j,k,1)-mint3d(i-1,j-1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,2)-mint3d(i,j-1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,3)-mint3d(i+1,j-1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,4)-mint3d(i-1,j,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,5)-mint3d(i,j,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,6)-mint3d(i+1,j,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,7)-mint3d(i-1,j+1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,8)-mint3d(i,j+1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,9)-mint3d(i+1,j+1,k-1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,10)-mint3d(i-1,j-1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,11)- mint3d(i,j-1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,12)-mint3d(i+1,j-1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,13)-mint3d(i-1,j,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,14)-mint3d(i,j,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,15)-mint3d(i+1,j,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,16)-mint3d(i-1,j+1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,17)-mint3d(i,j+1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,18)-mint3d(i+1,j+1,k)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,19)-mint3d(i-1,j-1,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,20)-mint3d(i,j-1,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,21)-mint3d(i+1,j-1,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,22)-mint3d(i-1,j,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,23)-mint3d(i,j,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,24)-mint3d(i+1,j,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,25)-mint3d(i-1,j+1,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,26)-mint3d(i,j+1,k+1)).gt.1.d-14.or. &
!               abs(alocont_nn3d(i,j,k,27)-mint3d(i+1,j+1,k+1)).gt.1.d-14) then
!
!               write(*,*) 'error of alo at'
!               write(*,*) 'indxx', indxx
!               write(*,*) 'indxy', indxy
!               write(*,*) 'indxz', indxz
!
!               write(*,*) 'level z(k-1)', z(k-1)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k-1), mint3d(i,j-1,k-1), mint3d(i+1,j-1,k-1)
!               write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k-1), mint3d(i,j,k-1), mint3d(i+1,j,k-1)
!               write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k-1), mint3d(i,j+1,k-1), mint3d(i+1,j+1,k-1)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1), alocont_nn3d(i,j,k,2), alocont_nn3d(i,j,k,3)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4), alocont_nn3d(i,j,k,5), alocont_nn3d(i,j,k,6)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7), alocont_nn3d(i,j,k,8), alocont_nn3d(i,j,k,9)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,1)-mint3d(i-1,j-1,k-1), alocont_nn3d(i,j,k,2)-mint3d(i,j-1,k-1), alocont_nn3d(i,j,k,3)-mint3d(i+1,j-1,k-1)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,4)-mint3d(i-1,j,k-1), alocont_nn3d(i,j,k,5)-mint3d(i,j,k-1), alocont_nn3d(i,j,k,6)-mint3d(i+1,j,k-1)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,7)-mint3d(i-1,j+1,k-1), alocont_nn3d(i,j,k,8)-mint3d(i,j+1,k-1), alocont_nn3d(i,j,k,9)-mint3d(i+1,j+1,k-1)
!               write(*,*)
!   
!               write(*,*) 'level z(k)', z(k)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k), mint3d(i,j-1,k), mint3d(i+1,j-1,k)
!               write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k), mint3d(i,j,k), mint3d(i+1,j,k)
!               write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k), mint3d(i,j+1,k), mint3d(i+1,j+1,k)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,10), alocont_nn3d(i,j,k,11), alocont_nn3d(i,j,k,12)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,13), alocont_nn3d(i,j,k,14), alocont_nn3d(i,j,k,15)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,16), alocont_nn3d(i,j,k,17), alocont_nn3d(i,j,k,18)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,10)-mint3d(i-1,j-1,k), alocont_nn3d(i,j,k,11)- mint3d(i,j-1,k), alocont_nn3d(i,j,k,12)-mint3d(i+1,j-1,k)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,13)-mint3d(i-1,j,k), alocont_nn3d(i,j,k,14)-mint3d(i,j,k), alocont_nn3d(i,j,k,15)-mint3d(i+1,j,k)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,16)-mint3d(i-1,j+1,k), alocont_nn3d(i,j,k,17)-mint3d(i,j+1,k), alocont_nn3d(i,j,k,18)-mint3d(i+1,j+1,k)
!               write(*,*)
!   
!               write(*,*) 'level z(k+1)', z(k+1)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), mint3d(i-1,j-1,k+1), mint3d(i,j-1,k+1), mint3d(i+1,j-1,k+1)
!               write(*,'(f8.4, 3es25.14)') y(j), mint3d(i-1,j,k+1), mint3d(i,j,k+1), mint3d(i+1,j,k+1)
!               write(*,'(f8.4, 3es25.14)') y(j+1), mint3d(i-1,j+1,k+1), mint3d(i,j+1,k+1), mint3d(i+1,j+1,k+1)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,19), alocont_nn3d(i,j,k,20), alocont_nn3d(i,j,k,21)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,22), alocont_nn3d(i,j,k,23), alocont_nn3d(i,j,k,24)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,25), alocont_nn3d(i,j,k,26), alocont_nn3d(i,j,k,27)
!               write(*,*)
!               write(*,'(a8, 3f25.4)') 'x/y', x(i-1), x(i), x(i+1)
!               write(*,'(f8.4, 3es25.14)') y(j-1), alocont_nn3d(i,j,k,19)-mint3d(i-1,j-1,k+1), alocont_nn3d(i,j,k,20)-mint3d(i,j-1,k+1), alocont_nn3d(i,j,k,21)-mint3d(i+1,j-1,k+1)
!               write(*,'(f8.4, 3es25.14)') y(j), alocont_nn3d(i,j,k,22)-mint3d(i-1,j,k+1), alocont_nn3d(i,j,k,23)-mint3d(i,j,k+1), alocont_nn3d(i,j,k,24)-mint3d(i+1,j,k+1)
!               write(*,'(f8.4, 3es25.14)') y(j+1), alocont_nn3d(i,j,k,25)-mint3d(i-1,j+1,k+1), alocont_nn3d(i,j,k,26)-mint3d(i,j+1,k+1), alocont_nn3d(i,j,k,27)-mint3d(i+1,j+1,k+1)
!               write(*,*)
!               stop
!            endif
!         endif
!      enddo
!   enddo
!enddo
!stop
!stop 'calculate first order integrals in coeff_source3 again (deactivated only for debugging)'
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
s4b=4
!
!-----------------------------------------------------------------------
!
!setting start mean intensities and start deviations to zero
mint3d=0.d0
fcontx3d=0.d0
fconty3d=0.d0
fcontz3d=0.d0
kcontxx3d=zero
kcontyy3d=zero
kcontzz3d=zero
kcontxy3d=zero
kcontxz3d=zero
kcontyz3d=zero
eps3d=0.d0
normalization3d=0.d0
alocont_nn3d=0.d0
epsmaxc_arr=0.d0
!
!setting index for threshold
indx_threshold=5
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxc
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

!   indxx=24
!   indxy=25
!   indxz=30
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 1), imask3d(indxx-1,indxy-1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 2), imask3d(indxx,indxy-1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 3), imask3d(indxx+1,indxy-1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 4), imask3d(indxx-1,indxy,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 5), imask3d(indxx,indxy,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 6), imask3d(indxx+1,indxy,indxz-1)   
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 7), imask3d(indxx-1,indxy+1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 8), imask3d(indxx,indxy+1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 9), imask3d(indxx+1,indxy+1,indxz-1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 10), imask3d(indxx-1,indxy-1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 11), imask3d(indxx,indxy-1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 12), imask3d(indxx+1,indxy-1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 13), imask3d(indxx-1,indxy,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 14), imask3d(indxx,indxy,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 15), imask3d(indxx+1,indxy,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 16), imask3d(indxx-1,indxy+1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 17), imask3d(indxx,indxy+1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 18), imask3d(indxx+1,indxy+1,indxz)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 19), imask3d(indxx-1,indxy-1,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 20), imask3d(indxx,indxy-1,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 21), imask3d(indxx+1,indxy-1,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 22), imask3d(indxx-1,indxy,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 23), imask3d(indxx,indxy,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 24), imask3d(indxx+1,indxy,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 25), imask3d(indxx-1,indxy+1,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 26), imask3d(indxx,indxy+1,indxz+1)
!   write(*,*) alocont_nn3d(indxx, indxy, indxz, 27), imask3d(indxx+1,indxy+1,indxz+1)
!   stop
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
      if(epsmaxc_arr(i).ge.epsmaxc_arr(i-1).and. &
         epsmaxc_arr(i-1).le.epsmaxc_arr(i-2).and. &
         epsmaxc_arr(i-2).ge.epsmaxc_arr(i-3).and. &
         epsmaxc_arr(i-3).le.epsmaxc_arr(i-4).and. &
         epsmaxc_arr(i-4).ge.epsmaxc_arr(i-5)) then
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
         if(abs(eps_max).gt.1.d-14) then
            if(opt_ng_cont) call ng_expol3d(scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
            if(opt_ait_cont) call ait_expol3d(scont3d_ng,ndxmax,ndymax,ndzmax,scont3d)
         endif
      endif
!
   endif
!
enddo
!
ttot_it_sc=ttot_it_sc/float(i)
!
!-----------------------------------------------------------------------
!
!store everything in benchmark12-arrays
epsmaxc_sc = epsmaxc_arr
mint3d_sc = mint3d
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               rad = x(i)**2+y(j)**2+z(k)**2
               call get_angles_spc(x(i), y(j), z(k), theta, phi)
               sint = sin(theta)
               cost = cos(theta)
               sinp = sin(phi)
               cosp = cos(phi)
               fcontr3d_sc(i,j,k) = fcontx3d(i,j,k)*sint*cosp + fconty3d(i,j,k)*sint*sinp + fcontz3d(i,j,k)*cost
               fcontth3d_sc(i,j,k) = fcontx3d(i,j,k)*cost*cosp + fconty3d(i,j,k)*cost*sinp - fcontz3d(i,j,k)*sint
               fcontphi3d_sc(i,j,k) = -fcontx3d(i,j,k)*sinp + fconty3d(i,j,k)*cosp
               a11 = kcontxx3d(i,j,k)*sint*cosp + kcontxy3d(i,j,k)*sint*sinp + kcontxz3d(i,j,k)*cost
               a12 = kcontxx3d(i,j,k)*cost*cosp + kcontxy3d(i,j,k)*cost*sinp - kcontxz3d(i,j,k)*sint
               a13 = -kcontxx3d(i,j,k)*sinp + kcontxy3d(i,j,k)*cosp
               a21 = kcontxy3d(i,j,k)*sint*cosp + kcontyy3d(i,j,k)*sint*sinp + kcontyz3d(i,j,k)*cost
               a22 = kcontxy3d(i,j,k)*cost*cosp + kcontyy3d(i,j,k)*cost*sinp - kcontyz3d(i,j,k)*sint
               a23 = -kcontxy3d(i,j,k)*sinp + kcontyy3d(i,j,k)*cosp
               a31 = kcontxz3d(i,j,k)*sint*cosp + kcontyz3d(i,j,k)*sint*sinp + kcontzz3d(i,j,k)*cost
               a32 = kcontxz3d(i,j,k)*cost*cosp + kcontyz3d(i,j,k)*cost*sinp - kcontzz3d(i,j,k)*sint
               a33 = -kcontxz3d(i,j,k)*sinp + kcontyz3d(i,j,k)*cosp
               kcontrr3d_sc(i,j,k) = a11*sint*cosp + a21*sint*sinp + a31*cost
               kcontrth3d_sc(i,j,k) = a12*sint*cosp + a22*sint*sinp + a32*cost
               kcontrphi3d_sc(i,j,k) = a13*sint*cosp + a23*sint*sinp + a33*cost
               kcontthth3d_sc(i,j,k) = a12*cost*cosp + a22*cost*sinp - a32*sint
               kcontthphi3d_sc(i,j,k) = a13*cost*cosp + a23*cost*sinp - a33*sint
               kcontphiphi3d_sc(i,j,k) = -a13*sinp + a23*cosp
            case default
         end select
      enddo
   enddo
enddo

j=ndymax/2+1
i=ndxmax/2+1
do k=1, ndzmax
   write(*,*) mint3d(i,j,k), kcontzz3d(i,j,k), kcontrr3d_sc(i,j,k), kcontzz3d(i,j,k)/mint3d(i,j,k), kcontrr3d_sc(i,j,k)/mint3d(i,j,k)
enddo
!stop
!
!*************start finite volume method********************************
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
scont3d=0.d0
mint3d=0.d0
fcontx3d=0.d0
fconty3d=0.d0
fcontz3d=0.d0
kcontxx3d=zero
kcontyy3d=zero
kcontzz3d=zero
kcontxy3d=zero
kcontxz3d=zero
kcontyz3d=zero
eps3d=0.d0
normalization3d=0.d0
alocont_nn3d=0.d0
epsmaxc_arr=0.d0
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
ttot_it_fvm=ttot_it_fvm/float(i)
write(*,'(a35,es20.8)') 'total computation time sc', ttot_it_sc
write(*,'(a35,es20.8)') 'total computation time fvm', ttot_it_fvm
!
!-----------------------------------------------------------------------
!
!store everything in benchmark12-arrays
epsmaxc_fvm = epsmaxc_arr
mint3d_fvm = mint3d
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               rad = x(i)**2+y(j)**2+z(k)**2
               call get_angles_spc(x(i), y(j), z(k), theta, phi)
               sint = sin(theta)
               cost = cos(theta)
               sinp = sin(phi)
               cosp = cos(phi)               
               fcontr3d_fvm(i,j,k) = fcontx3d(i,j,k)*sint*cosp + fconty3d(i,j,k)*sint*sinp + fcontz3d(i,j,k)*cost
               fcontth3d_fvm(i,j,k) = fcontx3d(i,j,k)*cost*cosp + fconty3d(i,j,k)*cost*sinp - fcontz3d(i,j,k)*sint
               fcontphi3d_fvm(i,j,k) = -fcontx3d(i,j,k)*sinp + fconty3d(i,j,k)*cosp
               a11 = kcontxx3d(i,j,k)*sint*cosp + kcontxy3d(i,j,k)*sint*sinp + kcontxz3d(i,j,k)*cost
               a12 = kcontxx3d(i,j,k)*cost*cosp + kcontxy3d(i,j,k)*cost*sinp - kcontxz3d(i,j,k)*sint
               a13 = -kcontxx3d(i,j,k)*sinp + kcontxy3d(i,j,k)*cosp
               a21 = kcontxy3d(i,j,k)*sint*cosp + kcontyy3d(i,j,k)*sint*sinp + kcontyz3d(i,j,k)*cost
               a22 = kcontxy3d(i,j,k)*cost*cosp + kcontyy3d(i,j,k)*cost*sinp - kcontyz3d(i,j,k)*sint
               a23 = -kcontxy3d(i,j,k)*sinp + kcontyy3d(i,j,k)*cosp
               a31 = kcontxz3d(i,j,k)*sint*cosp + kcontyz3d(i,j,k)*sint*sinp + kcontzz3d(i,j,k)*cost
               a32 = kcontxz3d(i,j,k)*cost*cosp + kcontyz3d(i,j,k)*cost*sinp - kcontzz3d(i,j,k)*sint
               a33 = -kcontxz3d(i,j,k)*sinp + kcontyz3d(i,j,k)*cosp
               kcontrr3d_fvm(i,j,k) = a11*sint*cosp + a21*sint*sinp + a31*cost
               kcontrth3d_fvm(i,j,k) = a12*sint*cosp + a22*sint*sinp + a32*cost
               kcontrphi3d_fvm(i,j,k) = a13*sint*cosp + a23*sint*sinp + a33*cost
               kcontthth3d_fvm(i,j,k) = a12*cost*cosp + a22*cost*sinp - a32*sint
               kcontthphi3d_fvm(i,j,k) = a13*cost*cosp + a23*cost*sinp - a33*sint
               kcontphiphi3d_fvm(i,j,k) = -a13*sinp + a23*cosp
            case default
         end select
      enddo
   enddo
enddo
!
end subroutine benchmark12_solution
