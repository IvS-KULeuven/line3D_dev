subroutine benchmark08_solution
!
!------calculating 2d continuum radiative transfer for spherically------
!----symmetric problems, with the theoretical solution given from-------
!---------------------------JOs 1d programs-----------------------------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb
use dimecr, only: n1d_cr, r1d_cr, norm1d_cr, epsl1d_cr, sline1d_cr, aloline1d_diag_cr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, t3d, mintbar3d, int3d, eps_cont3d, &
                  opac3d, opalbar3d, sline3d, ssobo3d, aloline_nn3d, aloline_on_nn3d, &
                  imask_innreg3d, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, normalization3d
use angles, only: dim_mu, weight_mu, nodes_mu
use bcondition, only: xic1, xic2
use mod_benchmark, only: epsmaxl_sc, epsmaxl_fvm, opalbar1d_cr, opalbar1d_jo, ssobo1d_jo, sline1d_jo, &
                         t1d_cr, t1d_jo, opac1d_cr, opac1d_jo, sline1d_sc, sline1d_fvm, ssobo1d_cr, &
                         mintbar1d_sc, mintbar1d_fvm, velr1d_cr
use freq, only: xnue0, nxobs, nodes_xobs, weight_xobs
use params_input, only: eps_line
use iter, only: itmaxl, devmaxl
use ng_extra, only: ng_const
use options, only: opt_ng_line, opt_ait_line, input_mod_dim
use params_input, only: teff, rstar, beta, vmin, vmax, hei, yhe, xmloss, kcont, eps_line, kline, &
                        kappa0, alpha, vth_fiducial
use params_stellar, only: sr
use options, only: opt_opal
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, xobsindx, muindx, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: s1, s2, s3, s4
integer(i4b) :: nr_jo, opt_opal_jo
real(dp) :: eps_cont, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xnue0_jo, xmloss_jo, vth_jo, &
            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo, eps_line_jo, kline_jo, kappa0_jo, alpha_jo
real(dp) :: eps_max, fdum, dummy1, dummy2, rad
real(dp) :: velx_p, vely_p, velz_p, vth_p, vel_p
real(dp) :: n_x, n_y, n_z, xobs, phinorm
real(dp) :: trad
!
! ... local arrays
real(dp), dimension(:), allocatable :: bnue1d_cr
real(dp), dimension(:,:), allocatable :: sline1d_ng
real(dp), dimension(:), allocatable :: r_jo, t_jo, rho_jo, opac_jo, tau_jo, opalbar_jo, &
                                       sline_jo, ssoboc_jo, ssobo_jo, scont_jo, velr_jo
!
! ... local logicals
logical :: check_fname
!
! ... local functions
real(dp) :: bnue
!
! ... for debugging
character(len=100) :: fname
!
call calc_mean3d(ndxmax, ndymax, ndzmax, imask_totreg3d, eps_cont3d, eps_cont)
!
!---------------------calculate sobolev solution------------------------
!
call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
!
call calc_startval_line
!
!---------------------check if input model is consistent----------------
!
inquire(file='./models/jo/line_jo.dat', exist=check_fname)
if(.not.check_fname) stop 'error in read_t1d: input file does not exist'
write(*,*) 'reading temperature structure from ./models/jo/line_jo.dat'
write(*,*)

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
if(abs(teff-teff_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'mdot', xmloss, xmloss_jo
if(abs(xmloss-xmloss_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmin', vmin, vmin_jo
if(abs(vmin*1.d5-vmin_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmax', vmax, vmax_jo
if(abs(vmax*1.d5-vmax_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'beta', beta, beta_jo
if(abs(beta-beta_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'yhe', yhe, yhe_jo
if(abs(yhe-yhe_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'hei', hei, hei_jo
if(abs(hei-hei_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'rstar', sr, sr_jo
if(abs(sr-sr_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'eps_cont', eps_cont, eps_cont_jo
if(abs(eps_cont-eps_cont_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'kcont', kcont, kcont_jo
if(abs(kcont-kcont_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
xnue0_jo=cgs_clight/lambda_jo/1.d-8
write(*,'(a20, 2es20.8)') 'xnue', xnue0, xnue0_jo
if(abs(xnue0-xnue0_jo)/xnue0.gt.1.d-6) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2es20.8)') 'xic1', xic1, xic1_jo
if(abs(xic1-xic1_jo).gt.1.d-5) then
   write(*,*) 'error in benchmark08: input data not consistent'
   trad = log(2.d0*cgs_planck*xnue0**3/cgs_clight**2/xic1_jo + 1.d0)
   trad = cgs_planck*xnue0/trad/cgs_kb
   write(*,*) 'use as radiation temperature on boundary:', trad
   stop
endif
write(*,'(a20, 2es20.8)') 'eps_line', eps_line, eps_line_jo
if(abs(eps_line-eps_line_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
write(*,'(a20, 2i20)') 'opt_opal', opt_opal, opt_opal_jo
if(opt_opal.ne.opt_opal_jo) stop 'error in benchmark08: input data not consistent'
if(opt_opal.eq.0) then
   write(*,'(a20, 2es20.8)') 'kline', kline, kline_jo
   if(abs(kline-kline_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
elseif(opt_opal.eq.1) then
   write(*,'(a20, 2es20.8)') 'kappa0', kappa0, kappa0_jo
   if(abs(kappa0-kappa0_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
   write(*,'(a20, 2es20.8)') 'alpha', alpha, alpha_jo
   if(abs(alpha-alpha_jo).gt.1.d-10) stop 'error in benchmark08: input data not consistent'
endif

write(*,*)
!
!interpolate everything onto own grid
allocate(sline1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(ssobo1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(ssobo1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(opac1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(opalbar1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(velr1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(opac1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(opalbar1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(t1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(t1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
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
enddo
!
!in cgs
opac1d_cr=opac1d_cr/sr
!
!-----------------------------------------------------------------------
!
allocate(mintbar1d_sc(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(mintbar1d_fvm(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
!
allocate(sline1d_sc(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(sline1d_fvm(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
!
allocate(epsmaxl_sc(itmaxl), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(epsmaxl_fvm(itmaxl), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
!
allocate(bnue1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
allocate(sline1d_ng(4,n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark08'
!
epsmaxl_sc=0.d0
epsmaxl_fvm=0.d0
!
!*************start short characteristics method************************
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
!
!
!-----------------------------------------------------------------------
!
!calculating planck function along z direction
do i=1, n1d_cr
   bnue1d_cr(i) = bnue(xnue0, t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i))
enddo
!
!setting start mean intensities and start deviations to zero
mintbar1d_sc=0.d0
mintbar1d_fvm=0.d0
!
epsl1d_cr=0.d0
norm1d_cr=0.d0
aloline1d_diag_cr=0.d0
mintbar1d_sc=0.d0
!
do j=1, n1d_cr
   sline1d_sc(j)=sline3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j)
enddo
!
!------------------testing alo------------------------------------------
!
!xic1=0.d0
!sline3d=0.d0
!sline3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr+1)=1.d0

!write(*,*) opalbar3d(ndxmax/2+1, ndymax/2+1,:)
!stop
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxl
!
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j_bar(sc)', 'sline(sc)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mintbar1d_sc(j), sline1d_sc(j), epsl1d_cr(j), norm1d_cr(j), aloline1d_diag_cr(j)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   epsl1d_cr=mintbar1d_sc
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
!
   norm1d_cr=0.d0
   mintbar1d_sc=0.d0
   mintbar3d=0.d0
   normalization3d=0.d0
   aloline_nn3d=0.d0
   do xobsindx=1, nxobs   
      write(*,'(a55,i4, a1, i4)') 'calculating all angles for (freq-point/nxobs)', xobsindx, '/', nxobs
      xobs=nodes_xobs(xobsindx)

      do muindx=1, dim_mu
         n_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
         n_y=0.d0
         n_z=nodes_mu(muindx)
!
         int3d=0.d0
         call fsc_line2d(muindx,xobsindx)
!
         do j=1, n1d_cr
            mintbar1d_sc(j) = mintbar3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+j)
            norm1d_cr(j) =  normalization3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+j)
            aloline1d_diag_cr(j) = aloline_nn3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j,14)
!            write(*,'(9es20.8)') xobs-vel_p, n_x, n_z, xobs, velx_p, vely_p, velz_p, vel_p, int3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j)
         enddo
!
      enddo
   enddo
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsl1d_cr, mintbar1d_sc, n1d_cr, eps_max)
   epsmaxl_sc(i)=eps_max
!
   if(abs(eps_max).lt.devmaxl) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      exit
   endif
!
!-------calculating alo-corrected source-functions on central ray-------
!
   write(*,*) '--------------calculating new iterate for source function (alo)----------------'
   write(*,*)
!
   dummy2=1.d0-eps_line
   do j=1, n1d_cr
      dummy1=1.d0-(1.d0-eps_line)*aloline1d_diag_cr(j)
      sline1d_sc(j)=(dummy2/dummy1)*mintbar1d_sc(j) - &
             (dummy2/dummy1)*aloline1d_diag_cr(j)*sline3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j) + &
              eps_line*bnue1d_cr(j)/dummy1
   enddo

!
!*****************************DEBUG START*******************************
!   write(fname,'(i3.3)') i
!   fname='outputFILES_DEBUG/'//trim(fname)
!   write(*,*) 'output to ', trim(fname)
!   open(1, file=trim(fname), form='formatted')
!      do j=1, n1d_cr
!         write(1,'(5es20.8)') z(ndzmax-n1d_cr-1+j), velz3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j)*vth_fiducial/vmax/1.d5, sline1d_sc(j), ssobo3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j), xic1
!      enddo
!   close(1)
!********************************DEBUG END******************************
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
         sline1d_ng(1,:)=sline1d_sc
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         sline1d_ng(2,:)=sline1d_sc
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         sline1d_ng(3,:)=sline1d_sc
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         sline1d_ng(4,:)=sline1d_sc
         s4=s4+ng_const
         if(opt_ng_line) call ng_expol1d(sline1d_ng, n1d_cr)
         if(opt_ait_line) call ait_expol1d(sline1d_ng, n1d_cr)
         sline1d_sc=sline1d_ng(1,:)
      endif
!
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
               sline3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), sline1d_sc(iim1), sline1d_sc(ii), rad)
            endif
         enddo
      enddo
   enddo
enddo
!
!*************start finite volume method********************************
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
!
s1=1
s2=2
s3=3
s4=4
!
!-----------------------------------------------------------------------
!
!start values for continuum source function
call calc_startval_line
do j=1, n1d_cr
   sline1d_fvm(j)=sline3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j)
enddo
!
mintbar3d=0.d0
!
epsl1d_cr=0.d0
aloline1d_diag_cr=0.d0
norm1d_cr=0.d0
mintbar1d_fvm=0.d0
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxl
!
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j_bar(fvm)', 'sline(fvm)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mintbar1d_fvm(j), sline1d_fvm(j), epsl1d_cr(j), norm1d_cr(j), aloline1d_diag_cr(j)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   epsl1d_cr=mintbar1d_fvm
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
!
   norm1d_cr=0.d0
   mintbar1d_fvm=0.d0
   mintbar3d=0.d0
   normalization3d=0.d0
   aloline_nn3d=0.d0
   do xobsindx=1, nxobs
      write(*,'(a55,i4, a1, i4)') 'calculating all angles for (freq-point/nxobs)', xobsindx, '/', nxobs
      do muindx=1, dim_mu
         call ffvm_line2d(muindx,xobsindx)
         do j=1, n1d_cr
            mintbar1d_fvm(j) = mintbar3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+j)
            norm1d_cr(j) =  normalization3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+j)
            aloline1d_diag_cr(j) = aloline_nn3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j,14)
         enddo
      enddo
   enddo
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsl1d_cr, mintbar1d_fvm, n1d_cr, eps_max)
   epsmaxl_fvm(i)=eps_max
!
   if(abs(eps_max).lt.devmaxl) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      exit
   endif
!
!-------calculating alo-corrected source-functions on central ray-------
!
   write(*,*) '--------------calculating new iterate for source function (alo)----------------'
   write(*,*)
!
   dummy2=1.d0-eps_line
   do j=1, n1d_cr
      dummy1=1.d0-(1.d0-eps_line)*aloline1d_diag_cr(j)
      sline1d_fvm(j)=(dummy2/dummy1)*mintbar1d_fvm(j) - &
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
         sline1d_ng(1,:)=sline1d_fvm
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         sline1d_ng(2,:)=sline1d_fvm
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         sline1d_ng(3,:)=sline1d_fvm
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         sline1d_ng(4,:)=sline1d_fvm
         s4=s4+ng_const
         if(opt_ng_line) call ng_expol1d(sline1d_ng, n1d_cr)
         if(opt_ait_line) call ait_expol1d(sline1d_ng, n1d_cr)
         sline1d_fvm=sline1d_ng(1,:)
      endif
!
   endif


!
!*****************************DEBUG START*******************************
!   write(fname,'(i3.3)') i
!   fname='outputFILES_DEBUG/'//trim(fname)
!   write(*,*) 'output to ', trim(fname)
!   open(1, file=trim(fname), form='formatted')
!      do j=1, n1d_cr
!         write(1,'(4es20.8)') z(ndzmax-n1d_cr-1+j), sline1d_fvm(j), ssobo3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j), xic1
!      enddo
!   close(1)
!********************************DEBUG END******************************

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
               sline3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), sline1d_fvm(iim1), sline1d_fvm(ii), rad)
            endif
         enddo
      enddo
   enddo
enddo
!
!
end subroutine benchmark08_solution
