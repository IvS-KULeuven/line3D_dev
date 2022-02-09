subroutine benchmark05_solution
!
!------calculating 2d continuum radiative transfer for spherically------
!----symmetric problems, with the theoretical solution given from-------
!---------------------------JOs 1d programs-----------------------------
!
use prog_type
use fund_const, only: pi, cgs_clight, cgs_planck, cgs_kb
use dimecr, only: n1d_cr, r1d_cr, norm1d_cr, alocont1d_diag_cr, epsc1d_cr, scont1d_cr
use dime3d, only: t3d, mint3d, int3d, opac3d, scont3d, ndxmax, ndymax, ndzmax, x, y, z, &
                  imask_innreg3d, alocont_nn3d
use angles, only: dim_mu, weight_mu, nodes_mu
use bcondition, only: xic1, xic2
use mod_benchmark, only: mint1d_sc, mint1d_fvm, mint1d_joray, mint1d_jomom, epsmaxc_sc, epsmaxc_fvm, &
                         t1d_cr, t1d_jo, opac1d_cr, opac1d_jo
use freq, only: xnue0
use iter, only: itmaxc, devmaxc
use ng_extra, only: ng_const
use options, only: opt_ng_cont, opt_ait_cont, input_mod_dim
use params_input, only: teff, rstar, beta, vmin, vmax, hei, yhe, xmloss, kcont
use params_stellar, only: sr
use mod_debug, only: indx1, indx2
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, muindx, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: nr_jo
real(dp) :: eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xnue0_jo, xmloss_jo, &
            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
real(dp) :: s1, s2, s3, s4, eps_max, fdum, dummy1, dummy2, rad, trad
real(dp) :: ts, te, eps_cont
!
! ... local arrays
real(dp), dimension(:), allocatable :: bnue1d_cr
real(dp), dimension(:,:), allocatable :: scont1d_ng
real(dp), dimension(:), allocatable :: r_jo, t_jo, rho_jo, opac_jo, tau_jo, mintray_jo, mintmom_jo, scontray_jo, scontmom_jo
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
if(input_mod_dim.ne.1) stop 'error in benchmark05_solution: input dimension ne 1'
!
check_fname=.false.
inquire(file='./models/jo/continuum_jo.dat', exist=check_fname)
if(.not.check_fname) stop 'error in benchmark05: 1d solution file does not exist'
!
!read 1d solution
open(1, file='./models/jo/continuum_jo.dat', form='formatted')
   read(1,*)
   read(1,'(i5, 13es20.8)') nr_jo, eps_cont_jo, kcont_jo, corrfc_jo, teff_jo, xic1_jo, lambda_jo, xmloss_jo, &
                            vmin_jo, vmax_jo, beta_jo, yhe_jo, hei_jo, sr_jo
   read(1,*)
   read(1,*)
   allocate(r_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(t_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(rho_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(opac_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(tau_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(mintray_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(mintmom_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(scontray_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
   allocate(scontmom_jo(nr_jo), stat=err)
      if(err.ne.0) stop 'allocation error in benchmark05'
!
   do i=1, nr_jo
      read(1,'(11es20.8)') r_jo(i), t_jo(i), rho_jo(i), opac_jo(i), tau_jo(i), mintray_jo(i), mintmom_jo(i), &
                           scontray_jo(i), scontmom_jo(i), fdum, fdum
   enddo
!
close(1)
!
write(*,'(3a20)') 'quantity', 'own input', 'JOs input'
write(*,'(a20, 2es20.8)') 'teff', teff, teff_jo
if(abs(teff-teff_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'mdot', xmloss, xmloss_jo
if(abs(xmloss-xmloss_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmin', vmin, vmin_jo
if(abs(vmin*1.d5-vmin_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'vmax', vmax, vmax_jo
if(abs(vmax*1.d5-vmax_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'beta', beta, beta_jo
if(abs(beta-beta_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'yhe', yhe, yhe_jo
if(abs(yhe-yhe_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'hei', hei, hei_jo
if(abs(hei-hei_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'rstar', sr, sr_jo
if(abs(sr-sr_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'eps_cont', eps_cont, eps_cont_jo
if(abs(eps_cont-eps_cont_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'kcont', kcont, kcont_jo
if(abs(kcont-kcont_jo).gt.1.d-10) stop 'error in benchmark05: input data not consistent'
xnue0_jo=cgs_clight/lambda_jo/1.d-8
write(*,'(a20, 2es20.8)') 'xnue', xnue0, xnue0_jo
if(abs(xnue0-xnue0_jo)/xnue0.gt.1.d-8) stop 'error in benchmark05: input data not consistent'
write(*,'(a20, 2es20.8)') 'xic1', xic1, xic1_jo
if(abs(xic1-xic1_jo).gt.1.d-5) then
   write(*,*) 'error in benchmark05: input data not consistent'
   trad = log(2.d0*cgs_planck*xnue0**3/cgs_clight**2/xic1_jo + 1.d0)
   trad = cgs_planck*xnue0/trad/cgs_kb
   write(*,*) 'use as radiation temperature on boundary:', trad
   stop
endif
!if(abs(xic1-xic1_jo).gt.1.d-5) stop 'error in benchmark05: input data not consistent'
write(*,*)
!
!interpolate everything onto own grid
allocate(mint1d_joray(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(mint1d_jomom(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(opac1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(opac1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(t1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(t1d_jo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
!
do i=1, n1d_cr
   t1d_cr(i) = t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   opac1d_cr(i) = opac3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i)
   rad = z(ndzmax-n1d_cr-1+i)
   call find_index(rad, r_jo, nr_jo, iim2, iim1, ii, iip1)
   opac1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), opac_jo(iim1), opac_jo(ii), rad)
   t1d_jo(i) = interpol_yp(r_jo(iim1), r_jo(ii), t_jo(iim1), t_jo(ii), rad)
   mint1d_joray(i) = interpol_yp(r_jo(iim1), r_jo(ii), mintray_jo(iim1), mintray_jo(ii), rad)
   mint1d_jomom(i) = interpol_yp(r_jo(iim1), r_jo(ii), mintmom_jo(iim1), mintmom_jo(ii), rad)
enddo
!
!in cgs
opac1d_cr=opac1d_cr/sr
!
!-----------------------------------------------------------------------
!
allocate(mint1d_sc(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(mint1d_fvm(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
!
allocate(epsmaxc_sc(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(epsmaxc_fvm(itmaxc), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'

allocate(bnue1d_cr(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
allocate(scont1d_ng(4,n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark05'
!
epsmaxc_sc=0.d0
epsmaxc_fvm=0.d0
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
!-----------------------------------------------------------------------
!
!start values for continuum source function
!scont3d=xic1!0.d0!xic1
mint3d=0.d0
!
!calculating planck function along z direction
do i=1, n1d_cr
   bnue1d_cr(i) = bnue(xnue0, t3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+i))
enddo
!
!setting start mean intensities and start deviations to zero
mint1d_sc=0.d0
mint1d_fvm=0.d0
!
epsc1d_cr=0.d0
alocont1d_diag_cr=0.d0
norm1d_cr=0.d0
scont1d_cr=0.d0
!
!test for alo
!scont1d_cr(1)=1.d0
!do j=1, ndxmax
!   do k=1, ndymax
!      do l=1, ndzmax
!         if(imask_innreg3d(j,k,l).ne.0) then
!            rad=sqrt(x(j)**2 + y(k)**2 + z(l)**2)
!            call find_index(rad, r1d_cr, n1d_cr, iim2, iim1, ii, iip1)
!            scont3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), scont1d_cr(iim1), scont1d_cr(ii), rad)
!         endif
!      enddo
!   enddo
!enddo
!xic1=0.d0

!do i=1, ndxmax
!   write(*,*) x(i), imask_innreg3d(i,ndymax/2+1,ndzmax/2+1), scont3d(i,ndymax/2+1,ndzmax/2+1), xic1
!enddo
!stop

!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxc
!
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j(sc)', 'scont(sc)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mint1d_sc(j), scont1d_cr(j), epsc1d_cr(j), norm1d_cr(j), alocont1d_diag_cr(j)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   epsc1d_cr=mint1d_sc
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
!
!
!***debug start
!   write(fname_debug,'(i3.3)') i
!   open(1, file='TRASH/'//fname_debug, form='formatted')
!   indx1=0
!   indx2=0
!***debug end

   norm1d_cr=0.d0
   mint1d_sc=0.d0
   alocont1d_diag_cr=0.d0

   do muindx=1, dim_mu
      call fsc_cont2d(muindx)
      do j=1, n1d_cr
!on positive z-axis
         mint1d_sc(j)=mint1d_sc(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
         alocont1d_diag_cr(j)=alocont1d_diag_cr(j) + alocont_nn3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j,14)*weight_mu(muindx)
!on negative z-axis
!         mint1d_sc(j)=mint1d_sc(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
!         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
      enddo
!      write(*,'(6es20.8)') nodes_mu(muindx), alocont_nn3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr, 14), int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr), weight_mu(muindx), alocont1d_diag_cr(1), mint1d_sc(1)
!      write(*,'(4es20.8)') int3d(ndxmax/2-1, ndymax/2+1, ndzmax-n1d_cr), int3d(ndxmax/2, ndymax/2+1, ndzmax-n1d_cr), int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr), int3d(ndxmax/2+2, ndymax/2+1, ndzmax-n1d_cr)
!      write(*,*)
!      write(1,'(6es20.8)') nodes_mu(muindx), int3d(ndxmax/2+1,ndymax/2+1,64), int3d(ndxmax/2+1,ndymax/2+1,63), int3d(ndxmax/2+1,ndymax/2+1,62), int3d(ndxmax/2+1,ndymax/2+1,61), xic1
   enddo

!***debug start
!   close(1)
!   write(*,*)
!   write(*,*) mint1d_sc
!    write(*,*) 'rays hit the star', indx1
!    write(*,*) 'rays that didnt hit the star', indx2
!    write(*,*) 'number directions', dim_mu
!***debug end
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsc1d_cr, mint1d_sc, n1d_cr, eps_max)
   epsmaxc_sc(i)=eps_max
!
   if(abs(eps_max).lt.devmaxc) then
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
   dummy2=1.d0-eps_cont
   do j=1, n1d_cr
      dummy1=1.d0-(1.d0-eps_cont)*alocont1d_diag_cr(j)
      scont1d_cr(j)=(dummy2/dummy1)*mint1d_sc(j) - &
             (dummy2/dummy1)*alocont1d_diag_cr(j)*scont3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j) + &
              eps_cont*bnue1d_cr(j)/dummy1
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
   do j=1, ndxmax
      do k=1, ndymax
         do l=1, ndzmax
            if(imask_innreg3d(j,k,l).eq.0) then
               rad=sqrt(x(j)**2 + y(k)**2 + z(l)**2)
               call find_index(rad, r1d_cr, n1d_cr, iim2, iim1, ii, iip1)
               scont3d(j,k,l) = interpol_yp(r1d_cr(iim1), r1d_cr(ii), scont1d_cr(iim1), scont1d_cr(ii), rad)
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
s1=1
s2=2
s3=3
s4=4
!
!-----------------------------------------------------------------------
!
!start values for continuum source function
scont3d=0.d0
mint3d=0.d0
!
epsc1d_cr=0.d0
alocont1d_diag_cr=0.d0
norm1d_cr=0.d0
scont1d_cr=0.d0
!
!---------------------start iteration scheme----------------------------
!
do i=1, itmaxc
!
!-----------------------------------------------------------------------
!
   write(*,fmt='(a5, 6(a20))') '#', 'radius', 'j(fvm)', 'scont(fvm)', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, n1d_cr
      write(*,fmt='(i5, 6(e20.8))') j, r1d_cr(j), mint1d_fvm(j), scont1d_cr(j), epsc1d_cr(j), norm1d_cr(j), alocont1d_diag_cr(j)
   end do
   write(*,*)
!
!-----------------------------------------------------------------------
!
   epsc1d_cr=mint1d_fvm
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities on cr--------------'
   write(*,*) 'step', i
   write(*,*)
!
   norm1d_cr=0.d0
   mint1d_fvm=0.d0
   alocont1d_diag_cr=0.d0
   do muindx=1, dim_mu
      call ffvm_cont2d(muindx)
      do j=1, n1d_cr
!on positive z-axis
         mint1d_fvm(j)=mint1d_fvm(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
         alocont1d_diag_cr(j)=alocont1d_diag_cr(j) + alocont_nn3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j,14)*weight_mu(muindx)
!on negative z-axis
!         mint1d_fvm(j)=mint1d_fvm(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
!         norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
      enddo
!      write(*,'(6es20.8)') nodes_mu(muindx), int3d(ndxmax/2+1,ndymax/2+1,64), int3d(ndxmax/2+1,ndymax/2+1,63), int3d(ndxmax/2+1,ndymax/2+1,62), int3d(ndxmax/2+1,ndymax/2+1,61), xic1
   enddo

!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)

   call calc_dev(epsc1d_cr, mint1d_fvm, n1d_cr, eps_max)
   epsmaxc_fvm(i)=eps_max
!
   if(abs(eps_max).lt.devmaxc) then
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
   dummy2=1.d0-eps_cont
   do j=1, n1d_cr
      dummy1=1.d0-(1.d0-eps_cont)*alocont1d_diag_cr(j)
      scont1d_cr(j)=(dummy2/dummy1)*mint1d_fvm(j) - &
             (dummy2/dummy1)*alocont1d_diag_cr(j)*scont3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr-1+j) + &
              eps_cont*bnue1d_cr(j)/dummy1
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
enddo
!
!
end subroutine benchmark05_solution
