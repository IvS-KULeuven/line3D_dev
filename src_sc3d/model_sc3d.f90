subroutine setup_mod1d
!
use prog_type
use fund_const, only: rsu, zero
use modext, only: nr_modext, r_modext, rho_modext1d, t_modext1d, &
                  velr_modext1d, vth_modext1d, eps_cont_modext1d
use params_input, only: rmax, vth_fiducial, tmin
use params_stellar, only: sr, smajorax_a, smajorax_b, smajorax_c
use inf_reg, only: rlim, rmin
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, &
                  imask_innreg3d, imask_bpoint3d, imask3d, &
                  rho3d, opac3d, opalbar3d, t3d, velx3d, vely3d, velz3d, &
                  vth3d, eps_cont3d
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2, iim1, ii, iip1
real(dp) :: rad, rad_log, rho, velr, vth, vth_mean, temp, eps_cont
!
! ... local logicals
logical :: linfo_min, linfo_max, linfo_boundary
!
! ... local functions
real(dp) :: getr
!
!-----------------------------------------------------------------------
!
!read 1d model
call read_mod1d
!
!do i=nr_modext, 1, -1
!   write(*,*) r_modext(i)/sr, rho_modext1d(i), t_modext1d(i)
!enddo
!
!convert to different units
r_modext=r_modext/sr
!
!check if minimum radius is okay
if(abs(rmin-r_modext(1)).gt.1.d-14) then
   write(*,*) 'error in setup_mod1d: r_min does not match'
   write(*,*) 'set rstar to ', r_modext(1)*sr/rsu
   stop
endif
if(r_modext(nr_modext).lt.rlim) then
   write(*,*) 'error in setup_mod1d: r_lim does not match'
   write(*,*) 'set rlim to ', r_modext(nr_modext)
   stop
endif
!
!--------calculate calculation volume and perform interpolations--------
!
!calculate a mean thermal velocity (required outside information region)
vth_mean=0.d0
do i=1, nr_modext
   vth_mean= vth_mean+vth_modext1d(i)
enddo
vth_mean=vth_mean/nr_modext
!
!interpolation in log-space
rho_modext1d=log10(rho_modext1d)
vth_modext1d=log10(vth_modext1d)
velr_modext1d=log10(velr_modext1d)
t_modext1d=log10(t_modext1d)
r_modext=log10(r_modext)
!
!default values of 3d-grid
imask_totreg3d=0
imask_innreg3d=0
imask_bpoint3d=0
opac3d=0.d0
rho3d=0.d0
opalbar3d=0.d0
t3d=0.d0
velx3d=0.d0
vely3d=0.d0
velz3d=0.d0
vth3d=vth_mean
eps_cont3d = zero
!
!inside calculation volume
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
!         call info_region(x(i), y(j), z(k), rmin, rlim, linfo_min, linfo_max, linfo_boundary)
         call info_region_rot(x(i), y(j), z(k), smajorax_a, smajorax_b, smajorax_c, rlim, linfo_min, linfo_max, linfo_boundary)
         if(linfo_boundary) imask_bpoint3d(i,j,k)=1
         if(.not.linfo_min) imask_innreg3d(i,j,k)=1
         if(linfo_min.and.linfo_max) then
            imask_totreg3d(i,j,k)=1
            rad=getr(x(i), y(j), z(k))
            rad_log=log10(rad)
            call find_index(rad_log, r_modext, nr_modext, iim2, iim1, ii, iip1)
            rho=interpol_yp(r_modext(iim1), r_modext(ii), rho_modext1d(iim1), rho_modext1d(ii), rad_log)
            rho=10.d0**rho
            rho3d(i,j,k)=rho
            velr=interpol_yp(r_modext(iim1), r_modext(ii), velr_modext1d(iim1), velr_modext1d(ii), rad_log)
            velr=10.d0**velr
            velx3d(i,j,k) = velr*x(i)/rad/vth_fiducial
            vely3d(i,j,k) = velr*y(j)/rad/vth_fiducial
            velz3d(i,j,k) = velr*z(k)/rad/vth_fiducial
            vth=interpol_yp(r_modext(iim1), r_modext(ii), vth_modext1d(iim1), vth_modext1d(ii), rad_log)
            vth3d(i,j,k) = 10.d0**vth
            temp=interpol_yp(r_modext(iim1), r_modext(ii), t_modext1d(iim1), t_modext1d(ii), rad_log)
            t3d(i,j,k) = 10.d0**temp
            eps_cont = interpol_yp(10.d0**r_modext(iim1), 10.d0**r_modext(ii), eps_cont_modext1d(iim1), eps_cont_modext1d(ii), rad)
            eps_cont3d(i,j,k) = eps_cont
         endif
      enddo
   enddo
enddo
!
!
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         rad=getr(x(i), y(j), z(k))
         if(imask_totreg3d(i,j,k).eq.0) then
!         if(imask_totreg3d(i,j,k).eq.0.and.&
!            imask_innreg3d(i,j,k).eq.0.and.rad.le.rmax) then
            rad_log=log10(rad)
            call find_index(rad_log, r_modext, nr_modext, iim2, iim1, ii, iip1)
            rho=interpol_yp(r_modext(iim1), r_modext(ii), rho_modext1d(iim1), rho_modext1d(ii), rad_log)
            rho=10.d0**rho
            if(rho.lt.0.d0) stop 'error in setup_mod1d: density extrapolated below 0'
            rho3d(i,j,k)=rho
            velr=interpol_yp(r_modext(iim1), r_modext(ii), velr_modext1d(iim1), velr_modext1d(ii), rad_log)
            velr=10.d0**velr
            if(velr.lt.0.d0) stop 'error in setup_mod1d: radial velocity extrapolated below 0'
            velx3d(i,j,k) = velr*x(i)/rad/vth_fiducial
            vely3d(i,j,k) = velr*y(j)/rad/vth_fiducial
            velz3d(i,j,k) = velr*z(k)/rad/vth_fiducial
            temp=interpol_yp(r_modext(iim1), r_modext(ii), t_modext1d(iim1), t_modext1d(ii), rad_log)
            temp=max(temp,log10(tmin))
!            write(*,'(10es20.8)') r_modext(iim1), r_modext(ii), t_modext1d(iim1), t_modext1d(ii), rad_log, temp
            if(temp.lt.log10(tmin)) stop 'error in setup_mod1d: temperature extrapolated below tmin'
            t3d(i,j,k) = 10.d0**temp
            eps_cont = interpol_yp(10.d0**r_modext(iim1), 10.d0**r_modext(ii), eps_cont_modext1d(iim1), eps_cont_modext1d(ii), rad)
            eps_cont3d(i,j,k) = eps_cont
         endif
      enddo
   enddo
enddo
!
!write(*,*) velr_modext1d

!deallocate external model structure
deallocate(r_modext)
deallocate(rho_modext1d)
deallocate(velr_modext1d)
deallocate(t_modext1d)
deallocate(vth_modext1d)
deallocate(eps_cont_modext1d)
!
!j=ndymax/2+1
!k=ndzmax/2+1
!do i=2, ndxmax-1
!   write(*,'(10es20.8)') x(i), velx3d(i,j,k)*vth_fiducial, velx3d(i,j,k)-velx3d(i-1,j,k), velx3d(i+1,j,k)-velx3d(i,j,k)
!enddo
!stop 'go on in setup_mod1d'
!
end subroutine setup_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_mod2d
!
use prog_type
use fund_const
use modext, only: nr_modext, ntheta_modext, r_modext, theta_modext, &
                  rho_modext2d, t_modext2d, velr_modext2d, velth_modext2d, &
                  velphi_modext2d, vth_modext2d, eps_cont_modext2d
use params_input, only: rmax, vth_fiducial
use params_stellar, only: sr, smajorax_a, smajorax_b, smajorax_c
use inf_reg, only: rlim, rmin
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, &
                  rho3d, opac3d, opalbar3d, t3d, velx3d, vely3d, velz3d, &
                  vth3d, imask_innreg3d, imask_bpoint3d, eps_cont3d
use mod_interp1d, only: find_index, interpol_yp
use mod_interp2d, only: get_xy_values1, get_xy_values2, bilin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2_r, iim1_r, ii_r, iip1_r, iim2_th, iim1_th, ii_th, iip1_th
real(dp) :: r1, r2, rp, th1, th2, thp, phip
real(dp) :: vala, valb, valc, vald, valp
real(dp) :: velr, velth, velphi, vth_mean
real(dp) :: sint, cost, sinp, cosp
!
! ... local logicals
logical :: linfo_min, linfo_max, linfo_boundary, llogr, llogth, llogf
!
! ... local functions
real(dp) :: getr
!
!-----------------------------------------------------------------------
!
!read 2d model
call read_mod2d
!
!convert to different units
r_modext=r_modext/sr
!
!check if minimum radius is okay
if(abs(rmin-r_modext(1)).gt.1.d-14) then
   write(*,*) 'error in setup_mod2d: r_min does not match'
   write(*,*) 'set rstar to ', r_modext(1)*sr/rsu
   stop
endif
!
if(r_modext(nr_modext).lt.rlim) then
   write(*,*) 'error in setup_mod2d: r_lim does not match'
   write(*,*) 'set rlim to ', r_modext(nr_modext)
   stop
endif
!
!--------calculate calculation volume and perform interpolations--------
!
!calculate a mean thermal velocity (required outside information region)
vth_mean=0.d0
do i=1, nr_modext
   do j=1, ntheta_modext
      vth_mean= vth_mean+vth_modext2d(i,j)
   enddo
enddo
vth_mean=vth_mean/nr_modext/ntheta_modext
!
!default values of 3d-grid
imask_totreg3d=0
imask_innreg3d=0
imask_bpoint3d=0
opac3d=0.d0
opalbar3d=0.d0
t3d=0.d0
velx3d=0.d0
vely3d=0.d0
velz3d=0.d0
vth3d=vth_mean
eps_cont3d = zero
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
!         call info_region(x(i), y(j), z(k), rmin, rlim, linfo_min, linfo_max, linfo_boundary)
         call info_region_rot(x(i), y(j), z(k), smajorax_a, smajorax_b, smajorax_c, rlim, linfo_min, linfo_max, linfo_boundary)
         if(linfo_boundary) imask_bpoint3d(i,j,k)=1
         if(.not.linfo_min) imask_innreg3d(i,j,k)=1
         if(linfo_min.and.linfo_max) then
            imask_totreg3d(i,j,k)=1
            rp=getr(x(i), y(j), z(k))
            call get_angles_spc(x(i), y(j), z(k), thp, phip)
!get coordinates of vertices around given point
            call find_index(rp, r_modext, nr_modext, iim2_r, iim1_r, ii_r, iip1_r)
            call find_index(thp, theta_modext, ntheta_modext, iim2_th, iim1_th, ii_th, iip1_th)
            call get_xy_values1(rp, thp, r_modext, theta_modext, nr_modext, ntheta_modext, &
                                iim1_r, ii_r, iim1_th, ii_th, &
                                r1, r2, th1, th2, llogr, llogth)
!interpolation of density
            call get_xy_values2(nr_modext, ntheta_modext, rho_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            rho3d(i,j,k) = valp
!interpolation of temperature
            call get_xy_values2(nr_modext, ntheta_modext, t_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            t3d(i,j,k) = valp
!interpolation of eps_cont
            call get_xy_values2(nr_modext, ntheta_modext, eps_cont_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            eps_cont3d(i,j,k) = valp
!interpolation of thermal velocities
            call get_xy_values2(nr_modext, ntheta_modext, vth_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            vth3d(i,j,k) = valp
!interpolation of velocity components
            call get_xy_values2(nr_modext, ntheta_modext, velr_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velr)
            call get_xy_values2(nr_modext, ntheta_modext, velth_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velth)
            call get_xy_values2(nr_modext, ntheta_modext, velphi_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velphi)
            sint=sin(thp)
            cost=cos(thp)
            sinp=sin(phip)
            cosp=cos(phip)
            velx3d(i,j,k)=velr*sint*cosp + velth*cost*cosp-velphi*sinp
            vely3d(i,j,k)=velr*sint*sinp + velth*cost*sinp+velphi*cosp
            velz3d(i,j,k)=velr*cost - velth*sint
         endif
      enddo
   enddo
enddo
!
!
!j=ndymax/2+1
!k=ndzmax/2+1
!do i=1, ndxmax
!   write(*,*) x(i), imask_totreg3d(i,j,k), imask_innreg3d(i,j,k)
!enddo
!stop
!
!outside calculation volume (required for sc interpolations!)
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         rp=getr(x(i), y(j), z(k))         
         if(imask_totreg3d(i,j,k).eq.0) then
            call get_angles_spc(x(i), y(j), z(k), thp, phip)
!get coordinates of vertices around given point
            call find_index(rp, r_modext, nr_modext, iim2_r, iim1_r, ii_r, iip1_r)
            call find_index(thp, theta_modext, ntheta_modext, iim2_th, iim1_th, ii_th, iip1_th)
            call get_xy_values1(rp, thp, r_modext, theta_modext, nr_modext, ntheta_modext, &
                                iim1_r, ii_r, iim1_th, ii_th, &
                                r1, r2, th1, th2, llogr, llogth)
!interpolation of density
            call get_xy_values2(nr_modext, ntheta_modext, rho_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            if(valp.lt.0.d0) stop 'error in setup_mod2d: density extrapolated below 0'
            rho3d(i,j,k) = valp
!interpolation of density
            call get_xy_values2(nr_modext, ntheta_modext, eps_cont_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            if(valp.lt.0.d0) stop 'error in setup_mod2d: eps_cont extrapolated below 0'
            eps_cont3d(i,j,k) = valp
!interpolation of temperature
            call get_xy_values2(nr_modext, ntheta_modext, t_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, valp)
            if(valp.lt.0.d0) stop 'error in setup_mod2d: temperature extrapolated below 0'
            t3d(i,j,k) = valp
!interpolation of velocity components
            call get_xy_values2(nr_modext, ntheta_modext, velr_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velr)
            call get_xy_values2(nr_modext, ntheta_modext, velth_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velth)
            call get_xy_values2(nr_modext, ntheta_modext, velphi_modext2d, iim1_r, ii_r, &
                                iim1_th, ii_th, vala, valb, valc, vald, llogf)
            call bilin(rp, thp, r1, r2, th1, th2, vala, valb, valc, vald, &
                       llogr, llogth, llogf, velphi)
            sint=sin(thp)
            cost=cos(thp)
            sinp=sin(phip)
            cosp=cos(phip)
            velx3d(i,j,k)=velr*sint*cosp + velth*cost*cosp-velphi*sinp
            vely3d(i,j,k)=velr*sint*sinp + velth*cost*sinp+velphi*cosp
            velz3d(i,j,k)=velr*cost - velth*sint
         endif
      enddo
   enddo
enddo
!write(*,'(10es20.8)') maxval(sqrt(velr_modext2d**2+velth_modext2d**2+velphi_modext2d**2)), maxval(sqrt(velx3d**2+vely3d**2+velz3d**2))
!stop
!
!velocities in fiducial thermal velocities
velx3d=velx3d/vth_fiducial
vely3d=vely3d/vth_fiducial
velz3d=velz3d/vth_fiducial
!
!
!
!deallocate external model structure
deallocate(r_modext)
deallocate(theta_modext)
deallocate(rho_modext2d)
deallocate(velr_modext2d)
deallocate(velth_modext2d)
deallocate(velphi_modext2d)
deallocate(t_modext2d)
deallocate(vth_modext2d)
deallocate(eps_cont_modext2d)
!
end subroutine setup_mod2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_mod3d
!
use prog_type
use fund_const
use modext, only: nr_modext, ntheta_modext, nphi_modext, r_modext, &
                  theta_modext, phi_modext, vth_modext3d, &
                  rho_modext3d, t_modext3d, velr_modext3d, velth_modext3d, &
                  velphi_modext3d, eps_cont_modext3d
use params_input, only: rmax, vth_fiducial
use params_stellar, only: sr, smajorax_a, smajorax_b, smajorax_c
use inf_reg, only: rlim, rmin
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, &
                  rho3d, opac3d, opalbar3d, t3d, velx3d, vely3d, velz3d, &
                  vth3d, imask_innreg3d, imask_bpoint3d, eps_cont3d
use mod_interp3d, only: get_rtp_indx, get_rtp_values1, get_rtp_values2, trilin_complete
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2
real(dp) :: r1, r2, rp, th1, th2, thp, phi1, phi2, phip
real(dp) :: rada, radb
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velr, velth, velphi, vth_mean
real(dp) :: sint, cost, sinp, cosp
!
! ... local logicals
logical :: linfo_min, linfo_max, linfo_boundary, llogr, llogth, llogphi, llogf, expol
!
! ... local functions
real(dp) :: getr
!
!-----------------------------------------------------------------------
!
!read 3d model
call read_mod3d
!
!convert to different units
r_modext=r_modext/sr
!
!check if minimum radius is okay
if(abs(rmin-r_modext(1)).gt.1.d-14) then
   write(*,*) 'error in setup_mod3d: r_min does not match'
   write(*,*) 'set rstar to ', r_modext(1)*sr/rsu
   stop
endif
if(r_modext(nr_modext).lt.rlim) then
   write(*,*) 'error in setup_mod3d: r_lim does not match'
   write(*,*) 'set rlim to ', r_modext(nr_modext)
   stop
endif
!
!--------calculate calculation volume and perform interpolations--------
!
!calculate a mean thermal velocity (required outside information region)
vth_mean=0.d0
do i=1, nr_modext
   do j=1, ntheta_modext
      do k=1, nphi_modext
         vth_mean= vth_mean+vth_modext3d(i,j,k)
      enddo
   enddo
enddo
vth_mean=vth_mean/nr_modext/ntheta_modext/nphi_modext
!
!default values of 3d-grid
imask_totreg3d=0
imask_innreg3d=0
imask_bpoint3d=0
rho3d=0.d0
opac3d=0.d0
opalbar3d=0.d0
t3d=0.d0
velx3d=0.d0
vely3d=0.d0
velz3d=0.d0
vth3d=vth_mean
eps_cont3d=0.d0

do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
!         call info_region(x(i), y(j), z(k), rmin, rlim, linfo_min, linfo_max, linfo_boundary)
         call info_region_rot(x(i), y(j), z(k), smajorax_a, smajorax_b, smajorax_c, rlim, linfo_min, linfo_max, linfo_boundary)
         if(linfo_boundary) imask_bpoint3d(i,j,k)=1
         if(.not.linfo_min) imask_innreg3d(i,j,k)=1
         if(linfo_min.and.linfo_max) then
            imask_totreg3d(i,j,k)=1
            rp=getr(x(i), y(j), z(k))
            call get_angles_spc(x(i), y(j), z(k), thp, phip)
!get coordinates of vertices around given point
            call get_rtp_indx(rp, thp, phip, r_modext, theta_modext, phi_modext, &
                              nr_modext, ntheta_modext, nphi_modext, &
                              indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                              expol, rmin, rlim)
            call get_rtp_values1(rp, thp, phip, r_modext, theta_modext, phi_modext, &
                                 nr_modext, ntheta_modext, nphi_modext, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 r1, r2, th1, th2, phi1, phi2, &
                                 llogr, llogth, llogphi)
            rada=r_modext(indx_r1)
            radb=r_modext(indx_r2)
!interpolation of density
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, rho_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)
            rho3d(i,j,k) = valp
!interpolation of eps_cont
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, eps_cont_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)
            eps_cont3d(i,j,k) = valp
!interpolation of temperature
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, t_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)

            t3d(i,j,k) = valp
!interpolation of thermal velocities
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, vth_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)

            if(valp.le.zero) then
               write(*,*) vala, valb, valc, vald, vale, valf, valg, valh
               write(*,*) r1, rp, r2
               write(*,*) th1, thp, th2
               write(*,*) phi1, phip, phi2
               write(*,*) valp, t3d(i,j,k), rho3d(i,j,k)
               write(*,*) llogr, llogth, llogphi, llogf, expol
               write(*,*) r1, rp, r2, valp
               write(*,*) th1*180./pi, thp*180./pi, th2*180./pi
               write(*,*) phi1*180./pi, phip*180./pi, phi2*180./pi
               write(*,*)
               stop 'error in setup_model3d: vth <= zero not allowed'
            endif

            vth3d(i,j,k) = valp
!interpolation of velocity components
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velr_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velr)
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velth_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velth)
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velphi_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velphi)
            sint=sin(thp)
            cost=cos(thp)
            sinp=sin(phip)
            cosp=cos(phip)
            velx3d(i,j,k)=velr*sint*cosp + velth*cost*cosp-velphi*sinp
            vely3d(i,j,k)=velr*sint*sinp + velth*cost*sinp+velphi*cosp
            velz3d(i,j,k)=velr*cost - velth*sint

!            if(j.eq.ndymax/2+1.and.k.eq.ndzmax/2+1) then
!               write(*,*) x(i), velx3d(i,j,k), velr
!            endif

         endif
      enddo
   enddo
enddo

!i=ndxmax/2+1
!j=ndymax/2+1
!k=ndzmax/2+1
!do k=1, ndzmax
!   write(*,*) z(k), velz3d(i,j,k)/1.d5, rho3d(i,j,k)
!enddo
!stop
!
!write(*,*) minval(vth3d)
!
!stop

!
!
!outside calculation volume (required for sc interpolations!)
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         rp=getr(x(i), y(j), z(k))
         if(imask_totreg3d(i,j,k).eq.0) then         
!         if(imask_totreg3d(i,j,k).eq.0.and.&
!            imask_innreg3d(i,j,k).eq.0.and.rp.le.rmax) then
            call get_angles_spc(x(i), y(j), z(k), thp, phip)
!get coordinates of vertices around given point
            call get_rtp_indx(rp, thp, phip, r_modext, theta_modext, phi_modext, &
                              nr_modext, ntheta_modext, nphi_modext, &
                              indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                              expol, rmin, rlim)
            call get_rtp_values1(rp, thp, phip, r_modext, theta_modext, phi_modext, &
                                 nr_modext, ntheta_modext, nphi_modext, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 r1, r2, th1, th2, phi1, phi2, &
                                 llogr, llogth, llogphi)
            rada=r_modext(indx_r1)
            radb=r_modext(indx_r2)
!interpolation of density
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, rho_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)
            if(valp.lt.0.) stop 'error in setup_mod3d: density extrapolated below 0'
            rho3d(i,j,k) = valp
!interpolation of eps_cont
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, eps_cont_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)
            if(valp.lt.0.) stop 'error in setup_mod3d: eps_cont extrapolated below 0'
            eps_cont3d(i,j,k) = valp
!interpolation of temperature
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, t_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., valp)
            if(valp.lt.0.) stop 'error in setup_mod3d: temperature extrapolated below 0'
            t3d(i,j,k) = valp
!interpolation of velocity components
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velr_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velr)
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velth_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velth)
            call get_rtp_values2(nr_modext, ntheta_modext, nphi_modext, velphi_modext3d, &
                                 indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, llogf)
            call trilin_complete(rp, thp, phip, r1, r2, th1, th2, phi1, phi2, &
                                 vala, valb, valc, vald, vale, valf, valg, valh, &                                 
                                 rada, radb, rada, radb, rada, radb, rada, radb, rp, &
                                 expol, .false., llogr, llogth, llogphi, llogf, .false., velphi)
            sint=sin(thp)
            cost=cos(thp)
            sinp=sin(phip)
            cosp=cos(phip)
            velx3d(i,j,k)=velr*sint*cosp + velth*cost*cosp-velphi*sinp
            vely3d(i,j,k)=velr*sint*sinp + velth*cost*sinp+velphi*cosp
            velz3d(i,j,k)=velr*cost - velth*sint
         endif
      enddo
   enddo
enddo
!
!velocities in fiducial thermal velocities
velx3d=velx3d/vth_fiducial
vely3d=vely3d/vth_fiducial
velz3d=velz3d/vth_fiducial
!
!deallocate external model structure
deallocate(r_modext)
deallocate(rho_modext3d)
deallocate(velr_modext3d)
deallocate(velth_modext3d)
deallocate(velphi_modext3d)
deallocate(t_modext3d)
deallocate(vth_modext3d)
deallocate(eps_cont_modext3d)
!
!write(*,*) minval(vth3d)
!stop 'go on in setup_mod3d'

!j=ndymax/2+1
!k=ndzmax/2+1
!do i=1, ndxmax
!   write(*,*) x(i), imask_totreg3d(i,j,k)
!enddo
!stop


end subroutine setup_mod3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_mod1d
!
use prog_type
use mod_directories, only: model1d_file, model_dir
use modext, only: nr_modext, r_modext, rho_modext1d, t_modext1d, &
                  velr_modext1d, vth_modext1d, eps_cont_modext1d
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, attr_id, dset_id
integer(hsize_t), dimension(1) :: dims_rad, dims_scalars
!
!----------------read in radius and theta grid from 2d file-------------
!
write(*,*) '------------------read model atmosphere from 1d input file---------------------'
write(*,*) 'file name: ', trim(model_dir)//'/'//model1d_file
write(*,*) 
!
if(allocated(r_modext)) deallocate(r_modext)
if(allocated(rho_modext1d)) deallocate(rho_modext1d)
if(allocated(t_modext1d)) deallocate(t_modext1d)
if(allocated(velr_modext1d)) deallocate(velr_modext1d)
if(allocated(vth_modext1d)) deallocate(vth_modext1d) 
if(allocated(eps_cont_modext1d)) deallocate(eps_cont_modext1d) 
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(model_dir)//'/'//model1d_file, h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
allocate(rho_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
allocate(t_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
allocate(velr_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
allocate(vth_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
allocate(eps_cont_modext1d(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod1d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext1d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr_modext1d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext1d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vth_modext1d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'eps_cont', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, eps_cont_modext1d, dims_rad, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!
end subroutine read_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_mod2d
!
use prog_type
use mod_directories, only: model2d_file, model_dir
use modext, only: nr_modext, ntheta_modext, r_modext, theta_modext, &
                  rho_modext2d, t_modext2d, vth_modext2d, &
                  velr_modext2d, velth_modext2d, velphi_modext2d, eps_cont_modext2d
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, attr_id, dset_id
integer(hsize_t), dimension(1) :: dims_rad, dims_theta, dims_scalars
integer(hsize_t), dimension(2) :: dims
!
!----------------read in radius and theta grid from 2d file-------------
!
write(*,*) '-----------------read model atmosphere from 2d input file----------------------'
write(*,*) 'file name: ', trim(model_dir)//'/'//model2d_file
write(*,*)
!
if(allocated(r_modext)) deallocate(r_modext)
if(allocated(theta_modext)) deallocate(theta_modext)
if(allocated(rho_modext2d)) deallocate(rho_modext2d)
if(allocated(t_modext2d)) deallocate(t_modext2d)
if(allocated(velr_modext2d)) deallocate(velr_modext2d)
if(allocated(velth_modext2d)) deallocate(velth_modext2d)
if(allocated(velphi_modext2d)) deallocate(velphi_modext2d)
if(allocated(vth_modext2d)) deallocate(vth_modext2d)    
if(allocated(eps_cont_modext2d)) deallocate(eps_cont_modext2d)    
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(model_dir)//'/'//model2d_file, h5f_acc_rdonly_f, file_id, err)
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
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims=(/ nr_modext, ntheta_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation r_modext'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation theta_modext'
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation, rho_modext2d'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation t_modext2d'
allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation velr_modext2d'
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation velth_modext2d'
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation velphi_modext2d'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation vth_modext2d'
allocate(eps_cont_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod2d: allocation eps_cont_modext2d'
!
!
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims, err)
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
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'eps_cont', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, eps_cont_modext2d, dims, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!
end subroutine read_mod2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_mod3d
!
use prog_type
use mod_directories, only: model3d_file, model_dir
use modext, only: nr_modext, ntheta_modext, nphi_modext, r_modext, &
                  theta_modext, phi_modext, rho_modext3d, t_modext3d, &
                  vth_modext3d, velr_modext3d, velth_modext3d, &
                  velphi_modext3d, eps_cont_modext3d
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, attr_id, dset_id
integer(hsize_t), dimension(1) :: dims_rad, dims_theta, dims_phi, dims_scalars
integer(hsize_t), dimension(3) :: dims
!
!----------------read in radius and theta grid from 3d file-------------
!
write(*,*) '-----------------read model atmosphere from 3d input file----------------------'
write(*,*) 'file name: ', trim(model_dir)//'/'//model3d_file
write(*,*)
!
if(allocated(r_modext)) deallocate(r_modext)
if(allocated(theta_modext)) deallocate(theta_modext)
if(allocated(phi_modext)) deallocate(phi_modext)
if(allocated(rho_modext3d)) deallocate(rho_modext3d)
if(allocated(t_modext3d)) deallocate(t_modext3d)
if(allocated(velr_modext3d)) deallocate(velr_modext3d)
if(allocated(velth_modext3d)) deallocate(velth_modext3d)
if(allocated(velphi_modext3d)) deallocate(velphi_modext3d)
if(allocated(vth_modext3d)) deallocate(vth_modext3d)
if(allocated(eps_cont_modext3d)) deallocate(eps_cont_modext3d)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(model_dir)//'/'//model3d_file, h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'nphi', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nphi_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)

call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
dims_phi=(/nphi_modext/)
dims=(/ nr_modext, ntheta_modext, nphi_modext /)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(phi_modext(nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error in read_mod3d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'phi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, phi_modext, dims_phi, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vth', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vth_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'eps_cont', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, eps_cont_modext3d, dims, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!
end subroutine read_mod3d

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_phot1d
!
use prog_type
use mod_directories, only: model1d_phot_file, model_dir
use modext, only: nr_modphot, r_modphot, rho_modphot1d, t_modphot1d
use params_input, only: teff, rstar, yhe, vmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: teff_modphot, rstar_modphot, yhe_modphot, vmin_modphot
!
!----------------read in radius and theta grid from 2d file-------------
!
write(*,*) '------------------read photospheric model from 1d input file-------------------'
write(*,*) 'file name: ', trim(model_dir)//'/'//model1d_phot_file
write(*,*) 
!
open(1, file=trim(model_dir)//'/'//trim(model1d_phot_file), form='formatted')
   read(1,*)
   read(1,'(i5, 4e20.8)') nr_modphot, teff_modphot, rstar_modphot, yhe_modphot, vmin_modphot
   read(1,*)
   allocate(r_modphot(nr_modphot))
   allocate(rho_modphot1d(nr_modphot))
   allocate(t_modphot1d(nr_modphot))
   do i=1, nr_modphot
      read(1,'(3e20.8)') r_modphot(i), t_modphot1d(i), rho_modphot1d(i)
   enddo
close(1)
!
if(abs(teff-teff_modphot).gt.1.d-14) write(*,*) 'error in read_phot1d: input teff not consistent'
if(abs(rstar-rstar_modphot).gt.1.d-14) write(*,*) 'error in read_phot1d: input rstar not consistent'
if(abs(yhe-yhe_modphot).gt.1.d-14) write(*,*) 'error in read_phot1d: input yhe not consistent'
if(abs(vmin-vmin_modphot).gt.1.d-14) write(*,*) 'error in read_phot1d: input vmin not consistent'
!
!
!
end subroutine read_phot1d

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_temperatures1
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, imask3d, &
                  t3d, scont3d, x
use params_input, only: tmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
!
! ... local logicals
!
! ... local functions
!
write(*,*) 'using LTE temperatures'
write(*,*)
!
t3d=tmin
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(0) !outside of calculation volume
               t3d(i,j,k)=tmin
            case(1,2,3,4) !standard points, near boundary points, boundary points
               t3d(i,j,k)=(pi*scont3d(i,j,k)/cgs_sb)**0.25d0
            case default
         end select
      enddo
   enddo
enddo
!
j=ndymax/2+1
k=ndzmax/2+1
do i=1, ndxmax
   write(*,'(es20.8, i5, 10es20.8)') x(i), imask3d(i,j,k), scont3d(i,j,k), t3d(i,j,k)
enddo
!stop 'go on in setup_temperatures'
!
end subroutine setup_temperatures1

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_opac1
!
use prog_type
use params_input, only: kcont, yhe, hei
use params_stellar, only: sr
use dime3d, only: ndxmax, ndymax, ndzmax, imask_totreg3d, &
     rho3d, opac3d
use mod_opacities, only: opac_thomson, opalbar_model_kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: rho
!
! ... local logicals
!
! ... local functions
!
write(*,*) 'using thomson-opacity parameterization'
write(*,*)
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
         if(imask_totreg3d(i,j,k).ne.0) then
            rho=rho3d(i,j,k)
            opac3d(i,j,k) = opac_thomson(yhe, hei, rho, kcont)*sr
         else
            opac3d(i,j,k) = 1.d-10      !set to low value (ne 0) to avoid NaN's
         endif
      enddo
   enddo
enddo
!
!in units 1/rstar
!opac3d=opac3d*sr
!opalbar3d=opalbar3d*sr
!
!j=ndymax/2+1
!k=ndzmax/2+1
!do i=1, ndxmax
!   write(*,'(10es20.8)') x(i), opac3d(i,j,k), opalbar3d(i,j,k), t3d(i,j,k)
!enddo
!stop 'go on in setup_opacities1'
!
end subroutine setup_opac1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_opalbar1
!
use prog_type
use params_input, only: kline, yhe, hei
use params_stellar, only: sr
use dime3d, only: ndxmax, ndymax, ndzmax, imask_totreg3d, &
                  rho3d, opalbar3d, t3d, x
use mod_opacities, only: opac_thomson, opalbar_model_kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: rho
!
! ... local logicals
!
! ... local functions
!
write(*,*) 'using line-strength parameterization'
write(*,*)
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
         if(imask_totreg3d(i,j,k).ne.0) then
            rho=rho3d(i,j,k)
            opalbar3d(i,j,k) = opalbar_model_kline(yhe, hei, rho, kline)*sr
            if(t3d(i,j,k).ge.1.d5) opalbar3d(i,j,k)=1.d-10
         else
            opalbar3d(i,j,k) = 1.d-10   !set to low value (ne 0) to avoid NaN's
         endif
      enddo
   enddo
enddo
!
!in units 1/rstar
!opac3d=opac3d*sr
!opalbar3d=opalbar3d*sr
!
!j=ndymax/2+1
!k=ndzmax/2+1
!do i=1, ndxmax
!   write(*,'(10es20.8)') x(i), opac3d(i,j,k), opalbar3d(i,j,k), t3d(i,j,k)
!enddo
!stop 'go on in setup_opacities1'
!
end subroutine setup_opalbar1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_opalbar2
!
use prog_type
use fund_const, only: xmsu
use params_input, only: yhe, hei, vmax, kappa0, xmloss, vth_fiducial, alpha
use params_stellar, only: sr
use dime3d, only: ndxmax, ndymax, ndzmax, imask_totreg3d, &
                  rho3d, opalbar3d, x, y, z, velx3d, vely3d, velz3d, t3d
use mod_opacities, only: opalbar_model_hamann
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: rho, rad, xmloss_cgs, vinf, dum
!
! ... local logicals
!
! ... local functions
!
write(*,*) 'using parameterization by Hamann with'
write(*,*) 'mdot ', xmloss
write(*,*) 'vinf ', vmax
write(*,*) 
!
!transform vmax and mass-loss rate to cgs
vinf=vmax*1.d5
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
         if(imask_totreg3d(i,j,k).ne.0) then
            rho=rho3d(i,j,k)
            rad=sqrt(x(i)**2 + y(j)**2 + z(k)**2)*sr
            opalbar3d(i,j,k) = opalbar_model_hamann(sr, vinf, xmloss_cgs, kappa0, alpha, vth_fiducial, rad, rho)*sr
            if(t3d(i,j,k).ge.1.d5) opalbar3d(i,j,k)=1.d-10
         else
            opalbar3d(i,j,k) = 1.d-14   !set to low value (ne 0) to avoid NaN's
         endif           
      enddo
   enddo
enddo
!
!do i=1, ndxmax
!   write(*,*) mask_totreg3d(i,i,i), opac3d(i,i,i), opalbar3d(i,i,i)
!enddo
!stop
!in units 1/rstar
!opac3d=opac3d*sr
!opalbar3d=opalbar3d*sr
!
end subroutine setup_opalbar2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_phot1d
!
!interpolation 1d photospheric model
!
use prog_type
use fund_const, only: rsu
use modext, only: nr_modphot, r_modphot, rho_modphot1d, t_modphot1d
use params_stellar, only: sr
use params_input, only: trad
use inf_reg, only: rmin, rlim
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, &
                  imask_innreg3d, imask_bpoint3d, &
                  rho3d, opac3d, opalbar3d, t3d, vth3d, velx3d, vely3d, velz3d
use mod_interp1d, only: find_index, interpol_yp
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l
integer(i4b) :: iip, iin, jjp, jjn, kkp, kkn, &
                iim2, iim1, ii, iip1
real(dp) :: dx, dxp, dxn, dy, dyp, dyn, dz, dzp, dzn, &
            dxy, dxz, dyz, qc_ijk, ql_ijk, &
            qcx1_ijk, qcx2_ijk, qcx_ijk, &
            qcy1_ijk, qcy2_ijk, qcy_ijk, &
            qcz1_ijk, qcz2_ijk, qcz_ijk, &
            qlx1_ijk, qlx2_ijk, qlx_ijk, &
            qly1_ijk, qly2_ijk, qly_ijk, &
            qlz1_ijk, qlz2_ijk, qlz_ijk
real(dp), parameter :: p=2.d0   !power for inverse distance weighting
real(dp) :: rad, rho, vth_surface, temp
!
! ... local logicals
!
! ... local functions
real(dp) :: getr
!
!-----------------------------------------------------------------------
!
!read 1d photospheric model
!not required anymore
!call read_phot1d
!
!-----------------------------------------------------------------------
!
!use a constant thermal velocity
do i=1, ndxmax
   if(x(i).ge.1.d0) then
      vth_surface=vth3d(i,ndymax/2+1,ndzmax/2+1)
      exit
   endif
enddo
!
!-----------------------------------------------------------------------
!
do i=1, ndxmax-1
   do j=1, ndymax-1
      do k=1, ndzmax-1
         if(imask_innreg3d(i,j,k).ne.0) then
!
!-----------------------------------------------------------------------
!
!find indices for extrapolation along x, y, z
!
!along positive x
            do l=i, ndxmax
               if(imask_totreg3d(l,j,k).ne.0) then
                  iip=l
                  exit
                endif
            enddo
!along negative x
            do l=i, 1, -1
               if(imask_totreg3d(l,j,k).ne.0) then
                  iin=l
                  exit
                endif
            enddo
!along positive y
            do l=j, ndymax
               if(imask_totreg3d(i,l,k).ne.0) then
                  jjp=l
                  exit
                endif
            enddo
!along negative y
            do l=j, 1, -1
               if(imask_totreg3d(i,l,k).ne.0) then
                  jjn=l
                  exit
                endif
            enddo
!along positive z
            do l=k, ndzmax
               if(imask_totreg3d(i,j,l).ne.0) then
                  kkp=l
                  exit
                endif
            enddo
!along negative z
            do l=k, 1, -1
               if(imask_totreg3d(i,j,l).ne.0) then
                  kkn=l
                  exit
                endif
            enddo
!
!-----------------------------------------------------------------------
!
!extrapolation of continuum and line opacities
!
!perform extrapolation along x
            dxp=x(iip)-x(i)
            dxn=x(i)-x(iin)
            if(dxp.lt.dxn) then
!use values on positive x-axis
               dx=dxp
!               qcx_ijk = interpol_yp(x(iip), x(iip+1), opac3d(iip,j,k)/rho3d(iip,j,k)/sr, &
!                                     opac3d(iip+1,j,k)/rho3d(iip+1,j,k)/sr, x(i))
!               qlx_ijk = interpol_yp(x(iip), x(iip+1), opalbar3d(iip,j,k)/rho3d(iip,j,k)/sr, &
!                                     opalbar3d(iip+1,j,k)/rho3d(iip+1,j,k)/sr, x(i))
               qcx_ijk = interpol_yp(x(iip), x(iip+1), opac3d(iip,j,k), opac3d(iip+1,j,k), x(i))
               qlx_ijk = interpol_yp(x(iip), x(iip+1), opalbar3d(iip,j,k), opalbar3d(iip+1,j,k), x(i))
            elseif(dxp.gt.dxn) then
!use values on negative x-axis
               dx=dxn
!               qcx_ijk = interpol_yp(x(iin), x(iin-1), opac3d(iin,j,k)/rho3d(iin,j,k)/sr, &
!                                      opac3d(iin-1,j,k)/rho3d(iin-1,j,k)/sr, x(i))
!               qlx_ijk = interpol_yp(x(iin), x(iin-1), opalbar3d(iin,j,k)/rho3d(iin,j,k)/sr, &
!                                      opalbar3d(iin-1,j,k)/rho3d(iin-1,j,k)/sr, x(i))
               qcx_ijk = interpol_yp(x(iin), x(iin-1), opac3d(iin,j,k), opac3d(iin-1,j,k), x(i))
               qlx_ijk = interpol_yp(x(iin), x(iin-1), opalbar3d(iin,j,k), opalbar3d(iin-1,j,k), x(i))
             else
!use mean of both values
               dx=dxp
!               qcx1_ijk = interpol_yp(x(iip), x(iip+1), opac3d(iip,j,k)/rho3d(iip,j,k)/sr, &
!                                      opac3d(iip+1,j,k)/rho3d(iip+1,j,k)/sr, x(i))
!               qcx2_ijk = interpol_yp(x(iin), x(iin-1), opac3d(iin,j,k)/rho3d(iin,j,k)/sr, &
!                                      opac3d(iin-1,j,k)/rho3d(iin-1,j,k)/sr, x(i))
!               qcx_ijk = (qcx1_ijk+qcx2_ijk)/2.d0
!               qlx1_ijk = interpol_yp(x(iip), x(iip+1), opalbar3d(iip,j,k)/rho3d(iip,j,k)/sr, &
!                                      opalbar3d(iip+1,j,k)/rho3d(iip+1,j,k)/sr, x(i))
!               qlx2_ijk = interpol_yp(x(iin), x(iin-1), opalbar3d(iin,j,k)/rho3d(iin,j,k)/sr, &
!                                      opalbar3d(iin-1,j,k)/rho3d(iin-1,j,k)/sr, x(i))
!               qlx_ijk = (qlx1_ijk+qlx2_ijk)/2.d0
               qcx1_ijk = interpol_yp(x(iip), x(iip+1), opac3d(iip,j,k), opac3d(iip+1,j,k), x(i))
               qcx2_ijk = interpol_yp(x(iin), x(iin-1), opac3d(iin,j,k), opac3d(iin-1,j,k), x(i))
               qcx_ijk = (qcx1_ijk+qcx2_ijk)/2.d0
               qlx1_ijk = interpol_yp(x(iip), x(iip+1), opalbar3d(iip,j,k), opalbar3d(iip+1,j,k), x(i))
               qlx2_ijk = interpol_yp(x(iin), x(iin-1), opalbar3d(iin,j,k), opalbar3d(iin-1,j,k), x(i))
               qlx_ijk = (qlx1_ijk+qlx2_ijk)/2.d0
            endif
!
!perform extrapolation along y
            dyp=y(jjp)-y(j)
            dyn=y(j)-y(jjn)
            if(dyp.lt.dyn) then
!use values on positive y-axis
               dy=dyp
!               qcy_ijk = interpol_yp(y(jjp), y(jjp+1), opac3d(i,jjp,k)/rho3d(i,jjp,k)/sr, &
!                                      opac3d(i,jjp+1,k)/rho3d(i,jjp+1,k)/sr, y(j))
!               qly_ijk = interpol_yp(y(jjp), y(jjp+1), opalbar3d(i,jjp,k)/rho3d(i,jjp,k)/sr, &
!                                      opalbar3d(i,jjp+1,k)/rho3d(i,jjp+1,k)/sr, y(j))
               qcy_ijk = interpol_yp(y(jjp), y(jjp+1), opac3d(i,jjp,k), opac3d(i,jjp+1,k), y(j))
               qly_ijk = interpol_yp(y(jjp), y(jjp+1), opalbar3d(i,jjp,k), opalbar3d(i,jjp+1,k), y(j))
            elseif(dyp.gt.dyn) then
!use values on negative y-axis
               dy=dyn
!               qcy_ijk = interpol_yp(y(jjn), y(jjn-1), opac3d(i,jjn,k)/rho3d(i,jjn,k)/sr, &
!                                      opac3d(i,jjn-1,k)/rho3d(i,jjn-1,k)/sr, y(j))
!               qly_ijk = interpol_yp(y(jjn), y(jjn-1), opalbar3d(i,jjn,k)/rho3d(i,jjn,k)/sr, &
!                                      opalbar3d(i,jjn-1,k)/rho3d(i,jjn-1,k)/sr, y(j))
               qcy_ijk = interpol_yp(y(jjn), y(jjn-1), opac3d(i,jjn,k), opac3d(i,jjn-1,k), y(j))
               qly_ijk = interpol_yp(y(jjn), y(jjn-1), opalbar3d(i,jjn,k), opalbar3d(i,jjn-1,k), y(j))
            else
!use mean of both values
               dy=dyp
!               qcy1_ijk = interpol_yp(y(jjp), y(jjp+1), opac3d(i,jjp,k)/rho3d(i,jjp,k)/sr, &
!                                      opac3d(i,jjp+1,k)/rho3d(i,jjp+1,k)/sr, y(j))
!               qcy2_ijk = interpol_yp(y(jjn), y(jjn-1), opac3d(i,jjn,k)/rho3d(i,jjn,k)/sr, &
!                                      opac3d(i,jjn-1,k)/rho3d(i,jjn-1,k)/sr, y(j))
!               qcy_ijk = (qcy1_ijk+qcy2_ijk)/2.d0
!               qly1_ijk = interpol_yp(y(jjp), y(jjp+1), opalbar3d(i,jjp,k)/rho3d(i,jjp,k)/sr, &
!                                      opalbar3d(i,jjp+1,k)/rho3d(i,jjp+1,k)/sr, y(j))
!               qly2_ijk = interpol_yp(y(jjn), y(jjn-1), opalbar3d(i,jjn,k)/rho3d(i,jjn,k)/sr, &
!                                      opalbar3d(i,jjn-1,k)/rho3d(i,jjn-1,k)/sr, y(j))
!               qly_ijk = (qly1_ijk+qly2_ijk)/2.d0
               qcy1_ijk = interpol_yp(y(jjp), y(jjp+1), opac3d(i,jjp,k), opac3d(i,jjp+1,k), y(j))
               qcy2_ijk = interpol_yp(y(jjn), y(jjn-1), opac3d(i,jjn,k), opac3d(i,jjn-1,k), y(j))
               qcy_ijk = (qcy1_ijk+qcy2_ijk)/2.d0
               qly1_ijk = interpol_yp(y(jjp), y(jjp+1), opalbar3d(i,jjp,k), opalbar3d(i,jjp+1,k), y(j))
               qly2_ijk = interpol_yp(y(jjn), y(jjn-1), opalbar3d(i,jjn,k), opalbar3d(i,jjn-1,k), y(j))
               qly_ijk = (qly1_ijk+qly2_ijk)/2.d0
            endif
!
!perform extrapolation along z
            dzp=z(kkp)-z(k)
            dzn=z(k)-z(kkn)
            if(dzp.lt.dzn) then
!use values on positive z-axis
               dz=dzp
!               qcz_ijk = interpol_yp(z(kkp), z(kkp+1), opac3d(i,j,kkp)/rho3d(i,j,kkp)/sr, &
!                                      opac3d(i,j,kkp+1)/rho3d(i,j,kkp+1)/sr, z(k))
!               qlz_ijk = interpol_yp(z(kkp), z(kkp+1), opalbar3d(i,j,kkp)/rho3d(i,j,kkp)/sr, &
!                                      opalbar3d(i,j,kkp+1)/rho3d(i,j,kkp+1)/sr, z(k))
               qcz_ijk = interpol_yp(z(kkp), z(kkp+1), opac3d(i,j,kkp), opac3d(i,j,kkp+1), z(k))
               qlz_ijk = interpol_yp(z(kkp), z(kkp+1), opalbar3d(i,j,kkp), opalbar3d(i,j,kkp+1), z(k))
            elseif(dzp.gt.dzn) then
!use values on negative z-axis
               dz=dzn
!               qcz_ijk = interpol_yp(z(kkn), z(kkn-1), opac3d(i,j,kkn)/rho3d(i,j,kkn)/sr, &
!                                      opac3d(i,j,kkn-1)/rho3d(i,j,kkn-1)/sr, z(k))
!               qlz_ijk = interpol_yp(z(kkn), z(kkn-1), opalbar3d(i,j,kkn)/rho3d(i,j,kkn)/sr, &
!                                      opalbar3d(i,j,kkn-1)/rho3d(i,j,kkn-1)/sr, z(k))
               qcz_ijk = interpol_yp(z(kkn), z(kkn-1), opac3d(i,j,kkn), opac3d(i,j,kkn-1), z(k))
               qlz_ijk = interpol_yp(z(kkn), z(kkn-1), opalbar3d(i,j,kkn), opalbar3d(i,j,kkn-1), z(k))
            else
!use mean of both values
!               dz=dzp
!               qcz1_ijk = interpol_yp(z(kkp), z(kkp+1), opac3d(i,j,kkp)/rho3d(i,j,kkp)/sr, &
!                                      opac3d(i,j,kkp+1)/rho3d(i,j,kkp+1)/sr, z(k))
!               qcz2_ijk = interpol_yp(z(kkn), z(kkn-1), opac3d(i,j,kkn)/rho3d(i,j,kkn)/sr, &
!                                      opac3d(i,j,kkn-1)/rho3d(i,j,kkn-1)/sr, z(k))
!               qcz_ijk = (qcz1_ijk+qcz2_ijk)/2.d0
!               qlz1_ijk = interpol_yp(z(kkp), z(kkp+1), opalbar3d(i,j,kkp)/rho3d(i,j,kkp)/sr, &
!                                      opalbar3d(i,j,kkp+1)/rho3d(i,j,kkp+1)/sr, z(k))
!               qlz2_ijk = interpol_yp(z(kkn), z(kkn-1), opalbar3d(i,j,kkn)/rho3d(i,j,kkn)/sr, &
!                                      opalbar3d(i,j,kkn-1)/rho3d(i,j,kkn-1)/sr, z(k))
!               qlz_ijk = (qlz1_ijk+qlz2_ijk)/2.d0
               qcz1_ijk = interpol_yp(z(kkp), z(kkp+1), opac3d(i,j,kkp), opac3d(i,j,kkp+1), z(k))
               qcz2_ijk = interpol_yp(z(kkn), z(kkn-1), opac3d(i,j,kkn), opac3d(i,j,kkn-1), z(k))
               qcz_ijk = (qcz1_ijk+qcz2_ijk)/2.d0
               qlz1_ijk = interpol_yp(z(kkp), z(kkp+1), opalbar3d(i,j,kkp), opalbar3d(i,j,kkp+1), z(k))
               qlz2_ijk = interpol_yp(z(kkn), z(kkn-1), opalbar3d(i,j,kkn), opalbar3d(i,j,kkn-1), z(k))
               qlz_ijk = (qlz1_ijk+qlz2_ijk)/2.d0
            endif
!
!-----------------------------------------------------------------------
!
!calculate a weighted mean of qx_ijk, qy_ijk, qz_ijk (inverse distance weights)
            dxy=(dx*dy)**p
            dxz=(dx*dz)**p
            dyz=(dy*dz)**p
            qc_ijk = (dyz*qcx_ijk +dxz*qcy_ijk + dxy*qcz_ijk)/(dxy+dxz+dyz)
            ql_ijk = (dyz*qlx_ijk +dxz*qly_ijk + dxy*qlz_ijk)/(dxy+dxz+dyz)
!
!-----------------------------------------------------------------------
!
!interpolate density and temperature from the 1d photospheric structure
!            rad = getr(x(i), y(j), z(k))
!            call find_index(rad, r_modphot, nr_modphot, iim2, iim1, ii, iip1)
!            if(rad.ne.0.d0) then
!!density
!               rho=interpol_yp(log10(r_modphot(iim1)), log10(r_modphot(ii)), &
!                               log10(rho_modphot1d(iim1)), log10(rho_modphot1d(ii)), log10(rad))
!               rho3d(i,j,k)=10.d0**rho
!!continuum and line opacity (in rstar)
!               opac3d(i,j,k)=sr*qc_ijk*10.d0**rho
!               opalbar3d(i,j,k)=sr*ql_ijk*10.d0**rho
!!temperature
!               temp=interpol_yp(log10(r_modphot(iim1)), log10(r_modphot(ii)), &
!                                log10(t_modphot1d(iim1)), log10(t_modphot1d(ii)), log10(rad))
!               t3d(i,j,k)=10.d0**temp
!            else
!               rho=interpol_yp(r_modphot(iim1), r_modphot(ii), &
!                               rho_modphot1d(iim1), rho_modphot1d(ii), rad)
!               rho3d(i,j,k)=rho
!!continuum and line opacity (in rstar)
!               opac3d(i,j,k)=sr*qc_ijk*rho
!               opalbar3d(i,j,k)=sr*ql_ijk*rho
!!temperature
!               temp=interpol_yp(r_modphot(iim1), r_modphot(ii), &
!                                t_modphot1d(iim1), t_modphot1d(ii), rad)
!               t3d(i,j,k)=temp
!            endif
!
!continuum and line opacity, avoiding extrapolation to negative values
            if(qc_ijk.le.0.d0) then
               opac3d(i,j,k)=1.d-10
            else
               opac3d(i,j,k)=qc_ijk
            endif
            if(ql_ijk.le.0.d0) then
               opalbar3d(i,j,k)=1.d-10
            else
               opalbar3d(i,j,k)=ql_ijk
            endif
!
!temperature
            t3d(i,j,k)=trad
!
!set velocities to zero
            velx3d(i,j,k) = 0.d0
            vely3d(i,j,k) = 0.d0
            velz3d(i,j,k) = 0.d0
!
!set thermal velocity to the surface value
            vth3d(i,j,k)=vth_surface
!
         endif
      enddo
   enddo
enddo
!
!***debug start
!
!open(1,file='TRASH/test3d.dat', form='formatted')
!   do i=1, ndxmax
!      do j=1, ndymax
!          do k=1, ndzmax
!             write(1,'(8es20.8)') getr(x(i),y(j),z(k)), opac3d(i,j,k), opalbar3d(i,j,k), t3d(i,j,k), vth3d(i,j,k), velx3d(i,j,k), vely3d(i,j,k), velz3d(i,j,k)
!          enddo
!      enddo
!   enddo
!close(1)
!
!open(1,file='TRASH/test1d.dat', form='formatted')
!   do i=1, nr_modphot
!      write(1,'(2es20.8)') r_modphot(i), rho_modphot1d(i)
!   enddo
!close(1)
!
!***debug end
!
!deallocate external model structure (not required anymore)
!deallocate(r_modphot)
!deallocate(rho_modphot1d)
!deallocate(t_modphot1d)
!
!
!
end subroutine setup_phot1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_mask3d
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, &
                  imask_innreg3d, imask_bpoint3d, imask3d, imask_totreg3d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, l, m, n
!
! ... local logicals
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!default value: imask3d=0
imask3d=0
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(imask_totreg3d(i,j,k).eq.1) then 
!standard interpolation scheme
            imask3d(i,j,k)=1
!overwrite values if in vicinity of inner boundary
            if(imask_innreg3d(i-1,j-1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j-1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j-1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j-1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j-1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j-1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j-1,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j-1,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j-1,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j  ,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j  ,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j  ,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j  ,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j  ,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j  ,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j  ,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j  ,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j+1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j+1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j+1,k+1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j+1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j+1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j+1,k).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i-1,j+1,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i,  j+1,k-1).eq.1) imask3d(i,j,k)=2
            if(imask_innreg3d(i+1,j+1,k-1).eq.1) imask3d(i,j,k)=2
!overwrite values if point is directly on boundary
            if(imask_bpoint3d(i,j,k).eq.1) imask3d(i,j,k)=3
         elseif(imask_innreg3d(i,j,k).eq.1) then
!points for which boundary intensities need to be specified
            do l=-2, 2
               do m=-2, 2
                  do n=-2, 2
                     if(imask_totreg3d(i+l,j+m,k+n).eq.1) imask3d(i,j,k)=4
                  enddo
               enddo
            enddo
         endif
      enddo
   enddo
enddo
!
!ensure that origin is not in not set (otherwise divisions by zero may occurr)
imask3d(ndxmax/2+1,ndymax/2+1,ndzmax/2+1)=0
!
!
!l=0
!m=0
!n=0
!j=ndymax/2+1
!do j=1, ndymax
!do i=1, ndxmax
!   do k=1, ndzmax
!      if(imask3d(i,j,k).eq.4) l=l+1
!      if(imask3d(i,j,k).eq.3) m=m+1
!      if(imask3d(i,j,k).eq.2) n=n+1
!   enddo
!enddo
!enddo
!write(*,*) 'ndxmax, ndymax, ndzmax', ndxmax, ndymax, ndzmax
!write(*,*) 'number of points inside star', l
!write(*,*) 'number of points on boundary', m
!write(*,*) 'number of points near surfac', n
!stop
!
!
!
end subroutine setup_mask3d
