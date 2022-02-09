!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, beta, vmin, vmax, vth_fiducial, xic, xnue0, epsl, ssobo3d)
!
!----------------calculates sobolev source function in 3d---------------
!without continuum
!neglecting multiple resonances
!note: velocity gradients are calculated with finite differences from the 3d cartesian grid, 
!      and are affected by numerical errors for low resolution grids
!
!units in:  radius in rstar
!           velocity in vth_fiducial
!           opalbar: integrated line opacity in 1/rstar, and in xobs-space, normalized to vth
!           temp    temperature in kelvin
!           xic     core-intensity
!
use prog_type
use fund_const, only: pi, zero
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
real(dp), dimension(ndxmax), intent(in) :: x
real(dp), dimension(ndymax), intent(in) :: y
real(dp), dimension(ndzmax), intent(in) :: z
real(dp), dimension(ndxmax,ndymax,ndzmax), intent(in) :: velx3d, vely3d, velz3d, opalbar3d, t3d, vth3d
integer(i1b), dimension(ndxmax,ndymax,ndzmax), intent(in) :: imask_totreg3d
real(dp), intent(in) :: xic, xnue0, epsl, beta, vmin, vmax, vth_fiducial
real(dp), dimension(ndxmax,ndymax,ndzmax), intent(inout) :: ssobo3d

!
! ... local scalars
integer(i4b) :: i, j, k, imu, iphi
integer(i4b), parameter :: nmu_c=31, nmu_nc=61, nphi=33
real(dp) :: del, mustar, tau, r_ijk, bval, vth
real(dp) :: beta_sobolev, betac_sobolev
real(dp) :: n_x, n_y, n_z
real(dp) :: theta_tilde, phi_tilde, cost, sint, cosp, sinp, costt, sintt, cospt, sinpt
real(dp) :: dvxdx, dvydy, dvzdz, dvydx, dvzdx, dvxdy, dvzdy, dvxdz, dvydz, dvsds
!
! ... for debugging
!real(dp) :: x_ii, y_jj, z_kk, vdum, b, vmin, vmax, beta, dvdum
!real(dp) :: dvxdx2, dvydy2, dvzdz2, dvydx2, dvzdx2, dvxdy2, dvzdy2, dvxdz2, dvydz2, dvsds2
real(dp) :: velr, gradv, bconst
!
! ... local arrays
real(dp), dimension(nmu_c) :: nodes_mu_c, weight_mu_c
real(dp), dimension(nmu_nc) :: nodes_mu_nc, weight_mu_nc
real(dp), dimension(nphi) :: nodes_phi, weight_phi
!
! ... local functions
real(dp) :: bnue, sobo1d
!
!
!******for the moment, calculate only 1d sobolev source function********
!
write(*,*) '--------------------calculating 1d sobolev source function---------------------'
write(*,*)
!
ssobo3d=zero
!
bconst=1.d0-(vmin/vmax)**(1.d0/beta)
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(imask_totreg3d(i,j,k).ne.0) then
!in vthermal
            velr=sqrt(velx3d(i,j,k)**2 + vely3d(i,j,k)**2 + velz3d(i,j,k)**2)*vth_fiducial/vth3d(i,j,k)
            r_ijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
            gradv=velr * bconst*beta/r_ijk**2/(1.d0-bconst/r_ijk)
            if(abs(gradv).lt.0.d0) stop 'error in sobo3d: gradv=0 is not allowed'
            ssobo3d(i,j,k)=sobo1d(r_ijk, velr, gradv, opalbar3d(i,j,k), t3d(i,j,k), xic, xnue0, epsl)
         endif
      enddo
   enddo
enddo


!i=ndxmax/2+1
!j=ndymax/2+1
!do k=1, ndzmax
!   if(imask_totreg3d(i,j,k).ne.0) then   
!      velr=sqrt(velx3d(i,j,k)**2 + vely3d(i,j,k)**2 + velz3d(i,j,k)**2)*vth_fiducial/vth3d(i,j,k)
!      r_ijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
!      gradv=velr * bconst*beta/r_ijk**2/(1.d0-bconst/r_ijk)
!      ssobo3d(i,j,k)=sobo1d(r_ijk, velr, gradv, opalbar3d(i,j,k), t3d(i,j,k), xic, xnue0, epsl)
!      write(*,*) r_ijk, ssobo3d(i,j,k), velr, gradv, beta, bconst!, opalbar3d(i,j,k), t3d(i,j,k)!, xic, xnue0, epsl
!   endif
!enddo
!stop 'go on in sobo3d'

return
!
!***********************************************************************
!
!
write(*,*) '--------------------calculating 3d sobolev source function---------------------'
write(*,*)
!
!-----------------------create phi grids--------------------------------
!
do i=1, nphi
   nodes_phi(i) = (i-1)*2.d0*pi/(nphi-1)
enddo
!
call precalc_weight_simps(nodes_phi, nphi, weight_phi)
!
!-----------debug: for analytic beta-velocity law-----------------------
!
!vmin=1.d6
!vmax=2.d8
!beta=1.d0
!vth=vth_fiducial
!b=1.d0-(vmin/vmax)**(1.d0/beta)
!
!-----------------------------------------------------------------------
!
ssobo3d=0.d0
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
!do i=2, ndxmax-1
!   do j=ndymax/2+1, ndymax/2+1
!      do k=ndzmax/2+10, ndzmax/2+10
!
         if(imask_totreg3d(i,j,k).ne.0) then
!--------debug: for analytic beta-velocity law--------------------------
!            x_ii=x(i)
!            y_jj=y(j)
!            z_kk=z(k)
!            r_ijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
!            vdum=vmax*(1.d0-b/r_ijk)**beta
!            dvdum=vmax*beta*b*(1.d0-b/r_ijk)**(beta-1)
!            vdum=vdum/vth
!            dvdum=dvdum/vth
!            dvxdx2 = vdum/r_ijk - x_ii**2*vdum/r_ijk**3 + x_ii**2*dvdum/r_ijk**4
!            dvydy2 = vdum/r_ijk - y_jj**2*vdum/r_ijk**3 + y_jj**2*dvdum/r_ijk**4
!            dvzdz2 = vdum/r_ijk - z_kk**2*vdum/r_ijk**3 + z_kk**2*dvdum/r_ijk**4
!            dvydx2 = -x_ii*y_jj*vdum/r_ijk**3 + x_ii*y_jj*dvdum/r_ijk**4
!            dvzdx2 = -x_ii*z_kk*vdum/r_ijk**3 + x_ii*z_kk*dvdum/r_ijk**4
!            dvxdy2 = dvydx2
!            dvzdy2 = -y_jj*z_kk*vdum/r_ijk**3 + y_jj*z_kk*dvdum/r_ijk**4
!            dvxdz2 = dvzdx2
!            dvydz2 = dvzdy2
!
!----calculate derivatives of velocity components in all directions-----
!
            if(imask_totreg3d(i-1,j,k).ne.0.and.&
               imask_totreg3d(i+1,j,k).ne.0) then
!central differences
               vth=(vth3d(i+1,j,k)+vth3d(i-1,j,k))/2.d0
               dvxdx = (velx3d(i+1,j,k)-velx3d(i-1,j,k))/(x(i+1)-x(i-1)) * vth_fiducial/vth
               dvydx = (vely3d(i+1,j,k)-vely3d(i-1,j,k))/(x(i+1)-x(i-1)) * vth_fiducial/vth
               dvzdx = (velz3d(i+1,j,k)-velz3d(i-1,j,k))/(x(i+1)-x(i-1)) * vth_fiducial/vth
            elseif(imask_totreg3d(i-1,j,k).ne.0) then
!backward differences
               vth=(vth3d(i,j,k)+vth3d(i-1,j,k))/2.d0
               dvxdx = (velx3d(i,j,k)-velx3d(i-1,j,k))/(x(i)-x(i-1)) * vth_fiducial/vth
               dvydx = (vely3d(i,j,k)-vely3d(i-1,j,k))/(x(i)-x(i-1)) * vth_fiducial/vth
               dvzdx = (velz3d(i,j,k)-velz3d(i-1,j,k))/(x(i)-x(i-1)) * vth_fiducial/vth
            else
!forward differences
               vth=(vth3d(i+1,j,k)+vth3d(i,j,k))/2.d0
               dvxdx = (velx3d(i+1,j,k)-velx3d(i,j,k))/(x(i+1)-x(i)) * vth_fiducial/vth
               dvydx = (vely3d(i+1,j,k)-vely3d(i,j,k))/(x(i+1)-x(i)) * vth_fiducial/vth
               dvzdx = (velz3d(i+1,j,k)-velz3d(i,j,k))/(x(i+1)-x(i)) * vth_fiducial/vth
            endif
!
            if(imask_totreg3d(i,j-1,k).ne.0.and.&
               imask_totreg3d(i,j+1,k).ne.0) then
!central differences
               vth=(vth3d(i,j+1,k)+vth3d(i,j-1,k))/2.d0
               dvxdy = (velx3d(i,j+1,k)-velx3d(i,j-1,k))/(y(j+1)-y(j-1)) * vth_fiducial/vth
               dvydy = (vely3d(i,j+1,k)-vely3d(i,j-1,k))/(y(j+1)-y(j-1)) * vth_fiducial/vth
               dvzdy = (velz3d(i,j+1,k)-velz3d(i,j-1,k))/(y(j+1)-y(j-1)) * vth_fiducial/vth
            elseif(imask_totreg3d(i,j-1,k).ne.0) then
!backward differences
               vth=(vth3d(i,j,k)+vth3d(i,j-1,k))/2.d0
               dvxdy = (velx3d(i,j,k)-velx3d(i,j-1,k))/(y(j)-y(j-1)) * vth_fiducial/vth
               dvydy = (vely3d(i,j,k)-vely3d(i,j-1,k))/(y(j)-y(j-1)) * vth_fiducial/vth
               dvzdy = (velz3d(i,j,k)-velz3d(i,j-1,k))/(y(j)-y(j-1)) * vth_fiducial/vth
            else
!forward differences
               vth=(vth3d(i,j+1,k)+vth3d(i,j,k))/2.d0
               dvxdy = (velx3d(i,j+1,k)-velx3d(i,j,k))/(y(j+1)-y(j)) * vth_fiducial/vth
               dvydy = (vely3d(i,j+1,k)-vely3d(i,j,k))/(y(j+1)-y(j)) * vth_fiducial/vth
               dvzdy = (velz3d(i,j+1,k)-velz3d(i,j,k))/(y(j+1)-y(j)) * vth_fiducial/vth
            endif
!
            if(imask_totreg3d(i,j,k-1).ne.0.and.&
               imask_totreg3d(i,j,k+1).ne.0) then
!central differences
               vth=(vth3d(i,j,k+1)+vth3d(i,j,k-1))/2.d0
               dvxdz = (velx3d(i,j,k+1)-velx3d(i,j,k-1))/(z(k+1)-z(k-1)) * vth_fiducial/vth
               dvydz = (vely3d(i,j,k+1)-vely3d(i,j,k-1))/(z(k+1)-z(k-1)) * vth_fiducial/vth
               dvzdz = (velz3d(i,j,k+1)-velz3d(i,j,k-1))/(z(k+1)-z(k-1)) * vth_fiducial/vth
            elseif(imask_totreg3d(i,j,k-1).ne.0) then
!backward differences
               vth=(vth3d(i,j,k)+vth3d(i,j,k-1))/2.d0
               dvxdz = (velx3d(i,j,k)-velx3d(i,j,k-1))/(z(k)-z(k-1)) * vth_fiducial/vth
               dvydz = (vely3d(i,j,k)-vely3d(i,j,k-1))/(z(k)-z(k-1)) * vth_fiducial/vth
               dvzdz = (velz3d(i,j,k)-velz3d(i,j,k-1))/(z(k)-z(k-1)) * vth_fiducial/vth
            else
!forward differences
               vth=(vth3d(i,j,k+1)+vth3d(i,j,k))/2.d0
               dvxdz = (velx3d(i,j,k+1)-velx3d(i,j,k))/(z(k+1)-z(k)) * vth_fiducial/vth
               dvydz = (vely3d(i,j,k+1)-vely3d(i,j,k))/(z(k+1)-z(k)) * vth_fiducial/vth
               dvzdz = (velz3d(i,j,k+1)-velz3d(i,j,k))/(z(k+1)-z(k)) * vth_fiducial/vth
            endif
!
!
!get angles in spherical coordinates
            call get_angles_spc(x(i),y(j),z(k), theta_tilde, phi_tilde)
            costt=cos(theta_tilde)
            cospt=cos(phi_tilde)
            sintt=sin(theta_tilde)
            sinpt=sin(phi_tilde)
!
!at each point, integrate over complete solid angle, defined w.r.t. radial direction
            beta_sobolev=0.d0
            betac_sobolev=0.d0
!
!prepare mu-integration weights
            r_ijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
            mustar=sqrt(1.d0 - 1.d0/r_ijk**2)
            del=(mustar+1.d0)/(nmu_nc-1)
            do imu=1, nmu_nc
               nodes_mu_nc(imu)=-1.d0 + (imu-1)*del
            enddo
            del=(1.d0-mustar)/(nmu_c-1)
            do imu=1, nmu_c
               nodes_mu_c(imu)=mustar + del*(imu-1)
            enddo
            call precalc_weight_simps(nodes_mu_c, nmu_c, weight_mu_c)
            call precalc_weight_simps(nodes_mu_nc, nmu_nc, weight_mu_nc)
!
!non-core rays
            do imu=1, nmu_nc
               cost=nodes_mu_nc(imu)
               sint=sqrt(1.d0-cost**2)
               do iphi=1, nphi
!transformation of direction-vectors to 3d cartesian grid
                  cosp=cos(nodes_phi(iphi))
                  sinp=sin(nodes_phi(iphi))
                  n_x = cost*sintt*cospt + sint*cosp*costt*cospt - sint*sinp*sinpt
                  n_y = cost*sintt*sinpt + sint*cosp*cospt*sinpt + sint*sinp*cospt
                  n_z = cost*costt - sint*cosp*sintt
!derivative along direction
                  dvsds = dvxdx*n_x**2 + dvydy*n_y**2 + dvzdz*n_z**2 + &
                          (dvydx + dvxdy)*n_x*n_y + (dvzdx+dvxdz)*n_x*n_z + (dvzdy+dvydz)*n_y*n_z
!                  dvsds2 = dvxdx2*n_x**2 + dvydy2*n_y**2 + dvzdz2*n_z**2 + &
!                          (dvydx2 + dvxdy2)*n_x*n_y + (dvzdx2+dvxdz2)*n_x*n_z + (dvzdy2+dvydz2)*n_y*n_z
                  dvsds = abs(dvsds)
                  if(dvsds.eq.0.) dvsds=1.d-6
!sobolev solution for each direction, and integration
                  tau=opalbar3d(i,j,k)/dvsds
!                  if(tau.le.1.d-5) then
!                     beta_sobolev=beta_sobolev + weight_mu_nc(imu)*weight_phi(iphi)
!                  elseif(tau.ge.1.d5) then
!                     beta_sobolev=beta_sobolev + weight_mu_nc(imu)*weight_phi(iphi)/tau
!                  else
!                     beta_sobolev=beta_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_nc(imu)*weight_phi(iphi)
!                  endif
                  beta_sobolev=beta_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_nc(imu)*weight_phi(iphi)
               enddo
            enddo
!core rays
            do imu=1, nmu_c
               cost=nodes_mu_c(imu)
               sint=sqrt(1.d0-cost**2)
               do iphi=1, nphi
!transformation of direction-vectors to 3d cartesian grid
                  cosp=cos(nodes_phi(iphi))
                  sinp=sin(nodes_phi(iphi))
                  n_x = cost*sintt*cospt + sint*cosp*costt*cospt - sint*sinp*sinpt
                  n_y = cost*sintt*sinpt + sint*cosp*cospt*sinpt + sint*sinp*cospt
                  n_z = cost*costt - sint*cosp*sintt
!derivative along direction
                  dvsds = dvxdx*n_x**2 + dvydy*n_y**2 + dvzdz*n_z**2 + &
                          (dvydx + dvxdy)*n_x*n_y + (dvzdx+dvxdz)*n_x*n_z + (dvzdy+dvydz)*n_y*n_z
!                  dvsds2 = dvxdx2*n_x**2 + dvydy2*n_y**2 + dvzdz2*n_z**2 + &
!                          (dvydx2 + dvxdy2)*n_x*n_y + (dvzdx2+dvxdz2)*n_x*n_z + (dvzdy2+dvydz2)*n_y*n_z
                  dvsds = abs(dvsds)
                  if(dvsds.eq.0.) dvsds=1.d-6
!sobolev solution for each direction, and integration
                  tau=opalbar3d(i,j,k)/dvsds
!                  if(tau.le.1.d-5) then
!                     beta_sobolev=beta_sobolev + weight_mu_c(imu)*weight_phi(iphi)
!                     betac_sobolev=betac_sobolev + weight_mu_c(imu)*weight_phi(iphi)
!                  elseif(tau.ge.1.d5) then
!                     beta_sobolev=beta_sobolev + weight_mu_c(imu)*weight_phi(iphi)/tau
!                     betac_sobolev=betac_sobolev + weight_mu_c(imu)*weight_phi(iphi)/tau
!                  else
!                     beta_sobolev=beta_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_c(imu)*weight_phi(iphi)
!                     betac_sobolev=betac_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_c(imu)*weight_phi(iphi)
!                  endif
                  beta_sobolev=beta_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_c(imu)*weight_phi(iphi)
                  betac_sobolev=betac_sobolev + (1.d0 - exp(-tau))/tau * weight_mu_c(imu)*weight_phi(iphi)
               enddo
            enddo
!
!normalization of betas and calculation of source function
            beta_sobolev=beta_sobolev/4.d0/pi
            betac_sobolev=betac_sobolev/4.d0/pi
            bval = bnue(xnue0, t3d(i,j,k))
            ssobo3d(i,j,k) = ((1.d0 - epsl) * betac_sobolev * xic + epsl * bval) / (epsl + (1.d0 - epsl) * beta_sobolev)
!
         endif
      enddo
   enddo
enddo
!
!
end subroutine sobo3d
