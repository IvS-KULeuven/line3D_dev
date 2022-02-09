subroutine benchmark07_solution
!
!------calculating 2d continuum radiative transfer for spherically------
!----symmetric problems, with the theoretical solution given from-------
!---------------------------JOs 1d programs-----------------------------
!
use prog_type
use fund_const, only: pi, xmsu, sigmae, cgs_mp
use dimecr, only: n1d_cr, r1d_cr
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, ssobo3d, opalbar3d, velx3d, vely3d, velz3d, vth3d, t3d, imask_totreg3d
use bcondition, only: xic1
use mod_benchmark, only: ssobo1d_crx, ssobo1d_cry, ssobo1d_crz, velr1d_cr, opalbar1d_cr, t1d_cr, &
                         ssobo1d
use freq, only: xnue0
use params_input, only: eps_line
use params_input, only: teff, rstar, beta, vmin, vmax, hei, yhe, xmloss, eps_line, kline, &
                        kappa0, alpha, vth_fiducial
use params_stellar, only: sr
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: bconst, c1, c2, xmloss_cgs, rric
real(dp) :: r, velr, velx, vely, velz, gradv, rho, ne, chi_thomson, opalbar
!
! ... local arrays
real(dp), dimension(:), allocatable :: gradv1d_cr, vth1d_cr
!
! ... local logicals
!
! ... local functions
real(dp) :: bnue, sobo1d
!
!
!-----------------------------------------------------------------------
!
allocate(gradv1d_cr(n1d_cr))
allocate(vth1d_cr(n1d_cr))
allocate(t1d_cr(n1d_cr))
allocate(opalbar1d_cr(n1d_cr))
allocate(velr1d_cr(n1d_cr))
!
if(ndxmax.ne.ndzmax) stop 'error in benchmark07: ndxmax, ndymax, ndzmax need to be the same'
if(ndymax.ne.ndzmax) stop 'error in benchmark07: ndxmax, ndymax, ndzmax need to be the same'
allocate(ssobo1d_crx(n1d_cr))
allocate(ssobo1d_cry(n1d_cr))
allocate(ssobo1d_crz(n1d_cr))
allocate(ssobo1d(n1d_cr))
!
vth1d_cr=vth_fiducial!/2.d0
vth3d=vth_fiducial!/2.d0
!
t1d_cr=teff
t3d=teff
!
!---------------------calculate model atmosphere------------------------
!
xic1=bnue(xnue0,teff)
!
!b-factor for beta-velocity-law
bconst=1.d0-(vmin/vmax)**(1.D0/beta)
!
c1=(1.d0+4.d0*yhe)*cgs_mp
c2=(1.d0+hei*yhe)/c1
!
!mass-loss rate in cgs
xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!1d model
do i=1, n1d_cr
   r=r1d_cr(i)
   velr = 1.d5*vmax*(1.d0-bconst/r)**beta
   gradv = velr * bconst*beta/r**2/(1.d0-bconst/r)
   velr1d_cr(i) = velr
   gradv1d_cr(i) = gradv
enddo
!
!3d model
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(imask_totreg3d(i,j,k).eq.1) then
            r=sqrt(x(i)**2+y(j)**2+z(k)**2)
            velr = 1.d5*vmax*(1.d0-bconst/r)**beta
            velx = x(i)*velr/r
            vely = y(j)*velr/r
            velz = z(k)*velr/r
            rho = xmloss_cgs/4.d0/pi/r**2/velr/sr**2
            ne=c2*rho
            chi_thomson=sigmae*ne
            opalbar = kline*chi_thomson
!
            velx3d(i,j,k)=velx
            vely3d(i,j,k)=vely
            velz3d(i,j,k)=velz
!
            opalbar3d(i,j,k)=opalbar
         endif
      enddo
   enddo
enddo
!
do i=1, n1d_cr
   opalbar1d_cr(i) = opalbar3d(ndxmax/2+1,ndymax/2+1, ndzmax-n1d_cr-1+i)
enddo
!
!in correct units: velocities in vth, opacities in x-space, and in 1/rstar
velr1d_cr=velr1d_cr/vth1d_cr
gradv1d_cr=gradv1d_cr/vth1d_cr
velx3d=velx3d/vth_fiducial
vely3d=vely3d/vth_fiducial
velz3d=velz3d/vth_fiducial
!
opalbar3d=opalbar3d*sr
opalbar1d_cr=opalbar1d_cr*sr
!
!-----------------------------------------------------------------------
!
!3d sobolev solution
call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
do i=1, n1d_cr
   j=ndxmax-n1d_cr-1+i
   ssobo1d_crx(i) = ssobo3d(j,ndymax/2+1,ndzmax/2+1)  
   ssobo1d_cry(i) = ssobo3d(ndxmax/2+1,j,ndzmax/2+1)  
   ssobo1d_crz(i) = ssobo3d(ndxmax/2+1,ndymax/2+1,j)
enddo
!
!1d sobolev solution
do i=1, n1d_cr
   write(*,*) r1d_cr(i), opalbar1d_cr(i)
   ssobo1d(i) = sobo1d(r1d_cr(i), velr1d_cr(i), gradv1d_cr(i), opalbar1d_cr(i), t1d_cr(i), xic1, xnue0, eps_line)
enddo
!
!
!
write(*,'(10a20)') 'r[rstar]', 'velr[vth]', 'opalbar[1/rstar]', 'vth[cm/s]', 'gradv', 't[k]', &
                   'ssobo3d_crx*r^2/ic', 'ssobo3d_cry*r^2/ic', 'ssobo3d_crz*r^2/ic', 'ssobo1d*r^2/ic'
do i=1, n1d_cr
   rric=r1d_cr(i)**2/xic1
   write(*,'(10es20.8)') r1d_cr(i), velr1d_cr(i), opalbar1d_cr(i), vth1d_cr(i), gradv1d_cr(i), t1d_cr(i), &
                        ssobo1d_crx(i)*rric, ssobo1d_cry(i)*rric, ssobo1d_crz(i)*rric, ssobo1d(i)*rric
enddo
write(*,*)
!
end subroutine benchmark07_solution
