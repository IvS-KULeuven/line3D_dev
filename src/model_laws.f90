function bvel(rad, vinf, bconst, beta)
!
! calculates velocities from beta velocity law in cgs
!           v(rad) = vinf*(1.-bconst/rad)^beta
!
!INPUT: vinf in cgs
!       rad  in r_star
!       beta, bconst
!OUTPUT: v(rad) in cgs
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rad, vinf, bconst, beta
real(dp) :: bvel
!
bvel=vinf*(1.d0-bconst/rad)**beta
!
end function bvel
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine bvel3d(vmin, vinf, beta, x, y, z, velx, vely, velz, gradv)
!
!calculates velocity components of beta velocity law
!for test (spherical symmetric) purposes
!
!INPUT: vinf, vmin in cgs
!       x, y, z coordinates in r_star
!       rad  in r_star
!       beta
!OUTPUT: velx, vely, velz in cgs
!        gradv in cgs
!
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vinf, beta, x, y, z
real(dp), intent(out) :: velx, vely, velz, gradv
!
! ... local scalars
real(dp) :: rad, velr, bconst
!
bconst=1.d0-(vmin/vinf)**(1.d0/beta)
rad=sqrt(x**2 + y**2 + z**2)
!
velr=vinf*(1.d0-bconst/rad)**beta
!
velx=velr*x/rad
vely=velr*y/rad
velz=velr*z/rad
!
gradv = vinf*beta*(1.d0-bconst/rad)**(beta-1)*bconst/rad**2
!
end subroutine bvel3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine vinf_bc(vmin, vrot, vesc, vcrit, beta, zeta, gamma, theta0, vinf, b)
!
!   calculates terminal velocity and b-factors for
!         bjorkman&casinelli (1993) model
!
!INPUT: vmin:   minimum velocity at stellar surface
!       vrot:   rotational velocity
!       vesc:   escape velocity
!       vcrit:  break up velocity
!       beta:   beta parameter (of beta velocity law)
!       zeta:   zeta parameter
!       gamma:  gamma parameter
!       theta0: co-latitude of plane
! 
!OUTPUT: vinf:   terminal velocity at theta0
!        b:      b-coefficient of radial velocity law at theta0
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vrot, vesc, vcrit, &
                        beta, zeta, gamma, theta0
real(dp), intent(out) :: vinf, b
!
! ... local scalars
!
vinf = zeta*vesc*(1.d0-sin(theta0)*vrot/vcrit)**gamma
b = 1.d0-(vmin/vinf)**(1.d0/beta)
!
end subroutine vinf_bc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function vthermal(vmicro, temp, matom)
!
!   calculates thermal velocity
!
!INPUT: vmicro: micro-turbulent velocity in cgs
!       temp: temperature in Kelvin
!       matom: atomic number
!

use prog_type
use fund_const, ONLY: cgs_kb, cgs_mp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmicro, temp
integer(i4b), intent(in) :: matom
real(dp) :: vthermal
!
vthermal = sqrt(2.d0*cgs_kb*temp/(matom*cgs_mp) + vmicro**2)
!
end function vthermal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function deldop(xnue0, vth)
!
!   calculates thermal doppler width
!
!INPUT: xnue0: transition frequency in 1/s
!       vth: thermal velocity in cm/s
!
use prog_type
use fund_const, ONLY: cgs_clight
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xnue0, vth
real(dp) :: deldop
!
deldop = xnue0 * vth / cgs_clight
!
end function deldop
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function calc_vmicro(vmin, vmax, vel)
!
!   calculates microt-turbulent velocity as a linear
!      function of absolute value of velocity
!
!INPUT:
!   vmin   minimum abs(velocity) = vmicro_min in cgs
!   vmax   maximum abs(velocity) = 10.*vmicro_max in cgs
!   vel    abs(velocity) at current position in cgs
!
USE prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vmax, vel
real(dp) :: calc_vmicro
!
! ... local scalars
real(dp) :: vmicro_min, vmicro_max
!
!set minimum and maximum microturbulent velocities
vmicro_min=vmin
vmicro_max=0.1d0*vmax
!
calc_vmicro = vmicro_min + (vmicro_min-vmicro_max)*(vel-vmin)/(vmin-vmax)
!
!
end function calc_vmicro
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function dilfac(r)  
!
!     dilution factor for a given radius r (measured in rstar)
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: r
!
! ... local scalars
real(dp) :: dilfac
!
dilfac = half*(one-sqrt(one-one/r**2))
!
!
return  
end function dilfac
!
!-----------------------------------------------------------------------
!----------------planck function in frequency space---------------------
!-----------------------------------------------------------------------
!
function bnue(xnue, t)  
!
!     planck function, nue in hertz, t in kelvin
!     bnue in erg per (cm**2 * sec * hertz)
!
use prog_type
use fund_const, only: cgs_clight, cgs_planck, cgs_kb
!
implicit none
!
! ... arguments
real(dp) :: xnue, t
!
! ... local scalars
real(dp) :: bnue

bnue = (2.d0*cgs_planck/(cgs_clight**2))*(xnue**3)/(exp(cgs_planck*xnue/(cgs_kb*t))-1.d0)
!
!
return  
end function bnue
!
!-----------------------------------------------------------------------
!
function bnue2(xnue, t, lint)  
!
!     planck function, nue in hertz, t in kelvin
!     bnue in erg per (cm**2 * sec * hertz)
!
use prog_type
use fund_const, only: cgs_clight, cgs_planck, cgs_kb, cgs_sb, pi, one
!
! ... arguments
logical, intent(in) :: lint
real(dp), intent(in) :: xnue, t
!
! ... local scalars
real(dp) :: bnue2

if(lint) then
!frequency integrated planck function
   bnue2=cgs_sb/pi * t**4
else
   bnue2 = (2.d0*cgs_planck/(cgs_clight**2))*(xnue**3)/(exp(cgs_planck*xnue/cgs_kb/t)-one)
endif
!
!
return  
end function bnue2

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine diffus(xnue, t_iim1, t_ii, t_iip1, r_iim1, r_ii, r_iip1, xic1, xic2)  
!
!   calculates planck function and its radial derivative db/dr
!                    for given input values at a point r_ii
!   note: calculation of temperature gradient depends crucially on used temperature points
!
! on input: xnue:  frequency at which everything is calculated
!           t_iim1:  temperature at point i-1 (corresponding to ghost point
!                    as diffus is calculated mainly on boundary)
!           t_ii:    temperature at point i     in K
!           t_iip1:  temperature at point i+1   in K
!           r_iim1:  radius at point i-1        in R_star
!           r_ii:    radius at point            in R_star
!           r_iip1:  radius at point i+1        in R_star
!
! on output: xic1:   planck function in cgs
!            xic2:   -db/dr in cgs/r_star

!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xnue, t_iim1, t_ii, t_iip1, &
                              r_iim1, r_ii, r_iip1
real(dp), intent(out) :: xic1, xic2
!
! ... local scalars
real(dp) :: dbdr, dbdt, dtdr
!
! ... local functions
real(dp) :: bnue
!
xic1 = bnue(xnue, t_ii)
dtdr = (t_iim1 - t_iip1)/(r_iip1-r_iim1)
dbdt = xic1 * cgs_planck*xnue / (1.d0-exp(-1.d0*cgs_planck*xnue/cgs_kb/t_ii)) / cgs_kb / t_ii**2
dbdr = dbdt * dtdr
xic2 = dbdr
!
end subroutine diffus
!
!-----------------------------------------------------------------------
!----------------mean molecular weight----------------------------------
!-----------------------------------------------------------------------
!
function mean_molecular_weight(ih,ihe,yhe)
!
!     calculates mean molecular weight (if contribution of metals besides helium
!              is negligible for both the electron density and the mass
!
!input: ih:   number free electrons per hydrogen atom
!       ihe:  number free electrons per helium atom
!       yhe:  helium abundance by number (N_He/N_H)
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: ih, ihe, yhe
real(dp) :: mean_molecular_weight
!
! ... local scalars
!
mean_molecular_weight = (one + four*yhe)/(one+ih + (one+ihe)*yhe)
!
!
return  
end function mean_molecular_weight
!
!-----------------------------------------------------------------------
!----------------isothermal sound speed---------------------------------
!-----------------------------------------------------------------------
!
function vsound(temp, mmw)
!
!     calculates isothermal sound speed (sqrt(p/rho)) for
!   given temperature (temp) and mean molecular weight (mmw)
!
!input: temp:   temperature in Kelvin
!       mmw:    mean molecular weight
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: temp, mmw
real(dp) :: vsound
!
! ... local scalars
!
vsound = sqrt(cgs_kb*temp/mmw/cgs_mp)
!
!
return  
end function vsound
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!y
subroutine calc_tlucy1d(sr, teff, nr, opac1d, r1d, t1d)
!
!calculate lucy temperature stratification
!on input: sr      stellar radius in cgs
!          teff    stellar effective temperature  
!          opac1d  (grey) opacity along a considered (radial) direction in 1/rstar
!          r1d     coordinates along the considered direction in rstar
!          nr      dimension of the 1d arrays  
!on output: t1d    temperature along r1d
use prog_type
use fund_const, only: zero, half, one, two, three, four
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nr
real(dp), intent(in) :: sr, teff
real(dp), dimension(nr), intent(in) :: opac1d, r1d
real(dp), dimension(nr), intent(out) :: t1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: tau, dtau, tau2, dtau2
!
! ... local arrays
real(dp), dimension(nr) :: x1d
!
! ... local functions
real(dp) :: dilfac
!
x1d = one/r1d
t1d = zero
!
tau=zero
tau2=zero
do i=nr-1, 1, -1
   dtau = half*(opac1d(i+1)+opac1d(i))*(x1d(i)-x1d(i+1))
   tau = tau+dtau
   dtau2 = half*(opac1d(i+1)+opac1d(i))*(r1d(i+1)-r1d(i))
   tau2 = tau2+dtau2
   t1d(i) = teff*(dilfac(r1d(i)) + three/four*tau)**0.25d0
!   write(*,*) r1d(i), x1d(i), dilfac(r1d(i)), opac1d(i), tau, dtau, tau2/3.d0
enddo
!outer-most point
t1d(nr)=t1d(nr-1)
!
!stop
!
!
!
end subroutine calc_tlucy1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_W(rad, mu, rho, temp, vel, vel_r, vel_theta, &
                       teff, v_esc, v_inf, rhoc_star)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                  FOR WIND UPFLOW COMPONENT
!
!INPUT: rad:   radius where values shall be calculated
!       mu:    co-latitude (w.r.t magnetic pole axis) where values shall be calculated
!       teff, v_esc, v_inf, rhoc_star: stellar parameter
!
!OUTPUT: density, temperature, velocity and velocity components
!        at those points
!
!NOTE: all quantities are defined in a carthesian coordinate system, where the z-axis
!      is aligned with the magnetic-pole-axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, teff, v_esc, v_inf, rhoc_star
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum1, dum2, dum3, dum4
!
!----------------------------density------------------------------------
!
dum1 = sqrt(rad-1.d0+mu*mu)*sqrt(1.d0+3.d0*mu*mu)
dum2 = (1.d0/rad)**(3.d0/2.d0)
dum3 = rad-1.d0
dum4 = 4.d0*rad-3.d0+3.d0*mu*mu
!
rho = 2.d0*(rhoc_star*v_esc/v_inf)*dum1*dum2/dum3/dum4
!
!--------------------------temperature----------------------------------
!
temp = TEFF
!
!----------------------------velocity-----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
if(mu.ge.0.d0) then
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
!
END SUBROUTINE COMPONENT_W
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_C(rad, mu, r_apex, rho, temp, vel, vel_r, vel_theta, &
                       teff, v_esc, rhoc_star, delta)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                  FOR COOL DOWNFLOW COMPONENT
!
!INPUT: rad:    radius where values shall be calculated
!       mu:     co-latitude where values shall be calculated
!       r_apex: radius of loop apex
!
!OUTPUT: rho:    density at the given point
!        temp:   temperature at the given pont
!        vel, vel_r, vel_theta: velocity components and absoulute velocity
!
!NOTE: all quantities are defined in carthensian coordinate system, 
!      where z-axis is aligned with magnetic pole axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, r_apex, teff, v_esc, rhoc_star, delta
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum1, dum2, dum3, dum4, mu_dum
!
!----------------------------density------------------------------------
!
mu_dum=mu
if(mu_dum.eq.0.d0) then
!avoid division by zero
   mu_dum=1.d-13
endif
!
dum1 = sqrt(rad-1.d0+mu_dum*mu_dum)*sqrt(1.d0+3.d0*mu_dum*mu_dum)
dum2 = 1.d0/rad/rad
dum3 = sqrt(mu_dum*mu_dum+delta*delta/rad/rad)
dum4 = 4.d0*rad-3.d0+3.d0*mu_dum*mu_dum
rho = 2.d0*rhoc_star*dum1*dum2/dum3/dum4
!
!--------------------------temperature----------------------------------
!
temp = TEFF
!
!----------------------------velocity-----------------------------------
!
!write(*,*) rad, r_apex
if(abs(rad-r_apex).lt.1.d-14) then
   vel = 0.d0
else
   vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
endif
!
if(mu.ge.0.d0) then
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
!
END SUBROUTINE COMPONENT_C
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_S(rad, mu, mu_star, mu_shock, rho, temp, vel, vel_r, vel_theta, &
                       teff, t_inf, v_esc, v_inf, rhoc_star, delta)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                FOR POST-SHOCK COMPONENT
!
!INPUT: rad: radius where values shall be calculated
!       mu:  co-latitude where values shall be calculated
!       mu_star: intersection point of photosphere with closed loop (attached to given point)
!       mu_shock: corresponding to the shock-radius on the closed loop
!
!OUTPUT: rho: density
!        temp: temperature
!        vel, vel_r, vel_theta: velocity-components and absolute velocity
!
!NOTE: all quantities are defined in carthensian coordinate system, 
!      where z-axis is aligned with magnetic pole axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, mu_star, mu_shock, teff, t_inf, v_esc, v_inf, rhoc_star, delta
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum0, dum1, dum2, dum3, dum4, dum_temp, rho_w
!
!--------------------------temperature----------------------------------
!
dum1 = 1.d0 - (1.d0-mu_star*mu_star)/(1.d0-mu_shock*mu_shock)
dum1 = dum1*dum1
dum2 = abs(mu - mu**3.d0 + 3.d0*mu**5.d0/5.d0 - mu**7.d0/7.d0)
dum3 = abs(mu_shock - mu_shock**3.d0 + 3.d0*mu_shock*5.d0/5.d0 - mu_shock**7.d0/7.d0)
dum4 = (dum2/dum3)**(1.d0/3.d0)
!
dum_temp= T_INF*dum1
temp = dum_temp*dum4
temp = max(temp, TEFF)
!
!----------------------------density------------------------------------
!
!calculate wind component at r_shock, mu_shock
dum0 = (1.d0-mu_shock*mu_shock)/(1.d0-mu_star*mu_star)
dum1 = sqrt(dum0-1.d0+mu_shock*mu_shock)*sqrt(1.d0+3.d0*mu_shock*mu_shock)
dum2 = (1.d0/dum0)**(3.d0/2.d0)
dum3 = dum0-1.d0
dum4 = 4.d0*dum0-3.d0+3.d0*mu*mu
rho_w = (rhoc_star*v_esc/v_inf)*dum1*dum2/dum3/dum4
!
!calculate density in post-shock region
rho = 4.d0 * rho_w * dum_temp / temp
!
!----------------------------velocity-----------------------------------
!
dum1 = (1.d0-(1.d0-mu_star*mu_star)/(1.d0-mu_shock*mu_shock))*V_INF/4.d0
dum2 = temp/dum_temp
dum3 = (1.d0-mu_shock*mu_shock)/(1.d0-mu*mu)
dum3 = dum3**3.d0
dum4 = sqrt(1.d0+3.d0*mu*mu)/sqrt(1.d0+3.d0*mu_shock*mu_shock)
vel = dum1*dum2*dum3*dum4
!
if(mu.ge.0.d0) then
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
END SUBROUTINE COMPONENT_S
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE VEL_ADM_WCS(rad, mu, v_inf, v_esc, ralfven, chi_inf, t_inf, teff, vel_r, vel_theta)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_inf    terminal velocity in cm/s
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!       chi_inf  cooling parameter
!       t_inf, teff
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: v_inf, v_esc, rad, mu, ralfven, chi_inf, teff, t_inf
REAL(DP) :: vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: mmu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
REAL(DP) :: dum1, dum2, vel, temp, temp_shock
!
!-----------------------------------------------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
!--------------------------wind upflow----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
if(mmu.ge.0.d0) then
   dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
   dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
else
   dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
   dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!------------------------cooled downflow--------------------------------
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
if(r_apex.lt.ralfven) then
     
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif

   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!-------------------------shock region----------------------------------
!
   mu_lower=0.d0
   mu_upper=mu_star
   CALL GET_MU_SHOCK(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
   r_shock=r_apex*(1.d0-mu_shock*mu_shock)

   if(rad.gt.r_shock) then

!temperature
      temp_shock = t_inf * (1.d0-1.d0/r_shock)**2
      temp = (abs(mmu - mmu**3 + 3.d0*mmu**5/5.d0 - mmu**7/7.d0) / &
              abs(mu_shock - mu_shock**3 + 3.d0*mu_shock*5/5.d0 - mu_shock**7/7.d0))**(1.d0/3.d0)
      temp = temp*temp_shock
      temp = max(temp, TEFF)

!velocity
      dum1 = (1.d0-1.d0/r_shock)*v_inf/4.d0 * temp/temp_shock
      vel = dum1 * (r_shock/rad)**3 * sqrt((1.d0+3.d0*mmu**2)/(1.d0+3.d0*mu_shock**2))
!
      if(mmu.ge.0.d0) then
         dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      else
         dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      endif
!
!      vel_r = vel*dum1
!      vel_theta = vel*dum2
   endif
!
endif


END SUBROUTINE VEL_ADM_WCS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE VEL_ADM_WCS_BVEL(rad, mu, v_inf, v_esc, ralfven, chi_inf, t_inf, teff, vel_r, vel_theta)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!---NOTE: IN WIND UPFLOW REGION, BETA-VELOCITY LAW
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_inf    terminal velocity in cm/s
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!       chi_inf  cooling parameter
!       t_inf, teff
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: v_inf, v_esc, rad, mu, ralfven, chi_inf, teff, t_inf
REAL(DP) :: vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: mmu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
REAL(DP) :: dum1, dum2, vel, temp, temp_shock
!
!-----------------------------------------------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
!--------------------------wind upflow----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
vel_r = vel
vel_theta = 0.d0
!
!------------------------cooled downflow--------------------------------
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
if(r_apex.lt.ralfven) then
     
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif

   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!-------------------------shock region----------------------------------
!
   mu_lower=0.d0
   mu_upper=mu_star
   CALL GET_MU_SHOCK(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
   r_shock=r_apex*(1.d0-mu_shock*mu_shock)

   if(rad.gt.r_shock) then

!temperature
      temp_shock = t_inf * (1.d0-1.d0/r_shock)**2
      temp = (abs(mmu - mmu**3 + 3.d0*mmu**5/5.d0 - mmu**7/7.d0) / &
              abs(mu_shock - mu_shock**3 + 3.d0*mu_shock*5/5.d0 - mu_shock**7/7.d0))**(1.d0/3.d0)
      temp = temp*temp_shock
      temp = max(temp, TEFF)

!velocity
      dum1 = (1.d0-1.d0/r_shock)*v_inf/4.d0 * temp/temp_shock
      vel = dum1 * (r_shock/rad)**3 * sqrt((1.d0+3.d0*mmu**2)/(1.d0+3.d0*mu_shock**2))
!
      if(mmu.ge.0.d0) then
         dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      else
         dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      endif
!
!      vel_r = vel*dum1
!      vel_theta = vel*dum2
   endif
!
endif


END SUBROUTINE VEL_ADM_WCS_BVEL
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE RHOVEL_ADM_C(rad, mu, v_esc, ralfven, rhoc_star, delta, vel_r, vel_theta, rho)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) ::  v_esc, rad, mu, ralfven, rhoc_star, delta
REAL(DP) :: vel_r, vel_theta, rho
!
! ... local scalars
REAL(DP) :: mmu, mu_star, r_apex
REAL(DP) :: dum1, dum2, dum3, dum4, vel
!
!------------------------cooled downflow--------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
!-----------------------------------------------------------------------
!
if(r_apex.ge.ralfven) then
   vel_r=0.d0
   vel_theta=0.d0
   rho=0.d0
   return
else
!
!----------------------velocity-----------------------------------------
!
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif
!
   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!----------------------density------------------------------------------
!
   dum1 = sqrt(rad-1.d0+mmu**2)*sqrt(1.d0+3.d0*mmu**2)
   dum2 = 1.d0/rad**2
   dum3 = sqrt(mmu**2+(delta/rad)**2)
   dum4 = 4.d0*rad-3.d0+3.d0*mmu**2
   rho = 2.d0*rhoc_star*dum1*dum2/dum3/dum4
!
endif


END SUBROUTINE RHOVEL_ADM_C
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE GET_MU_SHOCK(mu_star, chi_inf, xl, xu, mu_shock)
!
!-----------------------------------------------------------------------
!
!   CALCULATE MU_SHOCK (CORRESPONDING TO SHOCK-RADIUS) FOR GIVEN 
!                       r_apex, chi_inf
!
!   A TRANSCENDENTAL FUNCTION (SEE Owocki et al 2016, eq. 19)
!         NEEDS TO BE SOLVED
!
!   THIS IS DONE BY APPLYING THE REGULA-FALSI METHOD
!
!ON INPUT: mu_star: INTERSECTION POINT OF DIPOLE-LOOP WITH PHOTOSPHERE
!          chi_inf: COOLING PARAMETER
!          xl: LOWER GUESS OF MU_SHOCK (WILL BE DESTROYED)
!          xu: UPPER GUESS OF MU_SHOCK (WILL BE DESTROYED)
!
!ON OUTPUT: mu_shock: SOLUTION OF THE PROBLEM
!
!-----------------------------------------------------------------------
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: mu_star, chi_inf
REAL(DP) :: xl, xu, mu_shock
!
!
! ... local scalars
INTEGER(I4B), PARAMETER :: nit=100
INTEGER(I4B) :: i
REAL(DP), PARAMETER :: eps=1.d-14
REAL(DP) :: x_lower, x_upper, y_lower, y_upper, x_new, y_new
REAL(DP) :: swap
!
! ... local functions
REAL(DP) :: g_fct
!
!-----------------------------------------------------------------------
!
!swap lower and upper values, if necessary
if(xl.gt.xu) then
   swap=xl
   xl=xu
   xu=swap
endif
!
x_lower=xl
x_upper=xu
y_lower = g_fct(x_lower, mu_star, chi_inf)
y_upper = g_fct(x_upper, mu_star, chi_inf)
!
if(y_lower*y_upper .gt. 0.d0) then
!error if no interval is found
   STOP 'ERROR IN GET_MU_SHOCK: NO INTERVAL FOR WHICH f(xl) < f(xu) IS FOUND'
endif
!
!--now, when start values give unequal signs, make regula falsi method--
!
do i=1, nit
!
   x_new = x_lower - y_lower * (x_upper-x_lower)/(y_upper-y_lower)
   y_new = g_fct(x_new, mu_star, chi_inf)
!
   if(abs(y_new).le.eps) exit
!
   if(y_new*y_lower.lt.0.d0 ) then
      x_upper=x_new
      y_upper=y_new
   else
      x_lower=x_new
      y_lower=y_new
   endif
!
enddo
!
mu_shock=x_new
!
!
!
END SUBROUTINE GET_MU_SHOCK
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
FUNCTION G_FCT(x, mu_star, chi_inf)
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: x, mu_star, chi_inf
REAL(DP) :: g_fct
!
! ... local scalars
REAL(DP) :: term1, term2, term3, term4, dum1, dum2
!
dum1=1.d0-x*x
dum2=1.d0-mu_star*mu_star
!
term1 = x - x**3.d0 + 3.d0/5.d0 * x**5.d0 - x**7.d0 / 7.d0
term2 = chi_inf*(1.d0+3.d0*mu_star**2.d0)/6.d0/mu_star/dum2**2.d0
term3 = (1.d0-dum2/dum1)**4.d0
term4 = dum1**6.d0 / (1.d0+3.d0*x**2.d0)
!
g_fct = abs(term1) - term2*term3*term4
!
!
!
END FUNCTION G_FCT
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function calc_req(r_pole, m_star, v_rot)
!
!this function calculates the equatorial radius of a rotating star from
! which is approximated by an ellipsoid
!
!inputs:
!   r_pole: polar radius in r_sun
!   m_star: stellar mass in m_sun (without eddington factor!!!)
!   v_rot:  rotational velocity at equator in km/s
!
!outputs:
!   r_eq:   equatorial radius in units of r_pole
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu, xmsu, cgs_grav
!
implicit none
!
! ... arguments
real(dp), intent(in) :: r_pole, m_star, v_rot
real(dp) :: calc_req
!
! ... local scalars
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, r_eq_cgs, w0
!
r_pole_cgs = r_pole*rsu
m_star_cgs = m_star*xmsu
v_rot_cgs = v_rot*1.d5
!
w0 = v_rot_cgs**2*r_pole_cgs/2.d0/cgs_grav/m_star_cgs
!
r_eq_cgs = r_pole_cgs/(1.d0-w0)
calc_req = r_eq_cgs/r_pole_cgs
!
return
!
end function calc_req
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine surface_rot1(theta, r_pole, m_star, v_rot, r_eq, v_crit, omega, r_surface)
!
!       This function calculates the surface radius of a rotating star from
!       a roche model (e.g. Cranmer/Owocki 1995) for a given co-latitude
!
!input:
!       r_pole: polar radius in r_sun
!       m_star: stellar mass in m_sun (without eddington factor!!!)
!       v_rot:  rotational velocity at equator in km/s
!       theta:       co-latitude (for which radius is calculated)
!
!output:
!       r_eq:   equatorial radius in units of r_pole
!       v_crit: critical velocity in km/s
!       omega:  omega = omega_eq / omega_crit   (omega_eq, omega_crit the
!               angular velocities at equator and critical)
!       r_surface:   radius(theta) of stellar surface in units of r_pole
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu, xmsu, cgs_grav, pi
!
implicit none
!
! ... arguments
real(dp), intent(in) :: theta, r_pole, m_star, v_rot
real(dp), intent(out) :: r_eq, v_crit, omega, r_surface
!
! ... local scalars
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, r_eq_cgs, w0, v_crit_cgs
!
r_pole_cgs = r_pole*rsu
m_star_cgs = m_star*xmsu
v_rot_cgs = v_rot*1.d5
!
w0 = v_rot_cgs**2*r_pole_cgs/2.d0/cgs_grav/m_star_cgs
!
r_eq_cgs = r_pole_cgs/(1.d0-w0)
r_eq = r_eq_cgs/r_pole_cgs
!
v_crit_cgs = sqrt(2.d0*cgs_grav*m_star_cgs/3.d0/r_pole_cgs)
v_crit = v_crit_cgs/1.d5
!
omega = v_rot_cgs*1.5d0*r_pole_cgs/r_eq_cgs/v_crit_cgs
!
!calculate surface (in units of r_pole)
r_surface=3.d0/omega/sin(theta)*cos((pi+acos(omega*sin(theta)))/3.d0)
!
end subroutine surface_rot1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine model_adisc(rcyc, zcyc, rstar, mstar, alpha, mdot, mmw, temp, rho, velr, velphi)
!
!       This function calculates the structure of a standard alpha-disc model
!                following description by Frank-book
!Note two things here: 1. we include the mean-molecular weight
!                      2. there were some errors in the Frank-book (he has a smaller central temperature (also on wiki),
!                         maybe to account for forward and backward hemisphere. However then, the optical depth is not consistent.
!                         In my opinion, the Frank-formulation is simply wrong  
!
!input:
!       rcyc: cylindrical radius in rstar
!       zcyc: z-coordinate in rstar
!       rstar: in r_sun
!       mstar: in m_sun
!       mdot:  acrretion rate in msun/yr
!       alpha: alpha-prescription
!       mmw:   mean-molecular weight  
!
!output:
!       temp:   temperature
!       rho: density
!       velr,velphi:  velocity components
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rcyc, zcyc, rstar, mstar, alpha, mdot, mmw
real(dp), intent(out) :: temp, rho, velr, velphi
!
! ... local scalars
integer(i4b), parameter :: nz = 31
integer(i4b) :: i
real(dp) :: rstar_cgs, mstar_cgs, rcyc_10, rcyc_cgs, csound, teff
real(dp) :: kappaim1, kappai, rhoi, rhoim1, zi, zim1
real(dp) :: f_mid, rho_mid, velr_mid, temp_mid, tau_mid, scaleh_mid
real(dp) :: dtau, dz, tauz, zmin, zmax
real(dp) :: mdot_16, mdot_cgs

real(dp), parameter :: rho_exp1 = -7.d0/10.d0, rho_exp2 = 11.d0/20.d0, rho_exp3 = 5.d0/8.d0, rho_exp4 = -15.d0/8.d0, rho_exp5 = 11.d0/5.d0, rho_exp6=9.d0/8.d0, &
                       teff_exp1 = 0.d0, teff_exp2 = 1/4.d0, teff_exp3 = 1.d0/4.d0, teff_exp4 = -3.d0/4.d0, teff_exp5 = 1.d0, teff_exp6 = 0.d0, &
                       temp_exp1 = -1.d0/5.d0, temp_exp2 = 3.d0/10.d0, temp_exp3 = 1.d0/4.d0, temp_exp4 = -3.d0/4.d0, temp_exp5 = 6.d0/5.d0, temp_exp6 = 1.d0/4.d0, &
                       velr_exp1 = 4.d0/5.d0, velr_exp2 = 3.d0/10.d0, velr_exp3 = -1.d0/4.d0, velr_exp4 = -1.d0/4.d0, velr_exp5 = -14.d0/5.d0, &
                       tau_exp1 = -4.d0/5.d0, tau_exp2 = 1.d0/5.d0, tau_exp3 = zero, tau_exp4 = zero, tau_exp5 = 4.d0/5.d0, tau_exp6 = 1.d0, &
                       scaleh_exp1 = -1.d0/10.d0, scaleh_exp2 = 3.d0/20.d0, scaleh_exp3 = -3.d0/8.d0, scaleh_exp4 = 9.d0/8.d0, scaleh_exp5 = 3.d0/5.d0, scaleh_exp6=-3.d0/8.d0, &
                       rho_fac = 4.018d-8, temp_fac = 2.796d4, teff_fac=7.2707d3, tau_fac=291.d0, scaleh_fac=1.318d8
!_exp1 corresponding to alpha
!_exp2 corresponding to mdot
!_exp3 corresponding to mstar
!_exp4 corresponding to radius
!_exp5 corresponding to f-factor
!_exp6 corresponding to mean-molecular weight
!
! ... local functions
real(dp) :: vsound
!
if(rcyc.lt.one) stop 'error in model_adisc: rcyc < 1. not allowed'
if(zcyc.lt.zero) stop 'error in model_adisc: zcyc < 0. not allowed, use absolute value of zcyc'
!
!in normalized units
rstar_cgs = rstar*rsu
mstar_cgs = mstar*xmsu
rcyc_cgs = rcyc*rstar_cgs
rcyc_10 = rcyc_cgs/1.d10
mdot_cgs = mdot*xmsu/yr
mdot_16 = mdot_cgs/1.d16
!

!calculate midplane temperatures and densities
f_mid = (one-sqrt(one/rcyc))**0.25d0
scaleh_mid = scaleh_fac * alpha**scaleh_exp1 * mdot_16**scaleh_exp2 * mstar**scaleh_exp3 * rcyc_10**scaleh_exp4 * f_mid**scaleh_exp5 * mmw**scaleh_exp6
rho_mid = rho_fac * alpha**rho_exp1 * mdot_16**rho_exp2 * mstar**rho_exp3 * rcyc_10**rho_exp4 * f_mid**rho_exp5 * mmw**rho_exp6
temp_mid = temp_fac * alpha**temp_exp1 * mdot_16**temp_exp2 * mstar**temp_exp3 * rcyc_10**temp_exp4 * f_mid**temp_exp5 * mmw**temp_exp6
teff = teff_fac * alpha**teff_exp1 * mdot_16**teff_exp2 * mstar**teff_exp3 * rcyc_10**teff_exp4 * f_mid**teff_exp5 * mmw**teff_exp6
velr_mid = 2.7d4 * alpha**velr_exp1 * mdot_16**velr_exp2 * mstar**velr_exp3 * rcyc_10**velr_exp4 * f_mid**velr_exp5
tau_mid = tau_fac * alpha**tau_exp1 * mdot_16**tau_exp2 * mstar**tau_exp3 * rcyc_10**tau_exp4 * f_mid**tau_exp5 * mmw**tau_exp6
!
!calculate sound speed for mid-plane temperature
csound = vsound(temp_mid,mmw)
!
!calculate optical depth to given z-point
dtau = zero
zmax = max(2.d0*zcyc,10.d0*scaleh_mid/rstar_cgs)
zmin = zcyc
dz = (zmax-zmin)/float(nz-1)
zim1 = zmax
rhoim1 = rho_mid*exp(-(zim1*rstar_cgs/scaleh_mid)**2 / two)
kappaim1 = 5.d24*temp_mid**(-7.d0/2.d0) * rhoim1
!kappaim1 = 5.d24*temp_mid*rho_mid**2 * exp(-(zim1*rstar_cgs/scaleh_mid)**2)
tauz = 0.d0
do i = 2, nz
   zi = zim1-dz
   rhoi = rho_mid*exp(-(zi*rstar_cgs/scaleh_mid)**2 / two)
   kappai = 5.d24*temp_mid**(-7.d0/2.d0) * rhoi
!   kappai = 5.d24*temp_mid**(-7.d0/2.d0) *rho_mid**2 * exp(-(zi*rstar_cgs/scaleh_mid)**2)
   dtau = 0.5d0*(kappai*rhoi+kappaim1*rhoim1)*dz*rstar_cgs
!   dtau = 0.5d0*(kappai+kappaim1)*dz*rstar_cgs
   tauz = tauz+dtau
!    write(*,*) zi, rhoi, tauz, tau_mid, scaleh_mid, fdum, sqrt(csound**2 * rcyc_cgs**3/cgs_grav/mstar_cgs), mmw
   zim1=zi
   kappaim1=kappai
   rhoim1=rhoi
enddo
!
!since I am integrating here directly, and to be consistent with temperature estimates in center
!need to multiply tauz by 2/sqrt(pi), basically because:
!1. tau_mid = 5.d24*scaleh_mid*rho_mid^2*temp_mid^(-7/2) used to calculate the temperature
!2. when integrating directly (from Gaussian), we get:
!   tauz = 5.d24*rho_mid^2*temp_mid^(-7/2) * scaleh_mid * sqrt(pi)/2., where the last factor will give smaller optical depths than required to get correct temp_mid
tauz = tauz*two/sqrt(pi)
!
!------------------now, model at rcyc,zcyc--------------------
!
!effective temperature at rcyc
temp = (0.75d0*(tauz+two/three))**0.25d0 * teff
rho = rho_mid * exp(-0.5*(zcyc/scaleh_mid)**2)
velphi = sqrt(cgs_grav*mstar_cgs/rcyc_cgs)
velr = zero!velr_mid
!
!if(temp.ne.temp) then
!   write(*,*) temp, tauz, teff, rcyc, zcyc
!   stop 'error'
!endif
!
!write(*,*) rcyc*rstar_cgs/rsu, temp_mid, temp, teff, tauz, tau_mid
!stop
!
end subroutine model_adisc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function sline_depcoeff(xnue0, temp, b2, b3)
!
!calculated line source function for given departure coefficients
!of any given transition with frequency xnue0
!
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: b2, b3, temp, xnue0
real(dp) :: sline_depcoeff
!
! ... local scalars
real(dp) :: c1, c2
!
!2*h*nu^3/c^2
c1=two*cgs_planck*xnue0**3 / cgs_clight**2
!h*nu/k
c2=cgs_planck*xnue0/cgs_kb
!
sline_depcoeff = c1 / (b2*exp(c2/temp)/b3 - 1.d0)
!
end function sline_depcoeff
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function sline_petrenz(temp, b2, b3)
!
!-----------------according to petrenz & puls 1995----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: b2, b3, temp
real(dp) :: sline_petrenz
!
! ... local scalars
real(dp) :: c1, c2
!
!
!
!2*h*nu^3/c^2 (for h_alpha)
c1=1.404d-3
!h*nu/k (for h_alpha)
c2=2.192d4
!
sline_petrenz=c1 / (b2*exp(c2/temp)/b3 - 1.d0)
!
end function sline_petrenz  
