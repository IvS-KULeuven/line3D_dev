!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function sobo1d(r, velr, gradv, opalbar, t, xic, xnue0, epsl)
!
!---------calculate sobolev source function without continuum-----------
!-----------------for spherically symmetric outflows--------------------
!input: rad     radius in rstar
!       vel     velocity in vth
!       gradv   gradient of velocity along radial direction in vth/rstar
!       opalbar integrated line opacity in 1/rstar
!       temp    temperature in kelvin
!       xic     core-intensity
!       epsl    thermalization parameter
!output: sline
!
use prog_type
use mod_integ1d, only: precalc_weight_simps
!
implicit none
!
! ... arguments
real(dp), intent(in) :: r, velr, gradv, opalbar, t, xic, epsl, xnue0
real(dp) :: sobo1d
!
! ... local scalars
integer(i4b) :: i
integer(i4b), parameter :: nmu_c=31, nmu_nc=61
real(dp) :: del, mustar, tau, mu, qval, bval, beta, betac
!
! ... local arrays
real(dp), dimension(nmu_c) :: nodes_mu_c, weight_mu_c
real(dp), dimension(nmu_nc) :: nodes_mu_nc, weight_mu_nc
!
! ... local functions
real(dp) :: bnue
!
!------------------------create mu-grid---------------------------------
!
mustar=sqrt(1.d0 - 1.d0/r**2)
!
del=(mustar+1.d0)/(nmu_nc-1)
do i=1, nmu_nc
   nodes_mu_nc(i)=-1.d0 + del*(i-1)
enddo
!
del=(1.d0-mustar)/(nmu_c-1)
do i=1, nmu_c
   nodes_mu_c(i)=mustar + del*(i-1)
enddo
!
call precalc_weight_simps(nodes_mu_c, nmu_c, weight_mu_c)
call precalc_weight_simps(nodes_mu_nc, nmu_nc, weight_mu_nc)
!
!-----------------------------------------------------------------------
!
beta=0.d0
betac=0.d0
!
!write(*,*) r
do i=1, nmu_nc
   mu=nodes_mu_nc(i)
   qval=mu*mu*gradv + (1.d0 - mu*mu) * velr/r
   tau=opalbar/qval
   beta=beta + (1.d0 - exp(-tau))/tau * weight_mu_nc(i)
!   write(*,'(4es20.8)') mu, qval, tau, beta
enddo
!write(*,*)
!
do i=1, nmu_c
   mu=nodes_mu_c(i)
   qval=mu*mu*gradv + (1.d0 - mu*mu) * velr/r
   tau=opalbar/qval
   beta=beta + (1.d0 - exp(-tau))/tau * weight_mu_c(i)
   betac=betac + (1.d0 - exp(-tau))/tau * weight_mu_c(i)
enddo
!
beta=beta/2.d0
betac=betac/2.d0

!write(*,*) r, beta, betac

bval = bnue(xnue0, t)

sobo1d=((1.d0 - epsl) * betac * xic + epsl * bval) / (epsl + (1.d0 - epsl) * beta)
!
!
end function sobo1d
!
