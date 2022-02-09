!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark09
!
!calculate some timing properties, etc
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, &
                  sline3d, ssobo3d, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, &
                  opalbar3d, t3d
use bcondition, only: xic1
use freq, only: xnue0
use params_input, only: vmin, vmax, beta, eps_line, vth_fiducial
!
implicit none
!
! ... local scalars
integer(i4b) :: muindx, xobsindx
real(dp) :: ts_sobo, te_sobo, ts_linefvm, te_linefvm, ts_linesc, te_linesc
!
!-----------------------------------------------------------------------
!
ts_sobo=0.d0
te_sobo=0.d0
ts_linesc=0.d0
te_linesc=0.d0
ts_linefvm=0.d0
te_linefvm=0.d0
!
!-----------------timing of 3d sobolev approximation--------------------
!
call cpu_time(ts_sobo)
   call sobo3d(ndxmax, ndymax, ndzmax, x, y, z, imask_totreg3d, velx3d, vely3d, velz3d, vth3d, opalbar3d, t3d, &
            beta, vmin, vmax, vth_fiducial, xic1, xnue0, eps_line, ssobo3d)
call cpu_time(te_sobo)
!
!---------------timing of 2d line-calculation using sc and fvm----------
!-----------------(for a single direction and frequency)----------------
!
sline3d=ssobo3d
!
muindx=1
xobsindx=1
!
!
call cpu_time(ts_linesc)
   call fsc_line2d(muindx,xobsindx)
call cpu_time(te_linesc)
!
call cpu_time(ts_linefvm)
   call ffvm_line2d(muindx,xobsindx)
call cpu_time(te_linefvm)
!
!*****************************total timing******************************
!
write(*,'(a30, es20.8)') 'timing 1: -sobolev solution in 3d', te_sobo-ts_sobo
write(*,*)
write(*,'(a30)') 'timing 2: line solution in 2d for'
write(*,'(a30, es20.8)') '   a single direction and frequency'
write(*,'(a30, es20.8)') '   using sc method', te_linesc-ts_linesc
write(*,'(a30, es20.8)') '   using fvm method', te_linefvm-ts_linefvm
write(*,*)
!
!
end subroutine benchmark09
