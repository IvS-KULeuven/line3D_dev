!
!***********************************************************************
!***********************************************************************
!
!                     CONTINUUM ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for continuum transport
!
!             u------------p----------d
!                 dels_u      dels_d
!
!      bezier interpolations/integrations
!
!   input: upwind intensity at u:       int_u
!          continuum opacity at u:      opac_u
!          continuum source-fct at u:   scont_u
!          continuum opacity at p:      opac_p
!          continuum source-fct at p:   scont_p
!          continuum opacity at d:      opac_d
!          continuum source-fct at d:   scont_d
!          path-length from u to p:     dels_u
!          path-length from p to d:     dels_d
!          coordinates of point u:      x_u, z_u
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo_u, alo_p, alo_d:               alo-coefficients for points u, p, d
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use mod_integ1d, only: integ1d_tau_ud, coeff_source2
!
implicit none
!
! ... arguments
real(dp) :: int_u, opac_p, scont_p, &
                        opac_u, scont_u, dels_u, &
                        opac_d, scont_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
!
! ... local scalars
integer(i4b) :: i, nd
integer(i8b) :: alo_u_int, alo_p_int, alo_d_int, abs_sc_int, norm_int
real(sp) :: alo_u_sp, alo_p_sp, alo_d_sp, abs_sc_sp
real(dp) :: delt_u, delt_d, delt_ud, rerr_contr
real(dp) :: s1, s2, s3, s4, s5, r1, r2, r3, r4, r5, scont2, scont3, scont4, &
            opac2, opac3, opac4, grad, dum, h, f1, f2, f3, f4, f5
real(dp) :: s_iim1, s_ii, opac_iim1, opac_ii, scont_iim1, scont_ii, int_iim1, int_ii, &
            ru, rp, r_ii
real(dp) :: delt_u2, delt_u3, int_sc3, int_sc2
real(dp) :: velr, rho, c1, c2, opac_u2
real(dp) :: tau_u, tau_p, tau_d, tau_iim1, tau_ii, tau_iip1, scont_iip1
real(dp) :: e0, e1, e2, a, b, c
real(dp) :: norm
!
! ... for debug: analytic opacities
real(dp) :: bvel, opac_thomson2
real(dp) :: bconst, xmloss_cgs, alpha_hamann, kappa_hamann, &
            x_p, z_p, x_d, z_d, &
            r_u, r_p, r_d, &
            vel_u, vel_p, vel_d, &
            rho_u, rho_p, rho_d

! ... local functions
real(dp) :: integral0, integral1, integral2, integral2b, integral3, integral3b, &
            integral_source0, integral_source1, integral_source2, integral_source3
real(dp) :: interpol_ypl, interpol_typ_quad3
!
!calculate delta-tau steps via monotonic bezier interpolation
call integ1d_tau_ud(abs(dels_u), abs(dels_d), opac_u, opac_p, opac_d, delt_u, delt_d)
!
abs_sc = exp(-delt_u)
if(delt_u.gt.20.) abs_sc=zero

!
!call coeff_source1(delt_u, alo_u, alo_p)
!call coeff_source4(delt_u, alo_u, alo_p)
!alo_d=zero
!
call coeff_source2(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !bezier interpolation
!call coeff_source3(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !warning: might cause osciallations
!call coeff_source3b(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !warning: might cause osciallations


alo_p=one-alo_u-alo_d-abs_sc
norm=alo_u+alo_p+alo_d+abs_sc


!   write(*,*) 'testa', norm, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d, dels_u, dels_d
!
!if(norm.gt.one) then
!   alo_u_sp = alo_u
!   alo_p_sp = alo_p
!   alo_d_sp = alo_d
!   abs_sc_sp = abs_sc
!
!   alo_u_int = nint(alo_u_sp*1.d13, i8b)
!   alo_p_int = nint(alo_p_sp*1.d13, i8b)
!   alo_d_int = nint(alo_d_sp*1.d13, i8b)
!   abs_sc_int = nint(abs_sc_sp*1.d13, i8b)
!!
!   norm_int=alo_u_int+alo_p_int+alo_d_int+abs_sc_int
!!
!   alo_u = dble(alo_u_int)/dble(norm_int)
!   alo_p = dble(alo_p_int)/dble(norm_int)
!   alo_d = dble(alo_d_int)/dble(norm_int)
!   abs_sc = dble(abs_sc_int)/dble(norm_int)
!   norm=alo_u+alo_p+alo_d+abs_sc

!   write(*,*) 'testa', norm, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d
!   if(norm.gt.one) stop 'error in fsc_cont: alo coefficients too large'
!endif
!
!write(*,*) 'test', opac_u, opac_p, opac_d, dels_u, dels_d, delt_u, delt_d
!
contr_sc = alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
!
int_sc = int_u*abs_sc + contr_sc


!write(*,*) 'fsc_cont', delt_u, delt_d, opac_u, opac_p, opac_d, norm, contr_sc

return
!
end subroutine fsc_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                          dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for continuum transport
!
!     linear interpolations/integrations
!
!             u------------p---------d
!                 dels_u      dels_d
!
!   input: upwind intensity at u:       int_u
!          continuum opacity at u:      opac_u
!          continuum source-fct at u:   scont_u
!          continuum opacity at p:      opac_p
!          continuum source-fct at p:   scont_p
!          continuum opacity at d:      opac_d
!          continuum source-fct at d:   scont_d
!          path-length from p to d:     dels_d
!          path-length from u to p:     dels_u
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo_u, alo_p, alo_d:               alo-coefficients for points u, p, d
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use mod_integ1d, only: integ1d_tau_ud, coeff_source1
!
implicit none
!
! ... arguments
real(dp), intent(in) :: int_u, opac_p, scont_p, &
                        opac_u, scont_u, dels_u, &
                        opac_d, scont_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
!
! ... local scalars
integer(i4b) :: i, nd
integer(i8b) :: alo_u_int, alo_p_int, abs_sc_int, norm_int
real(sp) :: alo_u_sp, alo_p_sp, abs_sc_sp
real(dp) :: delt_u, delt_d, norm
!
! ... local functions
real(dp) :: integral0, integral1, integral2, integral2b, integral3, integral3b, &
            integral_source0, integral_source1, integral_source2, integral_source3
!
!
!write(*,'(5es20.8)') dels_u, dels_d, opac_u, opac_p, opac_d
!calculate delta-tau steps via monotonic bezier interpolation
call integ1d_tau_ud(abs(dels_u), abs(dels_d), opac_u, opac_p, opac_d, delt_u, delt_d)
!
!write(*,*) delt_u, delt_d
!write(*,*)
abs_sc = exp(-delt_u)
if(delt_u.gt.20.) abs_sc=zero
!
call coeff_source1(delt_u, alo_u, alo_p)
!
!norm=alo_u+alo_p+abs_sc
!
!if(norm.gt.one) then
!   alo_u_sp = alo_u
!   alo_p_sp = alo_p
!   abs_sc_sp = abs_sc
!
!   alo_u_int = nint(alo_u_sp*1.d13, i8b)
!   alo_p_int = nint(alo_p_sp*1.d13, i8b)
!   abs_sc_int = nint(abs_sc_sp*1.d13, i8b)
!!
!   norm_int=alo_u_int+alo_p_int+abs_sc_int
!!
!   alo_u = dble(alo_u_int)/dble(norm_int)
!   alo_p = dble(alo_p_int)/dble(norm_int)
!   abs_sc = dble(abs_sc_int)/dble(norm_int)
!   norm=alo_u+alo_p+abs_sc
!   write(*,*) 'testa', alo_u+alo_p+alo_d+abs_sc, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d
!   if(norm.gt.one) stop 'error in fsc_cont: alo coefficients too large'!
!endif


contr_sc = alo_u*scont_u + alo_p*scont_p
int_sc = int_u*abs_sc + contr_sc
!
!
end subroutine fsc_cont_lin

!
!***********************************************************************
!***********************************************************************
!
!                        LINE ROUTINES
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, &
                        int_u, opalbar_u, opalbar_p, &
                        sline_u, sline_p, &
                        vel_u, vel_p, &
                        vth_u, vth_p, &
                        dels_u, &
                        abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line transport
!
!             u------------p
!                 dels_u
!
!      linear interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_integ1d, only: coeff_source1
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, deltax, xcmf_min, xcmf_max, vth_fiducial, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                        opalbar_u, sline_u, vel_u, vth_u, dels_u
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
!
! ... local scalars
integer(i4b) :: i, nd
real(dp) :: m, q, s_p, s_u
real(dp) :: del_xcmf, xcmf_u, xcmf_p
real(dp) :: s_iim1, s_ii, sline_iim1, sline_ii, xcmf_iim1, xcmf_ii, &
            opal_iim1, opal_ii, vel_iim1, vel_ii
real(dp) :: delta_mean, vth_mean, phinorm
real(dp) :: delt_u, delt_ii, expdelt_ii
real(dp) :: a_ii, b_ii, a_iim1, b_iim1, a, b, c
!
! ... for debug: analytic structure
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
!
!use a mean delta
vth_mean=(vth_u+vth_p)/2.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!lrefine=.false.
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
!----------------------------start point--------------------------------
!
   s_iim1 = s_u
!
!linear interpolation coefficients
   b_iim1=s_iim1/s_p
   a_iim1=1.d0-b_iim1
!
   sline_iim1 = a_iim1*sline_u + b_iim1*sline_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iim1 = (a_iim1*opalbar_u + b_iim1*opalbar_p)*phinorm
!
   expdelt_ii=1.d0
   delt_u=0.d0
   alo_u=0.d0
   alo_p=0.d0
!
!-----------------------------------------------------------------------
!
   do i=2, nd
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_ii = xcmf_u + (i-1)*del_xcmf
      s_ii = s_u + (xcmf_ii-xcmf_u)/m
      b_ii=s_ii/s_p
      a_ii=1.d0-b_ii
      sline_ii = a_ii*sline_u + b_ii*sline_p
      vel_ii = a_ii*vel_u + b_ii*vel_p
      call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm)
      opal_ii = (a_ii*opalbar_u+b_ii*opalbar_p)*phinorm
!***debug start: analytic opacity
!      x_ii = xu_debug + s_ii*nnx_debug
!      y_ii = yu_debug + s_ii*nny_debug
!      z_ii = zu_debug + s_ii*nnz_debug
!      call model_debug(x_ii, y_ii, z_ii, velx_ii, vely_ii, velz_ii, opalbar_ii, opac_ii)
!      write(*,'(i5, 20es20.8)') i, a_ii*opalbar_u+b_ii*opalbar_p, opalbar_ii, x_ii, y_ii, z_ii, xu_debug, yu_debug, zu_debug, s_ii, nnx_debug, nny_debug, nnz_debug
!***debug end
!
!perform radiative transfer
      delt_ii = (opal_ii+opal_iim1)/2.d0*(s_ii-s_iim1)
      delt_u = delt_u + delt_ii
!
      call coeff_source1(delt_ii, a, b)
      expdelt_ii = exp(-delt_ii)
!
      alo_u = alo_u*expdelt_ii + (a*a_iim1 + b*a_ii)
      alo_p = alo_p*expdelt_ii + (a*b_iim1 + b*b_ii)
!
!referesh variables at iim1
      s_iim1=s_ii
      sline_iim1 = sline_ii
      opal_iim1 = opal_ii
      a_iim1 = a_ii
      b_iim1 = b_ii
   enddo

   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u+alo_p*sline_p
   int_sc = abs_sc*int_u+contr_sc
!
else
!
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iim1 = opalbar_u * phinorm
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm)
   opal_ii = opalbar_p*phinorm
!perform radiative transfer
   delt_u = (opal_ii+opal_iim1)/2.d0*(s_p-s_u)
   abs_sc = exp(-delt_u)
!
   call coeff_source1(delt_u, alo_u, alo_p)
!
   contr_sc = alo_u*sline_u + alo_p*sline_p
   int_sc = int_u*abs_sc + contr_sc
!
endif
!
!
end subroutine fsc_line_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line_lin_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                                opalbar_u, opalbar_p, &
                                sline_u, sline_p, &
                                vel_u, vel_p, &
                                vth_u, vth_p, &
                                dels_u, &
                                abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                                xcmf1d, vel1d, opalbar1d, opal1d, vth1d, sline1d, profile1d, int1d, s1d, tau1d, nd)
!
!-----------------------------------------------------------------------
!
!   debugging routine for fsc_line1d_lin: output all physical quantities
!               on a refined grid to test the routine
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line transport
!
!             u------------p
!                 dels_u
!
!      linear interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_integ1d, only: coeff_source1
!
implicit none
!
! ... arguments
logical, intent(in) :: luse_refine
integer, intent(out) :: nd
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                        opalbar_u, sline_u, vel_u, vth_u, dels_u
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, sline1d, &
                                                      vth1d, opal1d, int1d, profile1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: m, q, s_p, s_u
real(dp) :: del_xcmf, xcmf_u, xcmf_p
real(dp) :: s_iim1, s_ii, sline_iim1, sline_ii, xcmf_iim1, xcmf_ii, &
            opal_iim1, opal_ii, vel_iim1, vel_ii
real(dp) :: delta_mean, vth_mean, phinorm_iim1, phinorm_ii
real(dp) :: delt_u, delt_ii, expdelt_ii
real(dp) :: a_ii, b_ii, a_iim1, b_iim1, a, b, c
!
! ... for debug: analytic structure
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(allocated(s1d)) deallocate(s1d)
if(allocated(tau1d)) deallocate(tau1d)
if(allocated(vel1d)) deallocate(vel1d)
if(allocated(vth1d)) deallocate(vth1d)
if(allocated(xcmf1d)) deallocate(xcmf1d)
if(allocated(profile1d)) deallocate(profile1d)
if(allocated(opalbar1d)) deallocate(opalbar1d)
if(allocated(opal1d)) deallocate(opal1d)
if(allocated(sline1d)) deallocate(sline1d)
if(allocated(int1d)) deallocate(int1d)
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
!
!use a mean delta
vth_mean=(vth_u+vth_p)/2.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(luse_refine.and.lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(sline1d(nd))
   allocate(int1d(nd))
!
!----------------------------start point--------------------------------
!
   s_iim1 = s_u
!
!linear interpolation coefficients
   b_iim1=s_iim1/s_p
   a_iim1=1.d0-b_iim1
!
   sline_iim1 = a_iim1*sline_u + b_iim1*sline_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = (a_iim1*opalbar_u + b_iim1*opalbar_p)*phinorm_iim1
!
   delt_u=0.d0
   alo_u=0.d0
   alo_p=0.d0
!
   s1d(1) = s_iim1
   tau1d(1) = 0.d0
   vel1d(1) = vel_u
   vth1d(1) = vth_mean
   xcmf1d(1) = xcmf_u
   profile1d(1) = phinorm_iim1
   opalbar1d(1) = opalbar_u
   opal1d(1) = opal_iim1
   sline1d(1) = sline_iim1
   int1d(1) = int_u
!
!-----------------------------------------------------------------------
!
   do i=2, nd
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_ii = xcmf_u + (i-1)*del_xcmf
      s_ii = s_u + (xcmf_ii-xcmf_u)/m
      b_ii=s_ii/s_p
      a_ii=1.d0-b_ii
      sline_ii = a_ii*sline_u + b_ii*sline_p
      vel_ii = a_ii*vel_u + b_ii*vel_p
      call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
      opal_ii = (a_ii*opalbar_u+b_ii*opalbar_p)*phinorm_ii
!
!perform radiative transfer
      delt_ii = (opal_ii+opal_iim1)/2.d0*(s_ii-s_iim1)
      delt_u = delt_u + delt_ii
!
      call coeff_source1(delt_ii, a, b)
      expdelt_ii = exp(-delt_ii)
!
      int1d(i) = int1d(i-1)*expdelt_ii + a*sline_iim1 + b*sline_ii
!
      alo_u = alo_u*expdelt_ii + (a*a_iim1 + b*a_ii)
      alo_p = alo_p*expdelt_ii + (a*b_iim1 + b*b_ii)
!
!referesh variables at iim1
      s_iim1=s_ii
      sline_iim1 = sline_ii
      opal_iim1 = opal_ii
      a_iim1 = a_ii
      b_iim1 = b_ii
!
      s1d(i) = s_ii
      tau1d(i) = delt_u
      vel1d(i) = vel_ii
      vth1d(i) = vth_mean
      xcmf1d(i) = xcmf_ii
      profile1d(i) = phinorm_ii
      opalbar1d(i) = a_ii*opalbar_u + b_ii*opalbar_p
      opal1d(i) = opal_ii
      sline1d(i) = sline_ii

   enddo

   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u+alo_p*sline_p
   int_sc = abs_sc*int_u+contr_sc
!
else
!
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = opalbar_u * phinorm_iim1
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii = opalbar_p*phinorm_ii
!perform radiative transfer
   delt_u = (opal_ii+opal_iim1)/2.d0*(s_p-s_u)
   abs_sc = exp(-delt_u)
!
   call coeff_source1(delt_u, alo_u, alo_p)
!
   contr_sc = alo_u*sline_u + alo_p*sline_p
   int_sc = int_u*abs_sc + contr_sc
!
   nd=2
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(sline1d(nd))
   allocate(int1d(nd))
   s1d = (/ s_u, s_p /)
   tau1d = (/ 0.d0, delt_u /)
   vel1d = (/ vel_u, vel_p /)
   vth1d = (/ vth_mean, vth_mean /)
   xcmf1d = (/ xcmf_u, xcmf_p /)
   profile1d = (/ phinorm_iim1, phinorm_ii /)
   opalbar1d = (/ opalbar_u, opalbar_p /)
   opal1d = (/ opal_iim1, opal_ii /)
   sline1d = (/ sline_u, sline_p /)
   int1d = (/ int_u, int_sc /)
!
endif
!
!
end subroutine fsc_line_lin_debug
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, &
                    int_u, opalbar_u, opalbar_p, opalbar_d, &
                    sline_u, sline_p, sline_d, &
                    vel_u, vel_p, vel_d, &
                    vth_u, vth_p, vth_d, &
                    dels_u, dels_d, &
                    abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line transport
!
!             u------------p----------d
!                 dels_u      dels_d
!
!      bezier interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!          frequency integrated opacity at d:   opalbar_d
!          continuum source-fct at d:           sline_d
!          projected velocity at d:             vel_d   (in vth_fiducial)
!          thermal velocity at d:               vth_d 
!          path-length from u to p:             dels_u
!          path-length from p to d:             dels_d
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p, alo_d
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_interp2d, only: wp_integ1d, wp_interp1d
use fund_const, only: spi
use mod_integ1d, only: coeff_source2, coeff_source2b
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, deltax, xcmf_min, xcmf_max, vth_fiducial, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                        opalbar_u, sline_u, vel_u, vth_u, dels_u, &
                        opalbar_d, sline_d, vel_d, vth_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
!
! ... local scalars
integer(i4b) :: i, nd
real(dp) :: m, q, alpha_opalbar, alpha_sline, alpha_u, alpha_d, wu_sline, wp_sline, wd_sline, &
            wu_opalbar, wp_opalbar, wd_opalbar, opalbar_c, opalbar_cu, opalbar_cd, &
            sline_c, sline_cu, sline_cd
real(dp) :: a, b, c, au, bu, cu, ad, bd, cd, a_iim2, a_iim1, a_ii, a_iip1, &
            b_iim2, b_iim1, b_ii, b_iip1, c_ii, &
            at_iim1, bt_iim1, ct_iim1, at_ii, bt_ii, ct_ii, &
            at_iip1, bt_iip1, ct_iip1, at_iim2, bt_iim2, ct_iim2, &
            atline_iim1, btline_iim1, ctline_iim1, atline_ii, btline_ii, ctline_ii, &
            atline_iim2, btline_iim2, ctline_iim2, atline_iip1, btline_iip1, ctline_iip1
real(dp) :: s_u, s_p, s_d, s_iim2, s_iim1, s_ii, s_iip1, ts_iim2, ts_iim1, ts_ii, ts_iip1, &
            xcmf_u, xcmf_p, xcmf_d, del_xcmf, &
            xcmf_iim2, xcmf_iim1, xcmf_ii, xcmf_iip1, dxcmfs, &
            sline_iim2, sline_iim1, sline_ii, sline_iip1, &
            opalbar_iim2, opalbar_iim1, opalbar_ii, opalbar_iip1, vel_iim2, vel_iim1, vel_ii, vel_iip1, &
            phinorm_iim2, phinorm_iim1, phinorm_ii, phinorm_iip1, &
            opal_iim2, opal_iim1, opal_ii, opal_iip1, opal_cu, opal_cd, opalprime_ii, &
            expxcmf_u, expxcmf_p, expxcmf_d, expxcmf_iim2, expxcmf_iim1, expxcmf_ii, expxcmf_iip1, &
            errf_u, errf_p, errf_d, errf_iim2, errf_iim1, errf_ii, errf_iip1
real(dp) :: vth_mean, delta_mean
real(dp) :: expdtau, expdtau_ii, dtau_u, dtau_d, dtau_iim1, dtau_ii, dtau_iip1

!
! ... for debug: analytic structure
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: bconst, xmloss_cgs, alpha_hamann, kappa_hamann, &
            r_ii, r_iim1, r_iip1, &
            x_iim1, x_ii, x_iip1, &
            z_iim1, z_ii, z_iip1, &
            vr_iim1, vr_ii, vr_iip1, &
            vx_iim1, vx_ii, vx_iip1, &
            vz_iim1, vz_ii, vz_iip1, &
            rho_iim1, rho_ii, rho_iip1
real(dp) :: x_u, z_u, nx, nz
real(dp) :: bvel, opalbar_model_hamann, interpol_ypl
real(dp) :: opal_u, opal_p, opal_d, phinorm
!
! ... local arrays
real(dp), dimension(:), allocatable :: vel1d, s1d, opalbar1d, sline1d, vth1d, opal1d
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
s_d=s_p+abs(dels_d)
!
!use a mean delta
vth_mean=(vth_u+vth_p+vth_d)/3.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
xcmf_d=(xobs-vel_d)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!lrefine=.false.
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)
!
   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf=(xcmf_p-xcmf_u)/(nd-1)
!
!---------------calculate control points and coefficients---------------
!-----------------------for interpolation along s-----------------------
!
!derivative weights for opacity
   alpha_opalbar = dels_d/s_d
!derivative weights for source function
!   alpha_sline = alpha_opalbar
   alpha_sline = max(wp_interp1d,alpha_opalbar)
!control-point weights from weighted derivatives
   wu_sline = alpha_sline/2.d0
   wp_sline = ((2.d0-alpha_sline)*dels_d + (1.d0-alpha_sline)*dels_u)/2.d0/dels_d
   wd_sline = (alpha_sline-1.d0)*dels_u/2.d0/dels_d
   wu_opalbar = alpha_opalbar/2.d0
   wp_opalbar = ((2.d0-alpha_opalbar)*dels_d + (1.d0-alpha_opalbar)*dels_u)/2.d0/dels_d
   wd_opalbar = (alpha_opalbar-1.d0)*dels_u/2.d0/dels_d
!(non-monotonic) control points
   opalbar_c = wu_opalbar*opalbar_u + wp_opalbar*opalbar_p + wd_opalbar*opalbar_d
   sline_c = wu_sline*sline_u + wp_sline*sline_p + wd_sline*sline_d
!ensure monotonicity
   call pointcl1d_mbez (opalbar_c, opalbar_u, opalbar_p)
   call coeffcl1d_mbez(sline_c, sline_u, sline_p, wu_sline, wp_sline, wd_sline)
!
!----------------------------start point--------------------------------
!
   alo_u = 0.d0
   alo_p = 0.d0
   alo_d = 0.d0
   dtau_u = 0.d0
!
!at point u (<=> iim1)
   xcmf_iim1 = xcmf_u
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim1 = (s_iim1-s_u)/dels_u
   a_iim1 = 1.d0-b_iim1
!(bezier) interpolation coefficients
   ts_iim1 = (s_iim1-s_u)/dels_u
   at_iim1 = (1.d0-ts_iim1)**2
   bt_iim1 = 2.d0*ts_iim1*(1.d0-ts_iim1)
   ct_iim1 = ts_iim1**2
   atline_iim1 = at_iim1 + bt_iim1*wu_sline
   btline_iim1 = ct_iim1 + bt_iim1*wp_sline
   ctline_iim1 = bt_iim1*wd_sline
!actual interpolation
   sline_iim1 = atline_iim1*sline_u + btline_iim1*sline_p + ctline_iim1*sline_d
   opalbar_iim1 = at_iim1*opalbar_u + bt_iim1*opalbar_c + ct_iim1*opalbar_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = opalbar_iim1*phinorm_iim1
!
!
!at point ii
   xcmf_ii = xcmf_u + del_xcmf
   s_ii = s_u + (xcmf_ii-xcmf_u)/m
!(linear) interpolation coefficients
   b_ii = (s_ii-s_u)/dels_u
   a_ii = 1.d0-b_ii
!(bezier) interpolation coefficients
   ts_ii = (s_ii-s_u)/dels_u
   at_ii = (1.d0-ts_ii)**2
   bt_ii = 2.d0*ts_ii*(1.d0-ts_ii)
   ct_ii = ts_ii**2
   atline_ii = at_ii + bt_ii*wu_sline
   btline_ii = ct_ii + bt_ii*wp_sline
   ctline_ii = bt_ii*wd_sline
!actual interpolation
   sline_ii = atline_ii*sline_u + btline_ii*sline_p + ctline_ii*sline_d
   opalbar_ii = at_ii*opalbar_u + bt_ii*opalbar_c + ct_ii*opalbar_p
   vel_ii = a_ii*vel_u + b_ii*vel_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii = opalbar_ii*phinorm_ii
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   dtau_ii = (opal_iim1+opal_ii)*(s_ii-s_iim1)/2.d0
!
!---------------------------intergrid points----------------------------
!
   do i=2, nd-1
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_iip1 = xcmf_u + i*del_xcmf
      s_iip1 = s_u + (xcmf_iip1-xcmf_u)/m
!
!(linear) interpolation coefficients
      b_iip1 = (s_iip1-s_u)/dels_u
      a_iip1 = 1.d0-b_iip1
!(bezier) interpolation coefficients
      ts_iip1 = (s_iip1-s_u)/dels_u
      at_iip1 = (1.d0-ts_iip1)**2
      bt_iip1 = 2.d0*ts_iip1*(1.d0-ts_iip1)
      ct_iip1 = ts_iip1**2
      atline_iip1 = at_iip1 + bt_iip1*wu_sline
      btline_iip1 = ct_iip1 + bt_iip1*wp_sline
      ctline_iip1 = bt_iip1*wd_sline
!actual interpolation
      sline_iip1 = atline_iip1*sline_u + btline_iip1*sline_p + ctline_iip1*sline_d
      opalbar_iip1 = at_iip1*opalbar_u + bt_iip1*opalbar_c + ct_iip1*opalbar_p
      vel_iip1 = a_iip1*vel_u + b_iip1*vel_p
      call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
      opal_iip1 = opalbar_iip1*phinorm_iip1
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
      dtau_iip1 = (opal_ii+opal_iip1)*(s_iip1-s_ii)/2.d0
!
!calculate coefficients for source function integration
      call coeff_source2(dtau_ii, dtau_iip1, sline_iim1, sline_ii, sline_iip1, a_ii, b_ii, c_ii)
!      alpha_ii = a_ii*at_iim1 + b_ii*at_ii + c_ii*at_iip1
!      beta_ii = a_ii*bt_iim1 + b_ii*bt_ii + c_ii*bt_iip1
!      gamma_ii = a_ii*ct_iim1 + b_ii*ct_ii + c_ii*ct_iip1

      expdtau_ii=exp(-dtau_ii)
      dtau_u = dtau_u + dtau_ii
      alo_u = alo_u*expdtau_ii + a_ii*atline_iim1 + b_ii*atline_ii + c_ii*atline_iip1 !+ alpha_ii
      alo_p = alo_p*expdtau_ii + a_ii*btline_iim1 + b_ii*btline_ii + c_ii*btline_iip1 !+ beta_ii
      alo_d = alo_d*expdtau_ii + a_ii*ctline_iim1 + b_ii*ctline_ii + c_ii*ctline_iip1 !+ gamma_ii
!
!referesh variables at iim1
      sline_iim1 = sline_ii
      atline_iim1 = atline_ii
      btline_iim1 = btline_ii
      ctline_iim1 = ctline_ii
      opal_iim1 = opal_ii
      s_iim1 = s_ii
!
!refresh variables at ii
      sline_ii = sline_iip1
      atline_ii = atline_iip1
      btline_ii = btline_iip1
      ctline_ii = ctline_iip1
      opal_ii = opal_iip1
      dtau_ii = dtau_iip1
      s_ii = s_iip1
!
   enddo
!
!-------------------------interval [n-1,n]------------------------------
!
!recalculate everything at at n-2
   xcmf_iim2 = xcmf_p - 2*del_xcmf
   s_iim2 = s_u + (xcmf_iim2-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim2 = (s_iim2-s_u)/dels_u
   a_iim2 = 1.d0-b_iim2
!(bezier) interpolation coefficients
   ts_iim2 = (s_iim2-s_u)/dels_u
   at_iim2 = (1.d0-ts_iim2)**2
   bt_iim2 = 2.d0*ts_iim2*(1.d0-ts_iim2)
   ct_iim2 = ts_iim2**2
   atline_iim2 = at_iim2 + bt_iim2*wu_sline
   btline_iim2 = ct_iim2 + bt_iim2*wp_sline
   ctline_iim2 = bt_iim2*wd_sline
!actual interpolation
   sline_iim2 = atline_iim2*sline_u + btline_iim2*sline_p + ctline_iim2*sline_d
   opalbar_iim2 = at_iim2*opalbar_u + bt_iim2*opalbar_c + ct_iim2*opalbar_p
   vel_iim2 = a_iim2*vel_u + b_iim2*vel_p
   call calc_phinorm(vel_iim2, vth_mean, vth_fiducial, xobs, phinorm_iim2)
   opal_iim2 = opalbar_iim2*phinorm_iim2
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   dtau_iim1 = (opal_iim2+opal_iim1)*(s_iim1-s_iim2)/2.d0
!
!calculate coefficients for source function integration
   call coeff_source2b(dtau_iim1, dtau_ii, sline_iim2, sline_iim1, sline_ii, a_ii, b_ii, c_ii)

   expdtau_ii=exp(-dtau_ii)
   dtau_u = dtau_u + dtau_ii
   alo_u = alo_u*expdtau_ii + a_ii*atline_iim2 + b_ii*atline_iim1 + c_ii*atline_ii
   alo_p = alo_p*expdtau_ii + a_ii*btline_iim2 + b_ii*btline_iim1 + c_ii*btline_ii
   alo_d = alo_d*expdtau_ii + a_ii*ctline_iim2 + b_ii*ctline_iim1 + c_ii*ctline_ii
!
   abs_sc = exp(-dtau_u)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
   int_sc = int_u*abs_sc + contr_sc
!
!
else
!
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iim1 = opalbar_u * phinorm
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm)
   opal_ii = opalbar_p*phinorm
   call calc_phinorm(vel_d, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iip1 = opalbar_d*phinorm
!
!calculate delta-tau steps (assuming that opal is 2nd order bezier function in s)
   opalprime_ii = (opal_ii-opal_iim1)*dels_d/dels_u/(dels_d+dels_u) + (opal_iip1-opal_ii)*dels_u/dels_d/(dels_u+dels_d)
   opal_cu = opal_ii - dels_u*opalprime_ii/2.d0
   opal_cd = opal_ii + dels_d*opalprime_ii/2.d0
   call pointcl1d_mbez (opal_cu, opal_iim1, opal_ii)
   call pointcl1d_mbez (opal_cd, opal_ii, opal_iip1)
   dtau_u = (opal_iim1+opal_cu+opal_ii)*dels_u/3.d0
   dtau_d = (opal_ii+opal_cd+opal_iip1)*dels_d/3.d0
!   dtau_u = (opal_iim1+opal_ii)*dels_u/2.d0
!   dtau_d = (opal_iip1+opal_ii)*dels_d/2.d0
!WARNING: what happens if dtau_u->, dtau_d large, i.e. when resonance zone in downwind interval???
!
!perform radiative transfer
   abs_sc = exp(-dtau_u)
!
!   call coeff_source1(dtau_u, alo_u, alo_p)
!   alo_d=0.d0
   call coeff_source2(dtau_u, dtau_d, sline_u, sline_p, sline_d, alo_u, alo_p, alo_d)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
   int_sc = int_u*abs_sc + contr_sc
!
!
endif
!
!
!
end subroutine fsc_line
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                            opalbar_u, opalbar_p, opalbar_d, &
                            sline_u, sline_p, sline_d, &
                            vel_u, vel_p, vel_d, &
                            vth_u, vth_p, vth_d, &
                            dels_u, dels_d, &
                            abs_sc, contr_sc, int_sc, &
                            alo_u, alo_p, alo_d, &
                            luse_refine, xcmf1d, vel1d, opalbar1d, opal1d, vth1d, sline1d, profile1d, int1d, s1d, tau1d, nd)
!
!
!-----------------------------------------------------------------------

!   debugging routine for fsc_line1d: output all physical quantities
!              on a refined grid to test the routine
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line transport
!
!             u------------p----------d
!                 dels_u      dels_d
!
!      bezier interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!          frequency integrated opacity at d:   opalbar_d
!          continuum source-fct at d:           sline_d
!          projected velocity at d:             vel_d   (in vth_fiducial)
!          thermal velocity at d:               vth_d 
!          path-length from u to p:             dels_u
!          path-length from p to d:             dels_d
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p, alo_d
!
!additional output for debug routine: refined arrays + dimensions
!           nd
!           vel1d, opalbar1d, opal1d, vth1d, sline1d, int1d, s1d, profile1d, tau1d
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_interp2d, only: wp_integ1d, wp_interp1d
use fund_const, only: spi
use mod_integ1d, only: coeff_source2, coeff_source2b
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                        opalbar_u, sline_u, vel_u, vth_u, dels_u, &
                        opalbar_d, sline_d, vel_d, vth_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, sline1d, vth1d, opal1d, int1d, profile1d
logical, intent(in) :: luse_refine
!
! ... local scalars
integer(i4b) :: i, nd
real(dp) :: s_u, s_p, s_d, xcmf_u, xcmf_p, xcmf_d
real(dp) :: vth_mean, delta_mean, del_xcmf, m
real(dp) :: alpha_opalbar, wu_opalbar, wp_opalbar, wd_opalbar, opalbar_c, &   !inteprolation coefficients etc.
            alpha_sline, wu_sline, wp_sline, wd_sline, sline_c, &
            ts_iim2, ts_iim1, ts_ii, ts_iip1, &
            a_iim2, b_iim2, c_iim2, a_iim1, b_iim1, a_ii, b_ii, c_ii, a_iip1, b_iip1,  &
            at_iim2, bt_iim2, ct_iim2, at_iim1, bt_iim1, ct_iim1, &
            at_ii, bt_ii, ct_ii, at_iip1, bt_iip1, ct_iip1, &
            atline_iim2, btline_iim2, ctline_iim2, atline_iim1, btline_iim1, ctline_iim1, &
            atline_ii, btline_ii, ctline_ii, atline_iip1, btline_iip1, ctline_iip1
real(dp) :: s_iim2, s_iim1, s_ii, s_iip1, &
            vel_iim2, vel_iim1, vel_ii, vel_iip1, &
            xcmf_iim2, xcmf_iim1, xcmf_ii, xcmf_iip1, &
            opalbar_iim2, opalbar_iim1, opalbar_ii, opalbar_iip1, &
            opal_iim2, opal_iim1, opal_ii, opal_iip1, &
            phinorm_iim2, phinorm_iim1, phinorm_ii, phinorm_iip1, &
            sline_iim2, sline_iim1, sline_ii, sline_iip1, &
            delt_iim1, delt_ii, delt_iip1, delt_u, delt_d, expdtau_ii
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(allocated(xcmf1d)) deallocate(xcmf1d)
if(allocated(vel1d)) deallocate(vel1d)
if(allocated(s1d)) deallocate(s1d)
if(allocated(tau1d)) deallocate(tau1d)
if(allocated(opalbar1d)) deallocate(opalbar1d)
if(allocated(sline1d)) deallocate(sline1d)
if(allocated(vth1d)) deallocate(vth1d)
if(allocated(opal1d)) deallocate(opal1d)
if(allocated(profile1d)) deallocate(profile1d)
if(allocated(int1d)) deallocate(int1d)
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
s_d=s_p+abs(dels_d)
!
!use a mean delta
vth_mean=(vth_u+vth_p+vth_d)/3.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
xcmf_d=(xobs-vel_d)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(luse_refine.and.lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)
!
   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf=(xcmf_p-xcmf_u)/(nd-1)

   allocate(vel1d(nd))
   allocate(xcmf1d(nd))
   allocate(tau1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(sline1d(nd))
   allocate(vth1d(nd))
   allocate(profile1d(nd))
   allocate(int1d(nd))
   allocate(s1d(nd))
!
!---------------calculate control points and coefficients---------------
!-----------------------for interpolation along s-----------------------
!
!derivative weights for opacity
   alpha_opalbar = dels_d/s_d
!derivative weights for source function
   alpha_sline = max(wp_interp1d,alpha_opalbar)
!   alpha_sline=alpha_opalbar
!control-point weights from weighted derivatives
   wu_sline = alpha_sline/2.d0
   wp_sline = ((2.d0-alpha_sline)*dels_d + (1.d0-alpha_sline)*dels_u)/2.d0/dels_d
   wd_sline = (alpha_sline-1.d0)*dels_u/2.d0/dels_d
   wu_opalbar = alpha_opalbar/2.d0
   wp_opalbar = ((2.d0-alpha_opalbar)*dels_d + (1.d0-alpha_opalbar)*dels_u)/2.d0/dels_d
   wd_opalbar = (alpha_opalbar-1.d0)*dels_u/2.d0/dels_d
!(non-monotonic) control points
   opalbar_c = wu_opalbar*opalbar_u + wp_opalbar*opalbar_p + wd_opalbar*opalbar_d
   sline_c = wu_sline*sline_u + wp_sline*sline_p + wd_sline*sline_d
!ensure monotonicity
   call pointcl1d_mbez (opalbar_c, opalbar_u, opalbar_p)
   call coeffcl1d_mbez(sline_c, sline_u, sline_p, wu_sline, wp_sline, wd_sline)
!
!----------------------------start point--------------------------------
!
   alo_u = 0.d0
   alo_p = 0.d0
   alo_d = 0.d0
   delt_u = 0.d0
!
!at point u (<=> iim1)
   xcmf_iim1 = xcmf_u
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim1 = (s_iim1-s_u)/dels_u
   a_iim1 = 1.d0-b_iim1
!(bezier) interpolation coefficients
   ts_iim1 = (s_iim1-s_u)/dels_u
   at_iim1 = (1.d0-ts_iim1)**2
   bt_iim1 = 2.d0*ts_iim1*(1.d0-ts_iim1)
   ct_iim1 = ts_iim1**2
   atline_iim1 = at_iim1 + bt_iim1*wu_sline
   btline_iim1 = ct_iim1 + bt_iim1*wp_sline
   ctline_iim1 = bt_iim1*wd_sline
!actual interpolation
   sline_iim1 = atline_iim1*sline_u + btline_iim1*sline_p + ctline_iim1*sline_d
   opalbar_iim1 = at_iim1*opalbar_u + bt_iim1*opalbar_c + ct_iim1*opalbar_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = opalbar_iim1*phinorm_iim1
!
   s1d(1) = s_iim1
   tau1d(1) = 0.d0
   vel1d(1) = vel_iim1
   vth1d(1) = vth_mean
   xcmf1d(1) = xcmf_iim1
   profile1d(1) = phinorm_iim1
   opalbar1d(1) = opalbar_iim1
   opal1d(1) = opal_iim1
   sline1d(1) = sline_iim1
   int1d(1) = int_u
!
!
!
!at point ii
   xcmf_ii = xcmf_u + del_xcmf
   s_ii = s_u + (xcmf_ii-xcmf_u)/m
!(linear) interpolation coefficients
   b_ii = (s_ii-s_u)/dels_u
   a_ii = 1.d0-b_ii
!(bezier) interpolation coefficients
   ts_ii = (s_ii-s_u)/dels_u
   at_ii = (1.d0-ts_ii)**2
   bt_ii = 2.d0*ts_ii*(1.d0-ts_ii)
   ct_ii = ts_ii**2
   atline_ii = at_ii + bt_ii*wu_sline
   btline_ii = ct_ii + bt_ii*wp_sline
   ctline_ii = bt_ii*wd_sline
!actual interpolation
   sline_ii = atline_ii*sline_u + btline_ii*sline_p + ctline_ii*sline_d
   opalbar_ii = at_ii*opalbar_u + bt_ii*opalbar_c + ct_ii*opalbar_p
   vel_ii = a_ii*vel_u + b_ii*vel_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii = opalbar_ii*phinorm_ii
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   delt_ii = (opal_iim1+opal_ii)*(s_ii-s_iim1)/2.d0
!
!---------------------------intergrid points----------------------------
!
   do i=2, nd-1
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_iip1 = xcmf_u + i*del_xcmf
      s_iip1 = s_u + (xcmf_iip1-xcmf_u)/m
!
!(linear) interpolation coefficients
      b_iip1 = (s_iip1-s_u)/dels_u
      a_iip1 = 1.d0-b_iip1
!(bezier) interpolation coefficients
      ts_iip1 = (s_iip1-s_u)/dels_u
      at_iip1 = (1.d0-ts_iip1)**2
      bt_iip1 = 2.d0*ts_iip1*(1.d0-ts_iip1)
      ct_iip1 = ts_iip1**2
      atline_iip1 = at_iip1 + bt_iip1*wu_sline
      btline_iip1 = ct_iip1 + bt_iip1*wp_sline
      ctline_iip1 = bt_iip1*wd_sline
!actual interpolation
      sline_iip1 = atline_iip1*sline_u + btline_iip1*sline_p + ctline_iip1*sline_d
      opalbar_iip1 = at_iip1*opalbar_u + bt_iip1*opalbar_c + ct_iip1*opalbar_p
      vel_iip1 = a_iip1*vel_u + b_iip1*vel_p
      call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
      opal_iip1 = opalbar_iip1*phinorm_iip1
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
      delt_iip1 = (opal_ii+opal_iip1)*(s_iip1-s_ii)/2.d0
!
!calculate coefficients for source function integration
      call coeff_source2(delt_ii, delt_iip1, sline_iim1, sline_ii, sline_iip1, a_ii, b_ii, c_ii)
!      alpha_ii = a_ii*at_iim1 + b_ii*at_ii + c_ii*at_iip1
!      beta_ii = a_ii*bt_iim1 + b_ii*bt_ii + c_ii*bt_iip1
!      gamma_ii = a_ii*ct_iim1 + b_ii*ct_ii + c_ii*ct_iip1

      expdtau_ii=exp(-delt_ii)
      delt_u = delt_u + delt_ii
      alo_u = alo_u*expdtau_ii + a_ii*atline_iim1 + b_ii*atline_ii + c_ii*atline_iip1 !+ alpha_ii
      alo_p = alo_p*expdtau_ii + a_ii*btline_iim1 + b_ii*btline_ii + c_ii*btline_iip1 !+ beta_ii
      alo_d = alo_d*expdtau_ii + a_ii*ctline_iim1 + b_ii*ctline_ii + c_ii*ctline_iip1 !+ gamma_ii
!
      s1d(i) = s_ii
      tau1d(i) = delt_u
      vel1d(i) = vel_ii
      vth1d(i) = vth_mean
      xcmf1d(i) = xcmf_ii
      profile1d(i) = phinorm_ii
      opalbar1d(i) = opalbar_ii
      opal1d(i) = opal_ii
      sline1d(i) = sline_ii
      int1d(i) = int1d(i-1)*expdtau_ii + a_ii*sline_iim1 + b_ii*sline_ii + c_ii*sline_iip1
!
!referesh variables at iim1
      sline_iim1 = sline_ii
      atline_iim1 = atline_ii
      btline_iim1 = btline_ii
      ctline_iim1 = ctline_ii
      opal_iim1 = opal_ii
      s_iim1 = s_ii
!
!refresh variables at ii
      sline_ii = sline_iip1
      atline_ii = atline_iip1
      btline_ii = btline_iip1
      ctline_ii = ctline_iip1
      opal_ii = opal_iip1
      delt_ii = delt_iip1
      s_ii = s_iip1
      xcmf_ii = xcmf_iip1
      phinorm_ii = phinorm_iip1
      opalbar_ii = opalbar_iip1
      vel_ii = vel_iip1
!
   enddo
!
!-------------------------interval [n-1,n]------------------------------
!
!recalculate everything at at n-2
   xcmf_iim2 = xcmf_p - 2*del_xcmf
   s_iim2 = s_u + (xcmf_iim2-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim2 = (s_iim2-s_u)/dels_u
   a_iim2 = 1.d0-b_iim2
!(bezier) interpolation coefficients
   ts_iim2 = (s_iim2-s_u)/dels_u
   at_iim2 = (1.d0-ts_iim2)**2
   bt_iim2 = 2.d0*ts_iim2*(1.d0-ts_iim2)
   ct_iim2 = ts_iim2**2
   atline_iim2 = at_iim2 + bt_iim2*wu_sline
   btline_iim2 = ct_iim2 + bt_iim2*wp_sline
   ctline_iim2 = bt_iim2*wd_sline
!actual interpolation
   sline_iim2 = atline_iim2*sline_u + btline_iim2*sline_p + ctline_iim2*sline_d
   opalbar_iim2 = at_iim2*opalbar_u + bt_iim2*opalbar_c + ct_iim2*opalbar_p
   vel_iim2 = a_iim2*vel_u + b_iim2*vel_p
   call calc_phinorm(vel_iim2, vth_mean, vth_fiducial, xobs, phinorm_iim2)
   opal_iim2 = opalbar_iim2*phinorm_iim2
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   delt_iim1 = (opal_iim2+opal_iim1)*(s_iim1-s_iim2)/2.d0
!
!calculate coefficients for source function integration
   call coeff_source2b(delt_iim1, delt_ii, sline_iim2, sline_iim1, sline_ii, a_ii, b_ii, c_ii)

   expdtau_ii=exp(-delt_ii)
   delt_u = delt_u + delt_ii
   alo_u = alo_u*expdtau_ii + a_ii*atline_iim2 + b_ii*atline_iim1 + c_ii*atline_ii
   alo_p = alo_p*expdtau_ii + a_ii*btline_iim2 + b_ii*btline_iim1 + c_ii*btline_ii
   alo_d = alo_d*expdtau_ii + a_ii*ctline_iim2 + b_ii*ctline_iim1 + c_ii*ctline_ii
!
   i=nd
   s1d(i) = s_ii
   tau1d(i) = delt_u
   vel1d(i) = vel_ii
   vth1d(i) = vth_mean
   xcmf1d(i) = xcmf_ii
   profile1d(i) = phinorm_ii
   opalbar1d(i) = opalbar_ii
   opal1d(i) = opal_ii
   sline1d(i) = sline_ii
   int1d(i) = int1d(i-1)*expdtau_ii + a_ii*sline_iim2 + b_ii*sline_iim1 + c_ii*sline_ii
!
   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
   int_sc = int_u*abs_sc + contr_sc
!
!
else
!
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = opalbar_u * phinorm_iim1
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii = opalbar_p*phinorm_ii
   call calc_phinorm(vel_d, vth_mean, vth_fiducial, xobs, phinorm_iip1)
   opal_iip1 = opalbar_d*phinorm_iip1
!
!calculate delta-tau steps (assuming that opal is 2nd order bezier function in s)
!   opalprime_ii = (opal_ii-opal_iim1)*dels_d/dels_u/(dels_d+dels_u) + (opal_iip1-opal_ii)*dels_u/dels_d/(dels_u+dels_d)
!   opal_cu = opal_ii - dels_u*opalprime_ii/2.d0
!   opal_cd = opal_ii + dels_d*opalprime_ii/2.d0
!   call pointcl1d_mbez (opal_cu, opal_iim1, opal_ii)
!   call pointcl1d_mbez (opal_cd, opal_ii, opal_iip1)
!   delt_u = (opal_iim1+opal_cu+opal_ii)*dels_u/3.d0
!   delt_d = (opal_ii+opal_cd+opal_iip1)*dels_d/3.d0
!calculate delta-tau steps (assuming that opal is piecewise linear function in s)
   delt_u = (opal_iim1+opal_ii)*dels_u/2.d0
   delt_d = (opal_ii+opal_iip1)*dels_d/2.d0
!
!perform radiative transfer
   abs_sc = exp(-delt_u)
!
!   call coeff_source1(delt_u, alo_u, alo_p)
!   alo_d=0.d0
   call coeff_source2(delt_u, delt_d, sline_u, sline_p, sline_d, alo_u, alo_p, alo_d)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
   int_sc = int_u*abs_sc + contr_sc
!
   nd=2
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(sline1d(nd))
   allocate(int1d(nd))
   s1d = (/ s_u, s_p /)
   tau1d = (/ 0.d0, delt_u /)
   vel1d = (/ vel_u, vel_p /)
   vth1d = (/ vth_mean, vth_mean /)
   xcmf1d = (/ xcmf_u, xcmf_p /)
   profile1d = (/ phinorm_iim1, phinorm_ii /)
   opalbar1d = (/ opalbar_u, opalbar_p /)
   opal1d = (/ opal_iim1, opal_ii /)
   sline1d = (/ sline_u, sline_p /)
   int1d = (/ int_u, int_sc /)
!
!
endif
!
!
!
end subroutine fsc_line_debug
!
!***********************************************************************
!***********************************************************************
!
!                   LINE + CONTINUUM ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_linec_lin(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, &
                         opac_u, opac_p, &
                         opalbar_u, opalbar_p, &
                         scont_u, scont_p, &
                         sline_u, sline_p, &
                         vel_u, vel_p, &
                         vth_u, vth_p, &
                         dels_u, &
                         abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line+continuum transport
!
!             u------------p
!                 dels_u
!
!      linear interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          deltax:                              required xobs-spacing
!          xcmf_min, xcmf_max:                  range of xcmf  
!          upwind intensity at u:               int_u
!          continnum opacity at u:              opac_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           scont_u
!          line source-fct at u:                sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          continnum opacity at p:              opac_p
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at u:           scont_p
!          line source-fct at p:                sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!
!   output: absorption part from u to p:       abs_sc
!           (total) source contribution:       contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_integ1d, only: coeff_source1
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, dels_u, &
                        opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                        opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
!
! ... local scalars
integer(i4b) :: i, nd
real(dp) :: m, q, s_p, s_u
real(dp) :: del_xcmf, xcmf_u, xcmf_p
real(dp) :: s_iim1, s_ii, sline_iim1, sline_ii, xcmf_iim1, xcmf_ii, opac_iim1, opac_ii, &
            opal_iim1, opal_ii, opatot_iim1, opatot_ii, vel_iim1, vel_ii
real(dp) :: delta_mean, vth_mean, phinorm
real(dp) :: delt_u, delt_ii, expdelt_ii
real(dp) :: a_ii, b_ii, a_iim1, b_iim1, a, b, c, alo_ucont, alo_pcont
!
! ... for debug: analytic structure
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
!
!use a mean delta
vth_mean=(vth_u+vth_p)/2.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
!----------------------------start point--------------------------------
!
   s_iim1 = s_u
!
!linear interpolation coefficients
   b_iim1=s_iim1/s_p
   a_iim1=1.d0-b_iim1
!
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iim1 = (a_iim1*opalbar_u + b_iim1*opalbar_p)*phinorm
   opac_iim1 = a_iim1*opac_u + b_iim1*opac_p
   opatot_iim1 = opal_iim1+opac_iim1
!
   delt_u=0.d0
   alo_u=0.d0
   alo_p=0.d0
   alo_ucont=0.d0
   alo_pcont=0.d0
!
!-----------------------------------------------------------------------
!
   do i=2, nd
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_ii = xcmf_u + (i-1)*del_xcmf
      s_ii = s_u + (xcmf_ii-xcmf_u)/m
      b_ii=s_ii/s_p
      a_ii=1.d0-b_ii
      vel_ii = a_ii*vel_u + b_ii*vel_p
      call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm)
      opal_ii = (a_ii*opalbar_u+b_ii*opalbar_p)*phinorm
      opac_ii = a_ii*opac_u+b_ii*opac_p
      opatot_ii = opal_ii + opac_ii
!
!perform radiative transfer
      delt_ii = (opatot_ii+opatot_iim1)/2.d0*(s_ii-s_iim1)
      delt_u = delt_u + delt_ii
!
      call coeff_source1(delt_ii, a, b)
      expdelt_ii = exp(-delt_ii)
!
      alo_u = alo_u*expdelt_ii + (a*a_iim1*opal_iim1/opatot_iim1 + b*a_ii*opal_ii/opatot_ii)
      alo_p = alo_p*expdelt_ii + (a*b_iim1*opal_iim1/opatot_iim1 + b*b_ii*opal_ii/opatot_ii)
!
      alo_ucont = alo_ucont*expdelt_ii + (a*a_iim1*opac_iim1/opatot_iim1 + b*a_ii*opac_ii/opatot_ii)
      alo_pcont = alo_pcont*expdelt_ii + (a*b_iim1*opac_iim1/opatot_iim1 + b*b_ii*opac_ii/opatot_ii)
!
!referesh variables at iim1
      s_iim1=s_ii
      opatot_iim1 = opatot_ii
      opal_iim1 = opal_ii
      opac_iim1 = opac_ii
      a_iim1 = a_ii
      b_iim1 = b_ii
   enddo

   abs_sc = exp(-delt_u)                                 !absorption part
   contr_sc = alo_u*sline_u + alo_p*sline_p + &          !line contribution
              alo_ucont*scont_u + alo_pcont*scont_p      !continuum contribution
   int_sc = abs_sc*int_u + contr_sc
!
else
!
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm)
   opal_iim1 = opalbar_u * phinorm
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm)
   opal_ii = opalbar_p * phinorm
!
!calculate total delta-tau step
   opatot_iim1 = opal_iim1+opac_u
   opatot_ii   = opal_ii  +opac_p
   delt_u = (opatot_iim1+opatot_ii)*(s_p-s_u)/2.d0
!perform radiative transfer
   abs_sc = exp(-delt_u)
!source contribution coefficients (for both, continuum and line)
   call coeff_source1(delt_u, a, b)
   alo_u = a*opal_iim1/opatot_iim1
   alo_p = b*opal_ii/opatot_ii
!
   contr_sc = alo_u*sline_u + alo_p*sline_p + &                           !line contribution
              a*opac_u/opatot_iim1*scont_u + b*opac_p/opatot_ii*scont_p   !continuum contribution
!
   int_sc = int_u*abs_sc + contr_sc

!
endif
!
!
end subroutine fsc_linec_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_linec_lin_debug(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, &
                               opac_u, opac_p, &
                               opalbar_u, opalbar_p, &
                               scont_u, scont_p, &
                               sline_u, sline_p, &
                               vel_u, vel_p, &
                               vth_u, vth_p, &
                               dels_u, &
                               abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                               xcmf1d, vel1d, opac1d, opalbar1d, opal1d, opatot1d, &
                               vth1d, scont1d, sline1d, stot1d, profile1d, int1d, s1d, tau1d, nd)
!
!-----------------------------------------------------------------------
!
!   debugging routine for fsc_line1d: output all physical quantities
!              on a refined grid to test the routine
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line+continuum transport
!
!             u------------p
!                 dels_u
!
!      linear interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          continnum opacity at u:              opac_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           scont_u
!          line source-fct at u:                sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          continnum opacity at p:              opac_p
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at u:           scont_p
!          line source-fct at p:                sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!
!   output: absorption part from u to p:       abs_sc
!           (total) source contribution:       contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_integ1d, only: coeff_source1
!
implicit none
!
! ... arguments
logical, intent(in) :: luse_refine
integer, intent(out) :: nd
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, dels_u, &
                        opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                        opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, opac1d, opatot1d, opal1d, &
                                                      sline1d, scont1d, stot1d, vth1d, int1d, profile1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: m, q, s_p, s_u
real(dp) :: del_xcmf, xcmf_u, xcmf_p
real(dp) :: int_iim1, int_ii, opalbar_iim1, opalbar_ii, opac_iim1, opac_ii, &
            opal_iim1, opal_ii, opatot_iim1, opatot_ii, sline_iim1, sline_ii, &
            scont_iim1, scont_ii, stot_iim1, stot_ii, s_iim1, s_ii, xcmf_iim1, xcmf_ii, &
            vel_iim1, vel_ii, phinorm_iim1, phinorm_ii
real(dp) :: delta_mean, vth_mean
real(dp) :: delt_u, delt_ii, expdelt_ii
real(dp) :: a_ii, b_ii, a_iim1, b_iim1, a, b, c, alo_ucont, alo_pcont
!
! ... for debug: analytic structure
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(allocated(s1d)) deallocate(s1d)
if(allocated(tau1d)) deallocate(tau1d)
if(allocated(vel1d)) deallocate(vel1d)
if(allocated(vth1d)) deallocate(vth1d)
if(allocated(xcmf1d)) deallocate(xcmf1d)
if(allocated(profile1d)) deallocate(profile1d)
if(allocated(opalbar1d)) deallocate(opalbar1d)
if(allocated(opal1d)) deallocate(opal1d)
if(allocated(opac1d)) deallocate(opac1d)
if(allocated(opatot1d)) deallocate(opatot1d)
if(allocated(sline1d)) deallocate(sline1d)
if(allocated(scont1d)) deallocate(scont1d)
if(allocated(stot1d)) deallocate(stot1d)
if(allocated(int1d)) deallocate(int1d)
!
!-----------------------------------------------------------------------
!
!define s_p, s_u (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
!
!use a mean delta
vth_mean=(vth_u+vth_p)/2.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(luse_refine.and.lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(opac1d(nd))
   allocate(opatot1d(nd))
   allocate(sline1d(nd))
   allocate(scont1d(nd))
   allocate(stot1d(nd))
   allocate(int1d(nd))
!
!----------------------------start point--------------------------------
!
   s_iim1 = s_u
!
!linear interpolation coefficients
   b_iim1=s_iim1/s_p
   a_iim1=1.d0-b_iim1
!
   int_iim1 = int_u
   xcmf_iim1 = a_iim1*xcmf_u + b_iim1*xcmf_p
   sline_iim1 = a_iim1*sline_u + b_iim1*sline_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   opalbar_iim1 = a_iim1*opalbar_u + b_iim1*opalbar_p
   opac_iim1    = a_iim1*opac_u + b_iim1*opac_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1    = opalbar_iim1*phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   sline_iim1   = a_iim1*sline_u + b_iim1*sline_p
   scont_iim1   = a_iim1*scont_u + b_iim1*scont_p
   stot_iim1    = sline_iim1*opal_iim1/opatot_iim1 + scont_iim1*opac_iim1/opatot_iim1
!
!
   s1d(1) = s_iim1
   tau1d(1) = 0.d0
   vel1d(1) = vel_u
   vth1d(1) = vth_mean
   xcmf1d(1) = xcmf_u
   profile1d(1) = phinorm_iim1
   opalbar1d(1) = opalbar_u
   opal1d(1) = opal_iim1
   opac1d(1) = opac_iim1
   opatot1d(1) = opatot_iim1
   sline1d(1) = sline_iim1
   scont1d(1) = scont_iim1
   stot1d(1) = stot_iim1
   int1d(1) = int_u
!
   delt_u=0.d0
   alo_u=0.d0
   alo_p=0.d0
   alo_ucont=0.d0
   alo_pcont=0.d0
!
!-----------------------------------------------------------------------
!
   do i=2, nd
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_ii = xcmf_u + (i-1)*del_xcmf
      s_ii = s_u + (xcmf_ii-xcmf_u)/m
      b_ii=s_ii/s_p
      a_ii=1.d0-b_ii
      vel_ii = a_ii*vel_u + b_ii*vel_p
      opalbar_ii = a_ii*opalbar_u + b_ii*opalbar_p
      opac_ii    = a_ii*opac_u + b_ii*opac_p
      call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
      opal_ii    = opalbar_ii*phinorm_ii
      opatot_ii  = opal_ii + opac_ii
      sline_ii = a_ii*sline_u + b_ii*sline_p
      scont_ii   = a_ii*scont_u + b_ii*scont_p
      stot_ii    = sline_ii*opal_ii/opatot_ii + scont_ii*opac_ii/opatot_ii
!
!perform radiative transfer
      delt_ii = (opatot_ii+opatot_iim1)*(s_ii-s_iim1)/2.d0
      delt_u = delt_u + delt_ii
!
      call coeff_source1(delt_ii, a, b)
      expdelt_ii = exp(-delt_ii)
!
      int_ii = int_iim1*expdelt_ii + a*stot_iim1 + b*stot_ii
!
      alo_u = alo_u*expdelt_ii + (a*a_iim1*opal_iim1/opatot_iim1 + b*a_ii*opal_ii/opatot_ii)
      alo_p = alo_p*expdelt_ii + (a*b_iim1*opal_iim1/opatot_iim1 + b*b_ii*opal_ii/opatot_ii)
!
      alo_ucont = alo_ucont*expdelt_ii + (a*a_iim1*opac_iim1/opatot_iim1 + b*a_ii*opac_ii/opatot_ii)
      alo_pcont = alo_pcont*expdelt_ii + (a*b_iim1*opac_iim1/opatot_iim1 + b*b_ii*opac_ii/opatot_ii)
!
!referesh variables at iim1
      xcmf_iim1=xcmf_ii
      s_iim1=s_ii
      phinorm_iim1=phinorm_ii
      opalbar_iim1=opalbar_ii
      opatot_iim1 = opatot_ii
      opal_iim1 = opal_ii
      opac_iim1 = opac_ii
      sline_iim1 = sline_ii
      scont_iim1 = scont_ii
      stot_iim1 = stot_ii
      vel_iim1 = vel_ii
      a_iim1 = a_ii
      b_iim1 = b_ii
      int_iim1 = int_ii
!
      int1d(i) = int_ii
      s1d(i) = s_ii
      tau1d(i) = tau1d(i-1)+delt_ii
      vel1d(i) = vel_ii
      vth1d(i) = vth_mean
      xcmf1d(i) = xcmf_ii
      profile1d(i) = phinorm_ii
      opalbar1d(i) = opalbar_ii
      opac1d(i) = opac_ii
      opal1d(i) = opal_ii
      opatot1d(i) = opatot_ii
      sline1d(i) = sline_ii
      scont1d(i) = scont_ii
      stot1d(i) = stot_ii

   enddo

   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_ucont*scont_u + alo_pcont*scont_p
   int_sc = abs_sc*int_u+contr_sc
!
   if(abs(int_ii-int_sc).gt.1.d-14) then
      write(*,*) 'error in fsc_linec_lin_debug: alo-contribution and standard radiative transfer not consistent'
      stop
   endif
!
else
!
!for point iim1
   s_iim1       = s_u
   b_iim1       = (s_iim1-s_u)/(s_p-s_u)
   a_iim1       = 1.d0-b_iim1
   vel_iim1     = a_iim1*vel_u + b_iim1*vel_p
   opalbar_iim1 = a_iim1*opalbar_u + b_iim1*opalbar_p
   opac_iim1    = a_iim1*opac_u + b_iim1*opac_p
   call calc_phinorm(vel_u, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1 = opalbar_u * phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   sline_iim1   = a_iim1*sline_u + b_iim1*sline_p
   scont_iim1   = a_iim1*scont_u + b_iim1*scont_p
   stot_iim1    = sline_iim1*opal_iim1/opatot_iim1 + scont_iim1*opac_iim1/opatot_iim1

!for point ii
   s_ii       = s_p
   b_ii       = (s_ii-s_u)/(s_p-s_u)
   a_ii       = 1.d0-b_ii
   vel_ii     = a_ii*vel_u + b_ii*vel_p
   opalbar_ii = a_ii*opalbar_u + b_ii*opalbar_p
   opac_ii    = a_ii*opac_u + b_ii*opac_p
   call calc_phinorm(vel_p, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii = opalbar_p*phinorm_ii
   opatot_ii  = opal_ii + opac_ii
   sline_ii   = a_ii*sline_u + b_ii*sline_p
   scont_ii   = a_ii*scont_u + b_ii*scont_p
   stot_ii    = sline_ii*opal_ii/opatot_ii + scont_ii*opac_ii/opatot_ii
!
!calculate delta-tau step
   delt_u = (opatot_iim1+opatot_ii)*(s_p-s_u)/2.d0
   abs_sc = exp(-delt_u)
!
   call coeff_source1(delt_u, a, b)
   alo_u = a*opal_iim1/opatot_iim1
   alo_p = b*opal_ii/opatot_ii
!
   contr_sc = alo_u*sline_u + alo_p*sline_p + &                                   !line contribution
              a*opac_iim1/opatot_iim1*scont_iim1 + b*opac_ii/opatot_ii*scont_ii   !continuum contribution

   int_sc = int_u*abs_sc + contr_sc

!
   nd=2
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(opac1d(nd))
   allocate(opatot1d(nd))
   allocate(sline1d(nd))
   allocate(scont1d(nd))
   allocate(stot1d(nd))
   allocate(int1d(nd))
   s1d = (/ s_u, s_p /)
   tau1d = (/ 0.d0, delt_u /)
   vel1d = (/ vel_u, vel_p /)
   vth1d = (/ vth_mean, vth_mean /)
   xcmf1d = (/ xcmf_u, xcmf_p /)
   profile1d = (/ phinorm_iim1, phinorm_ii /)
   opalbar1d = (/ opalbar_iim1, opalbar_ii /)
   opal1d = (/ opal_iim1, opal_ii /)
   opac1d = (/ opac_iim1, opac_ii /)
   opatot1d = (/ opatot_iim1, opatot_ii /)
   sline1d = (/ sline_iim1, sline_ii /)
   scont1d = (/ scont_iim1, scont_ii /)
   stot1d = (/ stot_iim1, stot_ii /)
   int1d = (/ int_u, int_sc /)
!
endif
!
!
end subroutine fsc_linec_lin_debug
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_linec(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, &
                       opac_u, opac_p, opac_d, &
                       opalbar_u, opalbar_p, opalbar_d, &
                       scont_u, scont_p, scont_d, &
                       sline_u, sline_p, sline_d, &
                       vel_u, vel_p, vel_d, &
                       vth_u, vth_p, vth_d, &
                       dels_u, dels_d, &
                       abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line+continuum transport
!
!             u------------p----------d
!                 dels_u     dels_d
!
!      bezier interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          continnum opacity at u:              opac_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           scont_u
!          line source-fct at u:                sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          continnum opacity at p:              opac_p
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           scont_p
!          line source-fct at p:                sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!          continnum opacity at d:              opac_d
!          frequency integrated opacity at d:   opalbar_d
!          continuum source-fct at d:           scont_d
!          line source-fct at d:                sline_d
!          projected velocity at d:             vel_d   (in vth_fiducial)
!          thermal velocity at d:               vth_d 
!
!   output: absorption part from u to p:       abs_sc
!           (total) source contribution:       contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p, alo_d
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_interp2d, only: wp_interp1d
use mod_integ1d, only: coeff_source2, coeff_source2b
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, dels_u, dels_d, &
                        opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                        opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p, &
                        opac_d, opalbar_d, scont_d, sline_d, vel_d, vth_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
!
! ... local scalars
integer(i4b) :: i, nd
real(dp) :: m, q, s_u, s_p, s_d
real(dp) :: del_xcmf, xcmf_u, xcmf_p, xcmf_d
real(dp) :: alpha_opac, wu_opac, wp_opac, wd_opac, opac_c, &
            alpha_opalbar, wu_opalbar, wp_opalbar, wd_opalbar, opalbar_c, &
            alpha_scont, wu_scont, wp_scont, wd_scont, scont_c, &
            alpha_sline, wu_sline, wp_sline, wd_sline, sline_c
real(dp) :: a_iim2, a_iim1, a_ii, a_iip1, b_iim2, b_iim1, b_ii, b_iip1, &
            ts_iim2, ts_iim1, ts_ii, ts_iip1, &
            at_iim2, at_iim1, at_ii, at_iip1, &
            bt_iim2, bt_iim1, bt_ii, bt_iip1, &
            ct_iim2, ct_iim1, ct_ii, ct_iip1, &
            atline_iim2, atline_iim1, atline_ii, atline_iip1, &
            btline_iim2, btline_iim1, btline_ii, btline_iip1, &
            ctline_iim2, ctline_iim1, ctline_ii, ctline_iip1, &
            atcont_iim2, atcont_iim1, atcont_ii, atcont_iip1, &
            btcont_iim2, btcont_iim1, btcont_ii, btcont_iip1, &
            ctcont_iim2, ctcont_iim1, ctcont_ii, ctcont_iip1
real(dp) :: s_iim2, s_iim1, s_ii, s_iip1, &
            sline_iim2, sline_iim1, sline_ii, sline_iip1, &
            scont_iim2, scont_iim1, scont_ii, scont_iip1, &
            stot_iim2, stot_iim1, stot_ii, stot_iip1, &
            opac_iim2, opac_iim1, opac_ii, opac_iip1, &
            opalbar_iim2, opalbar_iim1, opalbar_ii, opalbar_iip1, &
            opal_iim2, opal_iim1, opal_ii, opal_iip1, &
            opatot_iim2, opatot_iim1, opatot_ii, opatot_iip1, &
            phinorm_iim2, phinorm_iim1, phinorm_ii, phinorm_iip1, &
            xcmf_iim2, xcmf_iim1, xcmf_ii, xcmf_iip1, &
            vel_iim2, vel_iim1, vel_ii, vel_iip1
real(dp) :: delta_mean, vth_mean, phinorm
real(dp) :: delt_u, delt_d, delt_iim1, delt_ii, delt_iip1, expdtau_ii, int_iim1, int_ii
real(dp) :: a, b, c, alo_ucont, alo_pcont, alo_dcont
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!define s_u, s_p, s_d (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
s_d=s_p+abs(dels_d)
!
!use a mean delta
vth_mean=(vth_u+vth_p+vth_d)/3.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
xcmf_d=(xobs-vel_d)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
!---------------calculate control points and coefficients---------------
!-----------------------for interpolation along s-----------------------
!
!derivative weights for opacity
   alpha_opalbar = dels_d/s_d
   alpha_opac = alpha_opalbar
!derivative weights for source function
   alpha_sline = max(wp_interp1d,alpha_opalbar)
   alpha_scont = alpha_opalbar
!control-point weights from weighted derivatives
   wu_sline = alpha_sline/2.d0
   wp_sline = ((2.d0-alpha_sline)*dels_d + (1.d0-alpha_sline)*dels_u)/2.d0/dels_d
   wd_sline = (alpha_sline-1.d0)*dels_u/2.d0/dels_d
   wu_scont = alpha_scont/2.d0
   wp_scont = ((2.d0-alpha_scont)*dels_d + (1.d0-alpha_scont)*dels_u)/2.d0/dels_d
   wd_scont = (alpha_scont-1.d0)*dels_u/2.d0/dels_d

   wu_opalbar = alpha_opalbar/2.d0
   wp_opalbar = ((2.d0-alpha_opalbar)*dels_d + (1.d0-alpha_opalbar)*dels_u)/2.d0/dels_d
   wd_opalbar = (alpha_opalbar-1.d0)*dels_u/2.d0/dels_d
   wu_opac = alpha_opac/2.d0
   wp_opac = ((2.d0-alpha_opac)*dels_d + (1.d0-alpha_opac)*dels_u)/2.d0/dels_d
   wd_opac = (alpha_opac-1.d0)*dels_u/2.d0/dels_d

!(non-monotonic) control points
   opalbar_c = wu_opalbar*opalbar_u + wp_opalbar*opalbar_p + wd_opalbar*opalbar_d
   opac_c = wu_opac*opac_u + wp_opac*opac_p + wd_opac*opac_d
   sline_c = wu_sline*sline_u + wp_sline*sline_p + wd_sline*sline_d
   scont_c = wu_scont*scont_u + wp_scont*scont_p + wd_scont*scont_d
!ensure monotonicity
   call pointcl1d_mbez(opalbar_c, opalbar_u, opalbar_p)
   call pointcl1d_mbez(opac_c, opac_u, opac_p)
   call coeffcl1d_mbez(scont_c, scont_u, scont_p, wu_scont, wp_scont, wd_scont)
   call coeffcl1d_mbez(sline_c, sline_u, sline_p, wu_sline, wp_sline, wd_sline)
!
!----------------------------start point--------------------------------
!
   alo_u = 0.d0
   alo_p = 0.d0
   alo_d = 0.d0
   alo_ucont = 0.d0
   alo_pcont = 0.d0
   alo_dcont = 0.d0
   delt_u = 0.d0
!
!at point u (<=> iim1)
   xcmf_iim1 = xcmf_u
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim1 = s_iim1/s_p
   a_iim1 = 1.d0-b_iim1
!(bezier) interpolation coefficients
   ts_iim1 = (s_iim1-s_u)/dels_u
   at_iim1 = (1.d0-ts_iim1)**2
   bt_iim1 = 2.d0*ts_iim1*(1.d0-ts_iim1)
   ct_iim1 = ts_iim1**2
   atline_iim1 = at_iim1 + bt_iim1*wu_sline
   btline_iim1 = ct_iim1 + bt_iim1*wp_sline
   ctline_iim1 = bt_iim1*wd_sline
   atcont_iim1 = at_iim1 + bt_iim1*wu_scont
   btcont_iim1 = ct_iim1 + bt_iim1*wp_scont
   ctcont_iim1 = bt_iim1*wd_scont
!actual interpolation
   sline_iim1 = atline_iim1*sline_u + btline_iim1*sline_p + ctline_iim1*sline_d
   scont_iim1 = atcont_iim1*scont_u + btcont_iim1*scont_p + ctcont_iim1*scont_p
   opalbar_iim1 = at_iim1*opalbar_u + bt_iim1*opalbar_c + ct_iim1*opalbar_p
   opac_iim1 = at_iim1*opac_u + bt_iim1*opac_c + ct_iim1*opac_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1    = opalbar_iim1*phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   stot_iim1 = opal_iim1*sline_iim1/opatot_iim1 + opac_iim1*scont_iim1/opatot_iim1

   int_iim1 = int_u
!
!at point ii
   xcmf_ii = xcmf_u+del_xcmf
   s_ii = s_u + (xcmf_ii-xcmf_u)/m
!(linear) interpolation coefficients
   b_ii = s_ii/s_p
   a_ii = 1.d0-b_ii
!(bezier) interpolation coefficients
   ts_ii = (s_ii-s_u)/dels_u
   at_ii = (1.d0-ts_ii)**2
   bt_ii = 2.d0*ts_ii*(1.d0-ts_ii)
   ct_ii = ts_ii**2
   atline_ii = at_ii + bt_ii*wu_sline
   btline_ii = ct_ii + bt_ii*wp_sline
   ctline_ii = bt_ii*wd_sline
   atcont_ii = at_ii + bt_ii*wu_scont
   btcont_ii = ct_ii + bt_ii*wp_scont
   ctcont_ii = bt_ii*wd_scont
!actual interpolation
   sline_ii = atline_ii*sline_u + btline_ii*sline_p + ctline_ii*sline_d
   scont_ii = atcont_ii*scont_u + btcont_ii*scont_p + ctcont_ii*scont_d
   opalbar_ii = at_ii*opalbar_u + bt_ii*opalbar_c + ct_ii*opalbar_p
   opac_ii = at_ii*opac_u + bt_ii*opac_c + ct_ii*opac_p
   vel_ii = a_ii*vel_u + b_ii*vel_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii    = opalbar_ii*phinorm_ii
   opatot_ii  = opal_ii + opac_ii
   stot_ii = opal_ii*sline_ii/opatot_ii + opac_ii*scont_ii/opatot_ii
!
    delt_ii = (opatot_iim1+opatot_ii)*(s_ii-s_iim1)/2.d0
!
!-----------------------------------------------------------------------
!
   do i=2, nd-1
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_iip1 = xcmf_u + i*del_xcmf
      s_iip1 = s_u + (xcmf_iip1-xcmf_u)/m

!(linear) interpolation coefficients
      b_iip1 = s_iip1/s_p
      a_iip1 = 1.d0-b_iip1
!(bezier) interpolation coefficients
      ts_iip1 = (s_iip1-s_u)/dels_u
      at_iip1 = (1.d0-ts_iip1)**2
      bt_iip1 = 2.d0*ts_iip1*(1.d0-ts_iip1)
      ct_iip1 = ts_iip1**2
      atline_iip1 = at_iip1 + bt_iip1*wu_sline
      btline_iip1 = ct_iip1 + bt_iip1*wp_sline
      ctline_iip1 = bt_iip1*wd_sline
      atcont_iip1 = at_iip1 + bt_iip1*wu_scont
      btcont_iip1 = ct_iip1 + bt_iip1*wp_scont
      ctcont_iip1 = bt_iip1*wd_scont
!actual interpolation
      sline_iip1 = atline_iip1*sline_u + btline_iip1*sline_p + ctline_iip1*sline_d
      scont_iip1 = atcont_iip1*scont_u + btcont_iip1*scont_p + ctcont_iip1*scont_d
      opalbar_iip1 = at_iip1*opalbar_u + bt_iip1*opalbar_c + ct_iip1*opalbar_p
      opac_iip1 = at_iip1*opac_u + bt_iip1*opac_c + ct_iip1*opac_p
      vel_iip1 = a_iip1*vel_u + b_iip1*vel_p
      call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
      opal_iip1    = opalbar_iip1*phinorm_iip1
      opatot_iip1  = opal_iip1 + opac_iip1
      stot_iip1 = opal_iip1*sline_iip1/opatot_iip1 + opac_iip1*scont_iip1/opatot_iip1
!
!calculate delta tau-step
      delt_iip1 = (opatot_ii+opatot_iip1)*(s_iip1-s_ii)/2.d0 
!
!calculate coefficients for source function integration
      call coeff_source2(delt_ii, delt_iip1, stot_iim1, stot_ii, stot_iip1, a, b, c)
!
      delt_u = delt_u + delt_ii
      expdtau_ii=exp(-delt_ii)
!      int_ii=int_iim1*expdtau_ii + a*stot_iim1 + b*stot_ii + c*stot_iip1
!
      alo_u = alo_u*expdtau_ii + a*atline_iim1*opal_iim1/opatot_iim1 + &
                                 b*atline_ii*opal_ii/opatot_ii + &
                                 c*atline_iip1*opal_iip1/opatot_iip1
      alo_p = alo_p*expdtau_ii + a*btline_iim1*opal_iim1/opatot_iim1 + &
                                 b*btline_ii*opal_ii/opatot_ii + &
                                 c*btline_iip1*opal_iip1/opatot_iip1
      alo_d = alo_d*expdtau_ii + a*ctline_iim1*opal_iim1/opatot_iim1 + &
                                 b*ctline_ii*opal_ii/opatot_ii + &
                                 c*ctline_iip1*opal_iip1/opatot_iip1
!
      alo_ucont = alo_ucont*expdtau_ii + a*atcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*atcont_ii*opac_ii/opatot_ii + &
                                         c*atcont_iip1*opac_iip1/opatot_iip1
      alo_pcont = alo_pcont*expdtau_ii + a*btcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*btcont_ii*opac_ii/opatot_ii + &
                                         c*btcont_iip1*opac_iip1/opatot_iip1
      alo_dcont = alo_dcont*expdtau_ii + a*ctcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*ctcont_ii*opac_ii/opatot_ii + &
                                         c*ctcont_iip1*opac_iip1/opatot_iip1

!referesh variables at iim1
!      xcmf_iim1 = xcmf_ii
!      s_iim1 = s_ii
!      phinorm_iim1 = phinorm_ii
!      opalbar_iim1 = opalbar_ii
!      sline_iim1 = sline_ii
!      scont_iim1 = scont_ii
!      vel_iim1 = vel_ii
      opac_iim1 = opac_ii
      opatot_iim1 = opatot_ii
      opal_iim1 = opal_ii
      stot_iim1 = stot_ii
!      int_iim1=int_ii
      atline_iim1 = atline_ii
      btline_iim1 = btline_ii
      ctline_iim1 = ctline_ii
      atcont_iim1 = atcont_ii
      btcont_iim1 = btcont_ii
      ctcont_iim1 = ctcont_ii
!
!referesh variables at ii
!      xcmf_ii = xcmf_iip1
!      phinorm_ii = phinorm_iip1
!      opalbar_ii = opalbar_iip1
!      sline_ii = sline_iip1
!      scont_ii = scont_iip1
!      vel_ii = vel_iip1
      s_ii = s_iip1
      opac_ii = opac_iip1
      opatot_ii = opatot_iip1
      opal_ii = opal_iip1
      stot_ii = stot_iip1
      delt_ii = delt_iip1
!
      atline_ii = atline_iip1
      btline_ii = btline_iip1
      ctline_ii = ctline_iip1
      atcont_ii = atcont_iip1
      btcont_ii = btcont_iip1
      ctcont_ii = ctcont_iip1


   enddo
!
!-------------------------interval [n-1,n]------------------------------
!
!recalculate everything at at n-2
   xcmf_iim2 = xcmf_p - 2*del_xcmf
   xcmf_iim1 = xcmf_p - del_xcmf
   s_iim2 = s_u + (xcmf_iim2-xcmf_u)/m
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim2 = (s_iim2-s_u)/dels_u
   a_iim2 = 1.d0-b_iim2
!(bezier) interpolation coefficients
   ts_iim2 = (s_iim2-s_u)/dels_u
   at_iim2 = (1.d0-ts_iim2)**2
   bt_iim2 = 2.d0*ts_iim2*(1.d0-ts_iim2)
   ct_iim2 = ts_iim2**2
   atline_iim2 = at_iim2 + bt_iim2*wu_sline
   btline_iim2 = ct_iim2 + bt_iim2*wp_sline
   ctline_iim2 = bt_iim2*wd_sline
   atcont_iim2 = at_iim2 + bt_iim2*wu_scont
   btcont_iim2 = ct_iim2 + bt_iim2*wp_scont
   ctcont_iim2 = bt_iim2*wd_scont
!actual interpolation
   sline_iim2 = atline_iim2*sline_u + btline_iim2*sline_p + ctline_iim2*sline_d
   scont_iim2 = atcont_iim2*scont_u + btcont_iim2*scont_p + ctcont_iim2*scont_d
   opalbar_iim2 = at_iim2*opalbar_u + bt_iim2*opalbar_c + ct_iim2*opalbar_p
   opac_iim2 = at_iim2*opac_u + bt_iim2*opac_c + ct_iim2*opac_p
   vel_iim2 = a_iim2*vel_u + b_iim2*vel_p
   call calc_phinorm(vel_iim2, vth_mean, vth_fiducial, xobs, phinorm_iim2)
   opal_iim2 = opalbar_iim2*phinorm_iim2
   opatot_iim2  = opal_iim2 + opac_iim2
   stot_iim2 = opal_iim2*sline_iim2/opatot_iim2 + opac_iim2*scont_iim2/opatot_iim2
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   delt_iim1 = (opatot_iim2+opatot_iim1)*(s_iim1-s_iim2)/2.d0
!
!calculate coefficients for source function integration
   call coeff_source2b(delt_iim1, delt_ii, stot_iim2, stot_iim1, stot_ii, a, b, c)

   expdtau_ii=exp(-delt_ii)
   delt_u = delt_u + delt_ii
!
!   int_ii=int_iim1*expdtau_ii + a*stot_iim2 + b*stot_iim1 + c*stot_ii
!
   alo_u = alo_u*expdtau_ii + a*atline_iim2*opal_iim2/opatot_iim2 + &
                              b*atline_iim1*opal_iim1/opatot_iim1 + &
                              c*atline_ii*opal_ii/opatot_ii
   alo_p = alo_p*expdtau_ii + a*btline_iim2*opal_iim2/opatot_iim2 + &
                              b*btline_iim1*opal_iim1/opatot_iim1 + &
                              c*btline_ii*opal_ii/opatot_ii
   alo_d = alo_d*expdtau_ii + a*ctline_iim2*opal_iim2/opatot_iim2 + &
                              b*ctline_iim1*opal_iim1/opatot_iim1 + &
                              c*ctline_ii*opal_ii/opatot_ii

   alo_ucont = alo_ucont*expdtau_ii + a*atcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*atcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*atcont_ii*opac_ii/opatot_ii
   alo_pcont = alo_pcont*expdtau_ii + a*btcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*btcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*btcont_ii*opac_ii/opatot_ii
   alo_dcont = alo_dcont*expdtau_ii + a*ctcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*ctcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*ctcont_ii*opac_ii/opatot_ii
!
   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d + &
              alo_ucont*scont_u + alo_pcont*scont_p + alo_dcont*scont_d
   int_sc = int_u*abs_sc + contr_sc
!
else
!
!for point iim1
   s_iim1       = s_u
   vel_iim1     = vel_u
   opalbar_iim1 = opalbar_u
   opac_iim1    = opac_u
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1    = opalbar_iim1*phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   sline_iim1   = sline_u
   scont_iim1    = scont_u
   stot_iim1    = sline_iim1*opal_iim1/opatot_iim1 + scont_iim1*opac_iim1/opatot_iim1
!
!for point ii
   s_ii       = s_p
   vel_ii     = vel_p
   opalbar_ii = opalbar_p
   opac_ii    = opac_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii    = opalbar_ii*phinorm_ii
   opatot_ii  = opal_ii + opac_ii
   sline_ii   = sline_p
   scont_ii   = scont_p
   stot_ii    = sline_ii*opal_ii/opatot_ii + scont_ii*opac_ii/opatot_ii
!for point ii
   s_iip1       = s_d
   vel_iip1     = vel_d
   opalbar_iip1 = opalbar_d
   opac_iip1    = opac_d
   call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
   opal_iip1    = opalbar_iip1*phinorm_iip1
   opatot_iip1  = opal_iip1 + opac_iip1
   sline_iip1   = sline_d
   scont_iip1   = scont_d
   stot_iip1    = sline_iip1*opal_iip1/opatot_iip1 + scont_iip1*opac_iip1/opatot_iip1
!
!calculate delta-tau
   delt_u = (opatot_iim1+opatot_ii)*(s_p-s_u)/2.d0
   delt_d = (opatot_iip1+opatot_ii)*(s_d-s_p)/2.d0
!
   abs_sc = exp(-delt_u)
   call coeff_source2(delt_u, delt_d, stot_iim1, stot_ii, stot_iip1, a, b, c)
!WARNING: what happens if delt_u->, delt_d large, i.e. when resonance zone in downwind interval???

   alo_u = a*opal_iim1/opatot_iim1
   alo_p = b*opal_ii/opatot_ii
   alo_d = c*opal_iip1/opatot_iip1

!   contr_sc = a *stot_iim1 + b*stot_ii + c*stot_iip1
    contr_sc = alo_u*sline_iim1 + alo_p*sline_ii + alo_d*sline_iip1 + &                                                  !line contribution
               a*opac_iim1/opatot_iim1*scont_iim1 + b*opac_ii/opatot_ii*scont_ii + c*opac_iip1/opatot_iip1*scont_iip1    !continuum contribution
!
   int_sc = int_u*abs_sc + contr_sc
!
!
!
endif
!
!
end subroutine fsc_linec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_linec_debug(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, &
                           opac_u, opac_p, opac_d, &
                           opalbar_u, opalbar_p, opalbar_d, &
                           scont_u, scont_p, scont_d, &
                           sline_u, sline_p, sline_d, &
                           vel_u, vel_p, vel_d, &
                           vth_u, vth_p, vth_d, &
                           dels_u, dels_d, &
                           abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                           xcmf1d, vel1d, opac1d, opalbar1d, opal1d, opatot1d, &
                           vth1d, scont1d, sline1d, stot1d, profile1d, int1d, s1d, tau1d, nd)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for line+continuum transport
!
!             u------------p----------d
!                 dels_u     dels_d
!
!      bezier interpolations/integrations
!
!   input: considered frequency shift:          xobs
!          fiducial thermal velocity:           vth_fiducial
!          upwind intensity at u:               int_u
!          continnum opacity at u:              opac_u
!          frequency integrated opacity at u:   opalbar_u
!          continuum source-fct at u:           scont_u
!          line source-fct at u:                sline_u
!          projected velocity at u:             vel_u   (in vth_fiducial)
!          thermal velocity at u:               vth_u 
!          continnum opacity at p:              opac_p
!          frequency integrated opacity at p:   opalbar_p
!          continuum source-fct at p:           scont_p
!          line source-fct at p:                sline_p
!          projected velocity at p:             vel_p   (in vth_fiducial)
!          thermal velocity at p:               vth_p 
!          continnum opacity at d:              opac_d
!          frequency integrated opacity at d:   opalbar_d
!          continuum source-fct at d:           scont_d
!          line source-fct at d:                sline_d
!          projected velocity at d:             vel_d   (in vth_fiducial)
!          thermal velocity at d:               vth_d 
!
!   output: absorption part from u to p:       abs_sc
!           (total) source contribution:       contr_sc
!           intensity at point p:              int_sc
!           alo-coefficients:                  alo_u, alo_p, alo_d
!
!-----------------------------------------------------------------------
!
use prog_type
use mod_interp2d, only: wp_interp1d
use mod_integ1d, only: coeff_source2, coeff_source2b
!
implicit none
!
! ... arguments
logical, intent(in) :: luse_refine
integer(i4b), intent(out) :: nd
real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, dels_u, dels_d, &
                        opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                        opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p, &
                        opac_d, opalbar_d, scont_d, sline_d, vel_d, vth_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, opac1d, opatot1d, opal1d, &
                                                      sline1d, scont1d, stot1d, vth1d, int1d, profile1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: m, q, s_u, s_p, s_d
real(dp) :: del_xcmf, xcmf_u, xcmf_p, xcmf_d
real(dp) :: alpha_opac, wu_opac, wp_opac, wd_opac, opac_c, &
            alpha_opalbar, wu_opalbar, wp_opalbar, wd_opalbar, opalbar_c, &
            alpha_scont, wu_scont, wp_scont, wd_scont, scont_c, &
            alpha_sline, wu_sline, wp_sline, wd_sline, sline_c
real(dp) :: a_iim2, a_iim1, a_ii, a_iip1, b_iim2, b_iim1, b_ii, b_iip1, &
            ts_iim2, ts_iim1, ts_ii, ts_iip1, &
            at_iim2, at_iim1, at_ii, at_iip1, &
            bt_iim2, bt_iim1, bt_ii, bt_iip1, &
            ct_iim2, ct_iim1, ct_ii, ct_iip1, &
            atline_iim2, atline_iim1, atline_ii, atline_iip1, &
            btline_iim2, btline_iim1, btline_ii, btline_iip1, &
            ctline_iim2, ctline_iim1, ctline_ii, ctline_iip1, &
            atcont_iim2, atcont_iim1, atcont_ii, atcont_iip1, &
            btcont_iim2, btcont_iim1, btcont_ii, btcont_iip1, &
            ctcont_iim2, ctcont_iim1, ctcont_ii, ctcont_iip1
real(dp) :: s_iim2, s_iim1, s_ii, s_iip1, &
            sline_iim2, sline_iim1, sline_ii, sline_iip1, &
            scont_iim2, scont_iim1, scont_ii, scont_iip1, &
            stot_iim2, stot_iim1, stot_ii, stot_iip1, &
            opac_iim2, opac_iim1, opac_ii, opac_iip1, &
            opalbar_iim2, opalbar_iim1, opalbar_ii, opalbar_iip1, &
            opal_iim2, opal_iim1, opal_ii, opal_iip1, &
            opatot_iim2, opatot_iim1, opatot_ii, opatot_iip1, &
            phinorm_iim2, phinorm_iim1, phinorm_ii, phinorm_iip1, &
            xcmf_iim2, xcmf_iim1, xcmf_ii, xcmf_iip1, &
            vel_iim2, vel_iim1, vel_ii, vel_iip1
real(dp) :: delta_mean, vth_mean, phinorm
real(dp) :: delt_u, delt_d, delt_iim1, delt_ii, delt_iip1, expdtau_ii, int_iim1, int_ii
real(dp) :: a, b, c, alo_ucont, alo_pcont, alo_dcont
!
! ... local arrays
!
! ... local logicals
logical :: lrefine
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(allocated(s1d)) deallocate(s1d)
if(allocated(tau1d)) deallocate(tau1d)
if(allocated(vel1d)) deallocate(vel1d)
if(allocated(vth1d)) deallocate(vth1d)
if(allocated(xcmf1d)) deallocate(xcmf1d)
if(allocated(profile1d)) deallocate(profile1d)
if(allocated(opalbar1d)) deallocate(opalbar1d)
if(allocated(opal1d)) deallocate(opal1d)
if(allocated(opac1d)) deallocate(opac1d)
if(allocated(opatot1d)) deallocate(opatot1d)
if(allocated(sline1d)) deallocate(sline1d)
if(allocated(scont1d)) deallocate(scont1d)
if(allocated(stot1d)) deallocate(stot1d)
if(allocated(int1d)) deallocate(int1d)
!
!-----------------------------------------------------------------------
!
!define s_u, s_p, s_d (not necessarily needed, but for the moment for debugging reasons)
s_u=0.d0
s_p=abs(dels_u)
s_d=s_p+abs(dels_d)
!
!use a mean delta
vth_mean=(vth_u+vth_p+vth_d)/3.d0
delta_mean=vth_mean/vth_fiducial
!
xcmf_u=(xobs-vel_u)/delta_mean
xcmf_p=(xobs-vel_p)/delta_mean
xcmf_d=(xobs-vel_d)/delta_mean
!
!check if refinement is required
lrefine=.false.
if(abs(xcmf_u-xcmf_p).gt.deltax) then
!delx-steps are too large
   if((xcmf_u*xcmf_p).le.0.d0) then
!resonance zone between point u and p
      lrefine=.true.
   elseif(abs(xcmf_u).le.xcmf_max) then
!resonance zone affects upwind point
      lrefine=.true.
   elseif(abs(xcmf_p).le.xcmf_max) then
!resonance zone affects point
      lrefine=.true.
   endif
endif
!
!-----------------------------------------------------------------------
!
!now calculate refinement
if(luse_refine.and.lrefine) then
!gradient and sign of xcmf
   m=(xcmf_p-xcmf_u)/(s_p-s_u)

   nd=ceiling(abs(xcmf_p-xcmf_u)/deltax)+1
   del_xcmf = (xcmf_p-xcmf_u)/(nd-1)
!
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(opac1d(nd))
   allocate(opatot1d(nd))
   allocate(sline1d(nd))
   allocate(scont1d(nd))
   allocate(stot1d(nd))
   allocate(int1d(nd))
!
!---------------calculate control points and coefficients---------------
!-----------------------for interpolation along s-----------------------
!
!derivative weights for opacity
   alpha_opalbar = dels_d/s_d
   alpha_opac = alpha_opalbar
!derivative weights for source function
   alpha_sline = max(wp_interp1d,alpha_opalbar)
   alpha_scont = alpha_opalbar
!control-point weights from weighted derivatives
   wu_sline = alpha_sline/2.d0
   wp_sline = ((2.d0-alpha_sline)*dels_d + (1.d0-alpha_sline)*dels_u)/2.d0/dels_d
   wd_sline = (alpha_sline-1.d0)*dels_u/2.d0/dels_d
   wu_scont = alpha_scont/2.d0
   wp_scont = ((2.d0-alpha_scont)*dels_d + (1.d0-alpha_scont)*dels_u)/2.d0/dels_d
   wd_scont = (alpha_scont-1.d0)*dels_u/2.d0/dels_d

   wu_opalbar = alpha_opalbar/2.d0
   wp_opalbar = ((2.d0-alpha_opalbar)*dels_d + (1.d0-alpha_opalbar)*dels_u)/2.d0/dels_d
   wd_opalbar = (alpha_opalbar-1.d0)*dels_u/2.d0/dels_d
   wu_opac = alpha_opac/2.d0
   wp_opac = ((2.d0-alpha_opac)*dels_d + (1.d0-alpha_opac)*dels_u)/2.d0/dels_d
   wd_opac = (alpha_opac-1.d0)*dels_u/2.d0/dels_d

!(non-monotonic) control points
   opalbar_c = wu_opalbar*opalbar_u + wp_opalbar*opalbar_p + wd_opalbar*opalbar_d
   opac_c = wu_opac*opac_u + wp_opac*opac_p + wd_opac*opac_d
   sline_c = wu_sline*sline_u + wp_sline*sline_p + wd_sline*sline_d
   scont_c = wu_scont*scont_u + wp_scont*scont_p + wd_scont*scont_d
!ensure monotonicity
   call pointcl1d_mbez(opalbar_c, opalbar_u, opalbar_p)
   call pointcl1d_mbez(opac_c, opac_u, opac_p)
   call coeffcl1d_mbez(scont_c, scont_u, scont_p, wu_scont, wp_scont, wd_scont)
   call coeffcl1d_mbez(sline_c, sline_u, sline_p, wu_sline, wp_sline, wd_sline)
!
!----------------------------start point--------------------------------
!
   alo_u = 0.d0
   alo_p = 0.d0
   alo_d = 0.d0
   alo_ucont = 0.d0
   alo_pcont = 0.d0
   alo_dcont = 0.d0
   delt_u = 0.d0
!
!at point u (<=> iim1)
   xcmf_iim1 = xcmf_u
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim1 = s_iim1/s_p
   a_iim1 = 1.d0-b_iim1
!(bezier) interpolation coefficients
   ts_iim1 = (s_iim1-s_u)/dels_u
   at_iim1 = (1.d0-ts_iim1)**2
   bt_iim1 = 2.d0*ts_iim1*(1.d0-ts_iim1)
   ct_iim1 = ts_iim1**2
   atline_iim1 = at_iim1 + bt_iim1*wu_sline
   btline_iim1 = ct_iim1 + bt_iim1*wp_sline
   ctline_iim1 = bt_iim1*wd_sline
   atcont_iim1 = at_iim1 + bt_iim1*wu_scont
   btcont_iim1 = ct_iim1 + bt_iim1*wp_scont
   ctcont_iim1 = bt_iim1*wd_scont
!actual interpolation
   sline_iim1 = atline_iim1*sline_u + btline_iim1*sline_p + ctline_iim1*sline_d
   scont_iim1 = atcont_iim1*scont_u + btcont_iim1*scont_p + ctcont_iim1*scont_p
   opalbar_iim1 = at_iim1*opalbar_u + bt_iim1*opalbar_c + ct_iim1*opalbar_p
   opac_iim1 = at_iim1*opac_u + bt_iim1*opac_c + ct_iim1*opac_p
   vel_iim1 = a_iim1*vel_u + b_iim1*vel_p
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1    = opalbar_iim1*phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   stot_iim1 = opal_iim1*sline_iim1/opatot_iim1 + opac_iim1*scont_iim1/opatot_iim1

   int_iim1 = int_u
!
   s1d(1) = s_iim1
   tau1d(1) = 0.d0
   vel1d(1) = vel_iim1
   vth1d(1) = vth_mean
   xcmf1d(1) = xcmf_iim1
   profile1d(1) = phinorm_iim1
   opalbar1d(1) = opalbar_iim1
   opal1d(1) = opal_iim1
   opac1d(1) = opac_iim1
   opatot1d(1) = opatot_iim1
   sline1d(1) = sline_iim1
   scont1d(1) = scont_iim1
   stot1d(1) = stot_iim1
   int1d(1) = int_u
!
!at point ii
   xcmf_ii = xcmf_u+del_xcmf
   s_ii = s_u + (xcmf_ii-xcmf_u)/m
!(linear) interpolation coefficients
   b_ii = s_ii/s_p
   a_ii = 1.d0-b_ii
!(bezier) interpolation coefficients
   ts_ii = (s_ii-s_u)/dels_u
   at_ii = (1.d0-ts_ii)**2
   bt_ii = 2.d0*ts_ii*(1.d0-ts_ii)
   ct_ii = ts_ii**2
   atline_ii = at_ii + bt_ii*wu_sline
   btline_ii = ct_ii + bt_ii*wp_sline
   ctline_ii = bt_ii*wd_sline
   atcont_ii = at_ii + bt_ii*wu_scont
   btcont_ii = ct_ii + bt_ii*wp_scont
   ctcont_ii = bt_ii*wd_scont
!actual interpolation
   sline_ii = atline_ii*sline_u + btline_ii*sline_p + ctline_ii*sline_d
   scont_ii = atcont_ii*scont_u + btcont_ii*scont_p + ctcont_ii*scont_d
   opalbar_ii = at_ii*opalbar_u + bt_ii*opalbar_c + ct_ii*opalbar_p
   opac_ii = at_ii*opac_u + bt_ii*opac_c + ct_ii*opac_p
   vel_ii = a_ii*vel_u + b_ii*vel_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii    = opalbar_ii*phinorm_ii
   opatot_ii  = opal_ii + opac_ii
   stot_ii = opal_ii*sline_ii/opatot_ii + opac_ii*scont_ii/opatot_ii
!
   delt_ii = (opatot_iim1+opatot_ii)*(s_ii-s_iim1)/2.d0
!
!-----------------------------------------------------------------------
!
   do i=2, nd-1
!calculate spatial grid (from s proportional to xcmf, i.e., assuming a liner velocity law)
!interpolate opacities, source functions, and velocities
      xcmf_iip1 = xcmf_u + i*del_xcmf
      s_iip1 = s_u + (xcmf_iip1-xcmf_u)/m

!(linear) interpolation coefficients
      b_iip1 = s_iip1/s_p
      a_iip1 = 1.d0-b_iip1
!(bezier) interpolation coefficients
      ts_iip1 = (s_iip1-s_u)/dels_u
      at_iip1 = (1.d0-ts_iip1)**2
      bt_iip1 = 2.d0*ts_iip1*(1.d0-ts_iip1)
      ct_iip1 = ts_iip1**2
      atline_iip1 = at_iip1 + bt_iip1*wu_sline
      btline_iip1 = ct_iip1 + bt_iip1*wp_sline
      ctline_iip1 = bt_iip1*wd_sline
      atcont_iip1 = at_iip1 + bt_iip1*wu_scont
      btcont_iip1 = ct_iip1 + bt_iip1*wp_scont
      ctcont_iip1 = bt_iip1*wd_scont
!actual interpolation
      sline_iip1 = atline_iip1*sline_u + btline_iip1*sline_p + ctline_iip1*sline_d
      scont_iip1 = atcont_iip1*scont_u + btcont_iip1*scont_p + ctcont_iip1*scont_d
      opalbar_iip1 = at_iip1*opalbar_u + bt_iip1*opalbar_c + ct_iip1*opalbar_p
      opac_iip1 = at_iip1*opac_u + bt_iip1*opac_c + ct_iip1*opac_p
      vel_iip1 = a_iip1*vel_u + b_iip1*vel_p
      call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
      opal_iip1    = opalbar_iip1*phinorm_iip1
      opatot_iip1  = opal_iip1 + opac_iip1
      stot_iip1 = opal_iip1*sline_iip1/opatot_iip1 + opac_iip1*scont_iip1/opatot_iip1
!
!calculate delta tau-step
      delt_iip1 = (opatot_ii+opatot_iip1)*(s_iip1-s_ii)/2.d0 
!
!calculate coefficients for source function integration
      call coeff_source2(delt_ii, delt_iip1, stot_iim1, stot_ii, stot_iip1, a, b, c)
!
      delt_u = delt_u + delt_ii
      expdtau_ii=exp(-delt_ii)
      int_ii=int_iim1*expdtau_ii + a*stot_iim1 + b*stot_ii + c*stot_iip1
!
      alo_u = alo_u*expdtau_ii + a*atline_iim1*opal_iim1/opatot_iim1 + &
                                 b*atline_ii*opal_ii/opatot_ii + &
                                 c*atline_iip1*opal_iip1/opatot_iip1
      alo_p = alo_p*expdtau_ii + a*btline_iim1*opal_iim1/opatot_iim1 + &
                                 b*btline_ii*opal_ii/opatot_ii + &
                                 c*btline_iip1*opal_iip1/opatot_iip1
      alo_d = alo_d*expdtau_ii + a*ctline_iim1*opal_iim1/opatot_iim1 + &
                                 b*ctline_ii*opal_ii/opatot_ii + &
                                 c*ctline_iip1*opal_iip1/opatot_iip1

      alo_ucont = alo_ucont*expdtau_ii + a*atcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*atcont_ii*opac_ii/opatot_ii + &
                                         c*atcont_iip1*opac_iip1/opatot_iip1
      alo_pcont = alo_pcont*expdtau_ii + a*btcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*btcont_ii*opac_ii/opatot_ii + &
                                         c*btcont_iip1*opac_iip1/opatot_iip1
      alo_dcont = alo_dcont*expdtau_ii + a*ctcont_iim1*opac_iim1/opatot_iim1 + &
                                         b*ctcont_ii*opac_ii/opatot_ii + &
                                         c*ctcont_iip1*opac_iip1/opatot_iip1
!
      int1d(i) = int_ii
      s1d(i) = s_ii
      tau1d(i) = tau1d(i-1)+delt_ii
      vel1d(i) = vel_ii
      vth1d(i) = vth_mean
      xcmf1d(i) = xcmf_ii
      profile1d(i) = phinorm_ii
      opalbar1d(i) = opalbar_ii
      opac1d(i) = opac_ii
      opal1d(i) = opal_ii
      opatot1d(i) = opatot_ii
      sline1d(i) = sline_ii
      scont1d(i) = scont_ii
      stot1d(i) = stot_ii

!referesh variables at iim1
!      xcmf_iim1 = xcmf_ii
!      s_iim1 = s_ii
!      phinorm_iim1 = phinorm_ii
!      opalbar_iim1 = opalbar_ii
!      sline_iim1 = sline_ii
!      scont_iim1 = scont_ii
!      vel_iim1 = vel_ii
      opac_iim1 = opac_ii
      opatot_iim1 = opatot_ii
      opal_iim1 = opal_ii
      stot_iim1 = stot_ii
      int_iim1=int_ii
      atline_iim1 = atline_ii
      btline_iim1 = btline_ii
      ctline_iim1 = ctline_ii
      atcont_iim1 = atcont_ii
      btcont_iim1 = btcont_ii
      ctcont_iim1 = ctcont_ii
!
!referesh variables at ii
      xcmf_ii = xcmf_iip1
      phinorm_ii = phinorm_iip1
      opalbar_ii = opalbar_iip1
      sline_ii = sline_iip1
      scont_ii = scont_iip1
      vel_ii = vel_iip1
      s_ii = s_iip1
      opac_ii = opac_iip1
      opatot_ii = opatot_iip1
      opal_ii = opal_iip1
      stot_ii = stot_iip1
      delt_ii = delt_iip1
!
      atline_ii = atline_iip1
      btline_ii = btline_iip1
      ctline_ii = ctline_iip1
      atcont_ii = atcont_iip1
      btcont_ii = btcont_iip1
      ctcont_ii = ctcont_iip1


   enddo
!
!-------------------------interval [n-1,n]------------------------------
!
!recalculate everything at at n-2
   xcmf_iim2 = xcmf_p - 2*del_xcmf
   xcmf_iim1 = xcmf_p - del_xcmf
   s_iim2 = s_u + (xcmf_iim2-xcmf_u)/m
   s_iim1 = s_u + (xcmf_iim1-xcmf_u)/m
!(linear) interpolation coefficients
   b_iim2 = (s_iim2-s_u)/dels_u
   a_iim2 = 1.d0-b_iim2
!(bezier) interpolation coefficients
   ts_iim2 = (s_iim2-s_u)/dels_u
   at_iim2 = (1.d0-ts_iim2)**2
   bt_iim2 = 2.d0*ts_iim2*(1.d0-ts_iim2)
   ct_iim2 = ts_iim2**2
   atline_iim2 = at_iim2 + bt_iim2*wu_sline
   btline_iim2 = ct_iim2 + bt_iim2*wp_sline
   ctline_iim2 = bt_iim2*wd_sline
   atcont_iim2 = at_iim2 + bt_iim2*wu_scont
   btcont_iim2 = ct_iim2 + bt_iim2*wp_scont
   ctcont_iim2 = bt_iim2*wd_scont
!actual interpolation
   sline_iim2 = atline_iim2*sline_u + btline_iim2*sline_p + ctline_iim2*sline_d
   scont_iim2 = atcont_iim2*scont_u + btcont_iim2*scont_p + ctcont_iim2*scont_d
   opalbar_iim2 = at_iim2*opalbar_u + bt_iim2*opalbar_c + ct_iim2*opalbar_p
   opac_iim2 = at_iim2*opac_u + bt_iim2*opac_c + ct_iim2*opac_p
   vel_iim2 = a_iim2*vel_u + b_iim2*vel_p
   call calc_phinorm(vel_iim2, vth_mean, vth_fiducial, xobs, phinorm_iim2)
   opal_iim2 = opalbar_iim2*phinorm_iim2
   opatot_iim2  = opal_iim2 + opac_iim2
   stot_iim2 = opal_iim2*sline_iim2/opatot_iim2 + opac_iim2*scont_iim2/opatot_iim2
!
!delta-tau step (assuming opal to be linear in (small) refined interval)
   delt_iim1 = (opatot_iim2+opatot_iim1)*(s_iim1-s_iim2)/2.d0
!
!calculate coefficients for source function integration
   call coeff_source2b(delt_iim1, delt_ii, stot_iim2, stot_iim1, stot_ii, a, b, c)

   expdtau_ii=exp(-delt_ii)
   delt_u = delt_u + delt_ii
!
   int_ii=int_iim1*expdtau_ii + a*stot_iim2 + b*stot_iim1 + c*stot_ii
!
   alo_u = alo_u*expdtau_ii + a*atline_iim2*opal_iim2/opatot_iim2 + &
                              b*atline_iim1*opal_iim1/opatot_iim1 + &
                              c*atline_ii*opal_ii/opatot_ii
   alo_p = alo_p*expdtau_ii + a*btline_iim2*opal_iim2/opatot_iim2 + &
                              b*btline_iim1*opal_iim1/opatot_iim1 + &
                              c*btline_ii*opal_ii/opatot_ii
   alo_d = alo_d*expdtau_ii + a*ctline_iim2*opal_iim2/opatot_iim2 + &
                              b*ctline_iim1*opal_iim1/opatot_iim1 + &
                              c*ctline_ii*opal_ii/opatot_ii

   alo_ucont = alo_ucont*expdtau_ii + a*atcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*atcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*atcont_ii*opac_ii/opatot_ii
   alo_pcont = alo_pcont*expdtau_ii + a*btcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*btcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*btcont_ii*opac_ii/opatot_ii
   alo_dcont = alo_dcont*expdtau_ii + a*ctcont_iim2*opac_iim2/opatot_iim2 + &
                                      b*ctcont_iim1*opac_iim1/opatot_iim1 + &
                                      c*ctcont_ii*opac_ii/opatot_ii
!
   abs_sc = exp(-delt_u)
   contr_sc = alo_u*sline_u + alo_p*sline_p + alo_d*sline_d + &
              alo_ucont*scont_u + alo_pcont*scont_p + alo_dcont*scont_d
   int_sc = int_u*abs_sc + contr_sc
!
!
!
   i=nd
   int1d(i) = int_ii
   s1d(i) = s_ii
   tau1d(i) = tau1d(i-1)+delt_ii
   vel1d(i) = vel_ii
   vth1d(i) = vth_mean
   xcmf1d(i) = xcmf_ii
   profile1d(i) = phinorm_ii
   opalbar1d(i) = opalbar_ii
   opac1d(i) = opac_ii
   opal1d(i) = opal_ii
   opatot1d(i) = opatot_ii
   sline1d(i) = sline_ii
   scont1d(i) = scont_ii
   stot1d(i) = stot_ii
!
else
!
!for point iim1
   s_iim1       = s_u
   vel_iim1     = vel_u
   opalbar_iim1 = opalbar_u
   opac_iim1    = opac_u
   call calc_phinorm(vel_iim1, vth_mean, vth_fiducial, xobs, phinorm_iim1)
   opal_iim1    = opalbar_iim1*phinorm_iim1
   opatot_iim1  = opal_iim1 + opac_iim1
   sline_iim1   = sline_u
   scont_iim1    = scont_u
   stot_iim1    = sline_iim1*opal_iim1/opatot_iim1 + scont_iim1*opac_iim1/opatot_iim1
!
!for point ii
   s_ii       = s_p
   vel_ii     = vel_p
   opalbar_ii = opalbar_p
   opac_ii    = opac_p
   call calc_phinorm(vel_ii, vth_mean, vth_fiducial, xobs, phinorm_ii)
   opal_ii    = opalbar_ii*phinorm_ii
   opatot_ii  = opal_ii + opac_ii
   sline_ii   = sline_p
   scont_ii   = scont_p
   stot_ii    = sline_ii*opal_ii/opatot_ii + scont_ii*opac_ii/opatot_ii
!for point ii
   s_iip1       = s_d
   vel_iip1     = vel_d
   opalbar_iip1 = opalbar_d
   opac_iip1    = opac_d
   call calc_phinorm(vel_iip1, vth_mean, vth_fiducial, xobs, phinorm_iip1)
   opal_iip1    = opalbar_iip1*phinorm_iip1
   opatot_iip1  = opal_iip1 + opac_iip1
   sline_iip1   = sline_d
   scont_iip1   = scont_d
   stot_iip1    = sline_iip1*opal_iip1/opatot_iip1 + scont_iip1*opac_iip1/opatot_iip1
!
!calculate delta-tau
   delt_u = (opatot_iim1+opatot_ii)*(s_p-s_u)/2.d0
   delt_d = (opatot_iip1+opatot_ii)*(s_d-s_p)/2.d0
!
   abs_sc = exp(-delt_u)
   call coeff_source2(delt_u, delt_d, stot_iim1, stot_ii, stot_iip1, a, b, c)

   alo_u = a*opal_iim1/opatot_iim1
   alo_p = b*opal_ii/opatot_ii
   alo_d = c*opal_iip1/opatot_iip1

!   contr_sc = a *stot_iim1 + b*stot_ii + c*stot_iip1
    contr_sc = alo_u*sline_iim1 + alo_p*sline_ii + alo_d*sline_iip1 + &                                                  !line contribution
               a*opac_iim1/opatot_iim1*scont_iim1 + b*opac_ii/opatot_ii*scont_ii + c*opac_iip1/opatot_iip1*scont_iip1    !continuum contribution
!
   int_sc = int_u*abs_sc + contr_sc
!
!
   nd=2
   allocate(s1d(nd))
   allocate(tau1d(nd))
   allocate(vel1d(nd))
   allocate(vth1d(nd))
   allocate(xcmf1d(nd))
   allocate(profile1d(nd))
   allocate(opalbar1d(nd))
   allocate(opal1d(nd))
   allocate(opac1d(nd))
   allocate(opatot1d(nd))
   allocate(sline1d(nd))
   allocate(scont1d(nd))
   allocate(stot1d(nd))
   allocate(int1d(nd))
   s1d = (/ s_u, s_p /)
   tau1d = (/ 0.d0, delt_u /)
   vel1d = (/ vel_u, vel_p /)
   vth1d = (/ vth_mean, vth_mean /)
   xcmf1d = (/ xcmf_u, xcmf_p /)
   profile1d = (/ phinorm_iim1, phinorm_ii /)
   opalbar1d = (/ opalbar_iim1, opalbar_ii /)
   opal1d = (/ opal_iim1, opal_ii /)
   opac1d = (/ opac_iim1, opac_ii /)
   opatot1d = (/ opatot_iim1, opatot_ii /)
   sline1d = (/ sline_iim1, sline_ii /)
   scont1d = (/ scont_iim1, scont_ii /)
   stot1d = (/ stot_iim1, stot_ii /)
   int1d = (/ int_u, int_sc /)
!
endif
!
!
end subroutine fsc_linec_debug
!
!***********************************************************************
!***********************************************************************
!
!                     MISCELLANEOUS
!
!***********************************************************************
!***********************************************************************
!
subroutine coeffcl1d_mbez(fc, f_im1, f_i, at, bt, ct)
!
!calculates new coefficients for 2d bezier interpolation to ensure monotonicity
!                 when interpolating in left interval
!
!on input: 
!   fc           - control point
!   f_im1, f_i   - function values that limit the interval
!   at, bt, ct   - coefficients that have been used for calculation of control point
!
!on output:
!   at, bt, ct - interpolation coefficients for control point sucht that monotonic interpolation can be used
!                   (e.g. at=0.,bt=1.,ct=0. if f_c > f_i > f_im1)
!
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: fc, f_im1, f_i
real(dp), intent(inout) :: at, bt, ct
!
if((f_i-f_im1)*(f_i-fc).le.0.) then
   at=0.
   bt=1.
   ct=0.
elseif((f_i-f_im1)*(fc-f_im1).le.0.) then
   at=1.
   bt=0.
   ct=0.
endif
!
end subroutine coeffcl1d_mbez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine pointcl1d_mbez (fc, f_im1, f_i)
!
!calculates new control point for bezier interpolation to ensure monotonicity
!                 when interpolating in left interval
!
!on input: 
!   fc           - current control point
!   f_im1, f_i   - function values that limit the interval
!
!on output:
!   fc - new control point such that monotonicity is ensured
!               (e.g. fc_new = f_i, if f_c > f_i > f_im1)
!
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: f_im1, f_i
real(dp), intent(inout) :: fc
!
!
!
if((f_i-f_im1)*(f_i-fc).le.0.) then
   fc=f_i
elseif((f_i-f_im1)*(fc-f_im1).le.0.) then
   fc=f_im1
endif
!
end subroutine pointcl1d_mbez
