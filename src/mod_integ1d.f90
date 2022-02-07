module mod_integ1d

  use prog_type
  use fund_const

  implicit none

  public :: integral_source0

contains
  !
  !***********************************************************************
  !
  !   for integration of source terms: weighted by del-tau
  !
  !***********************************************************************
  !
  function integral_source0(xim1, xi, fim1, fi)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1,xi
    !        with funct(x) approximated by average value: (fim1+fi)/2.
    !
    !
    !   input: function values:        fi, fim1
    !          at x-values:            xi
    !   output: evaluated integral:    integral_source0
    !
    !                    fim1--------fi
    !                    xim1--------xi
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, fim1, fi
    real(dp) :: integral_source0
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delt_i
    !
    delt_i=xi-xim1
    integral_source0 = (fim1+fi)/two * (one-exp(-delt_i))
    !
  end function integral_source0
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral_source1(xim1, xi, fim1, fi)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        with funct(x) approximated linearly by g(x): g(x)=a1*(x-xim1)+b1
    !                                                or   g(x)=a2*(x-xi)+b2
    !                                                 (integrals are the same)
    !
    !   input: function values:        fi, fim1
    !          at x-values:            xi
    !   output: evaluated integral:    integral_source1
    !
    !                    fim1--------fi
    !                    xim1--------xi
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, fim1, fi
    real(dp) :: integral_source1
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delt_i, dum
    !
    delt_i=xi-xim1
    !
    if(delt_i.le.1.d-6) then
       !thin limit
       integral_source1=delt_i*fim1
    elseif(delt_i.ge.20.) then
       !thick limit
       integral_source1=fi-(fi-fim1)/delt_i
    else
       dum=exp(-delt_i)
       integral_source1 = (one/delt_i-dum-dum/delt_i)*fim1 + (one-one/delt_i+dum/delt_i)*fi
    endif
    !
    !
  end function integral_source1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral_source2(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between 0., xi
    !        with funct(x) approximated by bezier curve g(x):
    !                g(x)=fim1*(1-q)^2 + fi*q^2 + 2.*fc*q*(1-q)
    !                  q=(x-xi)/(xi-xim1)
    !                  fc from weighted mean (see Hayek 2010, Auer 2003)
    !
    !
    !   input: function values:        fi, fim1
    !          at x-values:            xi
    !   output: evaluated integral:    integral1
    !
    !                    fim1--------fi-------fip1
    !                      0.--------xi-------xip1
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral_source2
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: fpi, fpb, fpf, delt_i, delt_ip1, delt_ii, u0, u1, u2, fc, fmin, fmax
    !
    !delta-taus
    delt_i=xi-xim1
    delt_ip1=xip1-xi
    delt_ii=xip1-xim1

    if(delt_i.le.1.d-6) then
    !thin limit
       integral_source2=delt_i*fim1
    elseif(delt_i.ge.20.) then
    !thick limit
       integral_source2=fi-(fi-fim1)/delt_i
    else
    !
    !derivative at point xi
       fpi = (fi-fim1)*delt_ip1/delt_i/delt_ii + (fip1-fi)*delt_i/delt_ip1/delt_ii
       fc = fi - delt_i*fpi/two
       fmin=min(fim1,fi)
       fmax=max(fim1,fi)
       if(fc.gt.fmax) then
          fc=fmax
       else if(fc.lt.fmin) then
          fc=fmin
       endif
    
       u0=one-exp(-delt_i)
       u1=delt_i - u0
       u2=delt_i**2 - two*u1
       integral_source2 = u0*fim1 + two*(fc-fim1)*u1/delt_i + (fi+fim1-two*fc)*u2/delt_i**2
    endif
    !
  end function integral_source2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral_source3(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        with funct(x) approximated quadratically by g(x): g(x)=a*(x-xi)**2 + b*(x-xi) + c
    !
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: evaluated integral:    integral_source3
    !
    !                    fim1--------fi---------fip1
    !                    xim1--------xi---------xip1
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral_source3
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delt_i, delt_ip1, delf_i, delf_ip1, a, b, c, dum
    !
    delt_i=xi-xim1
    delt_ip1=xip1-xi
    !
    delf_i=fi-fim1
    delf_ip1=fip1-fi
    !
    !write(*,'(a20, 3es20.8)') 'integral', xim1, xi, xip1
    !write(*,'(a20, 3es20.8)') 'integral', fim1, fi, fip1
    !write(*,'(a20, 2es20.8)') 'integral', delt_i, delt_ip1
    !
    if(delt_i.le.1.d-6) then
       !thin limit
       integral_source3=delt_i*fim1
    elseif(delt_i.ge.20.) then
       !thick limit
       integral_source3=fi-(fi-fim1)/delt_i
    else
       !calculate coefficients (numerically not the fastests method, but for the moment...)
       !   if(delt_ip1.eq.0..or.delt_i.eq.0..or.delt_ip1.eq.-delt_i) then
       !      write(*,*) xim1, xi, xip1
       !      write(*,*) delt_ip1, delt_i
       !   endif
       !avoid division by zero (e.g. if delt_i gt 1.d-6 and delt_ip1->0.)
       delt_ip1=max(delt_ip1,small_number)
       a=(delf_ip1/delt_ip1 - delf_i/delt_i)/(delt_ip1+delt_i)
       b=(delf_ip1*delt_i/delt_ip1 + delf_i*delt_ip1/delt_i)/(delt_ip1+delt_i)
       c=fi
       dum=exp(-delt_i)
       integral_source3 = c*(one-dum) + b*(delt_i*dum-one+dum) + a*(-delt_i**2*dum - two*delt_i*dum+two-two*dum)
    endif
    !
    !write(*,*) 'integral'
    !write(*,*) delt_i, delt_ip1
    !write(*,*) delf_ip1/delt_ip1, delf_i/delt_i, delt_ip1+delt_i
    !write(*,*) a, b, c
    !write(*,*) one-dum, (delt_i-one+dum), (delt_i**2 - two*delt_i+two-two*dum)
    !write(*,*) 'integral done'
    !stop
    !
  end function integral_source3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source1(dx, acoeff, bcoeff)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi
    !     funct(x) approximated linearly by g(x): g(x)=a*(x-xi) + b
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: coefficients for integration: acoeff, bcoeff
    !
    !                    fim1--------fi
    !                    xim1--------xi
    !                       dx=xi-xim1
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: acoeff, bcoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1

    if(dx.le.1.d-6) then
       acoeff = dx
       bcoeff = zero
    elseif(dx.ge.20.) then
       acoeff = one/dx
       bcoeff = one-acoeff
    else
       e0 = one-exp(-dx)
       e1 = dx-e0
       !first order integrals
       acoeff = e0 - e1/dx
       bcoeff = e1/dx
    endif
    !
  end subroutine coeff_source1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source2(dxi, dxip1, fim1, fi, fip1, acoeff, bcoeff, ccoeff)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi + ccoeff*fip1
    !     with funct(x) approximated by bezier curve g(x):
    !                g(x)=fim1*(1-q)^2 + fi*q^2 + 2.*fc*q*(1-q)
    !                  q=(x-xim1)/(xi-xim1)
    !                  fc from weighted mean (see Hayek 2010, Auer 2003)
    !
    !   input: function values:        fi, fim1
    !          at x-values:            dxi, dxip1
    !   output: coefficients           acoeff, bcoeff, ccoeff
    !
    !                    fim1-------------fi-------------fip1
    !                    xim1-------------xi-------------xip1
    !                         dxi=xi-xim1   dxip1=xip1-xi
    !
    !-----------------------------------------------------------------------
    !
    use mod_interp2d, only: wp_integ1d
    !
    ! ... arguments
    real(dp), intent(in) :: dxi, dxip1, fim1, fi, fip1
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1, e2, dx, fpi, fc, wi, wim1, wip1, mm, mp, norm, alpha
    real(dp), parameter :: fac=three/four, q=two
    !
    if(dxi.le.1.d-6) then
       acoeff = dxi
       bcoeff = zero
       ccoeff = zero
    elseif(dxip1.le.1.d-6) then
       !linear approach to avoid division by zero
       e0 = one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1 = dxi-e0
       acoeff = e0 - e1/dxi
       bcoeff = e1/dxi
       ccoeff = zero
       !elseif(dxi.ge.20.) then
       !   dx = dxi + dxip1
       !quadratic function
       !   acoeff = (two+dxip1)/dxi/dx
       !   bcoeff = one + (dxi-dxip1-two)/dxi/dxip1
       !   ccoeff = (two-dxi)/dxip1/dx
       !monotonic bezier (with predefined derivative weights)
       !   alpha = max(wp_integ1d,dxip1/dx)
       !   wim1 = alpha/two
       !   acoeff = -wim1
       !   bcoeff = one
       !   ccoeff = zero
       !linear approach
       !   acoeff = one/dxi
       !   bcoeff = one-acoeff
       !   ccoeff = zero
    else
       !
       dx=dxi+dxip1
       !
       !10: control point from derivative (monotonic)
       !20: control point from derivative (non-monotonic)
       !30: control point from predefined weights
       !40: control point from inverse distance weighting
       !50: control point from derivative (monotonic, with predefined weights for derivative)
       goto 20
       !
       !--------------control point from derivative (monotonic)----------------
       !
10     continue
       !derivative at point xi
       fpi = (fi-fim1)*dxip1/dxi/dx + (fip1-fi)*dxi/dxip1/dx
       fc = fi - dxi*fpi/two
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       !
       !standard interpolation
       acoeff = e0 + (e2-e1*(dxip1+two*dxi))/dxi/dx
       bcoeff = (e1*dx-e2)/dxi/dxip1
       ccoeff = (e2-e1*dxi)/dxip1/dx
       if(fi.ge.fim1) then
          if(fc.gt.fi) then
             !set fc to fi to ensure monotonicity
             bcoeff = two*e1/dxi - e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          elseif(fc.lt.fim1) then
             !set fc to fim1 to ensure monotonicity
             bcoeff = e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          endif
       elseif(fi.le.fim1) then
          if(fc.lt.fi) then
             !set fc to fi to ensure monotonicity
             bcoeff = two*e1/dxi - e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          elseif(fc.gt.fim1) then
             !set fc to fim1 to ensure monotonicity
             bcoeff = e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          endif
       endif
       return
       !
       !------------control point from derivative (non-monotonic)--------------
       !
20     continue
       !derivative at point xi
       fpi = (fi-fim1)*dxip1/dxi/dx + (fip1-fi)*dxi/dxip1/dx
       fc = fi - dxi*fpi/two
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       !
       !standard interpolation
       acoeff = e0 + (e2-e1*(dxip1+two*dxi))/dxi/dx
       bcoeff = (e1*dx-e2)/dxi/dxip1
       ccoeff = (e2-e1*dxi)/dxip1/dx
       !   write(*,*) 'coeff_source2a', acoeff, bcoeff, ccoeff
       !   goto 50
       !
       return
       !
       !-------------------control point from predefined weight----------------
       !
30     continue
       mm=abs((fi-fim1)/dxi)
       mp=abs((fip1-fi)/dxip1)
       if(mm.gt.mp) then
          wi=fac
          wim1=one-fac
       elseif(mm.lt.mp) then
          wi=one-fac
          wim1=fac
       else
          wim1=half
          wi=half
       endif
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       !
       bcoeff=two*e1*wi/dxi + (one-two*wi)*e2/dxi**2
       acoeff=e0-bcoeff
       ccoeff=zero
       return
       !
       !----------------control point from inverse distance weighting----------
       !
40     continue
       wim1=one/(half*dxi)**q
       wi=one/(half*dxi)**q
       wip1=one/(half*dxi+dxip1)**q
       norm=wim1+wi+wip1
       wim1=wim1/norm
       wi=wi/norm
       wip1=wip1/norm
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       !
       acoeff=e0+two*e1*(wim1-one)/dxi + e2*(one-two*wim1)/dxi**2
       bcoeff=two*wi*e1/dxi + e2*(one-two*wi)/dxi**2
       ccoeff=two*wip1*(e1/dxi-e2/dxi**2)
       return
       !
       !--------------control point from derivative (monotonic)----------------
       !----------------with predefined weight for derivative------------------
       !
50     continue
       !derivative at point xi
       alpha = max(wp_integ1d,dxip1/dx)
       !   write(*,*) wp_integ1d, dxip1, dx, dxip1/dx
       !   stop
       wim1 = alpha/two
       wi = ((two-alpha)*dxip1+(one-alpha)*dxi)/two/dxip1
       wip1 = (alpha-one)*dxi/two/dxip1
       !
       fc = wim1*fim1+wi*fi+wip1*fip1
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       !
       !standard interpolation
       acoeff=e0+two*e1*(wim1-one)/dxi + e2*(one-two*wim1)/dxi**2
       bcoeff=two*wi*e1/dxi + e2*(one-two*wi)/dxi**2
       ccoeff=two*wip1*(e1/dxi-e2/dxi**2)

       if(fi.ge.fim1) then
          if(fc.gt.fi) then
             !set fc to fi to ensure monotonicity
             bcoeff = two*e1/dxi - e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          elseif(fc.lt.fim1) then
             !set fc to fim1 to ensure monotonicity
             bcoeff = e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          endif
       elseif(fi.le.fim1) then
          if(fc.lt.fi) then
             !set fc to fi to ensure monotonicity
             bcoeff = two*e1/dxi - e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          elseif(fc.gt.fim1) then
             !set fc to fim1 to ensure monotonicity
             bcoeff = e2/dxi**2
             acoeff = e0-bcoeff
             ccoeff = zero
          endif
       endif

       !   stop 'debug thick problem, or when ration dxi/dxip1 tooo large/low'
       !   write(*,'(6es20.8)') acoeff, bcoeff, ccoeff, dxi, dxip1, dxi/dxip1!
       !    write(*,*) 'coeff_source2b', acoeff, bcoeff, ccoeff
       !   return

    endif
    !
    !
  end subroutine coeff_source2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source2b(dxi, dxip1, fim1, fi, fip1, acoeff, bcoeff, ccoeff)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xi, xip1
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi + ccoeff*fip1
    !     with funct(x) approximated by bezier curve g(x):
    !                g(x)=fi*(1-q)^2 + fip1*q^2 + 2.*fc*q*(1-q)
    !                  q=(x-xi)/(xip1-xi)
    !                  fc from weighted mean (see Hayek 2010, Auer 2003)
    !
    !   input: function values:        fi, fim1
    !          at x-values:            dxi, dxip1
    !   output: coefficients           acoeff, bcoeff, ccoeff
    !
    !                    fim1-------------fi-------------fip1
    !                    xim1-------------xi-------------xip1
    !                         dxi=xi-xim1   dxip1=xip1-xi
    !
    !-----------------------------------------------------------------------
    !
    use mod_interp2d, only: wp_integ1d
    !
    implicit none
    !
    ! ... arguments
    real(dp), intent(in) :: dxi, dxip1, fim1, fi, fip1
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1, e2, dx, fpi, fc, wi, wim1, wip1, mm, mp, norm, alpha
    !
    if(dxip1.le.1.d-6) then
       acoeff = zero
       bcoeff = dxip1
       ccoeff = zero
    elseif(dxi.le.1.d-6) then
       !linear approach to avoid division by zero
       e0 = one-exp(-dxip1)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1 = dxip1-e0
       acoeff = zero
       bcoeff = e0 - e1/dxip1
       ccoeff = e1/dxip1
    else
       !
       dx=dxi+dxip1
       !
       !derivative at point xi
       alpha = max(wp_integ1d,dxi/dx)
       wim1 = (alpha-one)*dxip1/two/dxi
       wi = ((two-alpha)*dxi+(one-alpha)*dxip1)/two/dxi
       wip1 = alpha/two
       !
       fc = wim1*fim1+wi*fi+wip1*fip1
       !
       e0=one-exp(-dxip1)
       if(dxip1.ge.20.) e0=one   !avoid floating overflow
       e1=dxip1 - e0
       e2=dxip1**2 - two*e1
!
       !standard interpolation
       acoeff = two*wim1*(e1/dxip1-e2/dxip1**2)
       bcoeff = e0 + two*e1*(wi-one)/dxip1 + e2*(one-two*wi)/dxip1**2
       ccoeff = two*wip1*e1/dxip1 + (one-two*wip1)*e2/dxip1**2

       if(fip1.ge.fi) then
          if(fc.gt.fip1) then
             !set fc to fip1 to ensure monotonicity
             acoeff = zero
             ccoeff = two*e1/dxip1 - e2/dxip1**2
             bcoeff = e0-ccoeff
          elseif(fc.lt.fi) then
             !set fc to fi to ensure monotonicity
             acoeff = zero
             ccoeff = e2/dxip1**2
             bcoeff = e0-bcoeff
          endif
       elseif(fip1.le.fi) then
          if(fc.lt.fip1) then
             !set fc to fip1 to ensure monotonicity
             acoeff = zero
             ccoeff = two*e1/dxip1 - e2/dxip1**2
             bcoeff = e0-ccoeff
          elseif(fc.gt.fi) then
             !set fc to fi to ensure monotonicity
             acoeff = zero
             ccoeff = e2/dxip1**2
             bcoeff = e0-bcoeff
          endif
       endif
       return

    endif
    !
    !
  end subroutine coeff_source2b
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source3(dxi, dxip1, fim1, fi, fip1, acoeff, bcoeff, ccoeff)
    !
    !-----------------------------------------------------------------------
    !

    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi + ccoeff*fip1
    !     funct(x) approximated quadratically by g(x): g(x)=a*(x-xim1)**2 + b*(x-xim1) + c
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: coefficients for integral: acoeff, bcoeff, ccoeff
    !
    !                    fim1------------fi-------------fip1
    !                    xim1------------xi-------------xip1
    !                        dxi=xi-xim1    dxip1=xip1-xi
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: dxi, dxip1, fim1, fi, fip1
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1, e2, dx, dfdx_im1, dfdx_i

    
    !stop 'coeff_source3 to be debugged'
    
    if(dxi.le.1.d-6) then
       acoeff = dxi
       bcoeff = zero
       ccoeff = zero
    elseif(dxi.ge.20.) then
       dx = dxi + dxip1
       acoeff = (two+dxip1)/dxi/dx
       bcoeff = one + (dxi-dxip1-two)/dxi/dxip1
       ccoeff = (two-dxi)/dxip1/dx
       !   acoeff = one/dxi
       !   bcoeff = one-acoeff
       !   ccoeff = zero
    else
       dx = dxi + dxip1
       dfdx_im1 = fim1*(-two*dxi-dxip1)/dxi/dx + fi*dx/dxi/dxip1 - fip1*dxi/dxip1/dx
       dfdx_i   = fim1*(-dxip1)/dxi/dx + fi*(dxip1-dxi)/dxi/dxip1 + fip1*dxi/dxip1/dx
       if(dfdx_im1*dfdx_i.lt.zero) then
          !first order integrals
          e0 = one-exp(-dxi)
          e1 = dxi-e0
          acoeff = e0 - e1/dxi
          bcoeff = e1/dxi
          ccoeff = zero
       else
          e0 = one-exp(-dxi)
          e1 = dxi-e0
          e2 = dxi**2 - two*e1
          acoeff = e0 + (e2-(two*dxi+dxip1)*e1)/dxi/dx
          bcoeff = (dx*e1-e2)/dxi/dxip1
          ccoeff = (e2-dxi*e1)/dxip1/dx
       endif
    endif
    !
  end subroutine coeff_source3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source3b(dxi, dxip1, fim1, fi, fip1, acoeff, bcoeff, ccoeff)
    !
    !-----------------------------------------------------------------------
    !
    
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi + ccoeff*fip1
    !     funct(x) approximated quadratically by g(x): g(x)=a*(x-xim1)**2 + b*(x-xim1) + c
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: coefficients for integral: acoeff, bcoeff, ccoeff
    !
    !                    fim1------------fi-------------fip1
    !                    xim1------------xi-------------xip1
    !                        dxi=xi-xim1    dxip1=xip1-xi
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: dxi, dxip1, fim1, fi, fip1
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1, e2, dx, dfdx_im1, dfdx_i
    

    !stop 'coeff_source3 to be debugged'

    if(dxi.le.1.d-6) then
       acoeff = dxi
       bcoeff = zero
       ccoeff = zero
    elseif(dxi.ge.20.) then
       dx = dxi + dxip1
       acoeff = (two+dxip1)/dxi/dx
       bcoeff = one + (dxi-dxip1-two)/dxi/dxip1
       ccoeff = (two-dxi)/dxip1/dx
       !   acoeff = one/dxi
       !   bcoeff = one-acoeff
       !   ccoeff = zero
    else
       dx = dxi + dxip1
       dfdx_im1 = fim1*(-two*dxi-dxip1)/dxi/dx + fi*dx/dxi/dxip1 - fip1*dxi/dxip1/dx
       dfdx_i   = fim1*(-dxip1)/dxi/dx + fi*(dxip1-dxi)/dxi/dxip1 + fip1*dxi/dxip1/dx
       e0 = one-exp(-dxi)
       e1 = dxi-e0
       e2 = dxi**2 - two*e1
       acoeff = e0 + (e2-(two*dxi+dxip1)*e1)/dxi/dx
       bcoeff = (dx*e1-e2)/dxi/dxip1
       ccoeff = (e2-dxi*e1)/dxip1/dx
    endif
    !
  end subroutine coeff_source3b
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff_source4(dxi, acoeff, bcoeff)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x)=funct(x)*exp(x-xi) between xim1, xi
    !        such that integral(f(x)) = acoeff*fim1 + bcoeff*fi
    !     with funct(x) approximated by 3rd order polynomial with zero derivatives at xim1, xi
    !
    !   input: function values:        fi, fim1
    !          at x-values:            dxi
    !   output: coefficients           acoeff, bcoeff
    !
    !                    fim1-------------fi
    !                    xim1-------------xi
    !                         dxi=xi-xim1
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: dxi
    real(dp), intent(out) :: acoeff, bcoeff
    !
    ! ... local scalars
    real(dp) :: e0, e1, e2, e3
    !
    if(dxi.le.1.d-6) then
       acoeff = dxi
       bcoeff = zero
    else
       !
       e0=one-exp(-dxi)
       if(dxi.ge.20.) e0=one   !avoid floating overflow
       e1=dxi - e0
       e2=dxi**2 - two*e1
       e3=dxi**3 - three*e2
       !
       !standard interpolation
       acoeff = e0 - three*e2/dxi**2 + two*e3/dxi**3
       bcoeff = e0 - acoeff
    endif
    !
    !
  end subroutine coeff_source4
  !
  !***********************************************************************
  !
  !      normal integration functions
  !
  !***********************************************************************
  !
  function integral0(xim1, xi, fim1, fi)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xim1,xi
    !     with f(x) approximated by average value: (fim1+fi)/2.
    !
    !
    !   input: function values:        fi, fim1
    !          at x-values:            xi
    !   output: evaluated integral:    integral0
    !
    !                    fim1--------fi
    !                    xim1--------xi
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, fim1, fi
    real(dp) :: integral0
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delt_i
    !
    integral0 = (fim1+fi)/two * (xi-xim1)
    !
  end function integral0
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral1(xim1, xi, fim1, fi)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xim1,xi
    !     with f(x) approximated linearly by: f(x)=m*x+t
    !
    !
    !   input: function values:        fi, fim1
    !          at x-values:            xi
    !   output: evaluated integral:    integral0
    !
    !                    fim1--------fi
    !                    xim1--------xi
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, fim1, fi
    real(dp) :: integral1
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    integral1 = fim1*xi - fi*xim1 + half*(fi-fim1)*(xi**2 - xim1**2)/(xi-xim1)
    !
  end function integral1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral2(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xim1,xi
    !     with f(x) approximated quadratic bezier spline:
    !            f(x)=fim1*(1-q)^2 + fi*q^2 + 2.*fc*q*(1-q)
    !               q=(x-xim1)/(xi-xim1)
    !               fc from weighted mean (see Hayek 2010, Auer 2003)
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: evaluated integral:    integral2
    !
    !                    fim1--------fi-------fip1
    !                    xim1--------xi-------xip1
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral2
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: dx_i, dx_ip1, dx, fpi, fc
    !
    !different path lengths
    dx_i=xi-xim1
    dx_ip1=xip1-xi
    dx=xip1-xim1
    !
    if(fi.eq.fim1) then
       !linear function to be integrated
       integral2 = (fi+fim1)*dx_i/two
    else
       !derivative at point xi
       fpi = (fi-fim1)*dx_ip1/dx_i/dx + (fip1-fi)*dx_i/dx_ip1/dx
       fc = fi - dx_i*fpi/two
       if((fc-fim1)/(fim1-fi).ge.zero) then
          !set fc to fim1 to ensure monotonicity
          integral2=(two*fim1+fi)*dx_i/three
       elseif((fc-fi)/(fi-fim1).ge.zero) then
          !set fc to fi to ensure monotonicity
          integral2=(fim1+two*fi)*dx_i/three
       else
          integral2 = (fim1+fc+fi)*dx_i/three
       endif
    endif
    !
  end function integral2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral2b(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xi,xip1
    !     with f(x) approximated quadratic bezier spline:
    !            f(x)=fi*(1-q)^2 + fip1*q^2 + 2.*fc*q*(1-q)
    !               q=(x-xim1)/(xi-xim1)
    !               fc from weighted mean (see Hayek 2010, Auer 2003)
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: evaluated integral:    integral2
    !
    !                    fim1--------fi-------fip1
    !                    xim1--------xi-------xip1
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral2b
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: dx_i, dx_ip1, dx, fpi, fc
    !
    !different path lengths
    dx_i=xi-xim1
    dx_ip1=xip1-xi
    dx=xip1-xim1
    !
    if(fi.eq.fip1) then
       !linear function to be integrated
       integral2b = (fi+fip1)*dx_ip1/two
    else
       !derivative at point xi
       fpi = (fi-fim1)*dx_ip1/dx_i/dx + (fip1-fi)*dx_i/dx_ip1/dx
       fc = fi + dx_ip1*fpi/two
       if((fc-fi)/(fi-fip1).ge.zero) then
          !set fc to fi to ensure monotonicity
          integral2b=(two*fi+fip1)*dx_ip1/three
       elseif((fc-fip1)/(fip1-fi).ge.zero) then
          !set fc to fip1 to ensure monotonicity
          integral2b=(fi+two*fip1)*dx_ip1/three
       else
          integral2b = (fi+fc+fip1)*dx_ip1/three
       endif
    endif
    !
  end function integral2b
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral3(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xim1,xi
    !     with f(x) approximated by quadratic function: 
    !            f(x)=a*(x-xi)^2 + b*(x-xi) + c
    !        or (if not monotonic) by linear function:
    !            f(x)=b*(x-xi)+c
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: evaluated integral:    integral2
    !
    !                    fim1--------fi-------fip1
    !                    xim1--------xi-------xip1
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral3
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delx_i, delx_ip1, delf_i, delf_ip1, a, b, c, dfdx_i, dfdx_im1
    !
    !different deltas
    delx_i=xi-xim1
    delx_ip1=xip1-xi
    !
    delf_i=fi-fim1
    delf_ip1=fip1-fi
    !
    !calculate coefficients (numerically not the fastests method, but for the moment...)
    a=(delf_ip1/delx_ip1 - delf_i/delx_i)/(delx_ip1+delx_i)
    b=(delf_ip1*delx_i/delx_ip1 + delf_i*delx_ip1/delx_i)/(delx_ip1+delx_i)
    c=fi
    !
    !ensure monotonic interpolation (take linear approximation if non-monotonic)
    dfdx_im1=-two*a*delx_i + b
    dfdx_i=b
    if(dfdx_im1*dfdx_i.lt.zero) then
       b=delf_i/delx_i
       integral3=c*delx_i - b*delx_i**2/two
    else
       integral3 = a*delx_i**3/three - b*delx_i**2/two + c*delx_i
    endif
    !
  end function integral3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function integral3b(xim1, xi, xip1, fim1, fi, fip1)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function f(x) between xi,xip1
    !     with f(x) approximated by quadratic function:
    !            f(x)=a*(x-xi)^2 + b*(x-xi) + c
    !
    !   input: function values:        fi, fim1, fip1
    !          at x-values:            xi, xim1, xip1
    !   output: evaluated integral:    integral2
    !
    !                    fim1--------fi-------fip1
    !                    xim1--------xi-------xip1
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, xip1, fim1, fi, fip1
    real(dp) :: integral3b
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: delx_i, delx_ip1, delf_i, delf_ip1, a, b, c, dfdx_i, dfdx_ip1
    !
    !different deltas
    delx_i=xi-xim1
    delx_ip1=xip1-xi
    !
    delf_i=fi-fim1
    delf_ip1=fip1-fi
    !
    !calculate coefficients (numerically not the fastests method, but for the moment...)
    a=(delf_ip1/delx_ip1 - delf_i/delx_i)/(delx_ip1+delx_i)
    b=(delf_ip1*delx_i/delx_ip1 + delf_i*delx_ip1/delx_i)/(delx_ip1+delx_i)
    c=fi
    !
    !ensure monotonic interpolation (take linear approximation if non-monotonic)
    dfdx_i=b
    dfdx_ip1=two*a*delx_ip1+b
    !
    if(dfdx_ip1*dfdx_i.lt.zero) then
       b=delf_ip1/delx_ip1
       integral3b=c*delx_ip1 + b*delx_ip1**2/two
    else
       integral3b = a*delx_ip1**3/three + b*delx_ip1**2/two + c*delx_ip1
    endif
    !
    !
  end function integral3b
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine integral4(xim1, xi, fim1, f2, f3, f4, fi, sum, rerr)
    !
    !calculates integral of a function f between [sim1,si]
    !
    !method: simpsons-rule corrected for error
    !
    !on input:
    !x-values at interval boundaries
    !   xim1-----------------------xi
    !f-values at 5 points, including the interval-boundaries
    !   AND   assuming an equidistant x-step size!!!
    !   fim1-----f2----f3----f4----fi
    !
    !
    ! ... arguments
    real(dp), intent(in) :: xim1, xi, fim1, f2, f3, f4, fi
    !
    ! ... local scalars
    real(dp), parameter :: tol=1.d-2
    real(dp) :: sum, sum1, sum2, h, rerr
    !
    !
    !
    h=(xi-xim1)/four
    sum1=fim1+four*f2+two*f3+four*f4+fi
    sum2=(-fim1+four*f2-six*f3+four*f4-fi)/fifteen
    !
    sum=h*(sum1+sum2)/three
    rerr=abs(sum2/(sum1+sum2))
    !
    !if(integral4.eq.zero) then
    !   rerr=zero
    !else
    !   rerr=abs(sum2/(sum1+sum2))
    !endif
    !if(rerr.gt.tol) then
    !   write(*,*) 'error in integral4, relative error of integration is:', rerr
    !   stop
    !endif
    !
  end subroutine integral4
  !
  !***********************************************************************
  !***********************************************************************
  !
  !                  integration using trapezoidal rule
  !                   integration using simpson's rule
  !               including error estiamtes in both cases
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine precalc_weight_trapez(x, n, weight)
    !
    !------calculates integration weights for weight function g(x)=1--------
    !
    !
    !   input: dimension n
    !          x-grid
    !
    !   output: w= weights for integration
    !
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: x, weight
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    ! ... local arrays
    real(dp), dimension(n) :: weightm, weightp
    !
    ! ... local characters
    
    do i=2, n
       weightm(i)=x(i)-x(i-1)
       weightp(i)=x(i)-x(i-1)
    end do
    !
    weightm=weightm/two
    weightp=weightp/two
    !
    !
    weight(1)=weightm(2)
    !
    do i=2, n-1
       weight(i)=weightp(i)+weightm(i+1)
    end do
    !
    weight(n)=weightp(n)
    !
    !
    !
  end subroutine precalc_weight_trapez
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_oweight_trapez(x, n, a, b, weight)
    !
    !---------calculates integration weights for trapezoidal rule-----------
    !------with open boundaries (values at a,b not included in nodes)-------
    !---------------assuming that f(a)=f(1) and f(b)=f(n)-------------------
    !
    !   input: dimension n
    !          x-grid
    !          boundary of integration: a, b
    !
    !   output: w= weights for integration
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n), intent(out) :: weight
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: dxip1, dxi
    !
    ! ... local arrays
    !
    ! ... local characters
    
    dxi=x(1)-a
    dxip1=x(2)-x(1)
    weight(1)=dxi+dxip1/two
    !
    do i=2, n-1
       dxi=dxip1
       dxip1=x(i+1)-x(i)
       weight(i) = (dxi+dxip1)/two
    end do
    !
    dxi=x(n)-x(n-1)
    dxip1=b-x(n)
    weight(n) = dxi/two + dxip1
    !
    !
    !
  end subroutine precalc_oweight_trapez
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_weight_trapez_err(x, n, neq, w1, w2, werr)
    !
    !------calculates integration weights for weight function g(x)=1--------
    !--------------including weights for error estimates--------------------
    !
    !   input: dimension n
    !          number of equidistant grid points neq (subintervals)
    !          x-grid
    !
    !   output: w1= weights for integration using complete grid
    !           w2= weights for integration using only every second grid point
    !           werr= error weights
    !
    ! ... arguments
    integer(i4b), intent(in) :: n, neq
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n) :: w1, w2, werr
    !
    ! ... local scalars
    integer(i4b) :: i,j
    real(dp) :: h1, h2
    !
    ! ... local arrays
    real(dp), dimension(n) :: weightm, weightp
    !
    ! ... local characters
    !
    !------------for trapezoidal rule: neq needs to be factor of 2----------
    !
    if(mod(neq,2).ne.0) stop 'error in precalc_weight_trapez: neq has to be factor of 2'
    if(neq.ne.2) then
       if(mod(neq/2,2).ne.0) stop 'error in precalc_weight_trapez: neq/2 has to be factor of 2'
    endif
    !
    !--------------check if neq following intervals are equidistant---------
    !
    do i=2, n, neq
       h1=x(i)-x(i-1)
       do j=1, neq
          h2=x(i+j-1)-x(i+j-2)
          !      write(*,'(2i5, 6e30.16)') i, j, x(i), x(i-1), x(i+j-1), x(i+j-2), h1, h2
          if(abs(h2-h1).gt.small_number) stop 'error in precalc_weight_trapez: grid not equidistant'
       enddo
    enddo
    !
    !----------------------half step weights--------------------------------
    !
    h1=(x(3)-x(1))/two
    w1(1)=h1/two
    w1(2)=h1
    werr(1)=-w1(1)/three
    werr(2)=w1(2)/three
    !
    do i=3, n-1, 2      
       h1=(x(i)-x(i-2))/two
       h2=(x(i+2)-x(i))/two
       w1(i) = h1/two + h2/two
       w1(i+1)= h2
       werr(i) = -w1(i)/three
       werr(i+1) = w1(i+1)/three
    enddo
    w1(n)=(x(n) - x(n-2))/four
    werr(n) = -w1(n)/three
    !
    !--------------------- full step weights--------------------------------
    !
    h1=x(3)-x(1)
    w2(1)=h1/two
    w2(2)=zero
    !
    do i=3, n-1, 2
       h1=(x(i)-x(i-2))/two
       h2=(x(i+2)-x(i))/two
       w2(i)=h1+h2
       w2(i+1)=zero
    enddo
    !
    w2(n)=(x(n)-x(n-2))/two
    !
    !
  end subroutine precalc_weight_trapez_err
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_weight_simps(x, n, weight)
    !
    !----------------calculates integration weights for---------------------
    !---------------simpson's 3-point integration formula-------------------
    !
    !
    !   input: dimension n
    !          x-grid
    !
    !   output: weight= weights for integration
    !
    !   note: x-array is overwritten if it is not equidistant yet
    !            for three subsequent integration-nodes
    !
    !
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: x, weight
    !
    ! ... local scalars
    integer(i4b) :: i, indx
    real(dp) :: del_m1, del_p1
    !
    ! ... local arrays
    !
    ! ... local characters
    !
    !-------------check if dimensions are okay (n has to be odd)------------
    !
    if(n.lt.3) stop 'error in precalc_weight_simps:: n has to be ge 3'
    if(mod(n,2).eq.0) stop 'error in precalc_weights_simps: n has to be odd'
    !
    !-------------make grid equidistant for three subsequent nodes-----------
    do i=2, n, 2
       x(i) = (x(i+1)+x(i-1))/two
    enddo
    !
    !------check if three subsequent x-points are really equidistant--------
    !
    do i=2, n-1, 2
       !has to round to 15 digits, otherwise machine-precision errors
       del_m1=x(i)-x(i-1)
       del_p1=x(i+1)-x(i)
       if(abs(del_m1-del_p1).gt.small_number) then
          stop 'error in precalc_weight_simps: x-grid not equidistant'
       endif
    enddo
    !
    !----------------------calculate weights--------------------------------
    !
    if(n.eq.3) then
       weight(1) = (x(3)-x(1))/six
       weight(2) = four*weight(1)
       weight(3) = weight(1)
    else
       weight(1) = (x(3)-x(1))/six
       do i=2, n-1
          if(mod(i,2).eq.0) then
             weight(i) = (x(i+1)-x(i-1))*four/six
          else
             weight(i) = (x(i+2)-x(i-2))/six
          endif
       enddo
       weight(n) = (x(n)-x(n-2))/six
    endif
    !
    !
  end subroutine precalc_weight_simps
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_weight_boole(x, n, w)
    !
    !----------------calculates integration weights for---------------------
    !---------------simpson's 3-point integration formula-------------------
    !--on grid x, and on grid with double resolution to increase accuracy---
    !-----------------and to calculate the error weights--------------------
    !-------------------(also known as boole's rule)------------------------
    !
    !   input: dimension n
    !          x-grid
    !
    !   output: w: weights for integration using complete grid
    !           werr: weights for error estimation
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(inout) :: x
    real(dp), dimension(n), intent(out) :: w!, werr
    !
    ! ... local scalars
    integer(i4b) :: i, j
    real(dp) :: h1, h2
    !
    ! ... local arrays
    !
    ! ... local characters
    !
    !-----------------only a number of 5+ni*4  grid points allowed-------
    !--------------------(with ni the number of intervals)---------------
    !
    if(mod(n-5,4).ne.0) stop 'error in precalc_weight_boole: n-5 has to be factor of 4'
    !
    !-------------make grid equidistant for five subsequent nodes-----------
    !
    do i=3, n, 4
       h1 = (x(i+2)-x(i-2))/four
       x(i-1) = x(i-2)+h1
       x(i) = x(i-1)+h1
       x(i+1) = x(i)+h1
    enddo
    !
    !--------------check if equidistant spacing within subintervals---------
    !
    do i=2, n, 4
       h1=x(i)-x(i-1)
       do j=1, 4
          h2=x(i+j-1)-x(i+j-2)
          !      write(*,'(2i5, 6e20.8)') i, j, x(i), x(i-1), x(i+j-1), x(i+j-2), h1, h2
          if(abs(h2-h1).gt.small_number) stop 'error in precalc_oweight_boole: grid not equidistant'
       enddo
    enddo
    !
    !------------------------complete grid----------------------------------
    !
    h1=(x(5)-x(1))/four
    w(1)=h1*fourteen/fortyfive
    w(2)=h1*sixtyfour/fortyfive
    w(3)=h1*twentyfour/fortyfive
    w(4)=h1*sixtyfour/fortyfive
    !
    do i=5, n-4, 4
       h1=(x(i)-x(i-4))/four
       h2=(x(i+4)-x(i))/four
       w(i) = (h1+h2)*fourteen/fortyfive
       w(i+1)=h2*sixtyfour/fortyfive
       w(i+2)=h2*twentyfour/fortyfive
       w(i+3)=h2*sixtyfour/fortyfive
    enddo
    h1=(x(n)-x(n-4))/four
    w(n)=h1*fourteen/fortyfive
    !
    !
    !
  end subroutine precalc_weight_boole
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_oweight_simps(x, n, a, b, weight)
    !
    !----------------calculates integration weights for---------------------
    !---------------simpson's 3-point integration formula-------------------
    !------with open boundaries (values at a,b not included in nodes)-------
    !---------------assuming linear extrapolation at the edges--------------
    !
    !   input: dimension n
    !          x-grid
    !          boundary of integration: a, b
    !
    !   output: w= weights for integration
    !
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), dimension(n), intent(inout) :: x
    real(dp), dimension(n), intent(out) :: weight
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: del_m1, del_p1, dxip1, dxi, hi
    !
    ! ... local arrays
    !
    ! ... local characters
    !
    !-------------check if dimensions are okay (n has to be odd)------------
    !
    if(n.lt.3) stop 'error in precalc_oweight_simps:: n has to be ge 3'
    if(mod(n,2).eq.0) stop 'error in precalc_oweights_simps: n has to be odd'
    !
    !-------------make grid equidistant for three subsequent nodes-----------
    do i=2, n, 2
       x(i) = (x(i+1)+x(i-1))/two
    enddo
    !
    !------check if three subsequent x-points are really equidistant--------
    !
    do i=2, n-1, 2
       !has to round to 15 digits, otherwise machine-precision errors
       del_m1=x(i)-x(i-1)
       del_p1=x(i+1)-x(i)
       if(abs(del_m1-del_p1).gt.small_number) then
          stop 'error in precalc_oweight_simps: x-grid not equidistant'
       endif
    enddo
    !
    !----------------------calculate weights--------------------------------
    !
    if(n.eq.3) then
       hi=x(3)-x(1)
       dxi=x(1)-a
       dxip1=b-x(3)
       weight(1) = hi/six + dxi + dxi**2/hi
       weight(2) = four*hi/six - dxi**2/hi - dxip1**2/hi
       weight(3) = hi/six + dxip1 + dxip1**2/hi
    else
       hi=x(3)-x(1)
       dxi=x(1)-a
       weight(1) = hi/six + dxi + dxi**2/hi
       weight(2) = four*hi/six - dxi**2/hi
       do i=3, n-2
          if(mod(i,2).eq.0) then
             weight(i) = (x(i+1)-x(i-1))*four/six
          else
             weight(i) = (x(i+1)-x(i-1))/three
          endif
       enddo
       hi=x(n)-x(n-2)
       dxi=b-x(n)
       weight(n-1) = four*hi/six - dxi**2/hi
       weight(n) = hi/six + dxi + dxi**2/hi
    endif
    !
    !
    !
  end subroutine precalc_oweight_simps
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_oweight_boole(x, n, a, b, w)
    !
    !----------------calculates integration weights for---------------------
    !---------------simpson's 3-point integration formula-------------------
    !--on grid x, and on grid with double resolution to increase accuracy---
    !-----------------and to calculate the error weights--------------------
    !------with open boundaries (values at a,b not included in nodes)-------
    !---------------assuming linear extrapolation at the edges--------------
    !
    !   input: dimension n
    !          x-grid
    !          integration limits [a,b]
    !
    !   output: w: weights for integration using complete grid
    !           werr: weights for error estimation
    !
    ! ... arguments
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), dimension(n), intent(inout) :: x
    real(dp), dimension(n), intent(out) :: w!, werr
    !
    ! ... local scalars
    integer(i4b) :: i, j
    real(dp) :: h1, h2, hp, hq
    !
    ! ... local arrays
    !
    ! ... local characters
    !
    !-----------------only a number of 5+ni*4  grid points allowed-------
    !--------------------(with ni the number of intervals)---------------
    !
    if(mod(n-5,4).ne.0) stop 'error in precalc_oweight_boole: n-5 has to be factor of 4'
    !
    !-------------make grid equidistant for five subsequent nodes-----------
    !
    do i=3, n, 4
       h1 = (x(i+2)-x(i-2))/four
       x(i-1) = x(i-2)+h1
       x(i) = x(i-1)+h1
       x(i+1) = x(i)+h1
    enddo
    !
    !--------------check if equidistant spacing within subintervals---------
    !
    do i=2, n, 4
       h1=x(i)-x(i-1)
       do j=1, 4
          h2=x(i+j-1)-x(i+j-2)
          !      write(*,'(2i5, 6e20.8)') i, j, x(i), x(i-1), x(i+j-1), x(i+j-2), h1, h2
          if(abs(h2-h1).gt.small_number) stop 'error in precalc_oweight_boole: grid not equidistant'
       enddo
    enddo
    !
    !------------------------complete grid----------------------------------
    !
    h1=(x(5)-x(1))/four
    hp=x(1)-a
    w(1)=h1*fourteen/fortyfive + hp+hp**2/two/h1
    w(2)=h1*sixtyfour/fortyfive - hp**2/two/h1
    w(3)=h1*twentyfour/fortyfive
    w(4)=h1*sixtyfour/fortyfive
    !
    do i=5, n-4, 4
       h1=(x(i)-x(i-4))/four
       h2=(x(i+4)-x(i))/four
       w(i) = (h1+h2)*fourteen/fortyfive
       w(i+1)=h2*sixtyfour/fortyfive
       w(i+2)=h2*twentyfour/fortyfive
       w(i+3)=h2*sixtyfour/fortyfive
    enddo
    h1=(x(n)-x(n-4))/four
    hq=b-x(n)
    w(n-1)=w(n-1)-hq**2/two/h1
    w(n)=h1*fourteen/fortyfive+hq+hq**2/two/h1
    !
    !
    !
  end subroutine precalc_oweight_boole
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_weight_simps_err(x, n, neq, w1, w2, werr)
    !
    !----------------calculates integration weights for---------------------
    !---------------simpson's 3-point integration formula-------------------
    !--------------including calculation of error weights-------------------
    !
    !   input: dimension n
    !          number of equidistant grid points neq (subintervals)
    !          x-grid
    !
    !   output: w1= weights for integration using complete grid
    !           w2= weights for integration using only every second grid point
    !           werr= weights for error estimation
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: n, neq
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n) :: w1, w2, werr
    !
    ! ... local scalars
    integer(i4b) :: i, j
    real(dp) :: h2, h1
    !
    ! ... local arrays
    !
    ! ... local characters
    !
    !----------------------neq has to be odd number-------------------------
    !
    if(mod(neq,4).ne.0) stop 'error in precalc_weight_simps: neq has to be factor of 4'
    !
    !--------------check if neq following intervals are equidistant---------
    !
    do i=2, n, neq
       h1=x(i)-x(i-1)
       do j=1, neq
          h2=x(i+j-1)-x(i+j-2)
          !      write(*,'(2i5, 6e20.8)') i, j, x(i), x(i-1), x(i+j-1), x(i+j-2), h1, h2
          if(abs(h2-h1).gt.small_number) stop 'error in precalc_weight_simps_err: grid not equidistant'
       enddo
    enddo
    !
    !------------------------complete grid----------------------------------
    !
    h1=(x(5)-x(1))/four
    w1(1)=h1/three
    w1(2)=four*h1/three
    w1(3)=two*h1/three
    w1(4)=four*h1/three
    !
    werr(1)=-w1(1)/fifteen
    werr(2)=w1(2)/fifteen
    werr(3)=-three*w1(3)/fifteen
    werr(4)=w1(4)/fifteen
    !
    do i=5, n-4, 4
       h1=(x(i)-x(i-4))/four
       h2=(x(i+4)-x(i))/four
       w1(i) = h1/three + h2/three
       w1(i+1)=4*h2/three
       w1(i+2)=2*h2/three
       w1(i+3)=4*h2/three
       !
       werr(i) = -w1(i)/fifteen
       werr(i+1) = w1(i+1)/fifteen
       werr(i+2) = -three*w1(i+2)/fifteen
       werr(i+3) = w1(i+3)/fifteen
       !   write(*,'(i5,  5e20.8)') i, x(i-4), x(i), x(i+4), h1, h2
    enddo
    w1(n)=(x(n) - x(n-4))/twelve
    werr(n) = -w1(n)/fifteen
    !
    !---------------------------half grid-----------------------------------
    !
    h1=(x(5)-x(1))/two
    w2(1)=h1/three
    w2(2)=zero
    w2(3)=four*h1/three
    w2(4)=zero

    do i=5, n-4, 4
       h1=(x(i)-x(i-4))/two
       h2=(x(i+4)-x(i))/two
       w2(i) = h1/three + h2/three
       w2(i+1) = zero
       w2(i+2) = four * h2/three
       w2(i+3) = zero
    enddo
    w2(n) = (x(n)-x(n-4))/six
    !
    !
    !
  end subroutine precalc_weight_simps_err
  !
  !***********************************************************************
  !***********************************************************************
  !
  !                     spline integration routines
  !                    legendre integration routines
  !                    chebyshev integration routines
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine precalc_weight_spline(x, n, weight, periodic)
    !
    !----------------calculates integration weights for---------------------
    !-----------cubic spline integration, where first derivatives-----------
    !------are approximated by central differences: catmull-rom-spline------
    !
    !   input: dimension n
    !          x-grid
    !          periodic: logical which is true, when periodic boundary conditions are used
    !                                    false, when boundary derivative is backward difference
    !
    !   output: weight= weights for integration
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: n
    logical, intent(in) :: periodic
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n) :: weight
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: dum1, dum2
    !
    ! ... local arrays
    real(dp), dimension(n) :: h
    !
    ! ... local characters
    !
    !
    if(n.lt.4) stop 'error in precalc_weight_spline: n has to be ge 4'
    !
    !-------------------store array of delta-x previously-------------------
    !
    h(1)=zero
    do i=2, n
       h(i) = x(i)-x(i-1)
    enddo
    !
    !----------------------calculate weights-on boundaries------------------
    !
    if(periodic) then
       !inner boundary condition periodic: central difference between grid point x(2) and x(n-1)
       !(x(n)=x(1))
       dum1 = twelve * (h(3) + h(2))
       weight(1) = h(2)/two + (h(2)*h(2)-h(3)*h(3))/dum1
       !
       dum1 = twelve * (h(2) + h(n))
       dum2 = twelve * (h(4) + h(3))
       weight(2) = (h(2)*h(2)-h(n)*h(n))/dum1 + (h(3)*h(3)-h(4)*h(4))/dum2 + (h(3)+h(4))/two
       !
       dum1 = twelve * (h(n-1)+h(n-2))
       dum2 = twelve * (h(2) + h(n))
       weight(n-1) = (h(n-1)*h(n-1)-h(n-2)*h(n-2))/dum1 + (h(n)*h(n)-h(2)*h(2))/dum2 + (h(n)+h(n-1))/two
       !
       dum1 = twelve * (h(n) + h(n-1))
       weight(n) = h(n)/two + (h(n)*h(n)-h(n-1)*h(n-1))/dum1
    else
       !make backward differences on boundaries
       dum1 = twelve * (h(3) + h(2))
       weight(1) = five*h(2)/twelve + h(2)*h(2)/dum1 - h(3)*h(3)/dum1
       !
       dum1 = twelve * (h(4) + h(3))
       weight(2) = seven*h(2)/twelve + h(3)/two + h(3)*h(3)/dum1 - h(4)*h(4)/dum1
       !
       dum1 = twelve * (h(n-1) + h(n-2))
       weight(n-1) = seven*h(n)/twelve + h(n-1)/two + h(n-1)*h(n-1)/dum1 - h(n-2)*h(n-2)/dum1
       !
       dum1 = twelve * (h(n) + h(n-1))
       weight(n) = five*h(n)/twelve + h(n)*h(n)/dum1 - h(n-1)*h(n-1)/dum1
       !
    endif
    !
    !-----------------calculate weights away from boundaries----------------
    !
    do i=3, n-2
       !
       dum1 = twelve * (h(i) + h(i-1))
       dum2 = twelve * (h(i+2) + h(i+1))
       weight(i) = (h(i)*h(i) - h(i-1)*h(i-1))/dum1 + (h(i)+h(i+1))/two + &
            (h(i+1)*h(i+1) - h(i+2)*h(i+2))/dum2             
       !
    enddo
    !
  end subroutine precalc_weight_spline
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine precalc_weight_legendre(nd, xlower, xupper, nodes,  weights)
    !
    !-----------------------------------------------------------------------
    !
    !         calculating nodes and weights of legendre polynomials 
    !                 for gauss-legendre-integration
    !
    !              nodes are eigenvalues of tridiagonal matrix
    !  weights can be calculated from first value of normalized eigenvector:
    !                       weight=2*ev(1)**2
    !
    !on input: nd: dimension of arrays
    !          xupper, xlower: integration bounds
    !
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xlower, xupper
    real(dp), dimension(nd) :: nodes, weights
    !
    ! ... local scalars
    integer(i4b) :: i, j, err
    real(dp), parameter :: upper_boundary=two*pi, lower_boundary=zero
    real(dp), dimension(nd) :: eigenvalues, subdiag
    real(dp), dimension(nd, nd) :: eigenvectors
    !
    !
    !-----------setting up tridiagonal-matrix / input for imtql2------------
    !
    eigenvalues = zero
    !
    subdiag(1)=zero
    do i=1, nd-1
       subdiag(i+1) = (one*i)/(sqrt(four*(one*i)**2 - one))
    end do
    !
    eigenvectors=zero
    !
    do i=1, nd
       eigenvectors(i,i) = one
    end do
    !
    !---------------calculating eigenvalues of matrix------------------------
    !---------------calculating eigenvectors of matrix-----------------------
    !
    call imtql2 (nd, eigenvalues, subdiag, eigenvectors, err)
    !
    if (err.gt.0) then
       write(*,*) 'error precalc_weight_legendre'
       stop
    end if
    !
    !---------------calculating nodes (=eigenvalues)-------------------------
    !---------------calculating weights from eigenvectors--------------------
    !
    do i=1, nd
       if(abs(eigenvalues(i)).lt.small_number) then
          nodes(i)=zero
       else
          nodes(i)=eigenvalues(i)
       endif
    enddo
    !
    !
    do i=1, nd
       weights(i)=two*eigenvectors(1,i)**2
    end do
    !
    !-------------variable transformation on range [xlower, xupper]---------
    !
    do i=1, nd
       nodes(i)=nodes(i)*(xupper-xlower)/two + (xupper+xlower)/two
       weights(i)=weights(i)*(xupper-xlower)/two
    enddo
    !
    !
    !
  end subroutine precalc_weight_legendre
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine  precalc_weight_chebyshev(nd, xlower, xupper, nodes,  weights)
    !
    !-----------------------------------------------------------------------
    !         calculating nodes and weights of chebycheff polynomials 
    !                 for gauss-chebyshev-integration
    !
    !on input: nd: dimension of arrays
    !          xupper, xlower: integration bounds
    !
    !note: although nodes are calculated in the first instance for an integral
    !         from [-1,1], the nodes are from [1.,-1]
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xlower, xupper
    real(dp), dimension(nd) :: nodes, weights
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    ! ... local arrays
    real(dp), dimension(nd) :: weightfct
    !
    !
    do i=1, nd
       nodes(i)=cos(((two * i - one)/(two * nd))*pi)
       if(abs(nodes(i)).lt.small_number) then
          nodes(i)=zero
       endif
    end do
    !
    !--------------------store weight-function at nodes---------------------
    !
    weightfct=sqrt(one-nodes*nodes)
    !
    !-----------------------calculate weights-------------------------------
    !
    weights=pi/(nd*one)
    !
    weights=weights*weightfct
    !
    !-------------variable transformation on range [xlower,xupper]----------
    !
    do i=1, nd
       nodes(i)=nodes(i)*(xupper-xlower)/two + (xupper+xlower)/two
       weights(i)=weights(i)*(xupper-xlower)/two
    enddo
    !
    !
    !
  end subroutine precalc_weight_chebyshev

  !***********************************************************************
  !***********************************************************************
  !
  !                  integration routines with weights
  !                           (pauldrach)
  
  !                   different weight functions:
  !               integral1: g(x)=1
  !               integral2: g(x)=x
  !               integral3: g(x)=sin(x)
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine integral_w1(x, y, n, sum)
    !
    !---------------------integrates a function f(x)=y----------------------
    !----------------------with weight-function g(x)=1----------------------
    !
    !   input: dimension n
    !          x-grid
    !          y-values
    !
    !   output: int= integral of y dx from x(1) to x(n)
    !
    !
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: x,y
    real(dp) :: sum
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    ! ... local arrays
    real(dp), dimension(n) :: weightm, weightp, weight
    !
    ! ... local characters
    !character(len=50) :: enter
    
    do i=2, n
       weightm(i)=x(i)-x(i-1)
       weightp(i)=x(i)-x(i-1)
    end do
    !
    weightm=weightm/two
    weightp=weightp/two
    !
    !
    weight(1)=weightm(2)
    !
    do i=2, n-1
       weight(i)=weightp(i)+weightm(i+1)
    end do
    !
    weight(n)=weightp(n)
    !
    !
    !
    sum=zero
    do i=1, n
       sum=sum+weight(i)*y(i)
    end do
    !
  end subroutine integral_w1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine integral_w2(x, y, n, sum)
    !
    !---------------------integrates a function f(x)=y----------------------
    !----------------------with weight-function g(x)=x----------------------
    !
    !   input: dimension n
    !          x-grid
    !          y-values
    !
    !   output: int= integral of y dx from x(1) to x(n)
    !
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: x,y
    real(dp) :: sum
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: del
    !
    ! ... local arrays
    real(dp), dimension(n) :: weightm, weightp, weight, m
    !
    ! ... local characters
    !character(len=50) :: enter
    !
    m=zero
    del=zero
    weightm=zero
    weightp=zero
    weight=zero

    
    do i=2, n
       !
       m(i)=(x(i)+x(i-1))/two
       del=x(i)-x(i-1)
       !
       weightm(i)=.5d0*del*(m(i)-del/six)
       weightp(i)=.5d0*del*(m(i)+del/six)
       !
    end do
    !
    !
    weight(1)=weightm(2)
    !
    do i=2, n-1
       weight(i)=weightp(i)+weightm(i+1)
    end do
    !
    weight(n)=weightp(n)
    !
    !
    !
    sum=zero
    do i=1, n
       sum=sum+weight(i)*y(i)/x(i)
    end do
    
    return
    
  end subroutine integral_w2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine integral_w3(x, y, n, sum)
    !
    !---------------------integrates a function f(x)=y----------------------
    !----------------------with weight-function g(x)=sin(x)----------------------
    !
    !   input: dimension n
    !          x-grid
    !          y-values
    !
    !   output: int= integral of y dx from x(1) to x(n)
    !
    !
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: x,y
    real(dp) :: sum
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: del
    !
    ! ... local arrays
    real(dp), dimension(n) :: weightm, weightp, weight
    !
    do i=2, n
       !
       del=x(i)-x(i-1)
       !
       weightm(i)=cos(x(i-1)) - (sin(x(i))-sin(x(i-1)))/del
       weightp(i)=-cos(x(i)) + (sin(x(i))-sin(x(i-1)))/del
       !
    end do
    !
    !
    weight(1)=weightm(2)
    !
    do i=2, n-1
       weight(i)=weightp(i)+weightm(i+1)
    end do
    !
    weight(n)=weightp(n)
    !
    !
    !
    sum=zero
    do i=1, n
       sum=sum+weight(i)*y(i)/sin(x(i))
    end do

    !
    !***********************************************************************
    !-----------------------debug-open--------------------------------------
    !***********************************************************************
    !
    do i=1, n
       if(abs(sin(x(i))).lt.small_number) then
          write(*,*) "integration at point", i, x(i), "not possible"
          !      stop
       end if
    end do

    write(*,*) "weight     ", "x     ", "y"
    do i=1, n
       write(*,*) weight(i), sin(x(i)), y(i)
    end do
    !
    !***********************************************************************
    !-----------------------debug-close-------------------------------------
    !***********************************************************************
    !

    
    return
    
  end subroutine integral_w3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine integ1d_tau_ud(dxm, dxp, fim1, fi, fip1, integm, integp)
    !
    !-----------------------------------------------------------------------
    !
    !     integrates function [xim1,xi] and [xi,xip1]
    !        with funct(x) approximated by monotonic quadratic bezier function
    !
    !
    !   input: function values:        fim1, fi, fip1
    !          at x-values:            dxm, dxp
    !   output: evaluated integrals:   integm, integp
    !
    !                    fim1------------fi-----------fip1
    !                    xim1------------xi-----------xip1
    !                    [------dxm------][------dxp-----]
    !                    [----integm-----][----integp----]
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp), intent(in) :: dxm, dxp, fim1, fi, fip1
    real(dp), intent(out) :: integm, integp
    !
    ! ... local scalars
    real(dp) :: dxi, dxip1, dx, fpi, fcm, fcp, fmin, fmax
    !
    !prepare x-distances
    dx=dxm+dxp
    !
    !calculate derivative at point i
    fpi=(fi-fim1)*dxp/dxm/dx + (fip1-fi)*dxm/dxp/dx
    !
    !calculate control point in interval [xim1,xi], and in [xi,xip1]
    fcm=fi-dxm*fpi/two
    fcp=fi+dxp*fpi/two
    !
    !ensure monotonicity
    fmin=min(fim1,fi)
    fmax=max(fim1,fi)
    fcm=max(fcm,fmin)
    fcm=min(fcm,fmax)
    !
    fmin=min(fi,fip1)
    fmax=max(fi,fip1)
    fcp=max(fcp,fmin)
    fcp=min(fcp,fmax)
    !
    !
    !perform integration
    integm=dxm*(fim1+fcm+fi)/three
    integp=dxp*(fi+fcp+fip1)/three

  end subroutine integ1d_tau_ud


  end module mod_integ1d
