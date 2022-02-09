module mod_interp1d
!
use prog_type
use fund_const
!
implicit none
!
contains 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine interpol_lin(nd_out, nd_in, xval_out, xval_in, yval_out, yval_in)
    !
    !--------------linear interpolation of yval_in onto yval_out--------
    !
    !-----input:  xval_in, xval_out: x-values
    !             yval_in: y-values
    !             nd_in, nd_out: dimensions of input/output-values
    !
    !-----output: yval_out---------------------------------------------------
    !
    !
    !
    !note: works only for decreasing vector1
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd_in, nd_out
    real(dp), dimension(nd_in), intent(in) :: xval_in, yval_in
    real(dp), dimension(nd_out), intent(in) :: xval_out
    real(dp), dimension(nd_out) :: yval_out
    !
    ! ... local scalars
    integer(i4b) :: i
    integer(i4b) :: iim2, iim1, ii, iip1
    !
    ! ... local functions
    !
    !------this routine works only for increasing xval_in-------------------
    !
    do i=1, nd_in-1
       if (xval_in(i).gt.xval_in(i+1)) then
          write(*,*) 'error in interpol: non monotonic increasing x-values'
          stop
       endif
    enddo
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nd_out
       !
       !find index for which will be interpolated
       !   call find_index(xval_out(i), xval_in, nd_in, ii)
       !   im1=ii-1
       call find_index(xval_out(i), xval_in, nd_in, iim2, iim1, ii, iip1)
       !
       !interpolate
       yval_out(i)=interpol_yp(xval_in(iim1), xval_in(ii), yval_in(iim1), yval_in(ii), xval_out(i))
       !
    enddo
    !
    !
  end subroutine interpol_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine interpol_cube(nd_out, nd_in, xval_out, xval_in, yval_out, yval_in)
    !
    !--------------cubic interpolation of yval_in onto yval_out-------------
    !
    !-----input:  xval_in, xval_out: x-values
    !             yval_in: y-values
    !             nd_in, nd_out: dimensions of input/output-values
    !
    !-----output: yval_out---------------------------------------------------
    !
    !
    !
    !note: works only for decreasing vector1
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd_in, nd_out
    real(dp), dimension(nd_in), intent(in) :: xval_in, yval_in
    real(dp), dimension(nd_out), intent(in) :: xval_out
    real(dp), dimension(nd_out) :: yval_out
    !
    ! ... local scalars
    integer(i4b) :: i
    integer(i4b) :: iim2, iim1, ii, iip1
    !
    ! ... local functions
    !
    !
    !------this routine works only for increasing xval_in-------------------
    !
    do i=1, nd_in-1
       if (xval_in(i).gt.xval_in(i+1)) then
          write(*,*) 'error in interpol: non monotonic increasing x-values'
          stop
       endif
    enddo
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nd_out
       !
       !find index for which will be interpolated
       !   call find_index(xval_out(i), xval_in, nd_in, ii)
       !special treatment of boundary derivatives at the limb
       !   im2=ii-2
       !   im1=ii-1
       !   ip1=ii+1
       !   if(ii.eq.2) im2=im1
       !   if(ii.eq.nd_in) ip1=ii
       call find_index(xval_out(i), xval_in, nd_in, iim2, iim1, ii, iip1)
       !interpolation
       yval_out(i)=interpol_fyp_cube(yval_in(iim2), yval_in(iim1), yval_in(ii), yval_in(iip1), &
                                     xval_in(iim2), xval_in(iim1), xval_in(ii), xval_in(iip1), &
                                     xval_out(i))
       !
    enddo
    !
    !
  end subroutine interpol_cube
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine interpol_cube_mono(nd_out, nd_in, xval_out, xval_in, yval_out, yval_in)
    !
    !--------------monotonic cubic interpolation of yval_in onto yval_out---
    !
    !-----input:  xval_in, xval_out: x-values
    !             yval_in: y-values
    !             nd_in, nd_out: dimensions of input/output-values
    !
    !-----output: yval_out---------------------------------------------------
    !
    !
    !
    !note: works only for decreasing vector1
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd_in, nd_out
    real(dp), dimension(nd_in), intent(in) :: xval_in, yval_in
    real(dp), dimension(nd_out), intent(in) :: xval_out
    real(dp), dimension(nd_out) :: yval_out
    !
    ! ... local scalars
    integer(i4b) :: i, iim2, iim1, ii, iip1
    !
    ! ... local arrays
    real(dp), dimension(nd_in) :: a_spline, b_spline, c_spline, d_spline
    !
    ! ... local functions
    !
    !
    !------this routine works only for increasing xval_in-------------------
    !
    do i=1, nd_in-1
       if (xval_in(i).gt.xval_in(i+1)) then
          write(*,*) 'error in interpol: non monotonic increasing x-values'
          stop
       endif
    enddo
    !
    !--------------------prepare spline coefficients------------------------
    !
    call cube_mono_coeff(nd_in, xval_in, yval_in, a_spline, b_spline, c_spline, d_spline)
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nd_out
       !
       !find index for which will be interpolated
       !   call find_index(xval_out(i), xval_in, nd_in, ii)
       call find_index(xval_out(i), xval_in, nd_in, iim2, iim1, ii, iip1)
       !interpolation
       yval_out(i)=interpol_yp_spline(a_spline(ii), b_spline(ii), c_spline(ii), &
                                      d_spline(ii), xval_in(ii), xval_out(i))
       !
    enddo
    !
    !
  end subroutine interpol_cube_mono
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine interpol_spline(nd_out, nd_in, xval_out, xval_in, yval_out, yval_in)
    !
    !--------------cubic spline interpolation of yval_in onto yval_out------
    !
    !-----input:  xval_in, xval_out: x-values
    !             yval_in: y-values
    !             nd_in, nd_out: dimensions of input/output-values
    !
    !-----output: yval_out---------------------------------------------------
    !
    !
    !
    !note: works only for decreasing vector1
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd_in, nd_out
    real(dp), dimension(nd_in), intent(in) :: xval_in, yval_in
    real(dp), dimension(nd_out), intent(in) :: xval_out
    real(dp), dimension(nd_out) :: yval_out
    !
    ! ... local scalars
    integer(i4b) :: i, ii, iim2, iim1, iip1
    !
    ! ... local arrays
    real(dp), dimension(nd_in) :: a_spline, b_spline, c_spline, d_spline
    !
    ! ... local functions
    !
    !------this routine works only for increasing xval_in-------------------
    !
    do i=1, nd_in-1
       if (xval_in(i).gt.xval_in(i+1)) then
          write(*,*) 'error in interpol: non monotonic increasing x-values'
          stop
       endif
    enddo
    !
    !--------------------prepare spline coefficients------------------------
    !
    call spline_coeff(nd_in, xval_in, yval_in, a_spline, b_spline, c_spline, d_spline)
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nd_out
       !
       !find index for which will be interpolated
       !   call find_index(xval_out(i), xval_in, nd_in, ii)
       call find_index(xval_out(i), xval_in, nd_in, iim2, iim1, ii, iip1)
       !interpolation
       yval_out(i)=interpol_yp_spline(a_spline(ii), b_spline(ii), c_spline(ii), &
                                      d_spline(ii), xval_in(ii), xval_out(i))

    enddo
    !
    !
  end subroutine interpol_spline
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine spline_coeff(nd, xin, yin, acoeff, bcoeff, ccoeff, dcoeff)
    !
    !       calculates coefficient for cubic spline interpolation:
    !
    !   s(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), dimension(nd), intent(in) :: xin, yin
    real(dp), dimension(nd) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: hi, hii
    !
    ! ... local arrays
    real(dp), dimension(nd) :: diag, sub_diag, sup_diag, spp_vec, d_vec
    real(dp), dimension(nd, nd) :: inverse
    !
    !------------------set up tri-diagonal matrix---------------------------
    !--------which determines second derivatives at the nodes---------------
    !
    diag(1) = one
    sub_diag(1) = zero
    sup_diag(1) = zero
    d_vec(1) = zero
    !
    do i=2, nd-1
       hi = xin(i) - xin(i-1)
       hii= xin(i+1) - xin(i)
       diag(i) = two*(hi+hii)
       sub_diag(i) = hi
       sup_diag(i) = hii
       d_vec(i) = six*(yin(i+1) - yin(i))/hii - six*(yin(i)-yin(i-1))/hi
    enddo
    !
    diag(nd) = one
    sub_diag(nd) = zero
    sup_diag(nd) = zero
    d_vec(nd) = zero
    !
    !do i=1, nd
    !   write(*,'(i5, 6e20.8)') i, xin(i), yin(i), sub_diag(i), diag(i), sup_diag(i), d_vec(i)
    !enddo
    !
    call invtri_lev(nd, sub_diag, diag, sup_diag, spp_vec, d_vec)
    !
    !------------------------calculate coefficients-------------------------
    !
    acoeff(1)=zero
    bcoeff(1)=zero
    ccoeff(1)=zero
    dcoeff(1)=zero
    !
    do i=2, nd
       hi=xin(i)-xin(i-1)
       acoeff(i) = (spp_vec(i)-spp_vec(i-1)) / six / hi
       bcoeff(i) = spp_vec(i) / two
       ccoeff(i) = (yin(i)-yin(i-1))/hi + (two*spp_vec(i) + spp_vec(i-1)) * hi / six
       dcoeff(i) = yin(i)
    enddo
    
  end subroutine spline_coeff
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine cube_mono_coeff(nd, xin, yin, acoeff, bcoeff, ccoeff, dcoeff)
    !
    !       calculates coefficient for piecewise cubic interpolation:
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are taken from the ansatz by steffen 1990
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), dimension(nd), intent(in) :: xin, yin
    real(dp), dimension(nd) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: hi, si, pi, him, sim, hip, sip
    real(dp), dimension(nd) :: y_prime
    !
    ! ... local arrays
    !
    !---------------------calculate left boundary---------------------------
    !
    hi=xin(2)-xin(1)
    hip=xin(3)-xin(2)
    si=(yin(2)-yin(1))/hi
    sip=(yin(3)-yin(2))/hip
    !
    pi=si*(one+hi/(hi+hip)) - sip*hi/(hip+hi)
    !
    y_prime(1)=(sign(one,pi)+sign(one,sip))*min(abs(sip), half*abs(pi))
    !
    !a, b, c, d are not needed at index 1, since interpolation begins
    !   earliest between i=2, i=1
    acoeff(1)=zero
    bcoeff(1)=zero
    ccoeff(1)=zero
    dcoeff(1)=zero
    !
    !----------------------calculate inner points---------------------------
    !
    do i=2, nd-1
       hi=xin(i)-xin(i-1)
       hip=xin(i+1)-xin(i)
       si=(yin(i)-yin(i-1))/hi
       sip=(yin(i+1)-yin(i))/hip
       pi=(si*hip + sip*hi)/(hi+hip)
       !
       y_prime(i)=(sign(one,si)+sign(one,sip))*min(abs(si),abs(sip),half*abs(pi))
       !
       acoeff(i) = -two*si/hi/hi + (y_prime(i-1)+y_prime(i))/hi/hi
       bcoeff(i) = -three*si/hi    + (y_prime(i-1)+two*y_prime(i))/hi
       ccoeff(i) = y_prime(i)
       dcoeff(i) = yin(i)
    enddo
    !
    !---------------------calculate right boundary--------------------------
    !
    hi=xin(nd)-xin(nd-1)
    him=xin(nd-1)-xin(nd-2)
    si=(yin(nd)-yin(nd-1))/hi
    sim=(yin(nd-1)-yin(nd-2))/him
    !
    pi=si*(one+hi/(hi+him)) - sim*hi/(hi+him)
    !
    y_prime(nd)=(sign(one,pi)+sign(one,sim))*min(abs(sim), half*abs(pi))
    !
    acoeff(nd) = -two*si/hi/hi + (y_prime(nd-1)+y_prime(nd))/hi/hi
    bcoeff(nd) = -three*si/hi    + (y_prime(nd-1)+two*y_prime(nd))/hi
    ccoeff(nd) = y_prime(nd)
    dcoeff(nd) = yin(nd)
    !
    !
    !
  end subroutine cube_mono_coeff
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine lin_coeff(nd, xin, yin, acoeff, bcoeff)
    !
    !       calculates coefficient for piecewise linear interpolation:
    !
    !   f(x) = acoeff * (x-xi) + bcoeff
    !
    !   first derivatives are taken from the ansatz by steffen 1990
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), dimension(nd), intent(in) :: xin, yin
    real(dp), dimension(nd) :: acoeff, bcoeff
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: hi, si
    !
    !a, b are not needed at index 1, since interpolation begins
    !   earliest between [x_2, x_1]
    acoeff(1)=zero
    bcoeff(1)=zero
    !
    do i=2, nd
       hi=xin(i)-xin(i-1)
       si=(yin(i)-yin(i-1))/hi
       acoeff(i) = si
       bcoeff(i) = yin(i)
    enddo
    !
    !
    !
  end subroutine lin_coeff
!
!***********************************************************************
!                            FUNCTIONS
!***********************************************************************
!
  function interpol_yp(xminus, xplus, yminus, yplus, xp)
    !
    !-----------------------------------------------------------------------
    !
    !     interpolates function values yminus, yplus
    !                               at xminus, xplus onto point xp
    !
    !   input: function values:                       yminus, yplus
    !          at x-values:                           xminus, xplus
    !          to be interpolated onto:               xp
    !   output: interpolated function value at xp:    interpol_yp
    !
    !   ansatz: f(x)=grad*x + yintercept
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp) :: interpol_yp
    real(dp), intent(in) :: xminus, xplus, yminus, yplus, xp
    !
    interpol_yp = yplus + (yplus-yminus)/(xplus-xminus) * (xp-xplus)
    !
  end function interpol_yp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine coeff_yp(f_im1, f_i, x_im1, x_i, x_p, acoeff, bcoeff)
    !         calculates coefficients for linear interpolation, such that
    !                 f(x) = acoeff*f_im1 + bcoeff*f_i
    !
    !   Ansatz: f(x) = a*(x-xi) + b
    !           (a, b  different from acoeff, bcoeff)
    !
    !on input:
    !                    f_im1----.------f_i
    !                    x_im1---x_p-----x_i
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: acoeff, bcoeff
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im1, f_i, x_im1, x_i, x_p
    real(dp), intent(out) :: acoeff, bcoeff
    !
    ! ... local scalars
    !
    ! ... local arrays
    !
    acoeff = -(x_p-x_i)/(x_i-x_im1)
    bcoeff = one-acoeff
    !
  end subroutine coeff_yp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_ypl(xim1, xi, fim1, fi, xp)
    !
    !-----------------------------------------------------------------------
    !
    !     interpolates function values fim1, fi
    !                               at xim1, xi onto point xp
    !         xim1------xp------xi
    !         fim1------.-------fi
    !
    !   input: function values:                       fim1, fi
    !          at x-values:                           xim1, xi
    !          to be interpolated onto:               xp
    !   output: interpolated function value at xp:    interpol_ypl
    !
    !   ansatz: log(f(x))=grad*log(x) + yintercept
    !               f(x)=fim1*(xp/xim1)**grad
    !
    !-----------------------------------------------------------------------
    !
    !
    ! ... arguments
    real(dp) :: interpol_ypl
    real(dp), intent(in) :: xim1, xi, fim1, fi, xp
    !
    ! ... local scalars
    real(dp) :: grad
    !
    ! ... local logicals
    logical :: llogf, llogx
    !
    llogx=.true.
    if(xim1.eq.zero) then
       llogx=.false.
    elseif(xp/xim1.le.zero) then
       llogx=.false.
    elseif(xi/xim1.le.zero) then
       llogx=.false.
    endif
    !
    llogf=.true.
    if(fim1.eq.zero) then
       llogf=.false.
    elseif(fi/fim1.le.zero) then
       llogf=.false.
    endif
    !
    if(llogx.and.llogf) then
       grad=log10(fi/fim1)/log10(xi/xim1)
       interpol_ypl = fim1*(xp/xim1)**grad
    else
       interpol_ypl = fim1 + (fi-fim1)/(xi-xim1) * (xp-xim1)
    endif
    !
    !
    !
  end function interpol_ypl
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_yp2(a, b, xi, xin)
    !
    !-----------------------------------------------------------------------
    !
    !     interpolates function values linearly with given (pre-calculated)
    !        a, b:
    !
    !      f(x) = a*(x-xi) + b
    !
    !      with x in [x(i-1), x(i)]
    !
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    real(dp) :: interpol_yp2
    real(dp), intent(in) :: a, b, xi, xin
    !
    interpol_yp2 = a*(xin-xi) + b
    !
  end function interpol_yp2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_2p_quad(f_im1, f_i, x_im1, x_i, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using quadratic bezier spline
    !
    !   f(t) = f_im1 * (1-t)^3 + f_i * t^2 + 2*f_c * t*(1-t)
    !      t = (x-x_i)/(x_i-x_im1)
    !
    !   control point is chosen from f_im1, f_i with weights w_lower=q*f_larger
    !
    !on input: 
    !
    !                    f_im1----.------f_i
    !                    x_im1---x_p-----x_i
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value at x_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_2p_quad
    real(dp), intent(in) :: f_im1, f_i, x_im1, x_i, x_p
    !
    ! ... local scalars
    real(dp), parameter :: fac = three
    real(dp) :: a, b, q, t
    !
    t = (x_p-x_im1)/(x_i-x_im1)
    !
    !if(f_im1.lt.f_i) then
    !   a = one-t*(one+t)/two
    !   b = t*(one+t)/two
    !elseif(f_i.lt.f_im1) then
    !   a = one-t*(three-t)two.
    !   b = t*(three-t)/two
    !else
    !   a = one-t
    !   b = t
    !endif
    

    if(f_im1.lt.f_i) then
       q=fac
    elseif(f_i.lt.f_im1) then
       q = one/fac
    else
       q = one
    endif
    !
    b=t*(two/(one+q)-t*(one-q)/(one+q))
    a=one-b
    !
    interpol_2p_quad = a*f_im1 + b*f_i
    !
  end function interpol_2p_quad
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_quad(f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using quadratic bezier spline
    !
    !   f(t) = f_im1 * (1-t)^3 + f_i * t^2 + 2*f_c * t*(1-t)
    !      t = (x-x_i)/(x_i-x_im1)
    !
    !   control point is calculated at the center of x(i) and x(i-1)
    !   first derivative is calculated from weighted mean
    !
    !on input:
    !
    !                    f_im1----.------f_i--------f_ip1
    !                    x_im1---x_p-----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_typ_quad
    real(dp), intent(in) :: f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: fp_f, fp_b, f_c
    real(dp) :: w, t
    !
    ! ... local arrays
    !
    !weight
    w = (x_ip1-x_i)/(x_ip1-x_im1)
    !
    !forward and backward differences, then weighted
    fp_f = (f_ip1-f_i)/(x_ip1-x_i)
    fp_b = (f_i-f_im1)/(x_i-x_im1)
    f_c = f_i - half*(x_i-x_im1)*(w*fp_b + (one-w)*fp_f)
    !
    !interpolation
    t = (x_p-x_im1)/(x_i-x_im1)
    interpol_typ_quad = f_im1*(one-t)**2 + f_i*t**2 + two*f_c*t*(one-t)
    !
    !
  end function interpol_typ_quad
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_quad2(f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using monotonic quadratic bezier spline
    !
    !   f(t) = f_im1 * (1-t)^2 + f_i * t^2 + 2*f_c * t*(1-t)
    !      t = (x-x_im1)/(x_i-x_im1)
    !    f_c = f_i - dfdx_i*(x_i-x_im1)/two
    !
    !   control point is calculated at the center of x(i) and x(i-1)
    !   first derivative is calculated from weighted mean (see Hayek 2010)
    !      (of backward and forward differences)
    !   three point interpolation
    !
    !on input: 
    !
    !                    f_im1----.-----f_i--------f_ip1
    !                    x_im1---x_p----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_typ_quad2
    real(dp), intent(in) :: f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: fp_f, fp_b, fp_i, f_c, f_min, f_max
    real(dp) :: w, t
    !
    ! ... local arrays
    !
    !derivative at point x_i from weighted mean
    fp_i = (f_i-f_im1)*(x_ip1-x_i)/(x_i-x_im1)/(x_ip1-x_im1) + (f_ip1-f_i)*(x_i-x_im1)/(x_ip1-x_i)/(x_ip1-x_im1)
    f_c = f_i - half*(x_i-x_im1)*fp_i
    f_min=min(f_im1, f_i)
    f_max=max(f_im1, f_i)
    if(f_c.gt.f_max) then
       f_c=f_max
    else if (f_c.lt.f_min) then
       f_c=f_min
    endif
    !
    !interpolation
    t = (x_p-x_im1)/(x_i-x_im1)
    interpol_typ_quad2 = f_im1*(one-t)**2 + f_i*t**2 + two*f_c*t*(one-t)
    !
    !
    !
  end function interpol_typ_quad2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_quad2b(f_im2, f_im1, f_i, x_im2, x_im1, x_i, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using monotonic quadratic bezier spline
    !
    !   f(t) = f_im1 * (1-t)^2 + f_i * t^2 + 2*f_c * t*(1-t)
    !      t = (x-x_im1)/(x_i-x_im1)
    !    f_c = f_im1 + dfdx_im1*(x_i-x_im1)/two
    !
    !   control point is calculated at the center of x(i) and x(i-1)
    !   first derivative is calculated from weighted mean (see Hayek 2010)
    !      (of both backward differences)
    !   three point interpolation
    !
    !on input: 
    !
    !                    f_im2----------f_im1----.-----f_i
    !                    x_im2----------x_im1---x_p----x_i
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_typ_quad2b
    real(dp), intent(in) :: f_im2, f_im1, f_i, &
                            x_im2, x_im1, x_i, x_p
    !
    ! ... local scalars
    real(dp) :: fp_f, fp_b, fp_im1, f_c, f_min, f_max
    real(dp) :: w, t
    !
    !real(dp) :: a, b, c, dx_i, dx_im1
    !
    ! ... local arrays
    !
    !
    !derivative at point x_im1 from weighted mean
    fp_im1 = (f_im1-f_im2)*(x_i-x_im1)/(x_im1-x_im2)/(x_i-x_im2) + (f_i-f_im1)*(x_im1-x_im2)/(x_i-x_im1)/(x_i-x_im2)
    f_c = f_im1 + half*(x_i-x_im1)*fp_im1
    f_min=min(f_im1, f_i)
    f_max=max(f_im1, f_i)
    if(f_c.gt.f_max) then
       f_c=f_max
    else if (f_c.lt.f_min) then
       f_c=f_min
    endif
    !
    !interpolation
    t = (x_p-x_im1)/(x_i-x_im1)
    interpol_typ_quad2b = f_im1*(one-t)**2 + f_i*t**2 + two*f_c*t*(one-t)
    !
    !
    !
  end function interpol_typ_quad2b
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_quad3(f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using quadratic function
    !
    !   f(x) = a*(x-xi)^2 + b*(x-xi) + c*(x-xi)
    !   ensure monotonicity in interval [x_im1,x_i]
    !
    !on input: 
    !
    !                    f_im1----.------f_i---------f_ip1
    !                    x_im1---x_p-----x_i---------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_typ_quad3
    real(dp), intent(in) :: f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: delx_i, delx_ip1, delf_i, delf_ip1, dfdx_im1, dfdx_i, dfdx_ip1, delx
    real(dp) :: a, b, c
    real(dp) :: adum, bdum, cdum, ddum
    !
    ! ... local arrays
    !
    delx_i=x_i-x_im1
    delx_ip1=x_ip1-x_i
    delf_i=f_i-f_im1
    delf_ip1=f_ip1-f_i
    !
    !calculate coefficients (numerically not the fastests method, but for the moment...)
    a=(delf_ip1/delx_ip1 - delf_i/delx_i)/(delx_ip1+delx_i)
    b=(delf_ip1*delx_i/delx_ip1 + delf_i*delx_ip1/delx_i)/(delx_ip1+delx_i)
    c=f_i
    !ensure monotonic interpolation (take linear approximation if non-monotonic)
    dfdx_im1=-two*a*delx_i + b
    dfdx_i=b
    if(dfdx_im1*dfdx_i.lt.zero) then
       b=delf_i/delx_i
       interpol_typ_quad3 = b*(x_p-x_i) + c
    else
       interpol_typ_quad3 = a*(x_p-x_i)**2 + b*(x_p-x_i) + c
    endif
    !
    !
  end function interpol_typ_quad3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_quad3b(f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using quadratic function
    !
    !   f(x) = a*(x-xi)^2 + b*(x-xi) + c*(x-xi)
    !   ensure monotonicity in interval [x_i,x_ip1]
    !
    !on input: 
    !
    !                    f_im1--------f_i-----.------f_ip1
    !                    x_im1--------x_i----x_p-----x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_typ_quad3b
    real(dp), intent(in) :: f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: delx_i, delx_ip1, delf_i, delf_ip1, dfdx_im1, dfdx_i, dfdx_ip1, delx
    real(dp) :: a, b, c
    real(dp) :: adum, bdum, cdum, ddum
    !
    ! ... local arrays
    !
    delx_i=x_i-x_im1
    delx_ip1=x_ip1-x_i
    !
    delf_i=f_i-f_im1
    delf_ip1=f_ip1-f_i
    !
    !calculate coefficients (numerically not the fastests method, but for the moment...)
    a=(delf_ip1/delx_ip1 - delf_i/delx_i)/(delx_ip1+delx_i)
    b=(delf_ip1*delx_i/delx_ip1 + delf_i*delx_ip1/delx_i)/(delx_ip1+delx_i)
    c=f_i
    !
    !ensure monotonic interpolation (take linear approximation if non-monotonic)
    !dfdx_im1=-two*a*delx_i + b
    dfdx_i=b
    dfdx_ip1=two*a*delx_ip1 + b
    if(dfdx_i*dfdx_ip1.lt.zero) then
       b=delf_ip1/delx_ip1
       interpol_typ_quad3b = b*(x_p-x_i) + c
    else
       interpol_typ_quad3b = a*(x_p-x_i)**2 + b*(x_p-x_i) + c
    endif
    !
    !
    !
  end function interpol_typ_quad3b
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine coeff_typ_quad(f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p, &
                            acoeff, bcoeff, ccoeff)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using standard quadratic interpolation
    !
    !on input:
    !
    !                    f_im1----.-----f_i--------f_ip1
    !                    x_im1---x_p----x_i--x_p---x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: coefficients acoeff, bcoeff, ccoeff
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: dx_i, dx_ip1, dx, dx_p
    !
    !
    !
    dx_i = x_i-x_im1
    dx_ip1 = x_ip1-x_i
    dx = dx_i+dx_ip1
    dx_p = x_p-x_i
    !
    acoeff = dx_p*(dx_p-dx_ip1)/dx_i/dx
    bcoeff = (dx_p+dx_i)*(dx_ip1-dx_p)/dx_i/dx_ip1
    ccoeff = dx_p*(dx_p+dx_i)/dx_ip1/dx
    !
    !
    !
  end subroutine coeff_typ_quad
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine coeff_typ_quad2(f_im1, f_i, f_ip1, &
                             x_im1, x_i, x_ip1, x_p, &
                             acoeff, bcoeff, ccoeff)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !                 using monotonic quadratic bezier spline
    !
    !   f(t) = f_im1 * (1-t)^2 + f_i * t^2 + 2*f_c * t*(1-t)
    !      t = (x-x_im1)/(x_i-x_im1)
    !    f_c = f_i - dfdx_i*(x_i-x_im1)/two
    !   such that f(x) = acoeff*f_im1 + bcoeff*f_i + ccoeff*f_ip1
    !
    !   control point is calculated at the center of x(i) and x(i-1)
    !   first derivative is calculated from weighted mean (see Hayek 2010)
    !      (of backward and forward differences)
    !   three point interpolation
    !
    !on input:
    !
    !                    f_im1----.-----f_i--------f_ip1
    !                    x_im1---x_p----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: coefficients acoeff, bcoeff, ccoeff
    !
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: dx_i, dx_ip1, dx, dx_p
    real(dp) :: p, q, r, t, f_c
    !
    !
    !
    dx_i = x_i-x_im1
    dx_ip1 = x_ip1-x_i
    dx = dx_i+dx_ip1
    dx_p = x_p-x_im1
    !
    if(f_i.eq.f_im1) then
       !linear interpolation
       bcoeff=dx_p/dx_i
       acoeff=one-bcoeff
       ccoeff=zero
    else 
       !
       p = dx_ip1/two/dx
       q = one+(dx_i**2-dx_ip1**2)/two/dx_ip1/dx
       r = -dx_i**2/two/dx_ip1/dx
       f_c = p*f_im1 + q*f_i + r*f_ip1
       !
       t=dx_p/dx_i
       if((f_c-f_im1)/(f_im1-f_i).ge.zero) then
          !set fc to fim1 to ensure monotonicity
          bcoeff=t**2
          acoeff=one-bcoeff
          ccoeff=zero
       elseif((f_c-f_i)/(f_i-f_im1).ge.zero) then
          !set fc to fi to ensure monotonicity
          acoeff=(one-t)**2
          bcoeff=two*t-t**2
          ccoeff=zero
       else
          acoeff = (one-t)**2+two*t*(one-t)*p
          bcoeff = t**2 + two*t*(one-t)*q
          ccoeff = two*t*(one-t)*r
       endif
    endif
    !
    !
    !
  end subroutine coeff_typ_quad2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine coeff_typ_quad3(f_im1, f_i, f_ip1, &
                             x_im1, x_i, x_ip1, x_p, &
                             acoeff, bcoeff, ccoeff)
    !
    !         calculates coefficients for quadratic interpolation, such that
    !                 f(x) = acoeff*f_im1 + bcoeff*f_i + ccoeff*f_ip1
    !
    !   Ansatz: f(x) = a*(x-xi)^2 + b*(x-xi) + c*(x-xi)
    !           (a, b, c different from acoeff, bcoeff, ccoeff)
    !   ensure monotonicity in interval [x_im1,x_i]
    !
    !on input:
    !                    f_im1----.------f_i--------f_ip1
    !                    x_im1---x_p-----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: acoeff, bcoeff, ccoeff
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: delx_i, delx_ip1, delx, delx_tot, dfdx_im1, dfdx_i
    real(dp) :: a, b, c
    real(dp) :: adum, bdum, cdum, ddum
    !
    ! ... local arrays
    !
    delx_i=x_i-x_im1
    delx_ip1=x_ip1-x_i
    delx=x_p-x_i
    delx_tot=x_ip1-x_im1
    !
    !derivative at point x_im1
    dfdx_im1 = -f_im1/delx_i - (delx_ip1*f_im1+delx_i*f_ip1)/delx_ip1/delx_tot + f_i*delx_tot/delx_i/delx_ip1
    !derivative at point x_i
    dfdx_i = ((f_i-f_im1)*delx_ip1**2 + (f_ip1-f_i)*delx_i**2)/delx_i/delx_ip1/delx_tot
    !
    if(dfdx_im1*dfdx_i.lt.0.) then
       !linear interpolation
       acoeff=-delx/delx_i
       bcoeff=one+delx/delx_i
       ccoeff=zero
    else
       acoeff=(delx**2 - delx*delx_ip1)/delx_i/delx_tot
       !   bcoeff=-delx**2/delx_i/delx_ip1 + delx*(delx_ip1**2-delx_i**2)/delx_i/delx_ip1/delx_tot + one
       !new version
       bcoeff=-delx**2/delx_i/delx_ip1 + delx*(delx_ip1-delx_i)/delx_i/delx_ip1 + one
       ccoeff=(delx**2 + delx*delx_i)/delx_ip1/delx_tot
    endif
    !
    !
  end subroutine coeff_typ_quad3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine coeff_typ_quad3b(f_im1, f_i, f_ip1, &
                              x_im1, x_i, x_ip1, x_p, &
                              acoeff, bcoeff, ccoeff)
    !
    !         calculates coefficients for quadratic interpolation, such that
    !                 f(x) = acoeff*f_im1 + bcoeff*f_i + ccoeff*f_ip1
    !
    !   Ansatz: f(x) = a*(x-xi)^2 + b*(x-xi) + c*(x-xi)
    !           (a, b, c different from acoeff, bcoeff, ccoeff)
    !   ensure monotonicity in interval [x_im1,x_i]
    !
    !on input:
    !                    f_im1-----------f_i----.-----f_ip1
    !                    x_im1-----------x_i---x_p----x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: acoeff, bcoeff, ccoeff
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im1, f_i, f_ip1, x_im1, x_i, x_ip1, x_p
    real(dp), intent(out) :: acoeff, bcoeff, ccoeff
    !
    ! ... local scalars
    real(dp) :: delx_i, delx_ip1, delx, delx_tot, dfdx_ip1, dfdx_i
    real(dp) :: a, b, c
    real(dp) :: adum, bdum, cdum, ddum
    !
    ! ... local arrays
    !
    delx_i=x_i-x_im1
    delx_ip1=x_ip1-x_i
    delx=x_p-x_i
    delx_tot=x_ip1-x_im1
    !
    !derivative at point x_i
    dfdx_i = -f_im1*delx_ip1/delx_i/delx_tot + f_i*(delx_ip1-delx_i)/delx_i/delx_ip1 + &
              f_ip1*delx_i/delx_ip1/delx_tot
    !derivative at point x_ip1
    dfdx_ip1 = f_im1*delx_ip1/delx_i/delx_tot - f_i*(delx_ip1+delx_i)/delx_i/delx_ip1 + &
               f_ip1*(two*delx_ip1+delx_i)/delx_ip1/delx_tot
    !
    if(dfdx_ip1*dfdx_i.lt.0.) then
       !linear interpolation
       acoeff=zero
       bcoeff=one-delx/delx_ip1
       ccoeff=delx/delx_ip1
    else
       acoeff=(delx**2 - delx*delx_ip1)/delx_i/delx_tot
       bcoeff=-delx**2/delx_i/delx_ip1 + delx*(delx_ip1-delx_i)/delx_i/delx_ip1 + one
       ccoeff=(delx**2 + delx*delx_i)/delx_ip1/delx_tot
    endif
    !
    !
  end subroutine coeff_typ_quad3b
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_yp_spline(a_spline, b_spline, c_spline, d_spline, xi, xin)
    !
    !   interpolates a function with given cubic spline coefficients
    !  a_spline, b_spline, c_spline and given point xi onto point xin
    !
    ! ... arguments
    real(dp) :: interpol_yp_spline
    real(dp), intent(in) :: a_spline, b_spline, c_spline, d_spline, xi, xin
    !
    !
    interpol_yp_spline = a_spline * (xin-xi)**3 + &
                         b_spline * (xin-xi)**2 + &
                         c_spline * (xin-xi) + &
                         d_spline
  end function interpol_yp_spline
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_fyp_cube(f_im2, f_im1, f_i, f_ip1, &
                             x_im2, x_im1, x_i, x_ip1, &
                             x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are calculated locally as secant:
    !                  dfdx_i = (f_ip1-fim1)/(x_ip1-x_im1)
    !                dfdx_im1 = (f_i-fim2)/(x_i-x_im2)
    !      => central differences: catmul-rom-spline
    !      => complete interpolation only accurate to second order
    !   first derivatives are calculated by forward / backward differences
    !      for the correct find_index procedure on boundaries
    !   four point interpolation
    !
    !on input:
    !
    !        f_im2-------f_im1----.-----f_i--------f_ip1
    !        x_im2-------x_im1---x_p----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    ! ... arguments
    real(dp) :: interpol_fyp_cube
    real(dp), intent(in) :: f_im2, f_im1, f_i, f_ip1, &
                            x_im2, x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: h_i, h_ip1, fp_i, fp_im1
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local arrays
    !
    !
    h_i=x_i-x_im1
    !
    !first derivative approximated by central differences
    fp_i = (f_ip1-f_im1)/(x_ip1 - x_im1)
    fp_im1 = (f_i - f_im2)/(x_i - x_im2)
    !
    acoeff = -two*(f_i-f_im1)/h_i**3 + (fp_i + fp_im1)/h_i**2
    bcoeff = -three*(f_i-f_im1)/h_i**2 + (fp_i*two + fp_im1)/h_i
    ccoeff = fp_i
    dcoeff=f_i
    !
    h_i=x_p-x_i
    interpol_fyp_cube = acoeff*h_i**3 + bcoeff*h_i**2 + ccoeff*h_i + dcoeff
    !
    !
    !
  end function interpol_fyp_cube
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_fyp_cubem(f_im2, f_im1, f_i, f_ip1, &
                              x_im2, x_im1, x_i, x_ip1, &
                              x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !               using monotonicity from Fritsch/Carlson 1980
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are calculated locally as secant:
    !                  dfdx_i = (f_ip1-fim1)/(x_ip1-x_im1)
    !                dfdx_im1 = (f_i-fim2)/(x_i-x_im2)
    !      => central differences: catmul-rom-spline
    !      => complete interpolation only accurate to second order
    !   four point interpolation
    !
    !on input: 
    !
    !        f_im2-------f_im1----.-----f_i--------f_ip1
    !        x_im2-------x_im1---x_p----x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp) :: interpol_fyp_cubem
    real(dp), intent(in) :: f_im2, f_im1, f_i, f_ip1, &
                            x_im2, x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: h_i, h_ip1, fp_i, fp_im1, s_i, rad, alpha, beta
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local arrays
    !
    ! ... local functions
    !
    !step size
    h_i=x_i-x_im1
    !
    !secant
    s_i=(f_i-f_im1)/h_i
    !
    if(s_i.eq.zero) then
       fp_i=zero
       fp_im1=zero
    else
       !first derivative approximated by central differences
       fp_i = (f_ip1-f_im1)/(x_ip1 - x_im1)
       fp_im1 = (f_i - f_im2)/(x_i - x_im2)
       if(sign(one,fp_i).ne.sign(one,fp_im1).and.sign(one,fp_i).ne.sign(one,s_i)) then
          !if different signs, no monotonic interpolation possible, use quadratic bezier approach
          interpol_fyp_cubem=interpol_typ_quad2b(f_im2,f_im1,f_i,x_im2,x_im1,x_i,x_p)
          return
       endif
       alpha = fp_i/s_i 
       beta = fp_im1/s_i
       rad = sqrt(alpha**2 + beta**2)
       if(rad.gt.three) then
          !manipulate derivatives such that monotonic interpolation is obtained
          fp_i=three*fp_i/rad
          fp_im1=three*fp_im1/rad
       endif
    endif
    !
    acoeff = -two*(f_i-f_im1)/h_i**3 + (fp_i + fp_im1)/h_i**2
    bcoeff = -three*(f_i-f_im1)/h_i**2 + (fp_i*two + fp_im1)/h_i
    ccoeff = fp_i
    dcoeff=f_i
    !
    h_i=x_p-x_i
    interpol_fyp_cubem = acoeff*h_i**3 + bcoeff*h_i**2 + ccoeff*h_i + dcoeff
    !
    !
    !
  end function interpol_fyp_cubem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_typ_cube(f_im1, f_i, f_ip1, &
                             x_im1, x_i, x_ip1, x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are
    !      central differences at point x(i)
    !      forward differences at point x(i-1)
    !   three point interpolation
    !
    !   note: point (x_ip1, f_ip1) not necessarily matched
    !
    !on input:
    !
    !        f_im1--------f_i--------f_ip1
    !        x_im1--------x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    ! ... arguments
    real(dp) :: interpol_typ_cube
    real(dp), intent(in) :: f_im1, f_i, f_ip1, &
                            x_im1, x_i, x_ip1, x_p
    !
    ! ... local scalars
    real(dp) :: h_i, fp_i, fp_im1
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local arrays
    !
    h_i=x_i-x_im1
    !
    !first derivative approximated by central differences
    fp_i = (f_ip1-f_im1)/(x_ip1 - x_im1)
    !first derivative approximated by forward differences
    fp_im1 = (f_i - f_im1)/(x_i - x_im1)
    !
    acoeff = -two*(f_i-f_im1)/h_i**3 + (fp_i + fp_im1)/h_i**2
    bcoeff = -three*(f_i-f_im1)/h_i**2    + (fp_i*two + fp_im1)/h_i
    ccoeff = fp_i
    dcoeff=f_i
    !
    h_i=x_p-x_i
    interpol_typ_cube = acoeff*h_i**3 + bcoeff*h_i**2 + ccoeff*h_i + dcoeff
    !
    !
  end function interpol_typ_cube
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_fyp_cube4(f_im2, f_im1, f_i, f_ip1, &
                              x_im2, x_im1, x_i, x_ip1, &
                              x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are calculated from weighted mean
    !   first derivatives are calculated by forward / backward differences
    !      for the correct find_index procedure on boundaries
    !   four point interpolation
    !
    !on input:
    !
    !        f_im2-------f_im1--------f_i--------f_ip1
    !        x_im2-------x_im1--------x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im2, f_im1, f_i, f_ip1, &
                            x_im2, x_im1, x_i, x_ip1, x_p
    real(dp) :: interpol_fyp_cube4
    !
    ! ... local scalars
    real(dp) :: h_i, fp_i, fp_im1, fp_f, fp_b, w
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local arrays
    !
    h_i=x_i-x_im1
    !
    !first derivative approximated by arithmetic mean
    fp_f = (f_ip1-f_i)/(x_ip1 - x_i)
    fp_b = (f_i-f_im1)/(x_i - x_im1)
    w = (x_ip1-x_i)/(x_ip1-x_im1)
    fp_i = w*fp_b + (one-w)*fp_f
    !
    fp_f = (f_i-f_im1)/(x_i - x_im1)
    fp_b = (f_im1-f_im2)/(x_im1 - x_im2)
    w = (x_i-x_im1)/(x_i-x_im2)
    fp_im1 = w*fp_b + (one-w)*fp_f
    !
    acoeff = -two*(f_i-f_im1)/h_i**3 + (fp_i + fp_im1)/h_i**2
    bcoeff = -three*(f_i-f_im1)/h_i**2    + (fp_i*two + fp_im1)/h_i
    ccoeff = fp_i
    dcoeff=f_i
    !
    h_i=x_p-x_i
    interpol_fyp_cube4 = acoeff*h_i**3 + bcoeff*h_i**2 + ccoeff*h_i + dcoeff
    !
    !
    !
  end function interpol_fyp_cube4
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_fyp_cube2(f_im2, f_im1, f_i, f_ip1, &
                              x_im2, x_im1, x_i, x_ip1, &
                              x_p)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are calculated by harmonic mean if they don't change direction
    !                     else, are set to zero (see Ibgui 2013)
    !      => monotonic interpolation
    !   works only if ghost zones are properly defined
    !   four point interpolation
    !
    !on input: 
    !
    !        f_im2-------f_im1--------f_i--------f_ip1
    !        x_im2-------x_im1--------x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    real(dp), intent(in) :: f_im2, f_im1, f_i, f_ip1, &
                            x_im2, x_im1, x_i, x_ip1, x_p
    real(dp) :: interpol_fyp_cube2
    !
    ! ... local scalars
    real(dp) :: h_i, fp_i, fp_im1
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, q
    real(dp) :: t1, t2, t3, alpha_im1, alpha_i, dum
    !
    ! ... local arrays
    !
    h_i=x_i-x_im1
    !
    t1 = (f_im1-f_im2)/(x_im1-x_im2)
    t2 = (f_i-f_im1)/(x_i-x_im1)
    t3 = (f_ip1-f_i)/(x_ip1-x_i)
    alpha_im1=(one+(x_i-x_im1)/(x_i-x_im2))/three
    alpha_i=(one+(x_ip1-x_i)/(x_ip1-x_im1))/three
    if(t2*t3.ge.zero) then
       dum=t2*t3
    else
       dum=zero
    endif
    if(abs(dum).lt.small_number) then
       fp_i = zero
    else
       fp_i = dum/((one-alpha_i)*t2 + alpha_i*t3)
    endif
    !
    if(t1*t2.ge.zero) then
       dum=t1*t2
    else
       dum=zero
    endif
    if(abs(dum).lt.small_number) then
       fp_im1 = zero
    else
       fp_im1 = dum/((one-alpha_im1)*t1 + alpha_im1*t2)
    endif
    !
    q=(x_p-x_i)/h_i
    acoeff = one - two*q**3 - three*q**2
    bcoeff = two*q**3 + three*q**2
    ccoeff = (q**3 + two*q**2 + q) * h_i
    dcoeff = (q**3 + q**2) * h_i
    !
    interpol_fyp_cube2 = acoeff*f_i + bcoeff*f_im1 + ccoeff*fp_i + dcoeff*fp_im1
    !
    !
    !
  end function interpol_fyp_cube2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function interpol_fyp_cube3(f_im2, f_im1, f_i, f_ip1, &
                              x_im2, x_im1, x_i, x_ip1, &
                              x_p, ii, nd)
    !
    !         interpolates values given on a 1d-axis onto point x_p
    !
    !   f(x) = acoeff * (x-xi)^3 + bcoeff * (x-xi)^2 + ccoeff * (x-xi) + dcoeff
    !
    !   first derivatives are calculated by harmonic mean if they don't change direction
    !                     else, are set to zero
    !      => monotonic interpolation according to Steffen 1990
    !   four point interpolation
    !
    !on input:
    !
    !        f_im2-------f_im1--------f_i--------f_ip1
    !        x_im2-------x_im1--------x_i--------x_ip1
    !
    !        x_p: coordinate of point onto which shall be interpolated
    !
    !on output: interpolated value f_p
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: ii, nd
    real(dp), intent(in) :: f_im2, f_im1, f_i, f_ip1, &
                            x_im2, x_im1, x_i, x_ip1, x_p
    real(dp) :: interpol_fyp_cube3
    !
    ! ... local scalars
    real(dp) :: h_im1, h_i, h_ip1, fp_im1, fp_i, &
                s_im1, s_i, s_ip1, p_im1, p_i
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    ! ... local arrays
    !
    !
    if(ii.eq.2) then
       h_i = x_i-x_im1
       h_ip1 = x_ip1-x_i
       s_i = (f_i-f_im1)/h_i
       s_ip1 = (f_ip1-f_i)/h_ip1
       p_i=s_i*(one+h_i/(h_i+h_ip1)) - s_ip1*h_i/(h_ip1+h_i)
       !
       fp_im1=(sign(one,p_i)+sign(one,s_ip1))*min(abs(s_ip1), half*abs(p_i))
       fp_i=(sign(one,s_i)+sign(one,s_ip1))*min(abs(s_i),abs(s_ip1),half*abs(p_i))
    else if(ii.eq.nd) then
       h_im1 = x_im1-x_im2
       h_i = x_i-x_im1
       s_im1 = (f_im1-f_im2)/h_im1
       s_i = (f_i-f_im1)/h_i
       p_im1=(s_im1*h_i + s_i*h_im1)/(h_im1+h_i)
       p_i=s_i*(one+h_i/(h_i+h_im1)) - s_im1*h_i/(h_i+h_im1)
       !
       fp_im1=(sign(one,s_im1)+sign(one,s_i))*min(abs(s_im1),abs(s_i),half*abs(p_im1))
       fp_i=(sign(one,p_i)+sign(one,s_im1))*min(abs(s_im1), half*abs(p_i))
    else
       h_im1 = x_im1-x_im2
       h_i = x_i-x_im1
       h_ip1 = x_ip1-x_i
       s_im1 = (f_im1-f_im2)/h_im1
       s_i = (f_i-f_im1)/h_i
       s_ip1 = (f_ip1-f_i)/h_ip1
       p_im1=(s_im1*h_i + s_i*h_im1)/(h_im1+h_i)
       p_i=(s_i*h_ip1 + s_ip1*h_i)/(h_i+h_ip1)
       !
       fp_im1=(sign(one,s_im1)+sign(one,s_i))*min(abs(s_im1),abs(s_i),half*abs(p_im1))
       fp_i=(sign(one,s_i)+sign(one,s_ip1))*min(abs(s_i),abs(s_ip1),half*abs(p_i))
    endif

    acoeff = -two*s_i/h_i/h_i + (fp_im1+fp_i)/h_i/h_i
    bcoeff = -three*s_i/h_i    + (fp_im1+two*fp_i)/h_i
    ccoeff = fp_i
    dcoeff = f_i
    !
    h_i=x_p-x_i
    interpol_fyp_cube3 = acoeff*h_i**3 + bcoeff*h_i**2 + ccoeff*h_i + dcoeff
    !
  end function interpol_fyp_cube3
!
!***********************************************************************
!                          INDEX FINDING
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine find_index(xp, x, nd, iim2, iim1, ii, iip1)
    !
    !find index for interpolation (binary search)
    !   xp <= x(indx-1) < x(indx) 
    !             or
    !   x(indx-1) < xp <= x(indx)
    !             or
    !   x(indx-1) < x(indx) <= xp
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xp
    real(dp), dimension(nd), intent(in) :: x
    integer(i4b) :: iim2, iim1, ii, iip1
    !
    ! ... local scalars
    integer(i4b) :: i, il, iu, im
    !
    !
    il=1
    iu=nd
    !
    do i=1, nd
       im=floor((il+iu)/two)
       if(x(im).gt.xp) then
          iu=im
       else
          il=im
       endif
       !
       if(iu-il.eq.1) exit
    enddo
    !
    iim2=max(1,iu-2)
    iim1=iu-1
    ii=iu
    iip1=min(nd,iu+1)
    !
  end subroutine find_index
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine find_index_linear(xp, x, nd, indx)
    !
    !find index for interpolation (linear search)
    !   xp <= x(indx-1) < x(indx)
    !             or
    !   x(indx-1) < xp <= x(indx)
    !             or
    !   x(indx-1) < x(indx) <= xp
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xp
    real(dp), dimension(nd), intent(in) :: x
    integer(i4b) :: indx
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    !start always with index 2
    indx=nd
    !
    do i=2, nd
       if(x(i).ge.xp) then
          indx=i
          exit
       endif
    enddo
    !
  end subroutine find_index_linear
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
end module mod_interp1d  
