!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_4p_lin(f_im1jm1, f_ijm1, f_im1j, f_ij, &
                           x_im1, x_i, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using biliniear interpolation
!
!   f(x,y) = f_jm1 + (f_j-f_jm1)*(y_p-y_im1)/(y_j-y_jm1)
!
!      with: f_j   = f_im1j   + (f_ij-f_im1j) *   (x_p-x_im1)/(x_i-x_im1)
!            f_jm1 = f_im1jm1 + (f_ijm1-f_im1jm1)*(x_p-x_im1)/(x_i-x_im1)
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im1jm1, f_ijm1, f_im1j, f_ij, &
                        x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_4p_lin
!
! ... local scalars
real(dp) :: dxi, dx, dyj, dy, rdx, rdy, rdxdy
real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
!
!define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
!
!define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdxdy=rdx*rdy
!
!bilinear interpolation
acoeff=1.d0-rdx-rdy+rdxdy
bcoeff=rdx-rdxdy
ccoeff=rdy-rdxdy
dcoeff=rdxdy
!
interpol2d_4p_lin = acoeff*f_im1jm1 + bcoeff*f_ijm1 + ccoeff*f_im1j + dcoeff*f_ij
!
end function interpol2d_4p_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_4p_log(f_im1jm1, f_ijm1, f_im1j, f_ij, &
                           x_im1, x_i, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!   using logarithmic interpolations in radius of points (im1,jm1), (i,j)
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im1jm1, f_ijm1, f_im1j, f_ij, &
                        x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_4p_log
!
! ... local scalars
real(dp) :: r_im1jm1, r_ijm1, r_im1j, r_ij, r_p, f1, f2, f3
!
!calculate radii
r_im1jm1=sqrt(x_im1**2+y_jm1**2)
r_ijm1=sqrt(x_i**2+y_jm1**2)
r_im1j=sqrt(x_im1**2+y_j**2)
r_ij=sqrt(x_i**2+y_j**2)
!
r_p=sqrt(x_p**2+y_p**2)
!
f1=f_im1jm1*10.d0**(log10(f_ij/f_im1jm1)*log10(r_p/r_ij)/log10(r_ij/r_im1jm1))
f2=f_ijm1*10.d0**(log10(f_ij/f_ijm1)*log10(r_p/r_ijm1)/log10(r_ij/r_ijm1))
f3=f_im1j*10.d0**(log10(f_ij/f_im1j)*log10(r_p/r_im1j)/log10(r_ij/r_im1j))
!
interpol2d_4p_log = (f1+f2+f3)/3.d0
!
end function interpol2d_4p_log
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_4p_lin(x_im1, x_i, y_jm1, y_j, x_p, y_p, &
                          acoeff, bcoeff, ccoeff, dcoeff)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using biliniear interpolation
!         returns coefficients, such that
!            f(x,y) = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: acoeff, bcoeff, ccoeff, dcoeff
!
! ... local scalars
real(dp) :: dxi, dx, dyj, dy, rdx, rdy, rdxdy
!
!define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
!
!define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdxdy=rdx*rdy
!
!bilinear interpolation
acoeff=1.d0-rdx-rdy+rdxdy
bcoeff=rdx-rdxdy
ccoeff=rdy-rdxdy
dcoeff=rdxdy
!
!
end subroutine coeff2d_4p_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_4p_bez(f_im1jm1, f_ijm1, f_im1j, f_ij, &
                           x_im1, x_i, y_jm1, y_j, x_p, y_p, &
                           a, b, c, d)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                    using 2d bezier interpolation
!
!   f(x,y) = (1-ty)^2*f_jm1 + 2*ty*(1-ty)*fc + ty^2*f_j
!
!      with: f_j =   (1-tx)^2*f_im1j   + 2*tx*(1-tx)*fc_j   + tx^2*f_ij
!            f_jm1 = (1-tx)^2*f_im1jm1 + 2*tx*(1-tx)*fc_jm1 + tx^2*f_ijm1
!
!       and: fc_j, fc_jm1, fc chosen from f_im1j, f_ij with weights w_lower=q*f_larger
!
!         returns coefficients, such that
!            f(x,y) = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im1jm1, f_ijm1, f_im1j, f_ij, &
                        x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d
!
! ... local scalars
real(dp), parameter :: fac=1.5d0
real(dp) :: axj, axjm1, bxj, bxjm1, qxj, qxjm1, ay, by, qy, tx, ty, f_j, f_jm1
!
!define deltax, deltay
tx = (x_p-x_im1)/(x_i-x_im1)
ty = (y_p-y_jm1)/(y_j-y_jm1)
!
!on level jm1
if(f_im1jm1.lt.f_ijm1) then
   qxjm1 = fac
elseif(f_ijm1.lt.f_im1jm1) then
   qxjm1 = 1./fac
else
   qxjm1 = 1.
endif
bxjm1 = tx*(2./(1.+qxjm1)-tx*(1.-qxjm1)/(1.+qxjm1))
axjm1 = 1.-bxjm1
f_jm1 = axjm1*f_im1jm1 + bxjm1*f_ijm1
!
!on level j
if(f_im1j.lt.f_ij) then
   qxj = fac
elseif(f_ij.lt.f_im1j) then
   qxj = 1./fac
else
   qxj = 1.
endif
bxj = tx*(2./(1.+qxj)-tx*(1.-qxj)/(1.+qxj))
axj = 1.-bxj
f_j = axj*f_im1j + bxj*f_ij
!
!along y
if(f_jm1.lt.f_j) then
   qy = fac
elseif(f_j.lt.f_jm1) then
   qy = 1./fac
else
   qy = 1.
endif
by = ty*(2./(1.+qy)-ty*(1.-qy)/(1.+qy))
ay = 1.-by
!
!coefficients
a = axjm1*ay
b = bxjm1*ay
c = axj*by
d = bxj*by
!
!
end subroutine coeff2d_4p_bez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_4p_idw(f_im1jm1, f_ijm1, f_im1j, f_ij, &
                           x_im1, x_i, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using inverse distance weighting with four neighbouring points
!
!   f(x,y) = (w_ij*f_ij + w_im1j*f_im1j + w_ijm1*f_ijm1 + w_im1jm1*f_im1jm1) /
!          (w_ij+w_im1j+w_ijm1+w_im1jm1)
!
!      with: w_ij = 1.d0/sqrt((x_p-x_i)^2+(y_p-y_j)^2)^p
!            p - weighting factor
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im1jm1, f_ijm1, f_im1j, f_ij, &
                        x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_4p_idw
!
! ... local scalars
real(dp) :: fac
real(dp) :: d_im1jm1, d_ijm1, d_im1j, d_ij
real(dp) :: w_im1jm1, w_ijm1, w_im1j, w_ij
!
!define weighting factor
fac=0.5d0
!
!calculate distance to each point
d_im1jm1 = (x_p-x_im1)**2+(y_p-y_jm1)**2
d_ijm1 = (x_p-x_i)**2+(y_p-y_jm1)**2
d_im1j = (x_p-x_im1)**2+(y_p-y_j)**2
d_ij = (x_p-x_i)**2+(y_p-y_j)**2
!
!avoid division by zero if d_ij=0, or d_im1j=0, ...
if(d_im1jm1.eq.0.) then 
   interpol2d_4p_idw = f_ijm1
elseif(d_ijm1.eq.0.) then
   interpol2d_4p_idw = f_ijm1
elseif(d_im1j.eq.0.) then
   interpol2d_4p_idw = f_im1j
elseif(d_ij.eq.0.) then
   interpol2d_4p_idw = f_ij
else
   w_im1jm1=1.d0/d_im1jm1**(fac/2.d0)
   w_ijm1=1.d0/d_ijm1**(fac/2.d0)
   w_im1j=1.d0/d_im1j**(fac/2.d0)
   w_ij=1.d0/d_ij**(fac/2.d0)
   interpol2d_4p_idw = (f_im1jm1*w_im1jm1 + f_ijm1*w_ijm1 + f_im1j*w_im1j + f_ij*w_ij) / &
                       (w_im1jm1 + w_ijm1 + w_im1j + w_ij)
endif
!
end function interpol2d_4p_idw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_9p_idw(f_im2jm2, f_im1jm2, f_ijm2, &
                           f_im2jm1, f_im1jm1, f_ijm1, &
                           f_im2j, f_im1j, f_ij, &
                           x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using inverse distance weighting with nine neighbouring points
!
!   f(x,y) = (w_im2jm2*f_im2jm2 + w_im1jm2*f_im1jm2 + w_ijm2*f_ijm2 + $
!            w_im2jm1*f_im2jm1 + w_im1jm1*f_im1jm1 + w_ijm1*f_ijm1 + $
!            w_im2j*f_im2j     + w_im1j*f_im1j     + w_ij*f_ij) / $
!              sum(w_ij)
!
!      with: w_ij = 1.d0/sqrt((x_p-x_i)^2+(y_p-y_j)^2)^p
!            p - weighting factor
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_9p_idw
!
! ... local scalars
real(dp) :: fac
real(dp) :: d_im2jm2, d_im1jm2, d_ijm2, &
            d_im2jm1, d_im1jm1, d_ijm1, &
            d_im2j, d_im1j, d_ij
real(dp) :: w_im2jm2, w_im1jm2, w_ijm2, &
            w_im2jm1, w_im1jm1, w_ijm1, &
            w_im2j, w_im1j, w_ij
!
!define weighting factor
fac=2.d0
!
!calculate distance to each point
d_im2jm2 = (x_p-x_im2)**2+(y_p-y_jm2)**2
d_im1jm2 = (x_p-x_im1)**2+(y_p-y_jm2)**2
d_ijm2 = (x_p-x_i)**2+(y_p-y_jm2)**2
d_im2jm1 = (x_p-x_im2)**2+(y_p-y_jm1)**2
d_im1jm1 = (x_p-x_im1)**2+(y_p-y_jm1)**2
d_ijm1 = (x_p-x_i)**2+(y_p-y_jm1)**2
d_im2j = (x_p-x_im2)**2+(y_p-y_j)**2
d_im1j = (x_p-x_im1)**2+(y_p-y_j)**2
d_ij = (x_p-x_i)**2+(y_p-y_j)**2
!
!avoid division by zero if d_ij=0, or d_im1j=0, ...
if(d_im2jm2.eq.0.) then
   interpol2d_9p_idw = f_im2jm2
elseif(d_im1jm2.eq.0.) then
   interpol2d_9p_idw = f_im1jm2
elseif(d_ijm2.eq.0.) then
   interpol2d_9p_idw = f_ijm2
elseif(d_im2jm1.eq.0.) then
   interpol2d_9p_idw = f_im2jm1
elseif(d_im1jm1.eq.0.) then
   interpol2d_9p_idw = f_im1jm1
elseif(d_ijm1.eq.0.) then
   interpol2d_9p_idw = f_ijm1
elseif(d_im2j.eq.0.) then
   interpol2d_9p_idw = f_im2j
elseif(d_im1j.eq.0.) then
   interpol2d_9p_idw = f_im1j
elseif(d_ij.eq.0.) then
   interpol2d_9p_idw = f_ij
else
!
   w_im2jm2=1.d0/d_im2jm2**(fac/2.d0)
   w_im1jm2=1.d0/d_im1jm2**(fac/2.d0)
   w_ijm2=1.d0/d_ijm2**(fac/2.d0)
   w_im2jm1=1.d0/d_im2jm1**(fac/2.d0)
   w_im1jm1=1.d0/d_im1jm1**(fac/2.d0)
   w_ijm1=1.d0/d_ijm1**(fac/2.d0)
   w_im2j=1.d0/d_im2j**(fac/2.d0)
   w_im1j=1.d0/d_im1j**(fac/2.d0)
   w_ij=1.d0/d_ij**(fac/2.d0)
!
   interpol2d_9p_idw = (f_im2jm2*w_im2jm2 + f_im1jm2*w_im1jm2 + f_ijm2*w_ijm2 + &
            f_im2jm1*w_im2jm1 + f_im1jm1*w_im1jm1 + f_ijm1*w_ijm1 + &
            f_im2j*w_im2j + f_im1j*w_im1j + f_ij*w_ij) / &
           (w_im2jm2 + w_im1jm2 + w_ijm2 + w_im2jm1 + w_im1jm1 + w_ijm1 + &
            w_im2j + w_im1j + w_ij)
endif
!
end function interpol2d_9p_idw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_9p_bez(f_im2jm2, f_im1jm2, f_ijm2, &
                          f_im2jm1, f_im1jm1, f_ijm1, &
                          f_im2j, f_im1j, f_ij, &
                          x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d bezier interpolation
!
!   f(x,y) = (1-ty)^2*f_jm1 + 2*ty*(1-ty)*fc + ty^2*f_j
!
!      with: f_j =   (1-tx)^2*f_im1j   + 2*tx*(1-tx)*fc_j   + tx^2*f_ij
!            f_jm1 = (1-tx)^2*f_im1jm1 + 2*tx*(1-tx)*fc_jm1 + tx^2*f_ijm1
!            f_jm2 = (1-tx)^2*f_im1jm2 + 2*tx*(1-tx)*fc_jm2 + tx^2*f_ijm2
!
!       and: fc_j, fc_jm1, fc_jm2, fc calculated from weighted mean of derivatives
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_9p_bez
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, &
            dyjm1, dyj, dy
real(dp) :: tx, ty, a, b, c, &
            f_jm2, f_jm1, f_j, &
            fc_jm2, fc_jm1, fc_j, fc, &
            fmin_j, fmin_jm1, fmin_jm2, fmin, &
            fmax_j, fmax_jm1, fmax_jm2, fmax
!
!define deltax, delty
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = (x_p-x_im1)/dxi
ty = (y_p-y_jm1)/dyj
!
!calculate control points fc_jm2, fc_jm1, fc_j
a = -dxi**2/2.d0/dxim1/dx
b = 1.d0+(dxi**2-dxim1**2)/2.d0/dxim1/dx
c = dxim1/2.d0/dx
fc_jm2 = f_im2jm2*a + f_im1jm2*b + f_ijm2*c
fc_jm1 = f_im2jm1*a + f_im1jm1*b + f_ijm1*c
fc_j   = f_im2j*a   + f_im1j*b   + f_ij*c
!
!ensure monotonicity on level j
fmin_j=min(f_im1j, f_ij)
fmax_j=max(f_im1j, f_ij)
if(fc_j.gt.fmax_j) then 
   fc_j=fmax_j
elseif(fc_j.lt.fmin_j) then 
   fc_j=fmin_j
endif
!
!ensure monotonicity on level j-1
fmin_jm1=min(f_im1jm1, f_ijm1)
fmax_jm1=max(f_im1jm1, f_ijm1)
if(fc_jm1.gt.fmax_jm1) then
   fc_jm1=fmax_jm1
elseif(fc_jm1.lt.fmin_jm1) then
   fc_jm1=fmin_jm1
endif
!
!ensure monotonicity on level j-2
fmin_jm2=min(f_im1jm2, f_ijm2)
fmax_jm2=max(f_im1jm2, f_ijm2)
if(fc_jm2.gt.fmax_jm2) then
   fc_jm2=fmax_jm2
elseif(fc_jm2.lt.fmin_jm2) then
   fc_jm2=fmin_jm2
endif
!
!calculate function values on each j-level
a=(1.d0-tx)**2
b=2.d0*tx*(1.d0-tx)
c=tx**2
f_jm2 = a*f_im1jm2 + b*fc_jm2 + c*f_ijm2
f_jm1 = a*f_im1jm1 + b*fc_jm1 + c*f_ijm1
f_j   = a*f_im1j   + b*fc_j   + c*f_ij
!
!calculate control point for interpolation along y
a = -dyj**2/2.d0/dyjm1/dy
b = 1.d0+(dyj**2-dyjm1**2)/2.d0/dyjm1/dy
c = dyjm1/2.d0/dy
fc = f_jm2*a + f_jm1*b + f_j*c
!
!write(*,'(a20,8es20.8)') 'function f_j', f_j, f_jm1, f_jm2, fc, a, b, c, fc_jm2
!
!ensure monotonicity along y
fmin=min(f_jm1, f_j)
fmax=max(f_jm1, f_j)
if(fc.gt.fmax) then
   fc=fmax
elseif(fc.lt.fmin) then
   fc=fmin
endif
!
interpol2d_9p_bez = (1.d0-ty)**2 * f_jm1 + 2.d0*ty*(1.d0-ty)*fc + ty**2*f_j
!
end function interpol2d_9p_bez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_9p_bez2(f_im2jm2, f_im1jm2, f_ijm2, &
                            f_im2jm1, f_im1jm1, f_ijm1, &
                            f_im2j, f_im1j, f_ij, &
                            x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using bezier surfaces
!
!   f(x,y) = (1-tx)^2 *    [(1-ty)^2*f_im1jm1 + 2*ty*(1-ty)*fc_01 + ty^2*f_im1j)]
!          + 2*tx*(1-tx) * [(1-ty)^2*fc_10    + 2*ty*(1-ty)*fc_11 + ty^2*fc_12)]
!          + tx^2 *        [(1-ty)^2*f_ijm1   + 2*ty*(1-ty)*fc_21 + ty^2*fc_ij)]
!
!      with: tx = (x_p-x_im1)/(x_i-x_im1)
!            ty = (y_p-y_jm1)/(y_j-y_jm1)
!
!       and: fc_01, fc_10, fc_11, fc_12, fc_21 calculated from weighted mean of derivatives
!
!on input: 
!
! y_j      f_im2j----------f_im1j-----fc_10----f_ij
!  |          |              |    x             |  
!  |          |              |(x_p,y_p)         |  
!  |          |            fc_01     fc_11    fc_21
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1----fc_12---f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_9p_bez2
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, &
            dyjm1, dyj, dy
real(dp) :: tx, ty, ax, bx, cx, ay, by, cy, a, b, c, &
            fc_01, fc_21, fc_10, fc_12, fc_11, &
            fmin, fmax
!
!define deltax, delty
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = (x_p-x_im1)/dxi
ty = (y_p-y_jm1)/dyj
!
!calculate all control points
ax = -dxi**2/2.d0/dxim1/dx
bx = 1.d0+(dxi**2-dxim1**2)/2.d0/dxim1/dx
cx = dxim1/2.d0/dx
!
ay = -dyj**2/2.d0/dyjm1/dy
by = 1.d0+(dyj**2-dyjm1**2)/2.d0/dyjm1/dy
cy = dyjm1/2.d0/dy
!
fc_01 = ay*f_im1jm2 + by*f_im1jm1 + cy*f_im1j
fc_21 = ay*f_ijm2   + by*f_ijm1   + cy*f_ij
fc_10 = ax*f_im2jm1 + bx*f_im1jm1 + cx*f_ijm1
fc_12 = ax*f_im2j   + bx*f_im1j   + cx*f_ij  
fc_11 = ay*f_im1jm2 + ax*f_im2jm1 + (bx+by-1.d0)*f_im1jm1 + cy*f_im1j + cx*f_ijm1
!
!ensure monotonicity
fmin = min(f_im1jm1, f_ijm1, f_im1j, f_ij)
fmax = max(f_im1jm1, f_ijm1, f_im1j, f_ij)
if(fc_01.gt.fmax) then 
   fc_01=fmax
elseif(fc_01.lt.fmin) then
   fc_01=fmin
endif
!
if(fc_21.gt.fmax) then 
   fc_21=fmax
elseif(fc_21.lt.fmin) then
   fc_21=fmin
endif
!
if(fc_10.gt.fmax) then 
   fc_10=fmax
elseif(fc_10.lt.fmin) then
   fc_10=fmin
endif
!
if(fc_12.gt.fmax) then 
   fc_12=fmax
elseif(fc_12.lt.fmin) then
   fc_12=fmin
endif
!
if(fc_11.gt.fmax) then
   fc_11=fmax
elseif(fc_11.lt.fmin) then
   fc_11=fmin
endif
!
a=(1.d0-ty)**2
b=2.d0*ty*(1.d0-ty)
c=ty**2
!
interpol2d_9p_bez2 = (1.d0-tx)**2 *       (a*f_im1jm1 + b*fc_01 + c*f_im1j) + &
                      2.d0*tx*(1.d0-tx) * (a*fc_10    + b*fc_11 + c*fc_12)  + &
                     tx**2*               (a*f_ijm1   + b*fc_21 + c*f_ij)
!
end function interpol2d_9p_bez2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_9p_quad(f_im2jm2, f_im1jm2, f_ijm2, &
                            f_im2jm1, f_im1jm1, f_ijm1, &
                            f_im2j, f_im1j, f_ij, &
                            x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d quadratic functions
!         ensure monotonicity via linear interpolation
!              if derivatives have opposite sign
!         ensure positive curvature via linear interpolation
!              (to avoid large overestimation of step-function interpolation)
!
!   f(x,y) = sum_i sum_j (a_ij * (x-x_i)^i * (y-y_j)^j)
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_9p_quad
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, &
            dyjm1, dyj, dy
real(dp) :: tx, ty, a, b, c, a2, b2, c2, a3, b3, c3, &
           fp_im1jm2, fp_im1jm1, fp_im1j, & 
           fp_ijm2, fp_ijm1, fp_ij, &
           fp_jm2, fp_jm1, fp_j, &
           f_jm2, f_jm1, f_j, &
           fp2_im1jm2, fp2_im1jm1, fp2_im1j, fp2_jm1
!
!define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = x_p-x_im1
ty = y_p-y_jm1
!
!derivatives at im1 on each j-level
a = -dxi/dxim1/dx
b = (dxi**2-dxim1**2)/dxi/dxim1/dx
c = dxim1/dxi/dx
fp_im1jm2 = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2
fp_im1jm1 = a*f_im2jm1 + b*f_im1jm1 + c*f_ijm1
fp_im1j   = a*f_im2j   + b*f_im1j   + c*f_ij  
!
!derivatives at i on each j-level
a2 = -a
b2 = b-2.d0/dxim1
c2 = c+2.d0/dx
fp_ijm2 = a2*f_im2jm2 + b2*f_im1jm2 + c2*f_ijm2
fp_ijm1 = a2*f_im2jm1 + b2*f_im1jm1 + c2*f_ijm1
fp_ij   = a2*f_im2j   + b2*f_im1j   + c2*f_ij
!
!second derivative at im1 on each j-level
a2 = 2.d0/dxim1/dx
b2 = -2.d0/dxi/dxim1
c2 = 2.d0/dxi/dx
fp2_im1jm2 = a2*f_im2jm2 + b2*f_im1jm2 + c2*f_ijm2
fp2_im1jm1 = a2*f_im2jm1 + b2*f_im1jm1 + c2*f_ijm1
fp2_im1j   = a2*f_im2j   + b2*f_im1j   + c2*f_ij
!
!interpolate on each j-level (and ensure monotonicity)
a3 = tx**2/dxim1/dx + a*tx
b3 = 1.d0-tx**2/dxi/dxim1 + b*tx
c3 = tx**2/dxi/dx + c*tx
if(fp_ij*fp_im1j.lt.0.) then
   f_j = f_im1j + (f_ij-f_im1j)*tx/dxi   !linear interpolation
else
!ensure positive curvature
!   if(fp2_im1j.lt.0.) then
!      f_j = f_im1j + (f_ij-f_im1j)*tx/dxi   !linear interpolation
!   else
!      f_j =a3*f_im2j + b3*f_im1j + c3*f_ij
!   endif
endif
!
if(fp_ijm1*fp_im1jm1.lt.0.) then
   f_jm1 = f_im1jm1 + (f_ijm1-f_im1jm1)*tx/dxi  !linear interpolation
else
!ensure positive curvature
!   if(fp2_im1jm1.lt.0.) then
!      f_jm1 = f_im1jm1 + (f_ijm1-f_im1jm1)*tx/dxi   !linear interpolation
!   else
      f_jm1 =a3*f_im2jm1 + b3*f_im1jm1 + c3*f_ijm1
!   endif
endif
!
if(fp_ijm2*fp_im1jm2.lt.0.) then
   f_jm2 = f_im1jm2 + (f_ijm2-f_im1jm2)*tx/dxi  !linear interpolation
else
!ensure positive curvature
!   if(fp2_im1jm2.lt.0.) then
!      f_jm2 = f_im1jm2 + (f_ijm2-f_im1jm2)*tx/dxi   !linear interpolation
!   else
      f_jm2 =a3*f_im2jm2 + b3*f_im1jm2 + c3*f_ijm2
!   endif
endif
!
!now, interpolation along y-axis
a = -dyj/dyjm1/dy
b = (dyj**2-dyjm1**2)/dyj/dyjm1/dy
c = dyjm1/dyj/dy
!derivatives at j and jm1
fp_jm1 =  a*f_jm2 + b*f_jm1 + c*f_j
fp_j   = -a*f_jm2 + (b-2.d0/dyjm1)*f_jm1 + (c+2.d0/dy)*f_j
!
!second derivative at jm1
fp2_jm1 = 2.d0*f_jm2/dyjm1/dy - 2.d0*f_jm1/dyj/dyjm1 + 2.d0*f_j/dyj/dy
!
if(fp_j*fp_jm1.lt.0.) then
   interpol2d_9p_quad = f_jm1 + (f_j-f_jm1)*ty/dyj   !linear interpolation
else
!ensure positive curvature
!   if(fp2_jm1.lt.0.) then
!      interpol2d_9p_quad = f_jm1 + (f_j-f_jm1)*ty/dyj   !linear interpolation
!   else
      interpol2d_9p_quad = (ty**2/dyjm1/dy + a*ty)*f_jm2 + (1.d0-ty**2/dyj/dyjm1 + b*ty)*f_jm1 + (ty**2/dyj/dy + c*ty)*f_j
!   endif
endif
!
!
!
end function interpol2d_9p_quad
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_9p_quad2(f_im2jm2, f_im1jm2, f_ijm2, &
                             f_im2jm1, f_im1jm1, f_ijm1, &
                             f_im2j, f_im1j, f_ij, &
                             x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d quadratic functions
!
!   f(x,y) = sum_i sum_j (a_ij * (x-x_i)^i * (y-y_j)^j)
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_9p_quad2
!
! ... local scalars
real(dp) :: dxim1, dxi, dxp, dx, &
            dyjm1, dyj, dyp, dy
real(dp) :: ax, bx, cx, ay, by, cy
!
!define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dxp = x_p-x_im1
dx = dxi+dxim1
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dyp = y_p-y_jm1
dy = dyj+dyjm1
!
!
!standard parabola
ax = dxp*(dxp-dxi)/dxim1/dx
bx = 1.d0+dxp*(dxi-dxim1-dxp)/dxi/dxim1
cx = dxp*(dxp+dxim1)/dxi/dx
ay = dyp*(dyp-dyj)/dyjm1/dy
by = 1.d0+dyp*(dyj-dyjm1-dyp)/dyj/dyjm1
cy = dyp*(dyp+dyjm1)/dyj/dy
!
interpol2d_9p_quad2 = ay*(ax*f_im2jm2 + bx*f_im1jm2 + cx*f_ijm2) + &
                      by*(ax*f_im2jm1 + bx*f_im1jm1 + cx*f_ijm1) + &
                      cy*(ax*f_im2j   + bx*f_im1j   + cx*f_ij)
!
!if fixed derivative at point i,j (zero)
!bx=2.d0*dxp/dxi-(dxp/dxi)**2
!ax=1.d0-bx
!by=2.d0*dyp/dyj-(dyp/dyj)**2
!ay=1.d0-by
!interpol2d_9p_quad2 = ay*(ax*f_im1jm1+bx*f_ijm1) + by*(ax*f_im1j+bx*f_ij)
!
!
!
!
end function interpol2d_9p_quad2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_9p_quad(f_im2jm2, f_im1jm2, f_ijm2, &
                          f_im2jm1, f_im1jm1, f_ijm1, &
                         f_im2j, f_im1j, f_ij, &
                         x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
                         a, b, c, d, e, f, g, h, i)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d quadratic functions
!         ensure monotonicity via linear interpolation
!              if derivatives have opposite sign
!         ensure positive curvature via linear interpolation
!              (to avoid large overestimation of step-function interpolation)
!
!   f(x,y) = sum_i sum_j (a_ij * (x-x_i)^i * (y-y_j)^j) = 
!          = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 + 
!            d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 + 
!            g*f_im2j   + h*f_im1j   + i*f_ij
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
use mod_interp2d, only: interpolation_threshold
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d, e, f, g, h, i
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, tx, &
            dyjm1, dyj, dy, ty
real(dp) :: fp_im1jm2, fp_im1jm1, fp_im1j, & 
            fp_ijm2, fp_ijm1, fp_ij, &
            fp_jm2, fp_jm1, fp_j, &
            f_jm2, f_jm1, f_j
real(dp) :: axjm2, bxjm2, cxjm2, &
            axjm1, bxjm1, cxjm1, &
            axj, bxj, cxj, &
            a1, b1, c1, a2, b2, c2, a3, b3, c3, &
            ay, by, cy
!
!define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = x_p-x_im1
ty = y_p-y_jm1
!
!derivatives at im1 on each j-level
a1 = -dxi/dxim1/dx
b1 = (dxi**2-dxim1**2)/dxi/dxim1/dx
c1 = dxim1/dxi/dx
fp_im1jm2 = a1*f_im2jm2 + b1*f_im1jm2 + c1*f_ijm2
fp_im1jm1 = a1*f_im2jm1 + b1*f_im1jm1 + c1*f_ijm1
fp_im1j   = a1*f_im2j   + b1*f_im1j   + c1*f_ij  
!
!derivatives at i on each j-level
a2 = -a1
b2 = b1-2.d0/dxim1
c2 = c1+2.d0/dx
fp_ijm2 = a2*f_im2jm2 + b2*f_im1jm2 + c2*f_ijm2
fp_ijm1 = a2*f_im2jm1 + b2*f_im1jm1 + c2*f_ijm1
fp_ij   = a2*f_im2j   + b2*f_im1j   + c2*f_ij
!
!interpolate on each j-level (and ensure monotonicity)
!note: use threshold of fp_im1/fp_i or fp_i/fp_im1 in order to avoid 
!      oscillations in iteration scheme
if(fp_ij*fp_im1j.le.0.) then
   axj=0.d0
   cxj=tx/dxi
   bxj=1.d0-cxj
elseif(fp_ij/fp_im1j.lt.interpolation_threshold.or.fp_im1j/fp_ij.lt.interpolation_threshold) then
   axj=0.d0
   cxj=tx/dxi
   bxj=1.d0-cxj
else
   axj = tx**2/dxim1/dx + a1*tx
   bxj = 1.d0-tx**2/dxi/dxim1 + b1*tx
   cxj = tx**2/dxi/dx + c1*tx
endif
f_j = axj*f_im2j + bxj*f_im1j + cxj*f_ij
!
if(fp_ijm1*fp_im1jm1.le.0.) then
   axjm1=0.d0
   cxjm1=tx/dxi
   bxjm1=1.d0-cxjm1
elseif(fp_ijm1/fp_im1jm1.lt.interpolation_threshold.or.fp_im1jm1/fp_ijm1.lt.interpolation_threshold) then
   axjm1=0.d0
   cxjm1=tx/dxi
   bxjm1=1.d0-cxjm1
else
   axjm1 = tx**2/dxim1/dx + a1*tx
   bxjm1 = 1.d0-tx**2/dxi/dxim1 + b1*tx
   cxjm1 = tx**2/dxi/dx + c1*tx
endif
f_jm1 = axjm1*f_im2jm1 + bxjm1*f_im1jm1 + cxjm1*f_ijm1
!
if(fp_ijm2*fp_im1jm2.le.0.) then
   axjm2=0.d0
   cxjm2=tx/dxi
   bxjm2=1.d0-cxjm2
elseif(fp_ijm2/fp_im1jm2.lt.interpolation_threshold.or.fp_im1jm2/fp_ijm2.lt.interpolation_threshold) then
   axjm2=0.d0
   cxjm2=tx/dxi
   bxjm2=1.d0-cxjm2
else
   axjm2 = tx**2/dxim1/dx + a1*tx
   bxjm2 = 1.d0-tx**2/dxi/dxim1 + b1*tx
   cxjm2 = tx**2/dxi/dx + c1*tx
endif
f_jm2 = axjm2*f_im2jm2 + bxjm2*f_im1jm2 + cxjm2*f_ijm2
!
!
!
!now, interpolation along y-axis
a3 = -dyj/dyjm1/dy
b3 = (dyj**2-dyjm1**2)/dyj/dyjm1/dy
c3 = dyjm1/dyj/dy
!derivatives at j and jm1
fp_jm1 =  a3*f_jm2 + b3*f_jm1 + c3*f_j
fp_j   = -a3*f_jm2 + (b3-2.d0/dyjm1)*f_jm1 + (c3+2.d0/dy)*f_j
!
if(fp_j*fp_jm1.le.0.) then
   ay=0.d0
   cy=ty/dyj
   by=1.d0-cy
elseif(fp_j/fp_jm1.lt.interpolation_threshold.or.fp_jm1/fp_j.lt.interpolation_threshold) then
   ay=0.d0
   cy=ty/dyj
   by=1.d0-cy
else
   ay = ty**2/dyjm1/dy + a3*ty
   by = 1.d0-ty**2/dyj/dyjm1 + b3*ty
   cy = ty**2/dyj/dy + c3*ty
endif
!
a=axjm2*ay
b=bxjm2*ay
c=cxjm2*ay
d=axjm1*by
e=bxjm1*by
f=cxjm1*by
g=axj*cy
h=bxj*cy
i=cxj*cy
!
!
!
end subroutine coeff2d_9p_quad
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_9p_quad2(f_im2jm2, f_im1jm2, f_ijm2, &
                            f_im2jm1, f_im1jm1, f_ijm1, &
                            f_im2j, f_im1j, f_ij, & 
                            x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
                            a, b, c, d, e, f, g, h, i)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d quadratic functions
!
!   f(x,y) = sum_i sum_j (a_ij * (x-x_i)^i * (y-y_j)^j) = 
!          = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 + 
!            d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 + 
!            g*f_im2j   + h*f_im1j   + i*f_ij
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments

real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d, e, f, g, h, i
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, &
            dyjm1, dyj, dy
real(dp) :: tx, ty, ax, bx, cx, ay, by, cy
!
!define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = x_p-x_im1
ty = y_p-y_jm1
!
a = -dxi/dxim1/dx
b = (dxi**2-dxim1**2)/dxi/dxim1/dx
c = dxim1/dxi/dx
ax = tx**2/dxim1/dx + a*tx
bx = 1.d0-tx**2/dxi/dxim1 + b*tx
cx = tx**2/dxi/dx + c*tx
!
a = -dyj/dyjm1/dy
b = (dyj**2-dyjm1**2)/dyj/dyjm1/dy
c = dyjm1/dyj/dy
ay = ty**2/dyjm1/dy + a*ty
by = 1.d0-ty**2/dyj/dyjm1 + b*ty
cy = ty**2/dyj/dy + c*ty
!
a=ax*ay
b=bx*ay
c=cx*ay
d=ax*by
e=bx*by
f=cx*by
g=ax*cy
h=bx*cy
i=cx*cy
!
!
!
end subroutine coeff2d_9p_quad2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_9p_bez(f_im2jm2, f_im1jm2, f_ijm2, &
                          f_im2jm1, f_im1jm1, f_ijm1, &
                          f_im2j, f_im1j, f_ij, &
                          x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
                          a, b, c, d, e, f, g, h, i)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using 2d bezier interpolation
!
!   f(x,y) = (1-ty)^2*f_jm1 + 2*ty*(1-ty)*fc + ty^2*f_j
!
!      with: f_j =   (1-tx)^2*f_im1j   + 2*tx*(1-tx)*fc_j   + tx^2*f_ij
!            f_jm1 = (1-tx)^2*f_im1jm1 + 2*tx*(1-tx)*fc_jm1 + tx^2*f_ijm1
!            f_jm2 = (1-tx)^2*f_im1jm2 + 2*tx*(1-tx)*fc_jm2 + tx^2*f_ijm2
!
!       and: fc_j, fc_jm1, fc_jm2, fc calculated from weighted mean of derivatives
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d, e, f, g, h, i
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, &
            dyjm1, dyj, dy
real(dp) :: tx, ty, ax, bx, cx, ay, by, cy, ay2 ,by2, cy2, &
            axt, bxt, cxt, ayt, byt, cyt, &
            axtjm2, bxtjm2, cxtjm2, axtjm1, bxtjm1, cxtjm1, &
            axtj, bxtj, cxtj, &
            axj, bxj, cxj, axjm1, bxjm1, cxjm1, axjm2, bxjm2, cxjm2, &
            f_jm2, f_jm1, f_j, &
            fc_jm2, fc_jm1, fc_j, fc, &
            fmin_j, fmin_jm1, fmin_jm2, fmin, &
            fmax_j, fmax_jm1, fmax_jm2, fmax
!
!define deltax, delty
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
!
tx = (x_p-x_im1)/dxi
ty = (y_p-y_jm1)/dyj
!
ax = (1.d0-tx)**2
bx = 2.d0*tx*(1.d0-tx)
cx = tx**2
!
!calculate control points fc_jm2, fc_jm1, fc_j
axt = -dxi**2/2.d0/dxim1/dx
!bxt = 1.d0+(dxi**2-dxim1**2)/2.d0/dxim1/dx
bxt = dx/2.d0/dxim1
cxt = dxim1/2.d0/dx
fc_jm2 = f_im2jm2*axt + f_im1jm2*bxt + f_ijm2*cxt
fc_jm1 = f_im2jm1*axt + f_im1jm1*bxt + f_ijm1*cxt
fc_j   = f_im2j*axt   + f_im1j*bxt   + f_ij*cxt
!
!ensure monotonicity on level j
if(f_ij.ge.f_im1j) then
   if(fc_j.gt.f_ij) then
      axtj=0.
      bxtj=0.
      cxtj=1.
   elseif(fc_j.lt.f_im1j) then
      axtj=0.
      bxtj=1.
      cxtj=0.
   else
      axtj=axt
      bxtj=bxt
      cxtj=cxt
   endif
elseif(f_ij.le.f_im1j) then
   if(fc_j.lt.f_ij) then
      axtj=0.
      bxtj=0.
      cxtj=1.
   elseif(fc_j.gt.f_im1j) then
      axtj=0.
      bxtj=1.
      cxtj=0.
   else
      axtj=axt
      bxtj=bxt
      cxtj=cxt
   endif
endif
axj = axtj*bx
bxj = ax+bxtj*bx
cxj = cx + cxtj*bx
f_j = axj*f_im2j + bxj*f_im1j + cxj*f_ij
!
!ensure monotonicity on level j-1
if(f_ijm1.ge.f_im1jm1) then
   if(fc_jm1.gt.f_ijm1) then
      axtjm1=0.
      bxtjm1=0.
      cxtjm1=1.
   elseif(fc_jm1.lt.f_im1jm1) then
      axtjm1=0.
      bxtjm1=1.
      cxtjm1=0.
   else
      axtjm1=axt
      bxtjm1=bxt
      cxtjm1=cxt
   endif
elseif(f_ijm1.le.f_im1jm1) then
   if(fc_jm1.lt.f_ijm1) then
      axtjm1=0.
      bxtjm1=0.
      cxtjm1=1.
   elseif(fc_jm1.gt.f_im1jm1) then
      axtjm1=0.
      bxtjm1=1.
      cxtjm1=0.
   else
      axtjm1=axt
      bxtjm1=bxt
      cxtjm1=cxt
   endif
endif
axjm1 = axtjm1*bx
bxjm1 = ax+bxtjm1*bx
cxjm1 = cx + cxtjm1*bx
f_jm1 = axjm1*f_im2jm1 + bxjm1*f_im1jm1 + cxjm1*f_ijm1
!
!ensure monotonicity on level j-2
if(f_ijm2.ge.f_im1jm2) then
   if(fc_jm2.gt.f_ijm2) then
      axtjm2=0.
      bxtjm2=0.
      cxtjm2=1.
!      write(*,*) 't1', f_ijm2, f_im1jm2, fc_jm2
   elseif(fc_jm2.lt.f_im1jm2) then
      axtjm2=0.
      bxtjm2=1.
      cxtjm2=0.
   else
      axtjm2=axt
      bxtjm2=bxt
      cxtjm2=cxt
   endif
elseif(f_ijm2.le.f_im1jm2) then
   if(fc_jm2.lt.f_ijm2) then
      axtjm2=0.
      bxtjm2=0.
      cxtjm2=1.
!      write(*,*) 't2', f_ijm2, f_im1jm2, fc_jm2
   elseif(fc_jm2.gt.f_im1jm2) then
      axtjm2=0.
      bxtjm2=1.
      cxtjm2=0.
   else
      axtjm2=axt
      bxtjm2=bxt
      cxtjm2=cxt
   endif
endif
axjm2 = axtjm2*bx
bxjm2 = ax + bxtjm2*bx
cxjm2 = cx + cxtjm2*bx
f_jm2 = axjm2*f_im2jm2 + bxjm2*f_im1jm2 + cxjm2*f_ijm2
!
!
!calculate control point for interpolation along y
ayt = -dyj**2/2.d0/dyjm1/dy
!byt = 1.d0+(dyj**2-dyjm1**2)/2.d0/dyjm1/dy
byt = dy/2.d0/dyjm1
cyt = dyjm1/2.d0/dy
fc = f_jm2*ayt + f_jm1*byt + f_j*cyt
!
!write(*,*) axjm2, bxjm2, cxjm2
!write(*,*) axtjm2, bxtjm2, cxtjm2
!write(*,*) f_im2jm2, f_im1jm2, f_ijm2
!write(*,'(a20,8es20.8)') 'routine f_j', f_j, f_jm1, f_jm2, fc, ayt, byt, cyt, fc_jm2
!
!ensure monotonicity along y
if(f_j.ge.f_jm1) then
   if(fc.gt.f_j) then
      ayt=0.
      byt=0.
      cyt=1.
   elseif(fc.lt.f_jm1) then
      ayt=0.
      byt=1.
      cyt=0.
   endif
elseif(f_j.le.f_jm1) then
   if(fc.lt.f_j) then
      ayt=0.
      byt=0.
      cyt=1.
   elseif(fc.gt.f_jm1) then
      ayt=0.
      byt=1.
      cyt=0.
   endif
endif
!
ay2 = (1.d0-ty)**2
by2 = 2.d0*ty*(1.d0-ty)
cy2 = ty**2
!
ay = ayt*by2
by = ay2+byt*by2
cy = cy2+cyt*by2
!
a = ay*axjm2
b = ay*bxjm2
c = ay*cxjm2
d = by*axjm1
e = by*bxjm1
f = by*cxjm1
g = cy*axj
h = cy*bxj
i = cy*cxj
!
end subroutine coeff2d_9p_bez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_9p_bez2(f_im2jm2, f_im1jm2, f_ijm2, &
                           f_im2jm1, f_im1jm1, f_ijm1, &
                           f_im2j, f_im1j, f_ij, &
                           x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
                           a, b, c, d, e, f, g, h, i)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using bezier surfaces
!
!   f(x,y) = (1-tx)^2 *    [(1-ty)^2*f_im1jm1 + 2*ty*(1-ty)*fc_01 + ty^2*f_im1j)]
!          + 2*tx*(1-tx) * [(1-ty)^2*fc_10    + 2*ty*(1-ty)*fc_11 + ty^2*fc_12)]
!          + tx^2 *        [(1-ty)^2*f_ijm1   + 2*ty*(1-ty)*fc_21 + ty^2*fc_ij)]
!
!      with: tx = (x_p-x_im1)/(x_i-x_im1)
!            ty = (y_p-y_jm1)/(y_j-y_jm1)
!
!       and: fc_01, fc_10, fc_11, fc_12, fc_21 calculated from inverse distance weighting
!
!on input: 
!
! y_j      f_im2j----------f_im1j-----fc_10----f_ij
!  |          |              |    x             |  
!  |          |              |(x_p,y_p)         |  
!  |          |            fc_01     fc_11    fc_21
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1----fc_12---f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d, e, f, g, h, i
!
! ... local scalars
real(dp) :: dx1, dx2, dx3, dx4, dx1s, dx2s, dx3s, dx4s, &
            dy1, dy2, dy3, dy4, dy1s, dy2s, dy3s, dy4s, &
            fac, norm01, norm10, norm12, norm21, norm11
real(dp) :: tx, ty, ax, bx, cx, ay, by, cy, &
            fc_01, fc_21, fc_10, fc_12, fc_11, &
            fmin, fmax
!
!define deltax, delty
stop 'to be debugged'
dx1 = x_i-x_im1
dx2 = x_i-x_im2
dx3 = (x_i+x_im1-2.d0*x_im2)/2.d0
dx4 = x_im1-x_im2
dx1s = dx1**2
dx2s = dx2**2
dx3s = dx3**2
dx4s = dx4**2
!
dy1 = y_j-y_jm1
dy2 = y_j-y_jm2
dy3 = (y_j+y_jm1-2.d0*y_jm2)/2.d0
dy4 = y_jm1-y_jm2
dy1s = dy1**2
dy2s = dy2**2
dy3s = dy3**2
dy4s = dy4**2
!
tx = (x_p-x_im1)/dx1
ty = (y_p-y_jm1)/dy1
!
ax=(1.d0-tx)**2
bx=2.d0*tx*(1.d0-tx)
cx=tx**2
!
ay=(1.d0-ty)**2
by=2.d0*ty*(1.d0-ty)
cy=ty**2
!
!calculate all control points via inverse distance weighting (weight = distance**fac)
fac = -4.d0
norm01 = (dx4s+dy3s)**fac + dy3s**fac/2. + (dx1s+dy3s)**fac + (dx4s+dy1s/4.)**fac + &
         (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dx4s+dy1s/4.)**fac + (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac
norm10 = (dx3s+dy4s)**fac + (dx1s/4.+dy4s)**fac + (dx1s/4. + dy4s)**fac + dx3s**fac + (dx1s/4.)**fac + &
         (dx1s/4.)**fac + (dx1s/4.)**fac + (dx3s+dy1s)**fac + (dx1s/4.+dy1s)**fac + (dx1s/4.+dy1s)**fac
norm12 = (dx3s+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx3s+dy1s)**fac + &
         (dx1s/4.+dy1s)**fac + (dx1s/2.+dy1s)**fac + dx3s**fac+2.*(dx1s/4.)**fac
norm21 = (dx2s+dx3s)**fac + (dx1s+dy3s)**fac + dy3s**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + &
         (dy1s/4.)**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dy1s/4.)**fac
norm11 = (dx3s+dy3s)**fac + 2.*(dx1s/4.+dy3s)**fac + 2.*(dx3s+dy1s/4.)**fac + 4.*(dx1s/4.+dy1s/4.)**fac
!
a = ax*by/norm01*(dx4s+dy3s)**fac + bx*ay/norm10*(dx3s+dy4s)**fac + bx*by/norm11*(dx3s+dy3s)**fac + &
    bx*cy/norm12*(dx3s+dy2s)**fac + cx*by/norm21*(dx2s+dy3s)**fac
b = ax*by/norm01*dy3s**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
    bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*(dx1s+dy3s)**fac
c = ax*by/norm01*(dx1s+dy3s)**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
    bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*dy3s**fac
d = ax*by/norm01*(dx4s+dy1s/4.)**fac + bx*ay/norm10*dx3s**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
    bx*cy/norm12*(dx3s+dy1s)**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
e = ay*ay + ay*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
    bx*cy/norm12*(dx1s/4.+dy1s)**fac + cy*bx/norm21*(dx1s+dy1s/4.)**fac
f = cx*ay + ax*by/norm01*(dx1s+dy1s/4.)**fac +bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
    bx*cy/norm12*(dx1s/4.+dy1s)**fac + cy*bx/norm21*(dy1s/4.)**fac
g = ax*by/norm01*(dx1s+dy1s/4.)**fac + bx*ay/norm10*(dx3s+dy1s)**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
    bx*cy/norm21*dx3s**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
h = ax*cy + ax*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
    bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dx1s+dy1s/4.)**fac
i = cx*cy + ax*by/norm01*(dx1s+dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
    bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dy1s/4.)**fac
!
end subroutine coeff2d_9p_bez2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_9p_idw(f_im2jm2, f_im1jm2, f_ijm2, &
                          f_im2jm1, f_im1jm1, f_ijm1, &
                          f_im2j, f_im1j, f_ij, &
                          x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
                          a, b, c, d, e, f, g, h, i)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using inverse distance weighting with nine neighbouring points
!
!   f(x,y) = (w_im2jm2*f_im2jm2 + w_im1jm2*f_im1jm2 + w_ijm2*f_ijm2 + $
!             w_im2jm1*f_im2jm1 + w_im1jm1*f_im1jm1 + w_ijm1*f_ijm1 + $
!             w_im2j*f_im2j     + w_im1j*f_im1j     + w_ij*f_ij) / $
!             sum(w_ij)
!
!      with: w_ij = 1.d0/sqrt((x_p-x_i)^2+(y_p-y_j)^2)^p
!            p - weighting factor
!
!on input: 
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |         x        |  
!  |          |              |     (x_p,y_p)    |  
!  |          |              |                  |  
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!  |          |              |                  |  
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im2jm2, f_im1jm2, f_ijm2, &
                        f_im2jm1, f_im1jm1, f_ijm1, &
                        f_im2j, f_im1j, f_ij, &
                        x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
real(dp), intent(out) :: a, b, c, d, e, f, g, h, i
!
! ... local scalars
real(dp) :: fac, norm
real(dp) :: d_im2jm2, d_im1jm2, d_ijm2, &
            d_im2jm1, d_im1jm1, d_ijm1, &
            d_im2j, d_im1j, d_ij
!
!define weighting factor
fac=1.5d0
!
!calculate distance to each point
d_im2jm2 = (x_p-x_im2)**2+(y_p-y_jm2)**2
d_im1jm2 = (x_p-x_im1)**2+(y_p-y_jm2)**2
d_ijm2 = (x_p-x_i)**2+(y_p-y_jm2)**2
d_im2jm1 = (x_p-x_im2)**2+(y_p-y_jm1)**2
d_im1jm1 = (x_p-x_im1)**2+(y_p-y_jm1)**2
d_ijm1 = (x_p-x_i)**2+(y_p-y_jm1)**2
d_im2j = (x_p-x_im2)**2+(y_p-y_j)**2
d_im1j = (x_p-x_im1)**2+(y_p-y_j)**2
d_ij = (x_p-x_i)**2+(y_p-y_j)**2
!
!avoid division by zero if d_ij=0, or d_im1j=0, ...
if(d_im2jm2.eq.0.) then
   a=1.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_im1jm2.eq.0.) then
   a=0.d0
   b=1.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_ijm2.eq.0.) then
   a=0.d0
   b=0.d0
   c=1.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_im2jm1.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=1.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_im1jm1.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=1.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_ijm1.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=1.d0
   g=0.d0
   h=0.d0
   i=0.d0
elseif(d_im2j.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=1.d0
   h=0.d0
   i=0.d0
elseif(d_im1j.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=1.d0
   i=0.d0
elseif(d_ij.eq.0.) then
   a=0.d0
   b=0.d0
   c=0.d0
   d=0.d0
   e=0.d0
   f=0.d0
   g=0.d0
   h=0.d0
   i=1.d0
else
!
   a=1.d0/d_im2jm2**(fac/2.d0)
   b=1.d0/d_im1jm2**(fac/2.d0)
   c=1.d0/d_ijm2**(fac/2.d0)
   d=1.d0/d_im2jm1**(fac/2.d0)
   e=1.d0/d_im1jm1**(fac/2.d0)
   f=1.d0/d_ijm1**(fac/2.d0)
   g=1.d0/d_im2j**(fac/2.d0)
   h=1.d0/d_im1j**(fac/2.d0)
   i=1.d0/d_ij**(fac/2.d0)

   norm = a+b+c+d+e+f+g+h+i

   a=a/norm
   b=b/norm
   c=c/norm
   d=d/norm
   e=e/norm
   f=f/norm
   g=g/norm
   h=h/norm
   i=i/norm
!
endif
!
end subroutine coeff2d_9p_idw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol2d_4p_cube(f_im1jm1, f_ijm1, f_im1j, f_ij, &
                            x_im1, x_i, y_jm1, y_j, x_p, y_p)
!
!         interpolates values given on a 2d grid onto point x_p, y_p
!                 using bicubic interpolation with predefined derivatives
!
!on input: 
!
!          f_im1j--------------f_ij          y_j
!            |                  |             |
!            |                  |             |
!            |         x        |             |
!            |     (x_p,y_p)    |             |
!            |                  |             |
!        f_im1jm1------------f_ijm1         x_im1--------------x_i
!                                           y_jm1
!                           
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p
!
use prog_type
!
implicit none
!
! ... argments
real(dp), intent(in) :: f_im1jm1, f_ijm1, f_im1j, f_ij, &
                        x_im1, x_i, y_jm1, y_j, x_p, y_p
real(dp) :: interpol2d_4p_cube
!
! ... local scalars
real(dp) :: dxi, dx, dyj, dy, rdx, rdy
real(dp) :: ax, bx, ay, by
!
rdx=(x_p-x_im1)/(x_i-x_im1)
rdy=(y_p-y_jm1)/(y_j-y_jm1)
!
!method 1: assuming zero derivatives at grid points
!ax = 3.d0*rdx**2 - 2.d0*rdx**3
!ay = 3.d0*rdy**2 - 2.d0*rdy**3
!
!method 2: assuming constant function towards 'ghost points' at i-2 and i+1
ax = 0.5d0*rdx + 1.5d0*rdx**2 - rdx**3
ay = 0.5d0*rdy + 1.5d0*rdy**2 - rdy**3
!
!
interpol2d_4p_cube = (1.d0-ax)*(1.d0-ay)*f_im1jm1 + ax*(1.d0-ay)*f_ijm1 + (1.d0-ax)*ay*f_im1j + ax*ay*f_ij
!
end function interpol2d_4p_cube
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
subroutine get_xy_indx(xp, yp, xarr, yarr, ndxmax, ndymax, iim1_x, ii_x, iim1_y, ii_y, expol, rmin, rmax)
!-----------------------------------------------------------------------
!
!  finds indices of grid for a square urrounding the point (xpn,yin,zin)
!
!  input:  coordinates of point:             xp, yp
!          dimension of grid:                ndxmax, ndymax
!          x,y-grid:                         xarr, yarr
!          boundaries of info-region:        rmin, rmax
!  output: indices of x-grid:                iim1_x, ii_x
!          indices of y-grid:                iim1_y, ii_y
!          flag if extrapolation is needed:  expol

use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax
real(dp), dimension(ndxmax), intent(in) :: xarr
real(dp), dimension(ndymax), intent(in) :: yarr
real(dp), intent(in) :: xp, yp
real(dp), intent(in) :: rmin, rmax
integer(i4b), intent(out) :: iim1_x, ii_x, iim1_y, ii_y
logical, intent(out) :: expol
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: alpha, beta
integer(i4b) :: startx, endx, starty, endy
real(dp) :: x_dum
logical :: linfo, linfo2, linfo_boundary
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(xp.ge.0.d0) then
   startx=ndxmax-1
   endx=1
   alpha=-1
else
   startx=2
   endx=ndxmax
   alpha=1
endif
!
if(yp.ge.0.d0) then
   starty=ndymax-1
   endy=1
   beta=-1
else
   starty=2
   endy=ndymax
   beta=1
endif
!
ii_x=startx
iim1_x=startx-alpha
do i=startx, endx, alpha
   if(alpha*xarr(i).ge.alpha*xp) then
      ii_x=i
      iim1_x=i-alpha
      exit
   endif
enddo
!
ii_y=starty
iim1_y=starty-beta
do i=starty, endy, beta
   if(beta*yarr(i).ge.beta*yp) then
      ii_y=i
      iim1_y=i-beta
      exit
   endif
enddo
!
!--------------extrapolation for grid-points near photosphere-----------
!
!check if extrapolation is needed at inner part of the star
call info_region(xarr(ii_x), yarr(ii_y), 0.d0, rmin, rmax, linfo, linfo2, linfo_boundary)
!
if(.not.linfo) then
   expol=.true.
   do i=1, 10
      ii_x=ii_x-alpha
      iim1_x=iim1_x-alpha
      ii_y=ii_y-beta
      iim1_y=iim1_y-beta
      call info_region(xarr(ii_x), yarr(ii_y), 0.d0, rmin, rmax, linfo, linfo2, linfo_boundary)
      if(linfo) exit
   enddo
   if(.not.linfo) stop 'error in get_xy_indx: linfo_phot eq false => extrapolation over more than 10 grid points'
   return
else
   expol=.false.
endif
!
!-------------extrapolation for grid points larger than rmax------------
!
!check if extrapolation is needed at outer part of the star
!call info_region(xarr(indx_x2), yarr(indx_y2), zarr(indx_z2), rmin, rmax, linfo, linfo2, linfo_boundary)
!
!if(.not.linfo2) then
 !  if(xin.eq.0.d0.and.yin.eq.0.d0) then
 !     write(*,*) 'point is on z-axis => no extrapolation needed'
 !     expol=.false.
 !  else if(xin.eq.0.d0.and.zin.eq.0.d0) then
 !     write(*,*) 'point is on y-axis => no extrapolation needed'
 !     expol=.false.
 !  else if(yin.eq.0.d0.and.zin.eq.0.d0) then
 !     write(*,*) 'point is on x-axis => no extrapolation needed'
 !     expol=.false.
 !  else
 !     expol=.true.
 !     do i=1, 10
 !        indx_x1=indx_x1+alpha
 !        indx_x2=indx_x2+alpha
 !        indx_y1=indx_y1+beta
 !        indx_y2=indx_y2+beta
 !        indx_z1=indx_z1+gamma
 !        indx_z2=indx_z2+gamma
 !        call info_region(xarr(indx_x2), yarr(indx_y2), zarr(indx_z2), rmin, rmax, linfo, linfo2, linfo_boundary)
 !        if(linfo2) exit
 !     enddo
 !     if(.not.linfo2) stop 'error in get_xyz_indx: linfo_max eq false => extrapolation over more than 10 grid points'
 !  endif
!else
!   expol=.false.
!endif
!
!
!
end subroutine get_xy_indx
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_xy_values1(x_coord, y_coord, xarr, yarr, ndxmax, ndymax, &
                          indx_x1, indx_x2, indx_y1, indx_y2, &
                          x1, x2, y1, y2, llogx, llogy)
!
!-----------------------------------------------------------------------
!
!   get coordinates of vertices of a plane,
!   with vertices given by indx_x1, ... indx_y2
!
!input: dimension of arrays x,y: ndxmax, ndymax
!       arrays x,y
!       coordinates of a point inside cube: x_coord, y_coord
!       indices of cube-vertices: indx_x1, ... indx_y2
!
!output: grid value: x1, x2, y1, y2
!        radii of vertices and of point: rada, ... radh, radp
!        flags to decide if logarithmic interpolation is allowed:
!               llogx, llogy
!
!-----------------------------------------------------------------------
!

use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax
integer(i4b), intent(in) :: indx_x1, indx_x2, indx_y1, indx_y2
real(dp), intent(in) :: x_coord, y_coord
real(dp), dimension(ndxmax), intent(in) :: xarr
real(dp), dimension(ndymax), intent(in) :: yarr
real(dp) :: x1, x2, y1, y2
logical :: llogx, llogy
!
! ... local scalars
!
! ... local functions
!
x1=xarr(indx_x1)
x2=xarr(indx_x2)
y1=yarr(indx_y1)
y2=yarr(indx_y2)
!
!all coordinates need the same sign
if(x1*x2.le.0.d0) then
   llogx=.false.
else
   if(x1*x_coord.le.0.d0) then
      llogx=.false.
   else
      llogx=.true.
   endif
endif
!
if(y1*y2.le.0.d0) then
   llogy=.false.
else
   if(y1*y_coord.le.0.d0) then
      llogy=.false.
   else
      llogy=.true.
   endif
endif
!
!
!
end subroutine get_xy_values1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_xy_values2(ndxmax, ndymax, yvalue2d, &
                          indx_x1, indx_x2, indx_y1, indx_y2, &
                          vala, valb, valc, vald, llogf)
!
!-----------------------------------------------------------------------
!
!   get physical quantities of a plane with vertices
!         given by indx_x1, ... indx_y2
!
!input: dimension of 2-d array: ndxmax, ndymax
!       physical value on grid: yvalue2d
!       indices of cube-vertices: indx_x1, ... indx_y2
!
!output: physical value on vertices: vala ... vald
!        flag to decide if logarithmic interpolation is allowed:
!               llogf
!
!examples:
!
!  in carthesian:     in spherical:
!
!                           c
!   c------d               / \
!   |      |              /   \
!   |      |             a     \
!   a------b            / \     \
!                      /   \     \
!                     /-----b-----d
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax
integer(i4b), intent(in) :: indx_x1, indx_x2, indx_y1, indx_y2
real(dp), dimension(ndxmax, ndymax), intent(in) :: yvalue2d
real(dp) :: vala, valb, valc, vald
logical :: llogf
!
! ... local scalars
!
! ... local functions
!
!-----------------------------------------------------------------------
!
vala = yvalue2d(indx_x1, indx_y1)
valb = yvalue2d(indx_x1, indx_y2)
valc = yvalue2d(indx_x2, indx_y1)
vald = yvalue2d(indx_x2, indx_y2)
!
!-------------check if logarithmic interpolation is allowed-------------
!
if(vala.le.0.d0.or.valb.le.0.d0.or.valc.le.0.d0.or.vald.le.0.d0) then
   llogf=.false.
else
   llogf=.true.
endif
!
end subroutine get_xy_values2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine bilin(xout, yout, &
                 x1, x2, y1, y2, &
                 vala, valb, valc, vald, &
                 llogx, llogy, llogf, yinterp)
!
use prog_type
use mod_interp1d
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xout, yout, x1, x2, y1, y2, vala, valb, valc, vald
real(dp) :: yinterp
logical, intent(in) :: llogx, llogy, llogf
!
! ... local scalars
real(dp) :: dum_vala, dum_valb, dum_valc, dum_vald
real(dp) :: dum_x1, dum_x2, dum_xout, dum_y1, dum_y2, dum_yout
real(dp) :: yvalue_n, yvalue_s
!
! ... local functions
!
!
!
if(llogf) then
!
!prepare input values for log-* interpolation
   dum_vala=log10(vala)
   dum_valb=log10(valb)
   dum_valc=log10(valc)
   dum_vald=log10(vald)
!
   if(llogx) then
!log-log interpolation in x-direction
      dum_x1=log10(abs(x1))
      dum_x2=log10(abs(x2))
      dum_xout=log10(abs(xout))
    else
!log-lin interpolation in x-direction
      dum_x1=x1
      dum_x2=x2
      dum_xout=xout
   endif
!
   if(llogy) then
!log-log-interpolation in y-direction
      dum_y1=log10(abs(y1))
      dum_y2=log10(abs(y2))
      dum_yout=log10(abs(yout))
   else
!log-lin-interpolation in y-direction
      dum_y1=y1
      dum_y2=y2
      dum_yout=yout
   endif

!perform interpolation in x-direction
   yvalue_n = interpol_yp(dum_x2, dum_x1, dum_valc, dum_vala, dum_xout)
   yvalue_s = interpol_yp(dum_x2, dum_x1, dum_vald, dum_valb, dum_xout)
!
!perform interpolation in y-direction
   yinterp = interpol_yp(dum_y2, dum_y1, yvalue_s, yvalue_n, dum_yout)
   yinterp = 10.d0**yinterp
!
else
!lin-lin-interpolation
!
   yvalue_n = interpol_yp(x2, x1, valc, vala, xout)
   yvalue_s = interpol_yp(x2, x1, vald, valb, xout)
!
   yinterp = interpol_yp(y2, y1, yvalue_s, yvalue_n, yout)
!
endif
!
!
end subroutine bilin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine bilin2(xp, yp, x1, x2, y1, y2, vala, valb, valc, vald, valp)
!
! ^
! |
! y2   valc---------vald
! |     |            |
! |     |            |
! |     |            |
! y1   vala---------valb
! |
! 0-----x1-----------x2--->
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xp, yp, x1, x2, y1, y2, vala, valb, valc, vald
real(dp) :: valp
!
! ... local scalars
real(dp) :: delx, dely, delxp, delyp
!
! ... local functions
!
!
delx=x2-x1
dely=y2-y1
!
delxp=xp-x2
delyp=yp-y2
!
valp = vald + (vald-valc)*delxp/delx + (vald-valb)*delyp/dely + &
       (vald-valc-valb+vala)*delxp*delyp/delx/dely

!write(*,*) vald, (vald-valc)/delx * (xp-x2), (vald-valb)/dely * (yp-y2),        (vald-valc-valb+vala)/delx/dely * (xp-x2)
!stop
!
!
end subroutine bilin2
