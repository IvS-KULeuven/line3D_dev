module mod_interp3d
!
use prog_type
use fund_const
use mod_interp1d, only: interpol_yp
!
implicit none
!
contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol3d_8p_lin(f_im1jm1km1, f_ijm1km1, f_im1jkm1, f_ijkm1, &
                           f_im1jm1k, f_ijm1k, f_im1jk, f_ijk, &
                           x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p)
!
!         interpolates values given on a 3d grid onto point x_p, y_p, z_p
!                      using triliniear interpolation
!
!   f(x,y,z) = f_km1 + (f_k-f_km1)*(z_p-z_km1)/(z_k-z_km1)
!
!      with:   f_km1 = f_jm1km1 + (f_jkm1-f_jm1km1)*(y_p-y_im1)/(y_j-y_jm1)
!              f_k   = f_jm1    + (f_jk  -f_jm1k  )*(y_p-y_im1)/(y_j-y_jm1)
!              f_jk     = f_im1jkm1   + (f_ijkm1-f_im1jkm1) *   (x_p-x_im1)/(x_i-x_im1)
!              f_jm1km1 = f_im1jm1km1 + (f_ijm1km1-f_im1jm1km1)*(x_p-x_im1)/(x_i-x_im1)
!
!on input: 
!
!          f_im1jk-------------f_ijk          z_k
!            /|                 /|             |  y_j 
!           / |                / |             |   /
!          /  |               /  |             |  /
!         /   |              /   |             | /
!        /    |             /    |             |/
!       / f_im1jm1---------/--f_ijkm1        x_im1--------------x_i
!      /     /  x(x,y,z)  /     /            y_jm1
! f_im1jm1k-----------f_ijm1k  /             z_km1
!     |    /             |    /
!     |   /              |   /
!     |  /               |  /
!     | /                | /
!     |/                 |/
!f_im1jm1km1---------f_ijm1km1
!                    
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p, z_p
!

!

!
! ... argments
real(dp), intent(in) :: f_im1jm1km1, f_ijm1km1, f_im1jkm1, f_ijkm1, &
                        f_im1jm1k, f_ijm1k, f_im1jk, f_ijk, &
                        x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p
real(dp) :: interpol3d_8p_lin
!
! ... local scalars
real(dp) :: dxi, dx, dyj, dy, dzk, dz, rdx, rdy, rdz, rdxdy
real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
real(dp) :: fdum1, fdum2, fdum3, fdum4
!
!define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
dzk=z_k-z_km1
dz=z_p-z_km1
!
!define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdz=dz/dzk
rdxdy=rdx*rdy
!
!trilinear interpolation
fdum1=1.d0-rdx-rdy+rdxdy
fdum2=rdx-rdxdy
fdum3=rdy-rdxdy
fdum4=rdxdy
!
ecoeff=fdum1*rdz
fcoeff=fdum2*rdz
gcoeff=fdum3*rdz
hcoeff=fdum4*rdz
acoeff=fdum1-ecoeff
bcoeff=fdum2-fcoeff
ccoeff=fdum3-gcoeff
dcoeff=fdum4-hcoeff
!
interpol3d_8p_lin = acoeff*f_im1jm1km1 + bcoeff*f_ijm1km1 + ccoeff*f_im1jkm1 + dcoeff*f_ijkm1 + &
                    ecoeff*f_im1jm1k   + fcoeff*f_ijm1k   + gcoeff*f_im1jk   + hcoeff*f_ijk
!
!
end function interpol3d_8p_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_8p_lin(x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p, &
                          acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff)
!
!         interpolates values given on a 3d grid onto point x_p, y_p, z_p
!                      using triliniear interpolation
!
!         returns coefficients, such that
!            f(x,y,z) = a*f_im1jm1km1 + b*f_ijm1km1 + c*f_im1jkm1 + d*f_ijkm1 + 
!                       e*f_im1jm1j   + f*f_ijm1k   + g*f_im1jk   + h*f_ijk
!
!   f(x,y,z) = f_km1 + (f_k-f_km1)*(z_p-z_km1)/(z_k-z_km1)
!
!      with:   f_km1 = f_jm1km1 + (f_jkm1-f_jm1km1)*(y_p-y_im1)/(y_j-y_jm1)
!              f_k   = f_jm1    + (f_jk  -f_jm1k  )*(y_p-y_im1)/(y_j-y_jm1)
!              f_jk     = f_im1jkm1   + (f_ijkm1-f_im1jkm1) *   (x_p-x_im1)/(x_i-x_im1)
!              f_jm1km1 = f_im1jm1km1 + (f_ijm1km1-f_im1jm1km1)*(x_p-x_im1)/(x_i-x_im1)
!
!on input: 
!
!          f_im1jk-------------f_ijk          z_k
!            /|                 /|             |  y_j 
!           / |                / |             |   /
!          /  |               /  |             |  /
!         /   |              /   |             | /
!        /    |             /    |             |/
!       / f_im1jm1---------/--f_ijkm1        x_im1--------------x_i
!      /     /  x(x,y,z)  /     /            y_jm1
! f_im1jm1k-----------f_ijm1k  /             z_km1
!     |    /             |    /
!     |   /              |   /
!     |  /               |  /
!     | /                | /
!     |/                 |/
!f_im1jm1km1---------f_ijm1km1
!                    
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p, z_p
!

!

!
! ... argments
real(dp), intent(in) :: x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p
real(dp), intent(out) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
real(dp) :: interpol3d_8p_lin
!
! ... local scalars
real(dp) :: dxi, dx, dyj, dy, dzk, dz, rdx, rdy, rdz, rdxdy
real(dp) :: fdum1, fdum2, fdum3, fdum4
!
!define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
dzk=z_k-z_km1
dz=z_p-z_km1
!
!define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdz=dz/dzk
rdxdy=rdx*rdy
!
!trilinear interpolation
fdum1=1.d0-rdx-rdy+rdxdy
fdum2=rdx-rdxdy
fdum3=rdy-rdxdy
fdum4=rdxdy
!
ecoeff=fdum1*rdz
fcoeff=fdum2*rdz
gcoeff=fdum3*rdz
hcoeff=fdum4*rdz
acoeff=fdum1-ecoeff
bcoeff=fdum2-fcoeff
ccoeff=fdum3-gcoeff
dcoeff=fdum4-hcoeff
!
!
end subroutine coeff3d_8p_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_27p_quad2(x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, &
                             z_km2, z_km1, z_k, x_p, y_p, z_p, &
                             c01, c02, c03, c04, c05, c06, c07, c08, c09, &
                             c10, c11, c12, c13, c14, c15, c16, c17, c18, &
                             c19, c20, c21, c22, c23, c24, c25, c26, c27)
!
!         interpolates values given on a 3d grid onto point x_p, y_p, z_p
!                 using 3d quadratic functions
!
!   f(x,y) = sum_i sum_j sum_k (a_ijk * (x-x_i)^i * (y-y_j)^j * (z-z_k)^k) = 
!          = c01*f_im2jm2km2 + c02*f_im1jm2km2 + c03*f_ijm2km2 + 
!            c04*f_im2jm1km2 + c05*f_im1jm1km2 + c06*f_ijm1km2 + 
!            c07*f_im2jkm2   + c08*f_im1jkm2   + c09*f_ijkm2 + 
!
!            c10*f_im2jm2km1 + c11*f_im1jm2km1 + c12*f_ijm2km1 + 
!            c13*f_im2jm1km1 + c14*f_im1jm1km1 + c15*f_ijm1km1 + 
!            c16*f_im2jkm1   + c17*f_im1jkm1   + c18*f_ijkm1 + 
!
!            c19*f_im2jm2k + c20*f_im1jm2k + c21*f_ijm2k + 
!            c22*f_im2jm1k + c23*f_im1jm1k + c24*f_ijm1k + 
!            c25*f_im2jk   + c26*f_im1jk   + c27*f_ijk
!
!on input: 
!
!                 f_im2jk---------------f_im1jk--------------f_ijk
!                    /|                    /|                 /|
!                   / |                   / |                / |
!                  /  |                  /  |               /  |
!                 /   |                 /   |              /   |
!                /    |                /    |             /    |
!               / f_im2jkm1-----------/-f_im1jkm1--------/--f_ijkm1
!              /     /|              /     /|           /     /|     
!       f_im2jm1k--------------f_im1jm1k-----------f_ijm1k   / |     
!            /|    /  |            /|    /  |         /|    /  |     
!           / |   /   |           / |   /   |        / |   /   |     
!          /  |  /    |          /  |  /    |       /  |  /    |     
!         /   | / f_im2jkm2-----/---|-/-f_im1jkm2--/---|-/--f_ijkm2
!        /    |/     /         /    |/     /      /    |/     /      
!       /f_im2jm1km1----------/-f_im1jm1km1------/--f_ijm1km1/
!      /     /|    /         /     /|    /      /     /|    /
!f_im2jm2k--------------f_im1jm2k------------f_ijm2k / |   / 
!     |    /  |  /          |    /  |  /       |    /  |  /  
!     |   /   | /           |   /   | /        |   /   | /   
!     |  /    |/            |  /    |/         |  /    |/
!     | / f_im2jm1km2-------|-/-f_im1jm1km2----|-/--f_ijm1km2
!     |/     /              |/     /           |/     /      
!f_im2jm2km1------------f_im1jm2km1---------f_ijm2km1/       
!     |    /                |    /             |    /
!     |   /                 |   /              |   /
!     |  /                  |  /               |  /
!     | /                   | /                | /
!     |/                    |/                 |/
!f_im2jm2km2----------f_im1jm2km2---------f_ijm2km2
!                    
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output: interpolated value at x_p, y_p, z_p
!

!

!
! ... argments
real(dp), intent(in) :: x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, &
                        z_km2, z_km1, z_k, x_p, y_p, z_p
real(dp), intent(out) :: c01, c02, c03, c04, c05, c06, c07, c08, c09, &
                         c10, c11, c12, c13, c14, c15, c16, c17, c18, &
                         c19, c20, c21, c22, c23, c24, c25, c26, c27
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, dyjm1, dyj, dy, dzkm1, dzk, dz
real(dp) :: tx, ty, tz, ax, bx, cx, ay, by, cy, az, bz, cz
!
!define deltax, deltay, deltaz
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
tx = x_p-x_im1
!
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
ty = y_p-y_jm1
!
dzkm1 = z_km1-z_km2
dzk = z_k-z_km1
dz = dzkm1+dzk
tz = z_p-z_km1
!
ax = (tx**2 - tx*dxi)/dxim1/dx
bx = 1.d0 - (tx**2+tx*(dxim1-dxi))/dxi/dxim1
cx = (tx**2 + tx*dxim1)/dxi/dx
!
ay = (ty**2 - ty*dyj)/dyjm1/dy
by = 1.d0 - (ty**2+ty*(dyjm1-dyj))/dyj/dyjm1
cy = (ty**2 + ty*dyjm1)/dyj/dy
!
az = (tz**2 - tz*dzk)/dzkm1/dz
bz = 1.d0 - (tz**2+tz*(dzkm1-dzk))/dzk/dzkm1
cz = (tz**2 + tz*dzkm1)/dzk/dz
!
!
c01=ax*ay*az
c02=bx*ay*az
c03=cx*ay*az
c04=ax*by*az
c05=bx*by*az
c06=cx*by*az
c07=ax*cy*az
c08=bx*cy*az
c09=cx*cy*az
!
c10=ax*ay*bz
c11=bx*ay*bz
c12=cx*ay*bz
c13=ax*by*bz
c14=bx*by*bz
c15=cx*by*bz
c16=ax*cy*bz
c17=bx*cy*bz
c18=cx*cy*bz
!
c19=ax*ay*cz
c20=bx*ay*cz
c21=cx*ay*cz
c22=ax*by*cz
c23=bx*by*cz
c24=cx*by*cz
c25=ax*cy*cz
c26=bx*cy*cz
c27=cx*cy*cz
!
!
end subroutine coeff3d_27p_quad2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function interpol3d_8p_idw(fa, fb, fc, fd, fe, ff, fg, fh, &
                           xa, ya, za, xb, yb, zb, xc, yc, zc, &
                           xd, yd, zd, xe, ye, ze, xf, yf, zf, &
                           xg, yg, zg, xh, yh, zh, xp, yp, zp)
!
!         interpolates values given on a 3d grid onto point xp, yp, zp
!                      using inverse distance weighting
!
!on input: 
!
!            fg-----------------fh
!            /|                 /|
!           / |                / |
!          /  |               /  |
!         /   |              /   |
!        /    |             /    |
!       /    fc------------/----fd
!      /     /  x(x,y,z)  /     / 
!    fe------------------ff    /  
!     |    /             |    /
!     |   /              |   /
!     |  /               |  /
!     | /                | /
!     |/                 |/
!    fa-----------------fb
!                    
!        xp, yp, zp: coordinates of point onto which shall be interpolated
!        xa, ya, za: coordinates of point A, B, ...
!
!on output: interpolated value at x_p, y_p, z_p
!


!

!
! ... argments
real(dp), intent(in) :: fa, fb, fc, fd, fe, ff, fg, fh, &
                        xa, xb, xc, xd, xe, xf, xg, xh, &
                        ya, yb, yc, yd, ye, yf, yg, yh, &
                        za, zb, zc, zd, ze, zf, zg, zh, &
                        xp, yp, zp
real(dp) :: interpol3d_8p_idw
!
! ... local scalars
real(dp) :: fac, norm
real(dp) :: da, db, dc, dd, de, df, dg, dh, &
            wa, wb, wc, wd, we, wf, wg, wh
real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
!
!
!define weighting factor
fac=two
!
!calculate distance to each point
da = (xp-xa)**2+(yp-ya)**2+(zp-za)**2
db = (xp-xb)**2+(yp-yb)**2+(zp-zb)**2
dc = (xp-xc)**2+(yp-yc)**2+(zp-zc)**2
dd = (xp-xd)**2+(yp-yd)**2+(zp-zd)**2
de = (xp-xe)**2+(yp-ye)**2+(zp-ze)**2
df = (xp-xf)**2+(yp-yf)**2+(zp-zf)**2
dg = (xp-xg)**2+(yp-yg)**2+(zp-zg)**2
dh = (xp-xh)**2+(yp-yh)**2+(zp-zh)**2
!
wa=zero
wb=zero
wc=zero
wd=zero
we=zero
wf=zero
wg=zero
wh=zero
!   
!avoid division by zero if da=0, ...
if(da.eq.zero) then
   wa=one
elseif(db.eq.zero) then
   wb=one
elseif(dc.eq.zero) then
   wc=one
elseif(dd.eq.zero) then
   wd=one
elseif(de.eq.zero) then
   we=one
elseif(df.eq.zero) then
   wf=one
elseif(dg.eq.zero) then
   wg=one
elseif(dh.eq.zero) then
   wh=one
else
!
   wa=1.d0/da**(fac/two)
   wb=1.d0/db**(fac/two)
   wc=1.d0/dc**(fac/two)
   wd=1.d0/dd**(fac/two)
   we=1.d0/de**(fac/two)
   wf=1.d0/df**(fac/two)
   wg=1.d0/dg**(fac/two)
   wh=1.d0/dh**(fac/two)      
!
endif
!
norm = wa + wb + wc + wd + we + wf + wg + wh
!
acoeff = wa/norm
bcoeff = wb/norm
ccoeff = wc/norm
dcoeff = wd/norm
ecoeff = we/norm
fcoeff = wf/norm
gcoeff = wg/norm
hcoeff = wh/norm
!
interpol3d_8p_idw = acoeff*fa + bcoeff*fb + ccoeff*fc + dcoeff*fd + &
                    ecoeff*fe + fcoeff*ff + gcoeff*fg + hcoeff*fh
!
end function interpol3d_8p_idw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_8p_idw(xa, ya, za, xb, yb, zb, xc, yc, zc, &
                          xd, yd, zd, xe, ye, ze, xf, yf, zf, &
                          xg, yg, zg, xh, yh, zh, xp, yp, zp, &
                          acoeff, bcoeff, ccoeff, dcoeff, &
                          ecoeff, fcoeff, gcoeff, hcoeff, fac)
!
!         interpolates values given on a 3d grid onto point xp, yp, zp
!                      using inverse distance weighting
!
!on input: 
!
!            fg-----------------fh
!            /|                 /|
!           / |                / |
!          /  |               /  |
!         /   |              /   |
!        /    |             /    |
!       /    fc------------/----fd
!      /     /  x(x,y,z)  /     / 
!    fe------------------ff    /  
!     |    /             |    /
!     |   /              |   /
!     |  /               |  /
!     | /                | /
!     |/                 |/
!    fa-----------------fb
!                    
!        xp, yp, zp: coordinates of point onto which shall be interpolated
!        xa, ya, za: coordinates of point A, B, ...
!
!on output: interpolated value at x_p, y_p, z_p
!


!

!
! ... argments
real(dp), intent(in) :: xa, xb, xc, xd, xe, xf, xg, xh, &
                        ya, yb, yc, yd, ye, yf, yg, yh, &
                        za, zb, zc, zd, ze, zf, zg, zh, &
                        xp, yp, zp
real(dp), intent(in), optional :: fac
real(dp), intent(out) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
!
! ... local scalars
real(dp) :: afac, norm
real(dp) :: da, db, dc, dd, de, df, dg, dh, &
            wa, wb, wc, wd, we, wf, wg, wh
!
!
!define weighting factor
if(.not.present(fac)) then
   afac=two
else
   afac=fac
endif
!
!calculate distance to each point
da = (xp-xa)**2+(yp-ya)**2+(zp-za)**2
db = (xp-xb)**2+(yp-yb)**2+(zp-zb)**2
dc = (xp-xc)**2+(yp-yc)**2+(zp-zc)**2
dd = (xp-xd)**2+(yp-yd)**2+(zp-zd)**2
de = (xp-xe)**2+(yp-ye)**2+(zp-ze)**2
df = (xp-xf)**2+(yp-yf)**2+(zp-zf)**2
dg = (xp-xg)**2+(yp-yg)**2+(zp-zg)**2
dh = (xp-xh)**2+(yp-yh)**2+(zp-zh)**2
!
wa=zero
wb=zero
wc=zero
wd=zero
we=zero
wf=zero
wg=zero
wh=zero
!   
!avoid division by zero if da=0, ...
if(da.eq.zero) then
   wa=one
elseif(db.eq.zero) then
   wb=one
elseif(dc.eq.zero) then
   wc=one
elseif(dd.eq.zero) then
   wd=one
elseif(de.eq.zero) then
   we=one
elseif(df.eq.zero) then
   wf=one
elseif(dg.eq.zero) then
   wg=one
elseif(dh.eq.zero) then
   wh=one
else
!
   wa=1.d0/da**(afac/two)
   wb=1.d0/db**(afac/two)
   wc=1.d0/dc**(afac/two)
   wd=1.d0/dd**(afac/two)
   we=1.d0/de**(afac/two)
   wf=1.d0/df**(afac/two)
   wg=1.d0/dg**(afac/two)
   wh=1.d0/dh**(afac/two)      
!
endif
!
norm = wa + wb + wc + wd + we + wf + wg + wh
!
acoeff = wa/norm
bcoeff = wb/norm
ccoeff = wc/norm
dcoeff = wd/norm
ecoeff = we/norm
fcoeff = wf/norm
gcoeff = wg/norm
hcoeff = wh/norm

!write(*,*) sqrt(da), sqrt(db), sqrt(dc), sqrt(dd), sqrt(de), sqrt(df), sqrt(dg), sqrt(dh)
!
end subroutine coeff3d_8p_idw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_xyz_indx(xin, yin, zin, xarr, yarr, zarr, ndxmax, ndymax, ndzmax, indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol, rmin, rmax)
!
!-----------------------------------------------------------------------
!
!  finds indices of grid for a cube surrounding the point (xin,yin,zin)
!
!  input:  coordinates of point:             xin, yin, zin
!          dimension of grid:                ndxmax, ndymax, ndzmax
!          x,y,z-grid:                       xarr, yarr, zarr
!          boundaries of info-region:        rmin, rmax
!  output: indices of x-grid:                indx_x1, indx_x2
!          indices of y-grid:                indx_y1, indx_y2
!          indices of z-grid:                indx_z1, indx_z2
!          flag if extrapolation is needed:  expol


!

!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
real(dp), dimension(ndxmax), intent(in) :: xarr
real(dp), dimension(ndymax), intent(in) :: yarr
real(dp), dimension(ndzmax), intent(in) :: zarr
real(dp), intent(in) :: xin, yin, zin
real(dp), intent(in) :: rmin, rmax
integer(i4b) :: indx_x1, indx_y1, indx_z1, indx_x2, indx_y2, indx_z2
logical :: expol
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: alpha, beta, gamma
integer(i4b) :: startx, endx, starty, endy, startz, endz
real(dp) :: x_dum
logical :: linfo, linfo2, linfo_boundary
!
! ... local functions
!
!-----------------------------------------------------------------------
!
if(xin.ge.0.d0) then
   startx=ndxmax-1
   endx=1
   alpha=-1
else
   startx=2
   endx=ndxmax
   alpha=1
endif
!
if(yin.ge.0.d0) then
   starty=ndymax-1
   endy=1
   beta=-1
else
   starty=2
   endy=ndymax
   beta=1
endif
!
if(zin.ge.0.d0) then
   startz=ndzmax-1
   endz=1
   gamma=-1
else
   startz=2
   endz=ndzmax
   gamma=1
endif
!
indx_x1=startx
indx_x2=startx-alpha
do i=startx, endx, alpha
   if(alpha*xarr(i).ge.alpha*xin) then
      indx_x1=i
      indx_x2=i-alpha
      exit
   endif
enddo
!
indx_y1=starty
indx_y2=starty-beta
do i=starty, endy, beta
   if(beta*yarr(i).ge.beta*yin) then
      indx_y1=i
      indx_y2=i-beta
      exit
   endif
enddo
!
indx_z1=startz
indx_z2=startz-gamma
do i=startz, endz, gamma
   if(gamma*zarr(i).ge.gamma*zin) then
      indx_z1=i
      indx_z2=i-gamma
      exit
   endif
enddo
!
!--------------extrapolation for grid-points near photosphere-----------
!
!check if extrapolation is needed at inner part of the star
call info_region(xarr(indx_x1), yarr(indx_y1), zarr(indx_z1), rmin, rmax, linfo, linfo2, linfo_boundary)
!
if(.not.linfo) then
   expol=.true.
   do i=1, 10
      indx_x1=indx_x1-alpha
      indx_x2=indx_x2-alpha
      indx_y1=indx_y1-beta
      indx_y2=indx_y2-beta
      indx_z1=indx_z1-gamma
      indx_z2=indx_z2-gamma
      call info_region(xarr(indx_x1), yarr(indx_y1), zarr(indx_z1), rmin, rmax, linfo, linfo2, linfo_boundary)
      if(linfo) exit
   enddo
   if(.not.linfo) stop 'error in get_xyz_indx: linfo_phot eq false => extrapolation over more than 10 grid points'
   return
else
   expol=.false.
endif
!
!-------------extrapolation for grid points larger than rmax------------
!
!check if extrapolation is needed at outer part of the star
call info_region(xarr(indx_x2), yarr(indx_y2), zarr(indx_z2), rmin, rmax, linfo, linfo2, linfo_boundary)
!
if(.not.linfo2) then
   if(xin.eq.0.d0.and.yin.eq.0.d0) then
      write(*,*) 'point is on z-axis => no extrapolation needed'
      expol=.false.
   else if(xin.eq.0.d0.and.zin.eq.0.d0) then
      write(*,*) 'point is on y-axis => no extrapolation needed'
      expol=.false.
   else if(yin.eq.0.d0.and.zin.eq.0.d0) then
      write(*,*) 'point is on x-axis => no extrapolation needed'
      expol=.false.
   else
      expol=.true.
      do i=1, 10
         indx_x1=indx_x1+alpha
         indx_x2=indx_x2+alpha
         indx_y1=indx_y1+beta
         indx_y2=indx_y2+beta
         indx_z1=indx_z1+gamma
         indx_z2=indx_z2+gamma
         call info_region(xarr(indx_x2), yarr(indx_y2), zarr(indx_z2), rmin, rmax, linfo, linfo2, linfo_boundary)
         if(linfo2) exit
      enddo
      if(.not.linfo2) stop 'error in get_xyz_indx: linfo_max eq false => extrapolation over more than 10 grid points'
   endif
else
   expol=.false.
endif
!
!
!
end subroutine get_xyz_indx

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_xyz_values1(x_coord, y_coord, z_coord, xarr, yarr, zarr, ndxmax, ndymax, ndzmax, &
                          indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                          x1, x2, y1, y2, z1, z2, &
                          rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                          llogx, llogy, llogz, expol)
!
!-----------------------------------------------------------------------
!
!   get coordinates and radii of a cube
!   with vertices given by indx_x1, ... indx_z2
!
!input: dimension of arrays x,y,z: ndxmax, ndymax, ndzmax
!       arrays x,y,z
!       coordinates of a point inside cube: x_coord, y_coord, z_coord
!       indices of cube-vertices: indx_x1, ... indx_z2
!       logical to decide if interpolation or extrapolation: expol
!
!output: grid value: x1, x2, y1, y2, z1, z2
!        radii of vertices and of point: rada, ... radh, radp
!        flags to decide if logarithmic interpolation is allowed:
!               llogx, llogy, llogz
!
!-----------------------------------------------------------------------
!


!

!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
integer(i4b), intent(in) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2
real(dp), intent(in) :: x_coord, y_coord, z_coord
real(dp), dimension(ndxmax), intent(in) :: xarr
real(dp), dimension(ndymax), intent(in) :: yarr
real(dp), dimension(ndzmax), intent(in) :: zarr
real(dp) :: x1, x2, y1, y2, z1, z2
real(dp) :: rada, radb, radc, radd, rade, radf, radg, radh, radp
logical :: llogx, llogy, llogz
logical, intent(in) :: expol
!
! ... local scalars
!
! ... local functions
real(dp) :: getr
!
!-----------------------------------------------------------------------
!
x1=xarr(indx_x1)
x2=xarr(indx_x2)
y1=yarr(indx_y1)
y2=yarr(indx_y2)
z1=zarr(indx_z1)
z2=zarr(indx_z2)
!
rada = getr(x1, y1, z1)
radb = getr(x2, y1, z1)
radc = getr(x1, y1, z2)
radd = getr(x2, y1, z2)
rade = getr(x1, y2, z1)
radf = getr(x2, y2, z1)
radg = getr(x1, y2, z2)
radh = getr(x2, y2, z2)
!
radp = getr(x_coord, y_coord, z_coord)
!
!-------------check if logarithmic interpolation is allowed-------------
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
if(z1*z2.le.0.d0) then
   llogz=.false.
else
   if(z1*z_coord.le.0.d0) then
      llogz=.false.
   else
      llogz=.true.
   endif
endif
!
!if extrapolation needs to be performed: logarithmic extrapolation is prohibited if
!   coordinates are less than 0.1 (otherwise, extrapolation from e.g. [1.d-2, 2.d-2] to 1.d-5
!                                  very inaccurate)
if(expol.and.(abs(x1).lt.0.1d0)) then
   llogx=.false.
endif
!
if(expol.and.(abs(y1).lt.0.1d0)) then
   llogy=.false.
endif
!
if(expol.and.(abs(z1).lt.0.1d0)) then
   llogz=.false.
endif
!
!
!
end subroutine get_xyz_values1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_xyz_values2(ndxmax, ndymax, ndzmax, yvalue3d, &
                           indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!
!-----------------------------------------------------------------------
!
!   get physical quantities of a cube 
!   with vertices given by indx_x1, ... indx_z2
!
!input: dimension of 3-d array: ndxmax, ndymax, ndzmax
!       physical value on grid: yvalue3d
!       indices of cube-vertices: indx_x1, ... indx_z2
!
!output: physical value on cube vertices: vala ... valh
!        flag to decide if logarithmic interpolation is allowed:
!               llogf
!
!-----------------------------------------------------------------------
!

!

!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
integer(i4b), intent(in) :: indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2
real(dp), dimension(ndxmax, ndymax, ndzmax), intent(in) :: yvalue3d
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh
logical :: llogf
!
! ... local scalars
!
! ... local functions
!
!-----------------------------------------------------------------------
!
vala = yvalue3d(indx_x1, indx_y1, indx_z1)
valb = yvalue3d(indx_x2, indx_y1, indx_z1)
valc = yvalue3d(indx_x1, indx_y1, indx_z2)
vald = yvalue3d(indx_x2, indx_y1, indx_z2)
vale = yvalue3d(indx_x1, indx_y2, indx_z1)
valf = yvalue3d(indx_x2, indx_y2, indx_z1)
valg = yvalue3d(indx_x1, indx_y2, indx_z2)
valh = yvalue3d(indx_x2, indx_y2, indx_z2)
!
if(vala.le.0.d0.or.valb.le.0.d0.or.valc.le.0.d0.or.vald.le.0.d0.or.vale.le.0.d0.or.valf.le.0.d0.or.valg.le.0.d0.or.valh.le.0.d0) then
   llogf=.false.
else
   llogf=.true.
endif
!
end subroutine get_xyz_values2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine trilin_complete(xout, yout, zout, &
                  x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  rada, radb, radc, radd, rade, radf, radg, radh, radp, &
                  expol, eexpol, llogx, llogy, llogz, llogf, lr2, yinterp)
!

!

!
! ... arguments
real(dp), intent(in) :: xout, yout, zout, x1, x2, y1, y2, z1, z2, &
                        rada, radb, radc, radd, rade, radf, radg, radh, radp
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, yinterp
logical, intent(in) :: expol, eexpol, llogx, llogy, llogz, llogf, lr2
!
!
if(lr2) then
!applying interpolation of function values * r^2
   vala=vala*rada*rada
   valb=valb*radb*radb
   valc=valc*radc*radc
   vald=vald*radd*radd
   vale=vale*rade*rade
   valf=valf*radf*radf
   valg=valg*radg*radg
   valh=valh*radh*radh
endif
!
!perform only interpolation (expol=.false.)
if(.not.expol) then
   call trilin(xout, yout, zout, &
               x1, x2, y1, y2, z1, z2, &
               vala, valb, valc, vald, vale, valf, valg, valh, &
               llogx, llogy, llogz, llogf, yinterp)
else
!set values to zero if extrapolation is needed (or perform extrapolation)
   if(eexpol) then
      call trilin(xout, yout, zout, &
                  x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  llogx, llogy, llogz, llogf, yinterp)
   else
      yinterp=0.d0
   endif
!
endif
!
if(lr2) yinterp=yinterp/radp/radp
!
!
!
end subroutine trilin_complete
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine trilin(xout, yout, zout, &
                  x1, x2, y1, y2, z1, z2, &
                  vala, valb, valc, vald, vale, valf, valg, valh, &
                  llogx, llogy, llogz, llogf, yinterp)        
!



!
! ... arguments
real(dp), intent(in) :: xout, yout, zout
real(dp), intent(in) :: vala, valb, valc, vald, vale, valf, valg, valh
real(dp), intent(in) :: x1, x2, y1, y2, z1, z2
real(dp) :: yinterp
logical, intent(in) :: llogx, llogy, llogz, llogf
!
! ... local scalars
real(dp) :: yvalue_s0, yvalue_s1, yvalue_n0, yvalue_n1, sols1, sols2, soln1, soln2
real(dp) :: yvalue_o0, yvalue_o1, yvalue_w0, yvalue_w1, solo1, solo2, solw1, solw2
real(dp) :: yvalue_t0, yvalue_t1, yvalue_b0, yvalue_b1, solt1, solt2, solb1, solb2
real(dp) :: yout1, yout2, yout3
real(dp) :: dum_vala, dum_valb, dum_valc, dum_vald, dum_vale, dum_valf, dum_valg, dum_valh
real(dp) :: dum_x1, dum_x2, dum_y1, dum_y2, dum_z1, dum_z2, dum_xout, dum_yout, dum_zout
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
!
!--------------------obtain values on all surfaces----------------------
!------------only to be convinced that solutions are all the same-------
!
!yvalue_s0 = interpol_yp(z2, z1, vald, valb, zout)
!yvalue_s1 = interpol_yp(z2, z1, valc, vala, zout)
!yvalue_n0 = interpol_yp(z2, z1, valh, valf, zout)
!yvalue_n1 = interpol_yp(z2, z1, valg, vale, zout)
!yvalue_o0 = interpol_yp(y2, y1, valh, vald, yout)
!yvalue_o1 = interpol_yp(y2, y1, valf, valb, yout)
!yvalue_w0 = interpol_yp(y2, y1, valg, valc, yout)
!yvalue_w1 = interpol_yp(y2, y1, vale, vala, yout)
!yvalue_t0 = interpol_yp(x2, x1, valh, valg, xout)
!yvalue_t1 = interpol_yp(x2, x1, vald, valc, xout)
!yvalue_b0 = interpol_yp(x2, x1, valf, vale, xout)
!yvalue_b1 = interpol_yp(x2, x1, valb, vala, xout)
!
!always two possible solutions on surfaces:
!sols1 = interpol_yp(x2, x1, yvalue_s0, yvalue_s1, xout)
!sols2 = interpol_yp(z2, z1, yvalue_t1, yvalue_b1, zout)
!soln1 = interpol_yp(x2, x1, yvalue_n0, yvalue_n1, xout)
!soln2 = interpol_yp(z2, z1, yvalue_t0, yvalue_b0, zout)
!solo1 = interpol_yp(y2, y1, yvalue_n0, yvalue_s0, yout)
!solo2 = interpol_yp(z2, z1, yvalue_o0, yvalue_o1, zout)
!solw1 = interpol_yp(y2, y1, yvalue_n1, yvalue_s1, yout)
!solw2 = interpol_yp(z2, z1, yvalue_w0, yvalue_w1, zout)
!solt1 = interpol_yp(y2, y1, yvalue_t0, yvalue_t1, yout)
!solt2 = interpol_yp(x2, x1, yvalue_o0, yvalue_w0, xout)
!solb1 = interpol_yp(y2, y1, yvalue_b0, yvalue_b1, yout)
!solb2 = interpol_yp(x2, x1, yvalue_o1, yvalue_w1, xout)
!
!top-bottom
!yout1 = interpol_yp(z2, z1, solt1, solb1, zout)
!north-south
!yout2 = interpol_yp(y2, y1, soln1, sols1, yout)
!east-west
!yout3 = interpol_yp(x2, x1, solo1, solw1, xout)
!
!   write(*,*) yout1, yout2, yout3
!
!
!********standard routine: interpolation in log-space if possible*******
!
if(llogf) then
!
!prepare input values for log-* interpolation
   dum_vala=log10(vala)
   dum_valb=log10(valb)
   dum_valc=log10(valc)
   dum_vald=log10(vald)
   dum_vale=log10(vale)
   dum_valf=log10(valf)
   dum_valg=log10(valg)
   dum_valh=log10(valh)
!
   if(llogz) then
!log-log interpolation in z-direction
      dum_z1=log10(abs(z1))
      dum_z2=log10(abs(z2))
      dum_zout=log10(abs(zout))
    else
!log-lin interpolation in z-direction
      dum_z1=z1
      dum_z2=z2
      dum_zout=zout
   endif
!
   if(llogx) then
!log-log-interpolation in x-direction
      dum_x1=log10(abs(x1))
      dum_x2=log10(abs(x2))
      dum_xout=log10(abs(xout))
   else
!log-lin-interpolation in x-direction
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

!perform interpolation in z-direction
   yvalue_s0 = interpol_yp(dum_z2, dum_z1, dum_vald, dum_valb, dum_zout)
   yvalue_s1 = interpol_yp(dum_z2, dum_z1, dum_valc, dum_vala, dum_zout)
   yvalue_n0 = interpol_yp(dum_z2, dum_z1, dum_valh, dum_valf, dum_zout)
   yvalue_n1 = interpol_yp(dum_z2, dum_z1, dum_valg, dum_vale, dum_zout)

!perform interpolation in x-direction
   soln1 = interpol_yp(dum_x2, dum_x1, yvalue_n0, yvalue_n1, dum_xout)
   sols1 = interpol_yp(dum_x2, dum_x1, yvalue_s0, yvalue_s1, dum_xout)
!
!perform interpolation in y-direction
   yinterp = interpol_yp(dum_y2, dum_y1, soln1, sols1, dum_yout)
   yinterp = 10.d0**yinterp
!
else
!lin-lin-interpolation
!
   yvalue_s0 = interpol_yp(z2, z1, vald, valb, zout)
   yvalue_s1 = interpol_yp(z2, z1, valc, vala, zout)
   yvalue_n0 = interpol_yp(z2, z1, valh, valf, zout)
   yvalue_n1 = interpol_yp(z2, z1, valg, vale, zout)
!
   soln1 = interpol_yp(x2, x1, yvalue_n0, yvalue_n1, xout)
   sols1 = interpol_yp(x2, x1, yvalue_s0, yvalue_s1, xout)
!
   yinterp = interpol_yp(y2, y1, soln1, sols1, yout)
!
endif
!
!
end subroutine trilin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_rtp_indx(rin, tin, pin, rarr, tarr, parr, nr, nt, np, indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, expol, rmin, rmax)
!
!-----------------------------------------------------------------------
!
!  finds indices of spherical grid for the grid-cell surrounding the point (rin,tin,pin)
!
!  input:  coordinates of point:             rin, tin, pin
!          dimension of grid:                nr, nt, np
!          x,y,z-grid:                       rarr, tarr, parr
!          boundaries of info-region:        rmin, rmax
!  output: indices of x-grid:                indx_r1, indx_r2
!          indices of y-grid:                indx_t1, indx_t2
!          indices of z-grid:                indx_p1, indx_p2
!          flag if extrapolation is needed:  expol


!

!
! ... arguments
integer(i4b), intent(in) :: nr, nt, np
real(dp), dimension(nr), intent(in) :: rarr
real(dp), dimension(nt), intent(in) :: tarr
real(dp), dimension(np), intent(in) :: parr
real(dp), intent(in) :: rin, tin, pin
real(dp), intent(in) :: rmin, rmax
integer(i4b) :: indx_r1, indx_t1, indx_p1, indx_r2, indx_t2, indx_p2
logical :: expol
!
! ... local scalars
integer(i4b) :: i
logical :: linfo, linfo2
!
! ... local functions
!
!
!
expol=.false.
!
!---------------------------hunt down theta-arr-------------------------
!
if(tarr(nt).lt.tin) then
   expol=.true.
   indx_t2=nt
   indx_t1=nt-1
elseif(tarr(nt).eq.tin) then
   indx_t2=nt
   indx_t1=nt-1   
elseif(tarr(1).gt.tin) then
   expol=.true.
   indx_t2=2
   indx_t1=1
elseif(tarr(1).eq.tin) then
   indx_t2=2
   indx_t1=1   
else   
   do i=2, nt
      if(tarr(i).ge.tin) then
         indx_t2=i
         indx_t1=i-1
         exit
      endif
   enddo
endif
!
!-------------------------hunt down phi-arr-----------------------------
!
if(parr(np).lt.pin) then
   expol=.true.
   indx_p2=np
   indx_p1=np-1
elseif(parr(np).eq.pin) then
   indx_p2=np
   indx_p1=np-1
elseif(parr(1).gt.pin) then
   expol=.true.
   indx_p2=2
   indx_p1=1
elseif(parr(1).eq.pin) then
   indx_p2=2
   indx_p1=1   
else
   do i=2, np
      if(parr(i).ge.pin) then
         indx_p2=i
         indx_p1=i-1
         exit
      endif
   enddo
endif
!
!------------------------hunt down radius-arr---------------------------
!
if(rarr(nr).lt.rin) then
   expol=.true.
   indx_r2=nr
   indx_r1=nr-1
elseif(rarr(nr).eq.rin) then
   indx_r2=nr
   indx_r1=nr-1   
elseif(rarr(1).gt.rin) then
   expol=.true.
   indx_r2=2
   indx_r1=1
elseif(rarr(1).eq.rin) then
   indx_r2=2
   indx_r1=1   
else
   do i=2, nr
      if(rarr(i).gt.rin) then
         indx_r2=i
         indx_r1=i-1
         exit
      endif
   enddo
endif
!
!
!
!
end subroutine get_rtp_indx


!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_rtp_values1(r_coord, t_coord, p_coord, rarr, tarr, parr, nr, nt, np, &
                          indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                          r1, r2, t1, t2, p1, p2, &
                          llogr, llogt, llogp)
!
!-----------------------------------------------------------------------
!
!   get coordinates of a spherical cell surrounding grid point (r_coord, t_coord, p_coord)
!   with vertices given by indx_r1, ... indx_p2
!
!input: dimension of arrays r,theta,phi: nr, nt, np
!       arrays r,theta,phi
!       coordinates of a point inside the cell: r_coord, t_coord, p_coord
!       indices of cube-vertices: indx_r1, ... indx_p2
!       logical to decide if interpolation or extrapolation: expol
!
!output: grid value: r1, r2, t1, t2, p1, p2
!        flags to decide if logarithmic interpolation is allowed:
!               llogr, llogt, llogp
!
!-----------------------------------------------------------------------
!


!

!
! ... arguments
integer(i4b), intent(in) :: nr, nt, np
integer(i4b), intent(in) :: indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2
real(dp), intent(in) :: r_coord, t_coord, p_coord
real(dp), dimension(nr), intent(in) :: rarr
real(dp), dimension(nt), intent(in) :: tarr
real(dp), dimension(np), intent(in) :: parr
real(dp) :: r1, r2, t1, t2, p1, p2
logical :: llogr, llogt, llogp
!
! ... local scalars
!
! ... local functions
!
!-----------------------------------------------------------------------
!
r1=rarr(indx_r1)
r2=rarr(indx_r2)
t1=tarr(indx_t1)
t2=tarr(indx_t2)
p1=parr(indx_p1)
p2=parr(indx_p2)
!
!-------------check if logarithmic interpolation is allowed-------------
!
!all coordinates need the same sign
!
!r always larger than 0
llogr=.true.
!
if(t1*t2*t_coord.eq.0.) then
   llogt=.false.
else
   llogt=.true.
endif
!
if(p1*p2*p_coord.eq.0.) then
   llogp=.false.
else
   llogp=.true.
endif
!
!
!
end subroutine get_rtp_values1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_rtp_values2(nr, nt, np, yvalue3d, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                           vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!
!-----------------------------------------------------------------------
!
!   get physical quantities of a spherical grid-cell
!   with vertices given by indx_r1, ... indx_p2
!
!input: dimension of 3-d array: nr, nt, np
!       physical value on grid: yvalue3d
!       indices of cube-vertices: indx_r1, ... indx_p2
!
!output: physical value on cell vertices: vala ... valh
!        flag to decide if logarithmic interpolation is allowed:
!               llogf
!
!-----------------------------------------------------------------------
!

!

!
! ... arguments
integer(i4b), intent(in) :: nr, nt, np
integer(i4b), intent(in) :: indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2
real(dp), dimension(nr, nt, np), intent(in) :: yvalue3d
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh
logical :: llogf
!
! ... local scalars
!
! ... local function
!
!-----------------------------------------------------------------------
!
vala = yvalue3d(indx_r1, indx_t1, indx_p1)
valb = yvalue3d(indx_r2, indx_t1, indx_p1)
valc = yvalue3d(indx_r1, indx_t1, indx_p2)
vald = yvalue3d(indx_r2, indx_t1, indx_p2)
vale = yvalue3d(indx_r1, indx_t2, indx_p1)
valf = yvalue3d(indx_r2, indx_t2, indx_p1)
valg = yvalue3d(indx_r1, indx_t2, indx_p2)
valh = yvalue3d(indx_r2, indx_t2, indx_p2)
!
if(vala.le.0.d0.or.valb.le.0.d0.or.valc.le.0.d0.or.vald.le.0.d0.or.vale.le.0.d0.or.valf.le.0.d0.or.valg.le.0.d0.or.valh.le.0.d0) then
   llogf=.false.
else
   llogf=.true.
endif
!
end subroutine get_rtp_values2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
end module mod_interp3d
