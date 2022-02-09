pro interpol_trilinear, value, x, y, z, pp, value_pp
;
;interpolates any physical value (value) on a 3-d grid (x, y, z) onto 
;                a given point r on that grid
;
;
;              f_im1jk----------------f_ijk
;                 /|                   /|
;                / |                  / |
;               /  |                 /  |
;              /   |                /   |
;        f_im1jm1k--------------f_ijm1k |
;             |    |               |    |
;             |    |               |    |
;             |f_im1jkm1-----------|f_ijkm1
;             |   /                |   /
;             |  /                 |  /
;             | /                  | /
;             |/                   |/ 
;       f_im1jm1km1------------f_ijm1km1
;
;
;input:    value: 3-d array, from which values at a given point shall be
;                 interpolated
;          x, y, z: axes of 3-d grid
;          pp: point p, at which value shall be calculated
;
;output:   value_pp: interpolated value at point p
;
;-------------------find dimensions of 3-d grid-------------------------
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
ndzmax=n_elements(z)
;
;--------------------check for dimensions of pp-------------------------
;
if(n_elements(pp) ne 3) then begin
   print, 'pp has number of elements ne 3'
   stop
endif
;
;-----------------------------------------------------------------------
;
find_indx, pp(0), x, ndxmax, i, im1
find_indx, pp(1), y, ndymax, j, jm1
find_indx, pp(2), z, ndzmax, k, km1
;
;-----------------------------------------------------------------------
;
f_im1jm1km1 = value(im1,jm1,km1)
f_ijm1km1 = value(i,jm1,km1)
f_im1jkm1 = value(im1,j,km1)
f_ijkm1 = value(i,j,km1)
f_im1jm1k = value(im1,jm1,k)
f_ijm1k = value(i,jm1,k)
f_im1jk = value(im1,j,k)
f_ijk = value(i,j,k)
;
delx=x(im1)-x(i)
dely=y(jm1)-y(j)
delz=z(km1)-z(k)
;
acoeff=(f_im1jk-f_ijk)/delx
bcoeff=(f_ijm1k-f_ijk)/dely
ccoeff=(f_ijkm1-f_ijk)/delz
dcoeff=(f_im1jm1k-f_im1jk-f_ijm1k+f_ijk)/delx/dely
ecoeff=(f_im1jkm1-f_im1jk-f_ijkm1+f_ijk)/delx/delz
fcoeff=(f_ijm1km1-f_ijm1k-f_ijkm1+f_ijk)/dely/delz
gcoeff=(f_im1jm1km1-f_im1jm1k-f_im1jkm1-f_ijm1km1+f_im1jk+f_ijm1k+f_ijkm1-f_ijk)/delx/dely/delz
hcoeff=f_ijk
;
delx=pp(0)-x(i)
dely=pp(1)-y(j)
delz=pp(2)-z(k)
;
value_pp = hcoeff + acoeff*delx + bcoeff*dely + ccoeff*delz + $
           dcoeff*delx*dely + ecoeff*delx*delz + fcoeff*dely*delz + gcoeff*delx*dely*delz
;
;
;
end
