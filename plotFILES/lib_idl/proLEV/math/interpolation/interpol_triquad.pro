function interp1d_lin, x_im1, x_i, f_im1, f_i, x_p
;
;1d linear interpolation:
;
;   x_im1-----x_p------x_i
;
;
f_p = f_im1 + (f_i-f_im1)*(x_p-x_im1)/(x_i-x_im1)
;
return, f_p
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
function interp1d_quad, x_im2, x_im1, x_i, f_im2, f_im1, f_i, x_p
;
;1d (monotonic) quadratic interpolation:
;
;   x_im2--------x_im1-----x_p------x_i
;
dxi=x_i-x_im1
dxim1=x_im1-x_im2
dfi=f_i-f_im1
dfim1=f_im1-f_im2
;
a = (dfi/dxi - dfim1/dxim1)/(dxi+dxim1)
b = dfim1/dxim1 + a*dxim1
c = f_im1
;
dx = x_p-x_im1
;
;ensure monotonicity in [x_im1,x_i]
;fp_im1=b
;fp_i=2.d0*a*dxi+b
;if(fp_im1*fp_i lt 0.) then begin
;   return,  f_im1 + dfi*dx/dxi
;endif else begin
   return, a*dx^2 + b*dx + c
;endelse
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
function interp1d_bez, x_im2, x_im1, x_i, f_im2, f_im1, f_i, x_p
;
;         interpolates values given on a 1d-axis onto point x_p
;                 using monotonic quadratic bezier spline
;
;   f(t) = f_im1 * (1-t)^2 + f_i * t^2 + 2*f_c * t*(1-t)
;      t = (x-x_im1)/(x_i-x_im1)
;    f_c = f_im1 + dfdx_im1*(x_i-x_im1)/2.d0
;
;   control point is calculated at the center of x(i) and x(i-1)
;   first derivative is calculated from weighted mean (see Hayek 2010)
;      (of both backward differences)
;   three point interpolation
;
;on input: 
;
;                    f_im2----------f_im1----.-----f_i
;                    x_im2----------x_im1---x_p----x_i
;
;        x_p: coordinate of point onto which shall be interpolated
;
;on output: interpolated value at x_p
;
;
;derivative at point x_im1 from weighted mean
fp_im1 = (f_im1-f_im2)*(x_i-x_im1)/(x_im1-x_im2)/(x_i-x_im2) + (f_i-f_im1)*(x_im1-x_im2)/(x_i-x_im1)/(x_i-x_im2)

;dx_i = x_i-x_im1
;dx_im1=x_im1-x_im2
;ddx = dx_i+dx_im1
;fc_dum = -dx_i^2/2.d0/dx_im1/ddx * f_im2 + (1.d0+(dx_i^2-dx_im1^2)/2.d0/dx_im1/ddx)*f_im1 + dx_im1/2.d0/ddx*f_i
f_c = f_im1 + 0.5d0*(x_i-x_im1)*fp_im1
;print, fc_dum, f_c

f_min=min([f_im1, f_i])
f_max=max([f_im1, f_i])
if(f_c gt f_max) then begin
   f_c=f_max
endif else begin
   if (f_c lt f_min) then begin
      f_c=f_min
   endif
endelse
;
;interpolation
t = (x_p-x_im1)/(x_i-x_im1)
return, f_im1*(1.d0-t)^2 + f_i*t^2 + 2.d0*f_c*t*(1.d0-t)
;
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro interpol_triquad, value, x, y, z, pp, value_pp
;
;interpolates any physical value (value) on a 3-d grid (x, y, z) onto 
;                a given point r on that grid
;
;
;                         f_im2jk-------------------f_im1jk------------------f_ijk 
;                           /|                        /|                      /|
;                          / |                       / |                     / |
;                         /  |                      /  |                    /  |
;                        /   |                     /   |                   /   |
;                   f_im2jm1k-----------------f_im1jm1k----------------f_ijm1k |
;                      /|    |                   /|    |      pp         /|    |
;                     / |    |                  / |    |       x        / |    |
;                    /  | f_im2jkm1------------/--|f_im1jkm1-----------/--|-f_ijkm1
;                   /   |   /|                /   |   /|              /   |   /|
;              f_im2jm2k-----------------f_im1jm2k----------------f_ijm2k |  / |
;                  |    | /  |               |    | /  |             |    | /  |
;                  |    |/   |               |    |/   |             |    |/   |
;                  |f_im2jm1km1--------------|f_im1jm1km1------------|f_ijm1km1|
;                  |   /|    |               |   /|    |             |   /|    |
;                  |  / |    |               |  / |    |             |  / |    |
;                  | /  |f_im2jkm2-----------|-/--|f_im1jkm2---------|-/--|-f_ijkm2
;                  |/   |   /                |/   |   /              |/   |   /
;            f_im2jm2km1---------------f_im1jm2km1---------------f_ijm2jm1|  /
;                  |    | /                  |    | /                |    | /
;                  |    |/                   |    |/                 |    |/
;                  |f_im2jm1km2--------------|f_im1jm1km2------------|f_ijm1km2
;                  |   /                     |   /                   |   /
;                  |  /                      |  /                    |  /
;                  | /                       | /                     | /
;                  |/                        |/                      |/
;             f_im2jm2km2--------------f_im1jm2km2----------------f_ijm2km2
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
find_indx, pp(0), x, ndxmax, im1, i
find_indx, pp(1), y, ndymax, jm1, j
find_indx, pp(2), z, ndzmax, km1, k
;
;----------------trilinear interpolation if i=1 or j=1 or k=1-----------
;
if(i eq 1 or j eq 1 or k eq 1) then begin
;
;between (i-1,j,k), (i,j,k)
   value_t1 = interp1d_lin(x(im1),x(i),value(im1,j,k),value(i,j,k),pp(0))
;between (i-1,j-1,k), (i,j-1,k)
   value_t2 = interp1d_lin(x(im1),x(i),value(im1,jm1,k),value(i,jm1,k),pp(0))
;between (i-1,j,k-1), (i,j,k-1)
   value_b1 = interp1d_lin(x(im1),x(i),value(im1,j,km1),value(i,j,km1),pp(0))
;between (i-1,j-1,k-1), (i,j-1,j-1)
   value_b2 = interp1d_lin(x(im1),x(i),value(im1,jm1,km1),value(i,jm1,km1),pp(0))
;
;between bottom and top
   value_n = interp1d_lin(z(km1),z(k),value_b1,value_t1,pp(2))
   value_s = interp1d_lin(z(km1),z(k),value_b2,value_t2,pp(2))
;
;between south and nord
   value_pp = interp1d_lin(y(jm1),y(j),value_s,value_n,pp(1))
;
endif else begin
;
;-----------------------tricubic interpolation--------------------------
;
   im2=im1-1
   jm2=jm1-1
   km2=km1-1
;
;on level k along x
;between (i-2,j,k), (i-1,j,k), (i,j,k)
   value_t1 = interp1d_quad(x(im2),x(im1),x(i),value(im2,j,k),value(im1,j,k),value(i,j,k),pp(0))
;between (i-2,j-1,k), (i-1,j-1,k), (i,j-1,k)
   value_t2 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm1,k),value(im1,jm1,k),value(i,jm1,k),pp(0))
;between (i-2,j-2,k), (i-1,j-2,k), (i,j-2,k)
   value_t3 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm2,k),value(im1,jm2,k),value(i,jm2,k),pp(0))

;on level k-1 along x
;between (i-2,j,k-1), (i-1,j,k-1), (i,j,k-1)
   value_b1 = interp1d_quad(x(im2),x(im1),x(i),value(im2,j,km1),value(im1,j,km1),value(i,j,km1),pp(0))
;between (i-2,j-1,k-1), (i-1,j-1,k-1), (i,j-1,k-1)
   value_b2 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm1,km1),value(im1,jm1,km1),value(i,jm1,km1),pp(0))
;between (i-2,j-2,k-1), (i-1,j-2,k-1), (i,j-2,k-1)
   value_b3 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm2,km1),value(im1,jm2,km1),value(i,jm2,km1),pp(0))
;between (i-2,j,k-1), (i-1,j,k-1), (i,j,k-1)

;on level k-2 along x
;between (i-2,j,k-2), (i-1,j,k-2), (i,j,k-2)
   value_q1 = interp1d_quad(x(im2),x(im1),x(i),value(im2,j,km2),value(im1,j,km2),value(i,j,km2),pp(0))
;between (i-2,j-1,k-2), (i-1,j-1,k-2), (i,j-1,k-2)
   value_q2 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm1,km2),value(im1,jm1,km2),value(i,jm1,km2),pp(0))
;between (i-2,j-2,k-2), (i-1,j-2,k-2), (i,j-2,k-2)
   value_q3 = interp1d_quad(x(im2),x(im1),x(i),value(im2,jm2,km2),value(im1,jm2,km2),value(i,jm2,km2),pp(0))
;
;on level k along y
   value_t = interp1d_quad(y(jm2),y(jm1),y(j),value_t3,value_t2,value_t1,pp(1))
;on level k-1 along y
   value_b = interp1d_quad(y(jm2),y(jm1),y(j),value_b3,value_b2,value_b1,pp(1))
;on level k-2 along y
   value_q = interp1d_quad(y(jm2),y(jm1),y(j),value_q3,value_q2,value_q1,pp(1))
;
;along x
   value_pp = interp1d_quad(z(km2),z(km1),z(k),value_q,value_b,value_t,pp(2))
;
endelse
;
;
;
end
