pro coeff3d_8p_lin, x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p, $
                    acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
;
;         interpolates values given on a 3d grid onto point x_p, y_p, z_p
;                      using triliniear interpolation
;
;         returns coefficients, such that
;            f(x,y,z) = a*f_im1jm1km1 + b*f_ijm1km1 + c*f_im1jkm1 + d*f_ijkm1 + 
;                       e*f_im1jm1j   + f*f_ijm1k   + g*f_im1jk   + h*f_ijk
;
;   f(x,y,z) = f_km1 + (f_k-f_km1)*(z_p-z_km1)/(z_k-z_km1)
;
;      with:   f_km1 = f_jm1km1 + (f_jkm1-f_jm1km1)*(y_p-y_im1)/(y_j-y_jm1)
;              f_k   = f_jm1    + (f_jk  -f_jm1k  )*(y_p-y_im1)/(y_j-y_jm1)
;              f_jk     = f_im1jkm1   + (f_ijkm1-f_im1jkm1) *   (x_p-x_im1)/(x_i-x_im1)
;              f_jm1km1 = f_im1jm1km1 + (f_ijm1km1-f_im1jm1km1)*(x_p-x_im1)/(x_i-x_im1)
;
;on input: 
;
;          f_im1jk-------------f_ijk          z_k
;            /|                 /|             |  y_j 
;           / |                / |             |   /
;          /  |               /  |             |  /
;         /   |              /   |             | /
;        /    |             /    |             |/
;       / f_im1jm1---------/--f_ijkm1        x_im1--------------x_i
;      /     /  x(x,y,z)  /     /            y_jm1
; f_im1jm1k-----------f_ijm1k  /             z_km1
;     |    /             |    /
;     |   /              |   /
;     |  /               |  /
;     | /                | /
;     |/                 |/
;f_im1jm1km1---------f_ijm1km1
;                    
;        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
;
;on output: coefficients
;
;define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
dzk=z_k-z_km1
dz=z_p-z_km1
;
;define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdz=dz/dzk
rdxdy=rdx*rdy
;
;trilinear interpolation
fdum1=1.d0-rdx-rdy+rdxdy
fdum2=rdx-rdxdy
fdum3=rdy-rdxdy
fdum4=rdxdy
;
ecoeff=fdum1*rdz
fcoeff=fdum2*rdz
gcoeff=fdum3*rdz
hcoeff=fdum4*rdz
acoeff=fdum1-ecoeff
bcoeff=fdum2-fcoeff
ccoeff=fdum3-gcoeff
dcoeff=fdum4-hcoeff
;
;
end
