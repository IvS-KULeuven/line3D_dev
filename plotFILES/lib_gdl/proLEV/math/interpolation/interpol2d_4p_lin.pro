function interpol2d_4p_lin, f_im1jm1, f_ijm1, f_im1j, f_ij, x_im1, x_i, y_jm1, y_j, x_p, y_p
;
;         interpolates values given on a 2d grid onto point x_p, y_p
;                 using biliniear interpolation
;
;   f(x,y) = f_jm1 + (f_j-f_jm1)*(y_p-y_im1)/(y_j-y_jm1)
;
;      with: f_j   = f_im1j   + (f_ij-f_im1j) *   (x_p-x_im1)/(x_i-x_im1)
;            f_jm1 = f_im1jm1 + (f_ijm1-f_im1jm1)*(x_p-x_im1)/(x_i-x_im1)
;
;on input: 
;
;          f_im1j--------------f_ij          y_j
;            |                  |             |
;            |                  |             |
;            |         x        |             |
;            |     (x_p,y_p)    |             |
;            |                  |             |
;        f_im1jm1------------f_ijm1         x_im1--------------x_i
;                                           y_jm1
;                           
;        x_p, y_p: coordinates of point onto which shall be interpolated
;
;on output: interpolated value at x_p, y_p
;
;define deltax, deltay
dxi=x_i-x_im1
dx=x_p-x_im1
dyj=y_j-y_jm1
dy=y_p-y_jm1
;
;define deltax, deltay-ratios
rdx=dx/dxi
rdy=dy/dyj
rdxdy=rdx*rdy
;
;bilinear interpolation
acoeff=1.d0-rdx-rdy+rdxdy
bcoeff=rdx-rdxdy
ccoeff=rdy-rdxdy
dcoeff=rdxdy
;
return, acoeff*f_im1jm1 + bcoeff*f_ijm1 + ccoeff*f_im1j + dcoeff*f_ij
;
end
