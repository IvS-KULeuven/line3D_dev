function interpol1d_2p_lin, f_im1, f_i, x_im1, x_i, x_p
;
;         interpolates values given on a 1d grid onto point x_p
;                 using linear interpolation
;
;   f(x) = f_im1 + (f_i-f_im1)*(x_p-x_im1)/(x_j-x_im1)
;
;on input: 
;
;          f_im1--------------f_i
;          x_im1--------------x_i
;                           
;        x_p: coordinates of point onto which shall be interpolated
;
;on output: interpolated value at x_p
;
;define deltax
dxi=x_i-x_im1
dx=x_p-x_im1
;
;define deltax-ratios
rdx=dx/dxi
;
;linear interpolation
acoeff=1.d0-rdx
bcoeff=rdx
;
return, acoeff*f_im1 + bcoeff*f_i
;
end
