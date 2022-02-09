function interpol2d_4p_idw, f_im1jm1, f_ijm1, f_im1j, f_ij, x_im1, x_i, y_jm1, y_j, x_p, y_p
;
;         interpolates values given on a 2d grid onto point x_p, y_p
;                 using inverse distance weighting with four neighbouring points
;
;   f(x,y) = (w_ij*f_ij + w_im1j*f_im1j + w_ijm1*f_ijm1 + w_im1jm1*f_im1jm1) /
;          (w_ij+w_im1j+w_ijm1+w_im1jm1)
;
;      with: w_ij = 1.d0/sqrt((x_p-x_i)^2+(y_p-y_j)^2)^p
;            p - weighting factor
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
;define weighting factor
fac=0.5d0
;
;calculate distance to each point
d_im1jm1 = (x_p-x_im1)^2+(y_p-y_jm1)^2
d_ijm1 = (x_p-x_i)^2+(y_p-y_jm1)^2
d_im1j = (x_p-x_im1)^2+(y_p-y_j)^2
d_ij = (x_p-x_i)^2+(y_p-y_j)^2
;
;avoid division by zero if d_ij=0, or d_im1j=0, ...
if(d_im1jm1 eq 0.) then return, f_im1jm1
if(d_ijm1 eq 0.) then return, f_ijm1
if(d_im1j eq 0.) then return, f_im1j
if(d_ij eq 0.) then return, f_ij
;
w_im1jm1=1.d0/d_im1jm1^(fac/2.d0)
w_ijm1=1.d0/d_ijm1^(fac/2.d0)
w_im1j=1.d0/d_im1j^(fac/2.d0)
w_ij=1.d0/d_ij^(fac/2.d0)
;
return, (f_im1jm1*w_im1jm1 + f_ijm1*w_ijm1 + f_im1j*w_im1j + f_ij*w_ij) / $
        (w_im1jm1 + w_ijm1 + w_im1j + w_ij)


end
