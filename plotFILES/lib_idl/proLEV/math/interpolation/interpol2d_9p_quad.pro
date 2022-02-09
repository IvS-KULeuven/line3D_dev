function interpol2d_9p_quad, f_im2jm2, f_im1jm2, f_ijm2, $
                             f_im2jm1, f_im1jm1, f_ijm1, $
                             f_im2j, f_im1j, f_ij, $
                             x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
;
;         interpolates values given on a 2d grid onto point x_p, y_p
;                 using 2d quadratic functions
;
;   f(x,y) = a*(y_p-y_jm1)^2 + b*(y_p-y_jm1) + c
;
;      with: a, b, c given from interpolations on each j-level
;
;on input: 
;
; y_j      f_im2j----------f_im1j--------------f_ij
;  |          |              |                  |  
;  |          |              |                  |  
;  |          |              |         x        |  
;  |          |              |     (x_p,y_p)    |  
;  |          |              |                  |  
;y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
;  |          |              |                  |  
;  |          |              |                  |  
;  |          |              |                  |  
;  |          |              |                  |  
;  |          |              |                  |  
;y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
;  |
;  --------x_im2-----------x_im1---------------x_i
;                           
;        x_p, y_p: coordinates of point onto which shall be interpolated
;
;on output: interpolated value at x_p, y_p
;
;define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
;
dyjm1 = y_jm1-y_jm2
dyj = y_j-y_jm1
dy = dyjm1+dyj
;
tx = x_p-x_im1
ty = (y_p-y_jm1)
;
;derivatives at im1 on each j-level
a = -dxi/dxim1/dx
b = (dxi^2-dxim1^2)/dxi/dxim1/dx
c = dxim1/dxi/dx
fp_im1jm2 = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2
fp_im1jm1 = a*f_im2jm1 + b*f_im1jm1 + c*f_ijm1
fp_im1j   = a*f_im2j   + b*f_im1j   + c*f_ij  
;
;derivatives at i on each j-level
a2 = -a
b2 = b-2.d0/dxim1
c2 = c+2.d0/dx
fp_ijm2 = a2*f_im2jm2 + b2*f_im1jm2 + c2*f_ijm2
fp_ijm1 = a2*f_im2jm1 + b2*f_im1jm1 + c2*f_ijm1
fp_ij   = a2*f_im2j   + b2*f_im1j   + c2*f_ij
;
;interpolate on each j-level (and ensure monotonicity)
a3 = tx^2/dxim1/dx + a*tx
b3 = 1.d0-tx^2/dxi/dxim1 + b*tx
c3 = tx^2/dxi/dx + c*tx
if(fp_ij*fp_im1j lt 0.) then begin
   f_j = f_im1j + (f_ij-f_im1j)*tx/dxi   ;linear interpolation
endif else begin
   f_j =a3*f_im2j + b3*f_im1j + c3*f_ij
endelse
;
if(fp_ijm1*fp_im1jm1 lt 0.) then begin
   f_jm1 = f_im1jm1 + (f_ijm1-f_im1jm1)*tx/dxi  ;linear interpolation
endif else begin
   f_jm1 =a3*f_im2jm1 + b3*f_im1jm1 + c3*f_ijm1
endelse
;
if(fp_ijm2*fp_im1jm2 lt 0.) then begin
   f_jm2 = f_im1jm2 + (f_ijm2-f_im1jm2)*tx/dxi  ;linear interpolation
endif else begin
   f_jm2 =a3*f_im2jm2 + b3*f_im1jm2 + c3*f_ijm2
endelse
;
;now, interpolation along y-axis
a = -dyj/dyjm1/dy
b = (dyj^2-dyjm1^2)/dyj/dyjm1/dy
c = dyjm1/dyj/dy
;derivatives at j and jm1
fp_jm1 =  a*f_jm2 + b*f_jm1 + c*f_j
fp_j   = -a*f_jm2 + (b-2.d0/dyjm1)*f_jm1 + (c+2.d0/dy)*f_j
if(fp_j*fp_jm1 lt 0.) then begin
   return, f_jm1 + (f_j-f_jm1)*ty/dyj   ;linear interpolation
endif else begin
   return, (ty^2/dyjm1/dy + a*ty)*f_jm2 + (1.d0-ty^2/dyj/dyjm1 + b*ty)*f_jm1 + (ty^2/dyj/dy + c*ty)*f_j
endelse



end
