pro interpol_biquad, value, x, y, pp, value_pp
;
;interpolates any physical value (value) on a 2-d grid (x, y) onto 
;                a given point pp
;
;input:    value: 2-d array, from which values at a given point shall be
;                 interpolated
;          x, y:  axes of 2-d grid
;          pp: point p, at which value shall be calculated
;
;output:   value_pp: interpolated value at point p
;
;-------------------find dimensions of 3-d grid-------------------------
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
;
;--------------------check for dimensions of pp-------------------------
;
if(n_elements(pp) ne 2) then begin
   print, 'pp has number of elements ne 2'
   stop
endif
;
;
find_indx, pp(0), x, ndxmax, i, im1
find_indx, pp(1), y, ndymax, j, jm1
;
if(i eq 1 or j eq 1) then begin
;linear interpolation
   value_pp = interpol2d_4p_lin(value(im1,jm1),value(i,jm1),value(im1,j),value(i,j), $
                                x(im1),x(i),y(jm1),y(j),pp(0),pp(1))
endif else begin
;quadratic interpolation
   im2=im1-1
   jm2=jm1-1
   value_pp = interpol2d_9p_quad(value(im2,jm2), value(im1,jm2), value(i,jm2), $
                                 value(im2,jm1), value(im1,jm1), value(i,jm1), $
                                 value(im2,j),   value(im1,j),   value(i,j), $
                                 x(im2), x(im1), x(i), y(jm2), y(jm1), y(j), $
                                 pp(0), pp(1))
endelse
;
end
