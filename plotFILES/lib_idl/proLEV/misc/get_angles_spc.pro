pro get_angles_spc, x, y, z, theta, phi
;
;--------------calculates angles in spherical coordinates---------------
;---------------for given carthesian coordinates x, y, z----------------
;----------------------for right-handed system--------------------------
;
;on input: coordinates x, y, z
;
;on output: angles theta, phi
;
;(from jon's formal solver)
;
;-----------------------calculation of theta----------------------------
;
rad = sqrt(x^2+y^2+z^2)
;
if(rad ne 0.) then begin
   theta = acos(z/rad)
endif else begin
   theta = 0.d0
endelse
;
;------------------------calculation of phi-----------------------------
;
if(z eq rad) then begin
   phi=0.d0
   return
endif
;
if(x eq 0.) then begin
   if(y eq 0.d0) then begin
      phi = 0.d0
   endif else begin
      if(y gt 0.) then begin
         phi = !pi/2.d0
      endif else begin
         phi = 3.d0*!pi/2.d0
      endelse
   endelse
endif else begin
   if(y gt 0.d and x gt 0.) then begin
;first quadrant
      phi = atan(y/x)
   endif else begin
      if (x gt 0.) then begin
;fourth quadrant
         phi = 2.d0*!pi + atan(y/x)
      endif else begin
;second and third quadrant
         phi = !pi + atan(y/x)
      endelse
   endelse
endelse

return
;
;-------------------------test if everything worked---------------------
;
x_test=rad*sin(theta)*cos(phi)
y_test=rad*sin(theta)*sin(phi)
z_test=rad*cos(theta)
;
if(abs(x_test-x) gt 1.d-10) then begin
   print, 'error in get_angles_spc', x_test, x
   stop
endif
;
if(abs(y_test-y) gt 1.d-10) then begin
   print, 'error in get_angles_spc', y_test, y
   stop
endif
;
if(abs(z_test-z) gt 1.d-10) then begin
   print, 'error in get_angles_spc', z_test, z
   stop
endif
;
;
end
