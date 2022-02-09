PRO angles_spc, x, y, z, theta, phi, help=print_help
;
;+
; NAME:
;	angles_spc
;
; PURPOSE:
;	This procedure calculates the latitude and polar angle of a vector
;       given in carthesian coordinates.
;
; CALLING SEQUENCE:
;	angles_spc, x, y, z, theta, phi
;
; INPUTS:
;	x:	x-coordinate of vector
;	y:	y-coordinate of vector
;	z:	z-coordinate of vector
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	theta:   Latitude-angle of vector
;       phi:     Azimuth-angle of vector
;
; EXAMPLE:
;	angles_spc, 1., 1., 1., theta, phi
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'angles_spc'
   return
endif
;
;-----------------------CALCULATION OF THETA----------------------------
;
rad = sqrt(x*x+y*y+z*z)
;
if(rad ne 0.d0) then begin
   theta = acos(z/rad)
endif else begin
   theta = 0.d0
endelse
;
;------------------------CALCULATION OF PHI-----------------------------
;
if(z eq rad) then begin
   phi=0.d0
   return
endif
;
if(x eq 0.d0) then begin
   if(y eq 0.d0) then begin
      phi = 0.d0
   endif else begin
      if(y gt 0.d0) then begin
         phi = !PI/2.d0
      endif else begin
         phi = 3.d0*!PI/2.d0
      endelse
   endelse
endif else begin
   if(y gt 0.d0 and x gt 0.d0) then begin
;first quadrant
      phi = atan(y/x)
   endif else begin
      if (x gt 0.d0) then begin
;fourth quadrant
         phi = 2.d0*!PI + atan(y/x)
      endif else begin
;second and third quadrant
         phi = !PI + atan(y/x)
      endelse
   endelse
endelse
;
;-------------------------TEST IF EVERYTHING WORKED---------------------
;
x_test=rad*sin(theta)*cos(phi)
y_test=rad*sin(theta)*sin(phi)
z_test=rad*cos(theta)
;
if(abs(x_test-x) gt 1.d-5) then begin
   print, 'ERROR IN GET_ANGLES_SPC', x_test, x
   stop
endif
;
if(abs(y_test-y) gt 1.d-5) then begin
   print, 'ERROR IN GET_ANGLES_SPC', y_test, y
   STOP
endif
;
if(abs(z_test-z) gt 1.d-5) then begin
   print, 'ERROR IN GET_ANGLES_SPC', z_test, z
   stop
endif
;
;
;
end







ex=fltarr(3)*0.d0
ey=fltarr(3)*0.d0
ez=fltarr(3)*0.d0
;
eex=fltarr(3)*0.d0
eey=fltarr(3)*0.d0
eez=fltarr(3)*0.d0
;
transmat=fltarr(3,3)*0.d0
;
PRINT, '----------CALCULATE TRANSFORMATION MATRIX--------------'
PRINT, 'alpha= ', alpha
PRINT, 'gamma= ', gamma
;
aalpha=alpha*!PI/180.d0
ggamma=gamma*!PI/180.d0
;
;calculate unit vectors of new coordinate system: ez=e_r, ex=e_phi, ey=e_theta
ex=[ -1.d0*sin(ggamma), cos(ggamma), 0.d0 ]
ey=[ cos(aalpha)*cos(ggamma), cos(aalpha)*sin(ggamma), -1.d0*sin(aalpha) ]
ez=[ sin(aalpha)*cos(ggamma), sin(aalpha)*sin(ggamma), cos(aalpha)]
;
if(abs(ex(0)) lt 1.d-14) then ex(0)=0.d0
if(abs(ex(1)) lt 1.d-14) then ex(1)=0.d0
if(abs(ex(2)) lt 1.d-14) then ex(2)=0.d0
if(abs(ey(0)) lt 1.d-14) then ey(0)=0.d0
if(abs(ey(1)) lt 1.d-14) then ey(1)=0.d0
if(abs(ey(2)) lt 1.d-14) then ey(2)=0.d0
if(abs(ez(0)) lt 1.d-14) then ez(0)=0.d0
if(abs(ez(1)) lt 1.d-14) then ez(1)=0.d0
if(abs(ez(2)) lt 1.d-14) then ez(2)=0.d0
;
;
;check for orthogonality
check1=dot_product(ex,ey)
check2=dot_product(ex,ez)
check3=dot_product(ey,ez)
;
if(abs(check1) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT2: ex, ey not orthogonal'
   stop
endif
if(abs(check2) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT2: ex, ez not orthogonal'
   stop
endif
if(abs(check3) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT2: ey, ez not orthogonal'
   stop
endif
;
;----TRANSFORMATION MATRIX FROM SYSTEM (ex,ey,ez) TO (eex, eey, eez)----
;
eex = [ 1.d0, 0.d0, 0.d0 ]
eey = [ 0.d0, 1.d0, 0.d0 ]
eez = [ 0.d0, 0.d0, 1.d0 ]
;
transmat = [ [dot_product(ex,eex), dot_product(ex,eey), dot_product(ex,eez)] , $
             [dot_product(ey,eex), dot_product(ey,eey), dot_product(ey,eez)] , $
             [dot_product(ez,eex), dot_product(ez,eey), dot_product(ez,eez)] ]
;

;
print, 'unit vectors:'
print, 'e_x_slice', ex
print, 'e_y_slice', ey
print, 'e_z_slice', ez
print, ' '
print, 'transformation matrix:'
print, transmat
;
;transmat#vec = vex_x*ex + vex_y*ey + vec_z*ez
;print, ' '
;print, nhat
;
END
