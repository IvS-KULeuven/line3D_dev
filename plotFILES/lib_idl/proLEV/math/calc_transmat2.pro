PRO CALC_TRANSMAT2, alpha, gamma, transmat, help=print_help
;
;+
; NAME:
;	CALC_TRANSMAT2
;
; PURPOSE:
;	This procedure calculates a transformation matrix to transform
;       a cylindrical coordinate system (e_x, e_y, e_z) to 
;       a carthesian coordinate system (ee_x, ee_y, ee_z).
;       The cylindrical coordinate system is specfied by viewing angles.
;       
;
; CALLING SEQUENCE:
;	CALC_TRANSMAT2, alpha, gamma, transmat
;
; INPUTS:
;	alpha:	viewing angle with respect to the z-axis of carthesian
;     	        coordinates in radiant
;	gamma:	viewing angle with respect to the x-axis of carthesian   
; 	        coordinates (in x-y-plane) in radiant
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	transformation-matrix for coordinate transformation, transmat
;
; EXAMPLE:
;	CALC_TRANSMAT2, !PI/2., 0., transmat
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'calc_transmat2'
   return
endif

if(alpha gt !pi) then begin
   print, 'error in calc_transmat2: alpha needs to be in radiant, range [0,pi]'
   stop
endif

if(gamma gt 2.d0*!pi) then begin
   print, 'error in calc_transmat2: gamma needs to be in radiant range [0,2*pi]'
   stop
endif
;
;---------------CALCULATION OF TRANSFORMATION MATRIX--------------------
;--------FOR SYSTEM (ex, ey, ez) TO SYSTEM (eex, eey, eez)--------------
;
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
aalpha=alpha;*!PI/180.d0
ggamma=gamma;*!PI/180.d0
;
;calculate unit vectors of new coordinate system:  ez=e_r, ex=e_phi, ey=e_theta
;or, to get correct azimuthal velocity directions: ez=e_r, ex_e_theta, ey=e_phi
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
