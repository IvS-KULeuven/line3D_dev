PRO CALC_TRANSMAT, nx, ny, nz, transmat, help=print_help
;
;+
; NAME:
;	CALC_TRANSMAT
;
; PURPOSE:
;	This procedure calculates a transformation matrix to transform
;       a cylindrical coordinate system (e_x, e_y, e_z) to 
;       a carthesian coordinate system (ee_x, ee_y, ee_z).
;       The cylindrical coordinate system is specfied by a direction
;       vector defining e_z (has to be normalized).
;
; CALLING SEQUENCE:
;	CALC_TRANSMAT, nx, ny, nz, transmat
;
; INPUTS:
;	nx:   	x-component in carthesian coordinates of input direction vector
;	ny:   	y-component in carthesian coordinates of input direction vector
;	nz:   	z-component in carthesian coordinates of input direction vector
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	transformation-matrix for coordinate transformation, transmat
;
; EXAMPLE:
;	CALC_TRANSMAT, 0.5, 0.5, 0.707, transmat
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'calc_transmat'
   return
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
nhat=[nx, ny, nz]
;
transmat=fltarr(3,3)*0.d0
;
PRINT, '----------CALCULATE TRANSFORMATION MATRIX--------------'
PRINT, 'n= ', nx, ny, nz
;
if(abs(nhat(0)) lt 1.d-14) then nhat(0)=0.d0
if(abs(nhat(1)) lt 1.d-14) then nhat(1)=0.d0
if(abs(nhat(2)) lt 1.d-14) then nhat(2)=0.d0
;
;CALCULATE CARTHESIAN COORDINATE SYSTEM THAT IS ALIGNED WITH 
;COORDINATE-SYSTEM ALONG Z_NEW= (nx, ny, nz)
;E.G. CYLINDRICAL COORDINATES (SET BY Z_CYLINDRICAL = (nx, ny, nz)
;
;ez predefined: direction to observer
ez=nhat
;
;calculate orthogonal ex from dot_product(ex,ez)=0
if(ez(0) eq 0.d0) then begin
   if(ez(1) eq 0.d0) then begin
      ex = [ 1.d0, 0.d0, 0.d0 ]
   endif else begin
      if(ez(2) eq 0.d0) then begin
         ex = [ 1.d0, 0.d0, 0.d0 ]
      endif else begin
         ex(0) = 1.d0
         ex(1) = 1.d0
         ex(2) = -ex(1)*ez(1)/ez(2)
      endelse
   endelse
endif else begin
   if (ez(1) eq 0.d0) then begin
      if(ez(0) eq 0.d0) then begin
         ex = [ 1.d0, 0.d0, 0.d0 ]
      endif else begin
         if(ez(2) eq 0.d0) then begin
            ex = [ 0.d0, 1.d0, 0.d0 ]
         endif else begin
            ex(2) = 1.d0
            ex(1) = 1.d0
            ex(0) = -ex(2)*ez(2)/ez(0)
         endelse
      endelse
   endif else begin
      if (ez(2) eq 0.d0) then begin
         if(ez(0) eq 0.d0) then begin
            ex = [ 1.d0, 0.d0, 0.d0 ]
         endif else begin
            if(ez(1) eq 0.d0) then begin
               ex = [ 0.d0, 1.d0, 0.d0 ]
            endif else begin
               ex(2) = 1.d0
               ex(1) = 1.d0
               ex(0) = -ex(1)*ez(1)/ez(0)
            endelse
         endelse
      endif else begin
         ex(0) = 1.d0
         ex(1) = 1.d0
         ex(2) = (-ex(0)*ez(0)-ex(1)*ez(1))/ez(2)
      endelse
   endelse
endelse
;
;calculate orthogonal ey from cross-product
ey = CROSSP(ez, ex)
;
;normalize unit vectors
ex = ex/sqrt(dot_product(ex,ex))
ey = ey/sqrt(dot_product(ey,ey))
ez = ez/sqrt(dot_product(ez,ez))
;
;check for orthogonality
check1=dot_product(ex,ey)
check2=dot_product(ex,ez)
check3=dot_product(ey,ez)
;
if(abs(check1) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT: ex, ey not orthogonal'
   stop
endif
if(abs(check2) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT: ex, ez not orthogonal'
   stop
endif
if(abs(check3) gt 1.d-14) then begin
   print, 'ERROR IN CALC_TRANSMAT: ey, ez not orthogonal'
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
print, ex
print, ey
print, ez
print, transmat
;
;transmat#vec = vex_x*ex + vex_y*ey + vec_z*ez
;print, ' '
;print, nhat
;
END
