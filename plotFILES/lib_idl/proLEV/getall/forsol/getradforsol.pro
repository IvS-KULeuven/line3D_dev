PRO getRADFORSOL, dir, FILE_EXIST, np, nzeta, nz, rad_ray
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/forsol_rad.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
rad_ray=FLTARR(np, nzeta, nz)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, rad_ray
CLOSE, 1
;
RETURN
;
END
