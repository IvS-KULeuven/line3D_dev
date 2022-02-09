PRO getzforsol, dir, FILE_EXIST, np, nzeta, nz, z_ray
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/forsol_zgrid.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
z_ray=FLTARR(np, nzeta, nz)*1.D0

OPENR, 1, fname1, /F77_UNFORMATTED
   READU, 1, z_ray
CLOSE, 1
;
;
RETURN
;
END
