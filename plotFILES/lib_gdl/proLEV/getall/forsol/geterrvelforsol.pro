PRO geterrvelforsol, dir, FILE_EXIST, np, nzeta, nz, errvx, errvy, errvz
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/forsol_errvx.dat'
fname2=dir+'/forsol_errvy.dat'
fname3=dir+'/forsol_errvz.dat'
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
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname3)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname3
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
errvx=FLTARR(np, nzeta, nz)*1.D0
errvy=FLTARR(np, nzeta, nz)*1.D0
errvz=FLTARR(np, nzeta, nz)*1.D0

OPENR, 1, fname1, /F77_UNFORMATTED
   READU, 1, errvx
CLOSE, 1
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, errvy
CLOSE, 1
;
OPENR, 1, fname3, /F77_UNFORMATTED
   READU, 1, errvz
CLOSE, 1
;
;
RETURN
;
END
