PRO getVEL, dir, FILE_EXIST, NDZMAX, NDXMAX, VELZ2D, VELX2D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/VELZ_IDL.dat'
fname2=dir+'/VELX_IDL.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
VELZ2D=FLTARR(NDZMAX, NDXMAX)*1.D0
VELX2D=FLTARR(NDZMAX, NDXMAX)*1.D0

OPENR, 1, fname1, /F77_UNFORMATTED
   READU, 1, VELZ2D
CLOSE, 1
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, VELX2D
CLOSE, 1
;
RETURN
;
END
